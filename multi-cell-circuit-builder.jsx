import { useState, useEffect, useCallback, useRef, useMemo } from "react";

/* ═══════════════════════════════════════════════
   CONSTANTS & CONFIG
   ═══════════════════════════════════════════════ */
const C = {
  bg: "#060a10", panel: "#0c1018", panel2: "#111922", border: "#1a2535", borderHi: "#2a3f5a",
  accent: "#00e8a2", accent2: "#ff5088", accent3: "#7060ff", accent4: "#ffa530", accent5: "#30c0ff",
  text: "#d0dae8", dim: "#556878", dimmer: "#354050",
  traces: ["#00e8a2","#ff5088","#7060ff","#ffa530","#30c0ff","#ff6666","#b0e830","#e050ff","#50ffb0","#ffcc00"],
  cellColors: ["#00e8a2","#30c0ff","#ffa530","#ff5088","#7060ff","#e050ff"],
};

const EDGE_UNITS = { K: "nM", n: "", w: "fold", alpha: "nM/min", beta: "1/min" };

/* ═══════════════════════════════════════════════
   MATH: Eigenvalue solver (QR iteration)
   ═══════════════════════════════════════════════ */
function matMul(A, B) {
  const n = A.length;
  return A.map((row, i) => B[0].map((_, j) => row.reduce((s, _, k) => s + A[i][k] * B[k][j], 0)));
}

function qrDecomp(A) {
  const n = A.length, Q = [], R = Array.from({ length: n }, () => Array(n).fill(0));
  for (let j = 0; j < n; j++) {
    let v = A.map(r => r[j]);
    for (let i = 0; i < j; i++) {
      const d = v.reduce((s, val, k) => s + val * Q[i][k], 0);
      R[i][j] = d;
      v = v.map((val, k) => val - d * Q[i][k]);
    }
    const norm = Math.sqrt(v.reduce((s, x) => s + x * x, 0));
    R[j][j] = norm;
    Q.push(norm > 1e-12 ? v.map(x => x / norm) : v);
  }
  const Qm = A.map((_, i) => Q.map(col => col[i]));
  return { Q: Qm, R };
}

function computeEigenvalues(M) {
  const n = M.length;
  if (n === 0) return [];
  if (n === 1) return [{ re: M[0][0], im: 0 }];
  let A = M.map(r => [...r]);
  for (let i = 0; i < 80; i++) {
    const { Q, R } = qrDecomp(A);
    A = matMul(R, Q);
  }
  const eigs = [];
  let i = 0;
  while (i < n) {
    if (i + 1 < n && Math.abs(A[i + 1][i]) > 1e-8) {
      const a = A[i][i], b = A[i][i + 1], c = A[i + 1][i], d = A[i + 1][i + 1];
      const tr = a + d, det = a * d - b * c, disc = tr * tr - 4 * det;
      if (disc < 0) {
        eigs.push({ re: tr / 2, im: Math.sqrt(-disc) / 2 });
        eigs.push({ re: tr / 2, im: -Math.sqrt(-disc) / 2 });
      } else {
        eigs.push({ re: (tr + Math.sqrt(disc)) / 2, im: 0 });
        eigs.push({ re: (tr - Math.sqrt(disc)) / 2, im: 0 });
      }
      i += 2;
    } else {
      eigs.push({ re: A[i][i], im: 0 });
      i++;
    }
  }
  return eigs;
}

function classifyStability(eigs) {
  if (eigs.length === 0) return { type: "unknown", label: "N/A", color: C.dim };
  const allNegReal = eigs.every(e => e.re < -1e-6);
  const anyPosReal = eigs.some(e => e.re > 1e-6);
  const hasImag = eigs.some(e => Math.abs(e.im) > 1e-6);
  if (allNegReal && hasImag) return { type: "stable_spiral", label: "Stable Spiral (Damped Oscillation)", color: C.accent };
  if (allNegReal && !hasImag) return { type: "stable_node", label: "Stable Node", color: C.accent };
  if (anyPosReal && hasImag) return { type: "unstable_spiral", label: "Unstable Spiral (Growing Oscillation)", color: C.accent2 };
  if (anyPosReal && !hasImag) return { type: "unstable_node", label: "Unstable Node", color: C.accent2 };
  if (eigs.some(e => e.re > 1e-6) && eigs.some(e => e.re < -1e-6)) return { type: "saddle", label: "Saddle Point", color: C.accent4 };
  if (eigs.every(e => Math.abs(e.re) < 1e-6) && hasImag) return { type: "center", label: "Center (Sustained Oscillation)", color: C.accent3 };
  return { type: "marginal", label: "Marginally Stable", color: C.accent4 };
}

/* ═══════════════════════════════════════════════
   SIMULATION ENGINE (Multi-cell)
   ═══════════════════════════════════════════════ */
function hillAct(x, K, n) { return Math.pow(x, n) / (Math.pow(K, n) + Math.pow(x, n)); }
function hillRep(x, K, n) { return Math.pow(K, n) / (Math.pow(K, n) + Math.pow(x, n)); }

function buildODESystem(cells, interEdges, template) {
  // Returns a function: (state, t, inputFn) => derivatives
  return (state, t, inputFn) => {
    const deriv = {};
    cells.forEach((cell, ci) => {
      template.nodes.forEach(nd => {
        const key = `${cell.id}_${nd.id}`;
        if (nd.type === "input") {
          deriv[key] = 0;
          return;
        }
        let prod = nd.alpha;
        // Internal edges
        template.edges.forEach(e => {
          if (e.to !== nd.id) return;
          const srcKey = `${cell.id}_${e.from}`;
          const srcVal = state[srcKey] || 0;
          if (e.type === "activate") prod += e.w * hillAct(srcVal, e.K, e.n);
          else prod *= hillRep(srcVal, e.K, e.n);
        });
        // Intercellular edges targeting this cell's node
        interEdges.forEach(ie => {
          if (ie.toCell !== cell.id || ie.toNode !== nd.id) return;
          const srcKey = `${ie.fromCell}_${ie.fromNode}`;
          const srcVal = state[srcKey] || 0;
          if (ie.type === "activate") prod += ie.w * hillAct(srcVal, ie.K, ie.n);
          else prod *= hillRep(srcVal, ie.K, ie.n);
        });
        deriv[key] = prod - nd.beta * (state[key] || 0);
      });
    });
    return deriv;
  };
}

function simulateSystem(cells, interEdges, template, tMax, dt, inputSignal, inputOn, inputOff) {
  const ode = buildODESystem(cells, interEdges, template);
  const steps = Math.floor(tMax / dt);
  const state = {};
  const result = { t: [] };

  // Init
  cells.forEach(cell => {
    template.nodes.forEach(nd => {
      const key = `${cell.id}_${nd.id}`;
      state[key] = nd.init || 0;
      result[key] = [];
    });
  });

  const inputFns = {
    step: (t) => t >= inputOn && t < inputOff ? 1 : 0,
    pulse: (t) => t >= inputOn && t < inputOn + (inputOff - inputOn) * 0.3 ? 1 : 0,
    ramp: (t) => t < inputOn ? 0 : Math.min(1, (t - inputOn) / ((inputOff - inputOn) * 0.5)),
    sine: (t) => 0.5 + 0.5 * Math.sin(t * 0.8),
    square: (t) => Math.sin(t * 0.5) > 0 ? 1 : 0,
  };
  const inputFn = inputFns[inputSignal] || inputFns.step;

  for (let i = 0; i <= steps; i++) {
    const t = i * dt;
    result.t.push(t);
    // Set input nodes
    cells.forEach(cell => {
      template.nodes.filter(n => n.type === "input").forEach(nd => {
        state[`${cell.id}_${nd.id}`] = inputFn(t);
      });
    });
    // Record
    Object.keys(state).forEach(k => result[k].push(state[k]));
    // Euler step
    const d = ode(state, t, inputFn);
    Object.keys(d).forEach(k => { state[k] = Math.max(0, (state[k] || 0) + d[k] * dt); });
  }
  return result;
}

function computeJacobian(cells, interEdges, template, steadyState) {
  const keys = [];
  cells.forEach(cell => {
    template.nodes.filter(n => n.type !== "input").forEach(nd => {
      keys.push(`${cell.id}_${nd.id}`);
    });
  });
  const n = keys.length;
  const J = Array.from({ length: n }, () => Array(n).fill(0));
  const eps = 1e-6;
  const ode = buildODESystem(cells, interEdges, template);
  const f0 = ode(steadyState, 0, () => 1);

  for (let j = 0; j < n; j++) {
    const perturbed = { ...steadyState };
    perturbed[keys[j]] = (perturbed[keys[j]] || 0) + eps;
    const f1 = ode(perturbed, 0, () => 1);
    for (let i = 0; i < n; i++) {
      J[i][j] = ((f1[keys[i]] || 0) - (f0[keys[i]] || 0)) / eps;
    }
  }
  return { J, keys };
}

/* ═══════════════════════════════════════════════
   INPUT SIGNALS
   ═══════════════════════════════════════════════ */
const INPUTS = {
  step: "Step", pulse: "Pulse", ramp: "Ramp", sine: "Sine", square: "Square",
};

/* ═══════════════════════════════════════════════
   PRESETS
   ═══════════════════════════════════════════════ */
const PRESETS = {
  negFB: {
    name: "Negative Feedback",
    nodes: [
      { id: "s", name: "Signal", type: "input", alpha: 0, beta: 0, init: 0 },
      { id: "x", name: "GeneX", type: "gene", alpha: 5, beta: 1, init: 0 },
    ],
    edges: [
      { from: "s", to: "x", type: "activate", K: 2, n: 2, w: 1 },
      { from: "x", to: "x", type: "repress", K: 2, n: 2, w: 1 },
    ],
  },
  toggle: {
    name: "Toggle Switch",
    nodes: [
      { id: "a", name: "GeneA", type: "gene", alpha: 6, beta: 1, init: 5 },
      { id: "b", name: "GeneB", type: "gene", alpha: 6, beta: 1, init: 0.1 },
    ],
    edges: [
      { from: "a", to: "b", type: "repress", K: 2, n: 3, w: 1 },
      { from: "b", to: "a", type: "repress", K: 2, n: 3, w: 1 },
    ],
  },
  repressilator: {
    name: "Repressilator",
    nodes: [
      { id: "a", name: "LacI", type: "gene", alpha: 10, beta: 1, init: 5 },
      { id: "b", name: "TetR", type: "gene", alpha: 10, beta: 1, init: 1 },
      { id: "c", name: "cI", type: "gene", alpha: 10, beta: 1, init: 0.5 },
    ],
    edges: [
      { from: "a", to: "b", type: "repress", K: 2, n: 3, w: 1 },
      { from: "b", to: "c", type: "repress", K: 2, n: 3, w: 1 },
      { from: "c", to: "a", type: "repress", K: 2, n: 3, w: 1 },
    ],
  },
  cascade: {
    name: "Signaling Cascade",
    nodes: [
      { id: "s", name: "Ligand", type: "input", alpha: 0, beta: 0, init: 0 },
      { id: "r", name: "Receptor", type: "protein", alpha: 3, beta: 0.5, init: 0 },
      { id: "k", name: "Kinase", type: "protein", alpha: 4, beta: 0.8, init: 0 },
      { id: "tf", name: "TF", type: "gene", alpha: 5, beta: 1, init: 0 },
    ],
    edges: [
      { from: "s", to: "r", type: "activate", K: 1, n: 2, w: 1 },
      { from: "r", to: "k", type: "activate", K: 1.5, n: 2, w: 1 },
      { from: "k", to: "tf", type: "activate", K: 2, n: 2, w: 1 },
    ],
  },
};

const NODE_META = {
  gene: { label: "Gene", icon: "🧬", color: C.accent },
  protein: { label: "Protein", icon: "⬡", color: C.accent5 },
  mrna: { label: "mRNA", icon: "〰", color: C.accent4 },
  input: { label: "Input", icon: "⚡", color: C.accent4 },
};

let _id = 0;
const uid = () => `n${++_id}`;
const cellUid = () => `c${++_id}`;

/* ═══════════════════════════════════════════════
   UI PRIMITIVES
   ═══════════════════════════════════════════════ */
const Btn = ({ children, onClick, active, color = C.accent, small, disabled, style: s }) => (
  <button onClick={onClick} disabled={disabled} style={{
    background: active ? `${color}20` : "transparent", border: `1px solid ${active ? color : C.border}`,
    color: active ? color : disabled ? C.dimmer : C.dim, borderRadius: 6,
    padding: small ? "3px 8px" : "5px 11px", fontSize: small ? 9 : 10,
    cursor: disabled ? "not-allowed" : "pointer", fontFamily: "inherit", transition: "all .15s", ...s,
  }}>{children}</button>
);
const Inp = ({ value, onChange, placeholder, type = "text", style: s, ...rest }) => (
  <input type={type} value={value} onChange={e => onChange(type === "number" ? parseFloat(e.target.value) || 0 : e.target.value)}
    placeholder={placeholder} {...rest}
    style={{ background: C.bg, border: `1px solid ${C.border}`, borderRadius: 5, color: C.text, padding: "4px 7px", fontSize: 11, fontFamily: "inherit", outline: "none", width: "100%", ...s }} />
);
const Lbl = ({ children, color = C.dim }) => (
  <div style={{ fontSize: 9, color, marginBottom: 2, textTransform: "uppercase", letterSpacing: 1, fontFamily: "inherit" }}>{children}</div>
);
const Sec = ({ title, color = C.accent, children, icon, extra }) => (
  <div style={{ background: C.panel, border: `1px solid ${C.border}`, borderRadius: 10, padding: 12, marginBottom: 10 }}>
    <div style={{ fontSize: 10, fontWeight: 700, color, letterSpacing: 2, textTransform: "uppercase", marginBottom: 10, display: "flex", alignItems: "center", gap: 6 }}>
      {icon && <span style={{ fontSize: 12 }}>{icon}</span>}{title}
      {extra && <span style={{ marginLeft: "auto" }}>{extra}</span>}
    </div>
    {children}
  </div>
);

/* ═══════════════════════════════════════════════
   EQUATION MODAL
   ═══════════════════════════════════════════════ */
function EquationModal({ data, onClose }) {
  if (!data) return null;
  return (
    <div onClick={onClose} style={{ position: "fixed", inset: 0, background: "rgba(0,0,0,.7)", zIndex: 999, display: "flex", alignItems: "center", justifyContent: "center", padding: 20 }}>
      <div onClick={e => e.stopPropagation()} style={{ background: C.panel, border: `1px solid ${C.borderHi}`, borderRadius: 12, padding: 20, maxWidth: 420, width: "100%" }}>
        <div style={{ fontSize: 13, fontWeight: 700, color: data.color || C.accent, marginBottom: 10 }}>{data.title}</div>
        <div style={{ background: C.bg, borderRadius: 8, padding: 14, fontFamily: "'JetBrains Mono', monospace", fontSize: 12, color: C.text, lineHeight: 2, marginBottom: 10, overflowX: "auto" }}>
          {data.equations.map((eq, i) => <div key={i}>{eq}</div>)}
        </div>
        {data.description && <div style={{ fontSize: 11, color: C.dim, lineHeight: 1.7, marginBottom: 10 }}>{data.description}</div>}
        {data.korean && <div style={{ fontSize: 10, color: C.dimmer, lineHeight: 1.6, padding: "8px 10px", background: `${C.accent}08`, borderRadius: 6 }}>🇰🇷 {data.korean}</div>}
        <button onClick={onClose} style={{ marginTop: 10, background: "none", border: `1px solid ${C.border}`, color: C.dim, borderRadius: 6, padding: "5px 16px", fontSize: 10, cursor: "pointer", fontFamily: "inherit" }}>Close</button>
      </div>
    </div>
  );
}

/* ═══════════════════════════════════════════════
   NETWORK DIAGRAM (with constants on arrows)
   ═══════════════════════════════════════════════ */
function CellNetworkDiagram({ template, cells, interEdges, onClickNode, onClickEdge, onClickInterEdge, selectedNode }) {
  const singleCell = cells.length <= 1;
  const cellW = singleCell ? 380 : 180, cellH = singleCell ? 200 : 150;
  const svgW = singleCell ? 400 : Math.min(800, cells.length * 200 + 50);
  const svgH = singleCell ? 220 : cells.length > 2 ? 360 : 200;

  const cellPositions = cells.map((c, i) => {
    if (singleCell) return { x: 10, y: 10 };
    const cols = Math.min(cells.length, 3);
    const row = Math.floor(i / cols), col = i % cols;
    return { x: 20 + col * (cellW + 20), y: 10 + row * (cellH + 30) };
  });

  const nodePositions = {};
  cells.forEach((cell, ci) => {
    const cp = cellPositions[ci];
    const ns = template.nodes;
    const n = ns.length;
    const cx = cp.x + cellW / 2, cy = cp.y + cellH / 2;
    const r = singleCell ? Math.min(70, 30 + n * 10) : Math.min(45, 20 + n * 8);
    ns.forEach((nd, ni) => {
      const a = -Math.PI / 2 + (2 * Math.PI * ni) / n;
      nodePositions[`${cell.id}_${nd.id}`] = { x: cx + r * Math.cos(a), y: cy + r * Math.sin(a) };
    });
  });

  const fontSize = singleCell ? 9 : 7;
  const nodeR = singleCell ? 17 : 12;

  return (
    <svg viewBox={`0 0 ${svgW} ${svgH}`} style={{ width: "100%", height: singleCell ? 220 : Math.min(360, svgH), display: "block" }}>
      <defs>
        <marker id="ma" markerWidth="7" markerHeight="5" refX="7" refY="2.5" orient="auto"><polygon points="0 0,7 2.5,0 5" fill={C.accent} /></marker>
        <marker id="mr" markerWidth="7" markerHeight="7" refX="4" refY="3.5" orient="auto"><rect x="3" y="0" width="2.5" height="7" fill={C.accent2} /></marker>
        <marker id="mia" markerWidth="7" markerHeight="5" refX="7" refY="2.5" orient="auto"><polygon points="0 0,7 2.5,0 5" fill={C.accent3} /></marker>
        <marker id="mir" markerWidth="7" markerHeight="7" refX="4" refY="3.5" orient="auto"><rect x="3" y="0" width="2.5" height="7" fill={C.accent4} /></marker>
      </defs>

      {/* Cell boundaries */}
      {cells.map((cell, ci) => {
        const cp = cellPositions[ci];
        return (
          <g key={cell.id}>
            <rect x={cp.x} y={cp.y} width={cellW} height={cellH} rx={12} fill={`${C.cellColors[ci % C.cellColors.length]}06`}
              stroke={C.cellColors[ci % C.cellColors.length]} strokeWidth={1.2} strokeDasharray="6,3" />
            <text x={cp.x + 8} y={cp.y + 14} fontSize={singleCell ? 10 : 8} fill={C.cellColors[ci % C.cellColors.length]} fontFamily="monospace" fontWeight="700">
              {cell.name}
            </text>
          </g>
        );
      })}

      {/* Internal edges with constants */}
      {cells.map((cell) =>
        template.edges.map((e, ei) => {
          const fk = `${cell.id}_${e.from}`, tk = `${cell.id}_${e.to}`;
          const fp = nodePositions[fk], tp = nodePositions[tk];
          if (!fp || !tp) return null;
          const isSelf = e.from === e.to;
          const isAct = e.type === "activate";
          const color = isAct ? C.accent : C.accent2;
          const marker = isAct ? "url(#ma)" : "url(#mr)";
          const labelText = `K=${e.K}${EDGE_UNITS.K} n=${e.n}`;

          if (isSelf) {
            const r2 = singleCell ? 30 : 22;
            return (
              <g key={`${cell.id}_ie_${ei}`} onClick={() => onClickEdge(e)} style={{ cursor: "pointer" }}>
                <path d={`M ${fp.x - 8} ${fp.y - nodeR} C ${fp.x - r2} ${fp.y - r2 * 1.6}, ${fp.x + r2} ${fp.y - r2 * 1.6}, ${fp.x + 8} ${fp.y - nodeR}`}
                  fill="none" stroke={color} strokeWidth={1.2} strokeDasharray={isAct ? "none" : "4,2"} markerEnd={marker} opacity={0.7} />
                <text x={fp.x} y={fp.y - r2 * 1.4 - 2} textAnchor="middle" fontSize={singleCell ? 7 : 5.5} fill={color} fontFamily="monospace" opacity={0.9}>{labelText}</text>
              </g>
            );
          }
          const dx = tp.x - fp.x, dy = tp.y - fp.y, len = Math.sqrt(dx * dx + dy * dy) || 1;
          const mx = (fp.x + tp.x) / 2, my = (fp.y + tp.y) / 2;
          const off = singleCell ? 14 : 10;
          const nx = -dy / len * off, ny = dx / len * off;
          return (
            <g key={`${cell.id}_ie_${ei}`} onClick={() => onClickEdge(e)} style={{ cursor: "pointer" }}>
              <path d={`M ${fp.x} ${fp.y} Q ${mx + nx} ${my + ny} ${tp.x} ${tp.y}`}
                fill="none" stroke={color} strokeWidth={1.2} strokeDasharray={isAct ? "none" : "4,2"} markerEnd={marker} opacity={0.7} />
              <rect x={mx + nx * 0.8 - 28} y={my + ny * 0.8 - 7} width={56} height={12} rx={3} fill={C.bg} fillOpacity={0.85} />
              <text x={mx + nx * 0.8} y={my + ny * 0.8 + 3} textAnchor="middle" fontSize={singleCell ? 7 : 5.5} fill={color} fontFamily="monospace">{labelText}</text>
            </g>
          );
        })
      )}

      {/* Intercellular edges */}
      {interEdges.map((ie, idx) => {
        const fk = `${ie.fromCell}_${ie.fromNode}`, tk = `${ie.toCell}_${ie.toNode}`;
        const fp = nodePositions[fk], tp = nodePositions[tk];
        if (!fp || !tp) return null;
        const isAct = ie.type === "activate";
        const color = isAct ? C.accent3 : C.accent4;
        const marker = isAct ? "url(#mia)" : "url(#mir)";
        const mx = (fp.x + tp.x) / 2, my = (fp.y + tp.y) / 2;
        return (
          <g key={`inter_${idx}`} onClick={() => onClickInterEdge(ie, idx)} style={{ cursor: "pointer" }}>
            <line x1={fp.x} y1={fp.y} x2={tp.x} y2={tp.y} stroke={color} strokeWidth={1.8} strokeDasharray={isAct ? "8,4" : "4,4"} markerEnd={marker} opacity={0.8} />
            <rect x={mx - 30} y={my - 8} width={60} height={14} rx={4} fill={C.bg} fillOpacity={0.9} stroke={color} strokeWidth={0.5} />
            <text x={mx} y={my + 3} textAnchor="middle" fontSize={7} fill={color} fontFamily="monospace" fontWeight="600">
              K={ie.K}{EDGE_UNITS.K} w={ie.w}
            </text>
          </g>
        );
      })}

      {/* Nodes */}
      {cells.map((cell) =>
        template.nodes.map(nd => {
          const key = `${cell.id}_${nd.id}`;
          const p = nodePositions[key];
          if (!p) return null;
          const meta = NODE_META[nd.type];
          const isSel = selectedNode === key;
          return (
            <g key={key} onClick={() => onClickNode(nd, cell)} style={{ cursor: "pointer" }}>
              <circle cx={p.x} cy={p.y} r={isSel ? nodeR + 3 : nodeR} fill={isSel ? `${meta.color}25` : `${meta.color}12`}
                stroke={isSel ? meta.color : `${meta.color}50`} strokeWidth={isSel ? 2 : 1} />
              <text x={p.x} y={p.y + 1} textAnchor="middle" fontSize={singleCell ? 13 : 10} fill="white" dominantBaseline="central">{meta.icon}</text>
              <text x={p.x} y={p.y + nodeR + (singleCell ? 10 : 8)} textAnchor="middle" fontSize={fontSize} fill={C.text} fontFamily="monospace">{nd.name}</text>
              <text x={p.x} y={p.y + nodeR + (singleCell ? 19 : 15)} textAnchor="middle" fontSize={fontSize - 1.5} fill={C.dimmer} fontFamily="monospace">
                α={nd.alpha} β={nd.beta}
              </text>
            </g>
          );
        })
      )}
    </svg>
  );
}

/* ═══════════════════════════════════════════════
   ELECTRONIC CIRCUIT ANALOG
   ═══════════════════════════════════════════════ */
function ElectronicAnalog({ template }) {
  const { nodes, edges } = template;
  if (nodes.length === 0) return <div style={{ color: C.dimmer, fontSize: 10, textAlign: "center", padding: 20 }}>Add nodes to see the electronic analog</div>;

  const w = 400, h = 180;
  const gap = w / (nodes.length + 1);

  const analogMap = {
    gene: { symbol: "Op-Amp", desc: "Nonlinear Amplifier" },
    protein: { symbol: "Buffer", desc: "Voltage Follower" },
    mrna: { symbol: "Amp", desc: "Linear Amplifier" },
    input: { symbol: "Vin", desc: "Signal Source" },
  };

  return (
    <svg viewBox={`0 0 ${w} ${h}`} style={{ width: "100%", height: 160, display: "block" }}>
      <text x={w / 2} y={12} textAnchor="middle" fontSize={8} fill={C.dim} fontFamily="monospace">≡ Equivalent Electronic Circuit</text>
      {/* Ground rail */}
      <line x1={20} y1={h - 20} x2={w - 20} y2={h - 20} stroke={C.dimmer} strokeWidth={0.8} />
      <text x={w - 15} y={h - 14} fontSize={7} fill={C.dimmer} fontFamily="monospace">GND</text>
      {/* Vdd rail */}
      <line x1={20} y1={25} x2={w - 20} y2={25} stroke={C.accent} strokeWidth={0.5} opacity={0.3} />
      <text x={w - 15} y={22} fontSize={7} fill={C.accent} fontFamily="monospace" opacity={0.5}>Vdd</text>

      {nodes.map((nd, i) => {
        const x = gap * (i + 1);
        const y = 80;
        const meta = NODE_META[nd.type];
        const am = analogMap[nd.type];

        return (
          <g key={nd.id}>
            {/* RC to ground (degradation) */}
            {nd.type !== "input" && (
              <>
                <line x1={x} y1={y + 22} x2={x} y2={h - 30} stroke={C.dim} strokeWidth={0.8} />
                <rect x={x - 6} y={h - 38} width={12} height={8} fill="none" stroke={C.accent4} strokeWidth={0.8} />
                <text x={x + 10} y={h - 32} fontSize={5.5} fill={C.accent4} fontFamily="monospace">τ=1/β={nd.beta > 0 ? (1/nd.beta).toFixed(1) : "∞"}min</text>
                <line x1={x} y1={h - 30} x2={x} y2={h - 20} stroke={C.dim} strokeWidth={0.8} />
              </>
            )}

            {/* Component body */}
            {nd.type === "input" ? (
              <>
                <circle cx={x} cy={y} r={14} fill="none" stroke={C.accent4} strokeWidth={1.2} />
                <text x={x} y={y - 2} textAnchor="middle" fontSize={8} fill={C.accent4} fontFamily="monospace">~</text>
                <text x={x} y={y + 5} textAnchor="middle" fontSize={5} fill={C.accent4} fontFamily="monospace">Vin</text>
              </>
            ) : (
              <>
                <polygon points={`${x - 16},${y - 16} ${x - 16},${y + 16} ${x + 18},${y}`}
                  fill={`${meta.color}10`} stroke={meta.color} strokeWidth={1.2} />
                <text x={x - 4} y={y - 3} fontSize={5} fill={meta.color} fontFamily="monospace">{nd.type === "gene" ? "+" : "1"}</text>
                <text x={x - 4} y={y + 6} fontSize={5} fill={meta.color} fontFamily="monospace">{nd.type === "gene" ? "−" : ""}</text>
              </>
            )}
            <text x={x} y={y + 30} textAnchor="middle" fontSize={6.5} fill={C.text} fontFamily="monospace">{nd.name}</text>
            <text x={x} y={y + 38} textAnchor="middle" fontSize={5} fill={C.dimmer} fontFamily="monospace">{am.desc}</text>
          </g>
        );
      })}

      {/* Connections */}
      {edges.map((e, ei) => {
        const fi = nodes.findIndex(n => n.id === e.from);
        const ti = nodes.findIndex(n => n.id === e.to);
        if (fi < 0 || ti < 0) return null;
        const fx = gap * (fi + 1) + 18, tx = gap * (ti + 1) - 16;
        const isSelf = fi === ti;
        const isAct = e.type === "activate";
        const color = isAct ? C.accent : C.accent2;

        if (isSelf) {
          const x = gap * (fi + 1);
          return (
            <g key={`ae${ei}`}>
              <path d={`M ${x + 18} ${80 - 5} L ${x + 32} ${80 - 5} L ${x + 32} ${48} L ${x - 20} ${48} L ${x - 20} ${80 - 10}`}
                fill="none" stroke={color} strokeWidth={0.8} strokeDasharray={isAct ? "none" : "3,2"} />
              <text x={x + 6} y={44} textAnchor="middle" fontSize={5} fill={color} fontFamily="monospace">
                {isAct ? "+" : "−"} feedback (R={e.K}kΩ)
              </text>
            </g>
          );
        }
        const y = fi < ti ? 80 - 3 : 80 + 3;
        return (
          <g key={`ae${ei}`}>
            <line x1={fx} y1={y} x2={tx} y2={y} stroke={color} strokeWidth={0.8} strokeDasharray={isAct ? "none" : "3,2"} />
            <text x={(fx + tx) / 2} y={y - 4} textAnchor="middle" fontSize={5} fill={color} fontFamily="monospace">
              {isAct ? "→" : "⊣"} gain={e.w}
            </text>
          </g>
        );
      })}
    </svg>
  );
}

/* ═══════════════════════════════════════════════
   TIME SERIES GRAPH
   ═══════════════════════════════════════════════ */
function Graph({ data, traceKeys, traceLabels, traceColors, zoomRange, setZoomRange }) {
  const canvasRef = useRef(null);
  const dragRef = useRef(null);
  const [hover, setHover] = useState(null);

  useEffect(() => {
    const cv = canvasRef.current;
    if (!cv || !data) return;
    const ctx = cv.getContext("2d");
    const dpr = window.devicePixelRatio || 1;
    const w = cv.clientWidth, h = cv.clientHeight;
    cv.width = w * dpr; cv.height = h * dpr;
    ctx.scale(dpr, dpr);

    const pad = { top: 16, right: 14, bottom: 30, left: 44 };
    const gw = w - pad.left - pad.right, gh = h - pad.top - pad.bottom;
    ctx.fillStyle = C.bg; ctx.fillRect(0, 0, w, h);

    const [tS, tE] = zoomRange || [data.t[0], data.t[data.t.length - 1]];
    const si = data.t.findIndex(t => t >= tS);
    const ei2 = data.t.findIndex(t => t > tE);
    const se = ei2 === -1 ? data.t.length : ei2;
    const ts = data.t.slice(si, se);

    let yMax = 0.1;
    traceKeys.forEach(k => {
      const sl = data[k]?.slice(si, se) || [];
      sl.forEach(v => { if (v > yMax) yMax = v; });
    });
    yMax *= 1.15;

    const toX = t => pad.left + ((t - tS) / (tE - tS)) * gw;
    const toY = v => pad.top + gh - (v / yMax) * gh;

    // grid
    ctx.strokeStyle = C.dimmer; ctx.lineWidth = 0.3;
    for (let i = 0; i <= 5; i++) { const y = pad.top + gh * i / 5; ctx.beginPath(); ctx.moveTo(pad.left, y); ctx.lineTo(w - pad.right, y); ctx.stroke(); }
    for (let i = 0; i <= 8; i++) { const x = pad.left + gw * i / 8; ctx.beginPath(); ctx.moveTo(x, pad.top); ctx.lineTo(x, h - pad.bottom); ctx.stroke(); }

    const step = Math.max(1, Math.floor(ts.length / (gw * 2)));
    traceKeys.forEach((k, idx) => {
      const yd = data[k]?.slice(si, se) || [];
      const col = traceColors[idx];
      ctx.strokeStyle = col; ctx.lineWidth = 1.8; ctx.lineJoin = "round"; ctx.beginPath();
      for (let i = 0; i < ts.length; i += step) { const px = toX(ts[i]), py = toY(yd[i]); i === 0 ? ctx.moveTo(px, py) : ctx.lineTo(px, py); }
      ctx.stroke();
      ctx.save(); ctx.shadowColor = col; ctx.shadowBlur = 4; ctx.strokeStyle = col; ctx.lineWidth = 0.6; ctx.beginPath();
      for (let i = 0; i < ts.length; i += step) { const px = toX(ts[i]), py = toY(yd[i]); i === 0 ? ctx.moveTo(px, py) : ctx.lineTo(px, py); }
      ctx.stroke(); ctx.restore();
    });

    // drag overlay
    if (dragRef.current?.dragging) {
      const { x1, x2 } = dragRef.current;
      ctx.fillStyle = `${C.accent}0a`; ctx.strokeStyle = C.accent; ctx.lineWidth = 1;
      const l = Math.min(x1, x2), r2 = Math.max(x1, x2);
      ctx.fillRect(l, pad.top, r2 - l, gh); ctx.strokeRect(l, pad.top, r2 - l, gh);
    }

    // axes
    ctx.strokeStyle = C.border; ctx.lineWidth = 1; ctx.beginPath();
    ctx.moveTo(pad.left, pad.top); ctx.lineTo(pad.left, h - pad.bottom); ctx.lineTo(w - pad.right, h - pad.bottom); ctx.stroke();
    ctx.fillStyle = C.dim; ctx.font = "9px 'JetBrains Mono',monospace"; ctx.textAlign = "center";
    for (let i = 0; i <= 4; i++) { const tv = tS + (tE - tS) * i / 4; ctx.fillText(tv.toFixed(1), toX(tv), h - pad.bottom + 13); }
    ctx.textAlign = "right";
    for (let i = 0; i <= 4; i++) { const yv = yMax * i / 4; ctx.fillText(yv.toFixed(1), pad.left - 5, toY(yv) + 3); }
    ctx.textAlign = "center"; ctx.fillText("time (min)", w / 2, h - 3);

    // hover
    if (hover) {
      const hx = hover.x;
      if (hx >= pad.left && hx <= w - pad.right) {
        ctx.strokeStyle = "rgba(255,255,255,.12)"; ctx.lineWidth = 1; ctx.setLineDash([3, 3]);
        ctx.beginPath(); ctx.moveTo(hx, pad.top); ctx.lineTo(hx, h - pad.bottom); ctx.stroke(); ctx.setLineDash([]);
        const tAtX = tS + ((hx - pad.left) / gw) * (tE - tS);
        const tIdx = data.t.findIndex(t => t >= tAtX);
        if (tIdx >= 0) {
          let ty = 28;
          ctx.fillStyle = "rgba(10,14,20,.92)"; ctx.fillRect(hx + 6, ty - 12, 110, traceKeys.length * 14 + 18); ctx.font = "8px monospace";
          ctx.textAlign = "left"; ctx.fillStyle = C.dim; ctx.fillText(`t = ${tAtX.toFixed(2)} min`, hx + 12, ty); ty += 12;
          traceKeys.forEach((k, i) => { ctx.fillStyle = traceColors[i]; ctx.fillText(`${traceLabels[i]}: ${(data[k]?.[tIdx] ?? 0).toFixed(3)} nM`, hx + 12, ty); ty += 13; });
        }
      }
    }
  }, [data, traceKeys, traceLabels, traceColors, zoomRange, hover]);

  const getX = e => e.clientX - canvasRef.current.getBoundingClientRect().left;
  return (
    <div>
      <canvas ref={canvasRef}
        onMouseDown={e => { dragRef.current = { dragging: true, x1: getX(e), x2: getX(e) }; }}
        onMouseMove={e => { setHover({ x: getX(e) }); if (dragRef.current?.dragging) dragRef.current.x2 = getX(e); }}
        onMouseUp={() => {
          if (!dragRef.current?.dragging || !data) { dragRef.current = null; return; }
          const { x1, x2 } = dragRef.current; dragRef.current = null;
          const cv = canvasRef.current, w = cv.clientWidth, pad = { left: 44, right: 14 }, gw = w - pad.left - pad.right;
          const [tS, tE] = zoomRange || [data.t[0], data.t[data.t.length - 1]];
          const pToT = px => tS + ((px - pad.left) / gw) * (tE - tS);
          if (Math.abs(x2 - x1) > 8) setZoomRange([Math.max(0, pToT(Math.min(x1, x2))), pToT(Math.max(x1, x2))]);
        }}
        onMouseLeave={() => { setHover(null); if (dragRef.current) dragRef.current = null; }}
        style={{ width: "100%", height: 260, display: "block", borderRadius: 8, cursor: "crosshair" }} />
      <div style={{ display: "flex", flexWrap: "wrap", gap: 8, padding: "6px 0", justifyContent: "center" }}>
        {traceLabels.map((l, i) => (
          <span key={i} style={{ fontSize: 9, color: traceColors[i], display: "flex", alignItems: "center", gap: 3 }}>
            <span style={{ width: 10, height: 2.5, background: traceColors[i], display: "inline-block", borderRadius: 1 }} />{l}
          </span>
        ))}
      </div>
    </div>
  );
}

/* ═══════════════════════════════════════════════
   MAIN APPLICATION
   ═══════════════════════════════════════════════ */
export default function BioCircuitWorkbench() {
  const [template, setTemplate] = useState({ nodes: [], edges: [] });
  const [cells, setCells] = useState([{ id: "c1", name: "Cell 1" }]);
  const [interEdges, setInterEdges] = useState([]);
  const [tab, setTab] = useState("design");
  const [selectedNode, setSelectedNode] = useState(null);
  const [eqModal, setEqModal] = useState(null);
  const [simData, setSimData] = useState(null);
  const [zoomRange, setZoomRange] = useState(null);
  const [stability, setStability] = useState(null);

  // sim params
  const [tMax, setTMax] = useState(50);
  const [dt] = useState(0.02);
  const [inputSig, setInputSig] = useState("step");
  const [inputOn, setInputOn] = useState(5);
  const [inputOff, setInputOff] = useState(45);

  // forms
  const [newName, setNewName] = useState("");
  const [newType, setNewType] = useState("gene");
  const [efrom, setEfrom] = useState("");
  const [eto, setEto] = useState("");
  const [etype, setEtype] = useState("activate");
  const [iefromCell, setIEfromCell] = useState("");
  const [iefromNode, setIEfromNode] = useState("");
  const [ietoCell, setIEtoCell] = useState("");
  const [ietoNode, setIEtoNode] = useState("");
  const [ieType, setIEtype] = useState("activate");

  const addNode = () => {
    if (!newName.trim()) return;
    const nd = { id: uid(), name: newName.trim(), type: newType, alpha: newType === "input" ? 0 : 5, beta: newType === "input" ? 0 : 1, init: 0 };
    setTemplate(t => ({ ...t, nodes: [...t.nodes, nd] }));
    setNewName("");
  };
  const removeNode = id => {
    setTemplate(t => ({ nodes: t.nodes.filter(n => n.id !== id), edges: t.edges.filter(e => e.from !== id && e.to !== id) }));
    setInterEdges(ie => ie.filter(e => e.fromNode !== id && e.toNode !== id));
  };
  const updateNode = (id, u) => setTemplate(t => ({ ...t, nodes: t.nodes.map(n => n.id === id ? { ...n, ...u } : n) }));
  const addEdge = () => {
    if (!efrom || !eto) return;
    setTemplate(t => ({ ...t, edges: [...t.edges, { from: efrom, to: eto, type: etype, K: 2, n: 2, w: 1 }] }));
  };
  const removeEdge = i => setTemplate(t => ({ ...t, edges: t.edges.filter((_, j) => j !== i) }));
  const updateEdge = (i, u) => setTemplate(t => ({ ...t, edges: t.edges.map((e, j) => j === i ? { ...e, ...u } : e) }));

  const addCell = () => {
    const c = { id: cellUid(), name: `Cell ${cells.length + 1}` };
    setCells(prev => [...prev, c]);
  };
  const removeCell = id => {
    setCells(prev => prev.filter(c => c.id !== id));
    setInterEdges(prev => prev.filter(e => e.fromCell !== id && e.toCell !== id));
  };
  const renameCell = (id, name) => setCells(prev => prev.map(c => c.id === id ? { ...c, name } : c));

  const addInterEdge = () => {
    if (!iefromCell || !iefromNode || !ietoCell || !ietoNode) return;
    setInterEdges(prev => [...prev, { fromCell: iefromCell, fromNode: iefromNode, toCell: ietoCell, toNode: ietoNode, type: ieType, K: 2, n: 2, w: 1 }]);
  };
  const removeInterEdge = i => setInterEdges(prev => prev.filter((_, j) => j !== i));
  const updateInterEdge = (i, u) => setInterEdges(prev => prev.map((e, j) => j === i ? { ...e, ...u } : e));

  const loadPreset = k => {
    _id = 0;
    setTemplate({ nodes: [...PRESETS[k].nodes], edges: [...PRESETS[k].edges] });
    setCells([{ id: "c1", name: "Cell 1" }]);
    setInterEdges([]);
    setSimData(null); setZoomRange(null); setStability(null);
  };

  const runSim = useCallback(() => {
    if (template.nodes.length === 0) return;
    const data = simulateSystem(cells, interEdges, template, tMax, dt, inputSig, inputOn, inputOff);
    setSimData(data);
    setZoomRange(null);
  }, [cells, interEdges, template, tMax, dt, inputSig, inputOn, inputOff]);

  const runStability = useCallback(() => {
    if (template.nodes.filter(n => n.type !== "input").length === 0) return;
    // Find steady state by simulating long
    const ssData = simulateSystem(cells, interEdges, template, 200, 0.05, inputSig, 0, 200);
    const lastIdx = ssData.t.length - 1;
    const ss = {};
    Object.keys(ssData).forEach(k => { if (k !== "t") ss[k] = ssData[k][lastIdx]; });
    const { J, keys } = computeJacobian(cells, interEdges, template, ss);
    const eigs = computeEigenvalues(J);
    const cls = classifyStability(eigs);
    setStability({ ss, J, keys, eigs, classification: cls });
  }, [cells, interEdges, template, inputSig]);

  const showEdgeEq = (e) => {
    const fromN = template.nodes.find(n => n.id === e.from);
    const toN = template.nodes.find(n => n.id === e.to);
    if (e.type === "activate") {
      setEqModal({
        title: `${fromN?.name} → ${toN?.name} (Activation)`,
        color: C.accent,
        equations: [
          `f(x) = w · xⁿ / (Kⁿ + xⁿ)`,
          `     = ${e.w} · [${fromN?.name}]^${e.n} / (${e.K}^${e.n} + [${fromN?.name}]^${e.n})`,
          ``,
          `Parameters:`,
          `  K  = ${e.K} nM  (half-max constant)`,
          `  n  = ${e.n}     (Hill coefficient)`,
          `  w  = ${e.w}     (fold activation)`,
        ],
        description: "Hill activation function: output increases sigmoidally as input rises. K is the input concentration giving half-maximal response, n controls the steepness (cooperativity).",
        korean: "힐 활성화 함수: 입력 농도가 증가하면 S자 형태로 출력이 증가합니다. K는 반최대 반응 농도, n은 협동성(시그모이드 기울기)을 결정합니다.",
      });
    } else {
      setEqModal({
        title: `${fromN?.name} ⊣ ${toN?.name} (Repression)`,
        color: C.accent2,
        equations: [
          `f(x) = Kⁿ / (Kⁿ + xⁿ)`,
          `     = ${e.K}^${e.n} / (${e.K}^${e.n} + [${fromN?.name}]^${e.n})`,
          ``,
          `Parameters:`,
          `  K  = ${e.K} nM  (repression threshold)`,
          `  n  = ${e.n}     (Hill coefficient)`,
          `  w  = ${e.w}     (repression weight)`,
        ],
        description: "Hill repression function: output decreases sigmoidally as repressor concentration rises. This acts as an inverting element — the biological analog of a NOT gate.",
        korean: "힐 억제 함수: 억제자 농도가 증가하면 출력이 감소합니다. 전자회로의 인버터(NOT gate)와 동일한 역할을 합니다.",
      });
    }
  };

  const showNodeEq = (nd) => {
    const incoming = template.edges.filter(e => e.to === nd.id);
    const acts = incoming.filter(e => e.type === "activate").map(e => {
      const src = template.nodes.find(n => n.id === e.from);
      return `+ ${e.w}·H⁺(${src?.name}, K=${e.K}, n=${e.n})`;
    });
    const reps = incoming.filter(e => e.type === "repress").map(e => {
      const src = template.nodes.find(n => n.id === e.from);
      return `× H⁻(${src?.name}, K=${e.K}, n=${e.n})`;
    });
    setEqModal({
      title: `${NODE_META[nd.type].icon} ${nd.name} — ODE`,
      color: NODE_META[nd.type].color,
      equations: [
        `d[${nd.name}]/dt = Production − Degradation`,
        ``,
        `Production = α ${acts.join(" ")}`,
        `           ${reps.join(" ")}`,
        `         = ${nd.alpha} nM/min ${acts.join(" ")} ${reps.join(" ")}`,
        ``,
        `Degradation = β · [${nd.name}]`,
        `            = ${nd.beta}/min · [${nd.name}]`,
        ``,
        `Half-life = ln(2)/β = ${nd.beta > 0 ? (Math.log(2) / nd.beta).toFixed(2) : "∞"} min`,
        `Time constant τ = 1/β = ${nd.beta > 0 ? (1 / nd.beta).toFixed(2) : "∞"} min`,
      ],
      description: `This node has ${incoming.length} incoming connection(s). H⁺ = Hill activation, H⁻ = Hill repression. The degradation follows first-order kinetics (equivalent to RC discharge in electronics).`,
      korean: `이 노드에는 ${incoming.length}개의 입력이 있습니다. 1차 분해 반응은 전자회로에서 RC 방전과 동치입니다. τ = 1/β가 시정수(time constant)예요.`,
    });
  };

  // Build trace info for graph
  const traceKeys = [], traceLabels = [], traceColors2 = [];
  cells.forEach((cell, ci) => {
    template.nodes.filter(n => n.type !== "input").forEach((nd, ni) => {
      const key = `${cell.id}_${nd.id}`;
      traceKeys.push(key);
      traceLabels.push(cells.length > 1 ? `${cell.name}:${nd.name}` : nd.name);
      traceColors2.push(C.traces[(ci * template.nodes.length + ni) % C.traces.length]);
    });
  });

  const tabs = [
    { key: "design", label: "🔧 Design", c: C.accent },
    { key: "tissue", label: "🧫 Tissue", c: C.accent5 },
    { key: "simulate", label: "▶ Simulate", c: C.accent4 },
    { key: "analysis", label: "📐 Analysis", c: C.accent3 },
  ];

  return (
    <div style={{ minHeight: "100vh", background: C.bg, color: C.text, fontFamily: "'JetBrains Mono','Fira Code',monospace", padding: 14 }}>
      <link href="https://fonts.googleapis.com/css2?family=JetBrains+Mono:wght@300;400;600;700&display=swap" rel="stylesheet" />

      <div style={{ textAlign: "center", marginBottom: 12, paddingTop: 4 }}>
        <div style={{ fontSize: 8, letterSpacing: 5, color: C.dim, textTransform: "uppercase" }}>Systems Biology Workbench</div>
        <h1 style={{ fontSize: 18, fontWeight: 700, margin: "3px 0", background: `linear-gradient(135deg,${C.accent},${C.accent5})`, WebkitBackgroundClip: "text", WebkitTextFillColor: "transparent" }}>
          Multi-Cell Bio-Circuit Builder
        </h1>
        <div style={{ fontSize: 9, color: C.dimmer }}>Design · Duplicate · Simulate · Analyze</div>
      </div>

      {/* Presets */}
      <div style={{ display: "flex", flexWrap: "wrap", gap: 4, marginBottom: 10, justifyContent: "center" }}>
        {Object.entries(PRESETS).map(([k, v]) => <Btn key={k} onClick={() => loadPreset(k)} small>{v.name}</Btn>)}
      </div>

      {/* Tabs */}
      <div style={{ display: "flex", gap: 3, marginBottom: 10 }}>
        {tabs.map(t => (
          <button key={t.key} onClick={() => setTab(t.key)} style={{
            flex: 1, padding: "7px 0", fontSize: 10, fontWeight: 600, fontFamily: "inherit",
            background: tab === t.key ? `${t.c}15` : C.panel, border: `1.5px solid ${tab === t.key ? t.c : C.border}`,
            color: tab === t.key ? t.c : C.dim, borderRadius: 7, cursor: "pointer",
          }}>{t.label}</button>
        ))}
      </div>

      {/* Network Diagram */}
      <Sec title="Network" icon="⬡" color={C.accent5} extra={
        <span style={{ fontSize: 8, color: C.dimmer, fontWeight: 400 }}>Click nodes/arrows for equations</span>
      }>
        <CellNetworkDiagram template={template} cells={cells} interEdges={interEdges}
          onClickNode={(nd) => { setSelectedNode(nd.id); showNodeEq(nd); }}
          onClickEdge={showEdgeEq}
          onClickInterEdge={(ie, idx) => {
            const fn = template.nodes.find(n => n.id === ie.fromNode);
            const tn = template.nodes.find(n => n.id === ie.toNode);
            const fc = cells.find(c => c.id === ie.fromCell);
            const tc = cells.find(c => c.id === ie.toCell);
            setEqModal({
              title: `${fc?.name}:${fn?.name} ${ie.type === "activate" ? "→" : "⊣"} ${tc?.name}:${tn?.name} (Intercellular)`,
              color: ie.type === "activate" ? C.accent3 : C.accent4,
              equations: [
                ie.type === "activate"
                  ? `f(x) = w · xⁿ / (Kⁿ + xⁿ)`
                  : `f(x) = Kⁿ / (Kⁿ + xⁿ)`,
                `K = ${ie.K} nM, n = ${ie.n}, w = ${ie.w}`,
              ],
              description: `Intercellular signal from ${fc?.name} to ${tc?.name}. This models diffusible signaling molecules (e.g., AHL in quorum sensing, morphogens).`,
              korean: "세포 간 신호 전달: 쿼럼 센싱의 AHL이나 모르포겐 같은 확산성 신호분자를 모델링합니다.",
            });
          }}
          selectedNode={selectedNode}
        />
      </Sec>

      {/* ═══ DESIGN TAB ═══ */}
      {tab === "design" && (
        <>
          <Sec title="Add Node" icon="➕" color={C.accent}>
            <div style={{ display: "flex", gap: 4, flexWrap: "wrap", marginBottom: 8 }}>
              {Object.entries(NODE_META).map(([k, v]) => <Btn key={k} onClick={() => setNewType(k)} active={newType === k} color={v.color} small>{v.icon} {v.label}</Btn>)}
            </div>
            <div style={{ display: "flex", gap: 5 }}>
              <Inp value={newName} onChange={setNewName} placeholder="Name (e.g. LacI, GFP...)" style={{ flex: 1 }} />
              <Btn onClick={addNode}>Add</Btn>
            </div>
          </Sec>

          <Sec title="Add Connection" icon="🔀" color={C.accent2}>
            <div style={{ display: "flex", gap: 4, marginBottom: 8 }}>
              <Btn onClick={() => setEtype("activate")} active={etype === "activate"} color={C.accent} small>→ Activate</Btn>
              <Btn onClick={() => setEtype("repress")} active={etype === "repress"} color={C.accent2} small>⊣ Repress</Btn>
            </div>
            <div style={{ display: "flex", gap: 4, alignItems: "center" }}>
              <select value={efrom} onChange={e => setEfrom(e.target.value)} style={{ flex: 1, background: C.bg, border: `1px solid ${C.border}`, borderRadius: 5, color: C.text, padding: 5, fontSize: 10, fontFamily: "inherit" }}>
                <option value="">From...</option>
                {template.nodes.map(n => <option key={n.id} value={n.id}>{n.name}</option>)}
              </select>
              <span style={{ color: etype === "activate" ? C.accent : C.accent2, fontSize: 14 }}>{etype === "activate" ? "→" : "⊣"}</span>
              <select value={eto} onChange={e => setEto(e.target.value)} style={{ flex: 1, background: C.bg, border: `1px solid ${C.border}`, borderRadius: 5, color: C.text, padding: 5, fontSize: 10, fontFamily: "inherit" }}>
                <option value="">To...</option>
                {template.nodes.map(n => <option key={n.id} value={n.id}>{n.name}</option>)}
              </select>
              <Btn onClick={addEdge} color={C.accent2}>Add</Btn>
            </div>
          </Sec>

          <Sec title={`Nodes (${template.nodes.length})`} icon="📋" color={C.dim}>
            {template.nodes.map(n => (
              <div key={n.id} style={{ display: "flex", alignItems: "center", gap: 6, padding: "5px 6px", marginBottom: 3, background: selectedNode === n.id ? `${NODE_META[n.type].color}0a` : "transparent", borderRadius: 5, cursor: "pointer" }}
                onClick={() => { setSelectedNode(n.id); showNodeEq(n); }}>
                <span style={{ fontSize: 12 }}>{NODE_META[n.type].icon}</span>
                <span style={{ flex: 1, fontSize: 10, color: NODE_META[n.type].color }}>{n.name}</span>
                <span style={{ fontSize: 8, color: C.dimmer }}>α={n.alpha} β={n.beta}</span>
                <button onClick={e => { e.stopPropagation(); removeNode(n.id); }} style={{ background: "none", border: "none", color: C.dimmer, cursor: "pointer", fontSize: 11 }}>✕</button>
              </div>
            ))}
            {template.nodes.length > 0 && selectedNode && (() => {
              const nd = template.nodes.find(n => n.id === selectedNode);
              if (!nd || nd.type === "input") return null;
              return (
                <div style={{ background: C.bg, borderRadius: 6, padding: 10, marginTop: 8 }}>
                  <div style={{ fontSize: 10, color: NODE_META[nd.type].color, marginBottom: 6 }}>{nd.name} — Parameters</div>
                  <div style={{ display: "grid", gridTemplateColumns: "1fr 1fr 1fr", gap: 6 }}>
                    <div><Lbl>α (nM/min)</Lbl><Inp type="number" value={nd.alpha} onChange={v => updateNode(nd.id, { alpha: v })} step="0.5" min="0" /></div>
                    <div><Lbl>β (1/min)</Lbl><Inp type="number" value={nd.beta} onChange={v => updateNode(nd.id, { beta: v })} step="0.1" min="0" /></div>
                    <div><Lbl>Init (nM)</Lbl><Inp type="number" value={nd.init} onChange={v => updateNode(nd.id, { init: v })} step="0.5" min="0" /></div>
                  </div>
                </div>
              );
            })()}
          </Sec>

          <Sec title={`Edges (${template.edges.length})`} icon="🔗" color={C.dim}>
            {template.edges.map((e, i) => {
              const fn = template.nodes.find(n => n.id === e.from), tn = template.nodes.find(n => n.id === e.to);
              return (
                <div key={i} style={{ background: C.bg, borderRadius: 5, padding: 6, marginBottom: 4 }}>
                  <div style={{ display: "flex", alignItems: "center", gap: 4, fontSize: 10, marginBottom: 4 }}>
                    <span>{fn?.name}</span>
                    <span style={{ color: e.type === "activate" ? C.accent : C.accent2 }}>{e.type === "activate" ? "→" : "⊣"}</span>
                    <span style={{ flex: 1 }}>{tn?.name}</span>
                    <button onClick={() => showEdgeEq(e)} style={{ background: "none", border: `1px solid ${C.border}`, color: C.dim, cursor: "pointer", fontSize: 8, borderRadius: 4, padding: "1px 6px", fontFamily: "inherit" }}>f(x)</button>
                    <button onClick={() => removeEdge(i)} style={{ background: "none", border: "none", color: C.dimmer, cursor: "pointer", fontSize: 11 }}>✕</button>
                  </div>
                  <div style={{ display: "grid", gridTemplateColumns: "1fr 1fr 1fr", gap: 4 }}>
                    <div><Lbl>K (nM)</Lbl><Inp type="number" value={e.K} onChange={v => updateEdge(i, { K: v })} step="0.25" min="0.1" /></div>
                    <div><Lbl>n (Hill)</Lbl><Inp type="number" value={e.n} onChange={v => updateEdge(i, { n: v })} step="0.5" min="0.5" /></div>
                    <div><Lbl>w (fold)</Lbl><Inp type="number" value={e.w} onChange={v => updateEdge(i, { w: v })} step="0.1" min="0" /></div>
                  </div>
                </div>
              );
            })}
          </Sec>
        </>
      )}

      {/* ═══ TISSUE TAB ═══ */}
      {tab === "tissue" && (
        <>
          <Sec title="Cell Instances" icon="🧫" color={C.accent5}>
            <div style={{ fontSize: 10, color: C.dim, marginBottom: 8, lineHeight: 1.6 }}>
              Each cell contains the same internal circuit. Duplicate cells and add intercellular signals.
            </div>
            {cells.map((c, i) => (
              <div key={c.id} style={{ display: "flex", alignItems: "center", gap: 6, padding: "6px 8px", marginBottom: 4, background: C.bg, borderRadius: 6, border: `1px solid ${C.cellColors[i % C.cellColors.length]}30` }}>
                <span style={{ width: 8, height: 8, borderRadius: "50%", background: C.cellColors[i % C.cellColors.length] }} />
                <Inp value={c.name} onChange={v => renameCell(c.id, v)} style={{ flex: 1, background: "transparent", border: "none", fontSize: 11 }} />
                {cells.length > 1 && (
                  <button onClick={() => removeCell(c.id)} style={{ background: "none", border: "none", color: C.dimmer, cursor: "pointer", fontSize: 11 }}>✕</button>
                )}
              </div>
            ))}
            <Btn onClick={addCell} color={C.accent5} style={{ marginTop: 6 }}>+ Duplicate Cell</Btn>
          </Sec>

          {cells.length > 1 && (
            <Sec title="Intercellular Connections" icon="📡" color={C.accent3}>
              <div style={{ fontSize: 9, color: C.dim, marginBottom: 8 }}>Connect outputs of one cell to inputs of another (e.g., quorum sensing, morphogens).</div>
              <div style={{ display: "flex", gap: 4, marginBottom: 6 }}>
                <Btn onClick={() => setIEtype("activate")} active={ieType === "activate"} color={C.accent3} small>→ Activate</Btn>
                <Btn onClick={() => setIEtype("repress")} active={ieType === "repress"} color={C.accent4} small>⊣ Repress</Btn>
              </div>
              <div style={{ display: "grid", gridTemplateColumns: "1fr 1fr", gap: 4, marginBottom: 6 }}>
                <div>
                  <Lbl>From Cell</Lbl>
                  <select value={iefromCell} onChange={e => setIEfromCell(e.target.value)} style={{ width: "100%", background: C.bg, border: `1px solid ${C.border}`, borderRadius: 5, color: C.text, padding: 4, fontSize: 10, fontFamily: "inherit" }}>
                    <option value="">Cell...</option>
                    {cells.map(c => <option key={c.id} value={c.id}>{c.name}</option>)}
                  </select>
                </div>
                <div>
                  <Lbl>From Node</Lbl>
                  <select value={iefromNode} onChange={e => setIEfromNode(e.target.value)} style={{ width: "100%", background: C.bg, border: `1px solid ${C.border}`, borderRadius: 5, color: C.text, padding: 4, fontSize: 10, fontFamily: "inherit" }}>
                    <option value="">Node...</option>
                    {template.nodes.map(n => <option key={n.id} value={n.id}>{n.name}</option>)}
                  </select>
                </div>
                <div>
                  <Lbl>To Cell</Lbl>
                  <select value={ietoCell} onChange={e => setIEtoCell(e.target.value)} style={{ width: "100%", background: C.bg, border: `1px solid ${C.border}`, borderRadius: 5, color: C.text, padding: 4, fontSize: 10, fontFamily: "inherit" }}>
                    <option value="">Cell...</option>
                    {cells.map(c => <option key={c.id} value={c.id}>{c.name}</option>)}
                  </select>
                </div>
                <div>
                  <Lbl>To Node</Lbl>
                  <select value={ietoNode} onChange={e => setIEtoNode(e.target.value)} style={{ width: "100%", background: C.bg, border: `1px solid ${C.border}`, borderRadius: 5, color: C.text, padding: 4, fontSize: 10, fontFamily: "inherit" }}>
                    <option value="">Node...</option>
                    {template.nodes.map(n => <option key={n.id} value={n.id}>{n.name}</option>)}
                  </select>
                </div>
              </div>
              <Btn onClick={addInterEdge} color={C.accent3}>Add Intercellular Link</Btn>

              {interEdges.length > 0 && (
                <div style={{ marginTop: 10 }}>
                  {interEdges.map((ie, i) => {
                    const fc = cells.find(c => c.id === ie.fromCell), tc = cells.find(c => c.id === ie.toCell);
                    const fn = template.nodes.find(n => n.id === ie.fromNode), tn = template.nodes.find(n => n.id === ie.toNode);
                    return (
                      <div key={i} style={{ background: C.bg, borderRadius: 5, padding: 6, marginBottom: 4 }}>
                        <div style={{ display: "flex", alignItems: "center", gap: 4, fontSize: 10, marginBottom: 4 }}>
                          <span style={{ color: C.accent3 }}>{fc?.name}:{fn?.name}</span>
                          <span style={{ color: ie.type === "activate" ? C.accent3 : C.accent4 }}>{ie.type === "activate" ? "→" : "⊣"}</span>
                          <span style={{ flex: 1, color: C.accent3 }}>{tc?.name}:{tn?.name}</span>
                          <button onClick={() => removeInterEdge(i)} style={{ background: "none", border: "none", color: C.dimmer, cursor: "pointer" }}>✕</button>
                        </div>
                        <div style={{ display: "grid", gridTemplateColumns: "1fr 1fr 1fr", gap: 4 }}>
                          <div><Lbl>K (nM)</Lbl><Inp type="number" value={ie.K} onChange={v => updateInterEdge(i, { K: v })} step="0.25" min="0.1" /></div>
                          <div><Lbl>n (Hill)</Lbl><Inp type="number" value={ie.n} onChange={v => updateInterEdge(i, { n: v })} step="0.5" min="0.5" /></div>
                          <div><Lbl>w (fold)</Lbl><Inp type="number" value={ie.w} onChange={v => updateInterEdge(i, { w: v })} step="0.1" min="0" /></div>
                        </div>
                      </div>
                    );
                  })}
                </div>
              )}
            </Sec>
          )}
        </>
      )}

      {/* ═══ SIMULATE TAB ═══ */}
      {tab === "simulate" && (
        <>
          <Sec title="Input Signal" icon="⚡" color={C.accent4}>
            <div style={{ display: "flex", flexWrap: "wrap", gap: 4, marginBottom: 8 }}>
              {Object.entries(INPUTS).map(([k, v]) => <Btn key={k} onClick={() => setInputSig(k)} active={inputSig === k} color={C.accent4} small>{v}</Btn>)}
            </div>
            <div style={{ display: "grid", gridTemplateColumns: "1fr 1fr 1fr", gap: 6 }}>
              <div><Lbl>Start (min)</Lbl><Inp type="number" value={inputOn} onChange={setInputOn} step="1" min="0" /></div>
              <div><Lbl>End (min)</Lbl><Inp type="number" value={inputOff} onChange={setInputOff} step="1" min="0" /></div>
              <div><Lbl>Total (min)</Lbl><Inp type="number" value={tMax} onChange={setTMax} step="5" min="10" /></div>
            </div>
          </Sec>

          <div style={{ display: "flex", gap: 5, marginBottom: 10 }}>
            <button onClick={runSim} disabled={template.nodes.length === 0} style={{
              flex: 1, padding: 11, fontSize: 12, fontWeight: 700, fontFamily: "inherit",
              background: `linear-gradient(135deg,${C.accent}25,${C.accent5}25)`, border: `2px solid ${C.accent}`,
              borderRadius: 9, color: C.accent, cursor: template.nodes.length === 0 ? "not-allowed" : "pointer",
            }}>▶ Run Simulation</button>
            {zoomRange && <Btn onClick={() => setZoomRange(null)} color={C.accent5} style={{ padding: "11px 14px" }}>🔍 Reset Zoom</Btn>}
          </div>

          {simData && (
            <Sec title={zoomRange ? `Zoomed: ${zoomRange[0].toFixed(1)} – ${zoomRange[1].toFixed(1)} min` : "Time Response"} icon="📈" color={C.accent}>
              <Graph data={simData} traceKeys={traceKeys} traceLabels={traceLabels} traceColors={traceColors2} zoomRange={zoomRange} setZoomRange={setZoomRange} />
              <div style={{ fontSize: 8, color: C.dimmer, textAlign: "center", marginTop: 4 }}>Drag to select a region · Hover for values</div>
            </Sec>
          )}
        </>
      )}

      {/* ═══ ANALYSIS TAB ═══ */}
      {tab === "analysis" && (
        <>
          <div style={{ display: "flex", gap: 5, marginBottom: 10 }}>
            <button onClick={runStability} disabled={template.nodes.filter(n => n.type !== "input").length === 0} style={{
              flex: 1, padding: 11, fontSize: 12, fontWeight: 700, fontFamily: "inherit",
              background: `linear-gradient(135deg,${C.accent3}25,${C.accent}25)`, border: `2px solid ${C.accent3}`,
              borderRadius: 9, color: C.accent3, cursor: "pointer",
            }}>📐 Run Stability Analysis</button>
          </div>

          {stability && (
            <>
              <Sec title="Stability Classification" icon="🎯" color={stability.classification.color}>
                <div style={{ fontSize: 16, fontWeight: 700, color: stability.classification.color, marginBottom: 6 }}>
                  {stability.classification.label}
                </div>
                <div style={{ fontSize: 10, color: C.dim, lineHeight: 1.7, marginBottom: 10 }}>
                  {stability.classification.type === "stable_node" && "All eigenvalues have negative real parts with no imaginary components. The system returns to equilibrium monotonically after perturbation."}
                  {stability.classification.type === "stable_spiral" && "Eigenvalues have negative real parts with nonzero imaginary parts. The system returns to equilibrium with damped oscillations."}
                  {stability.classification.type === "unstable_node" && "At least one eigenvalue has a positive real part. The system diverges from this equilibrium point."}
                  {stability.classification.type === "unstable_spiral" && "Complex eigenvalues with positive real parts. The system spirals away from equilibrium with growing oscillations."}
                  {stability.classification.type === "center" && "Purely imaginary eigenvalues indicate sustained oscillations (limit cycle behavior)."}
                  {stability.classification.type === "saddle" && "Mixed positive and negative real eigenvalues indicate a saddle point — stable in some directions, unstable in others."}
                </div>
                <div style={{ background: `${C.accent}08`, borderRadius: 6, padding: 8, fontSize: 9, color: C.dimmer, lineHeight: 1.6 }}>
                  🇰🇷 {stability.classification.type === "stable_node" && "모든 고유값의 실수부가 음수이고 허수부가 0입니다. 섭동 후 단조롭게 평형으로 복귀합니다."}
                  {stability.classification.type === "stable_spiral" && "고유값이 음의 실수부 + 허수부를 가집니다. 감쇠 진동하며 평형으로 돌아갑니다."}
                  {stability.classification.type === "unstable_spiral" && "고유값의 실수부가 양수이고 허수부가 있어 진동이 점점 커지며 발산합니다."}
                  {stability.classification.type === "center" && "순허수 고유값은 지속적인 진동(리밋 사이클)을 의미합니다."}
                  {stability.classification.type === "saddle" && "양수/음수 고유값이 섞여 있어 안장점입니다 — 일부 방향으로만 안정적입니다."}
                  {stability.classification.type === "unstable_node" && "양의 실수부 고유값이 존재하여 평형점에서 발산합니다."}
                </div>
              </Sec>

              <Sec title="Eigenvalues" icon="λ" color={C.accent3}>
                <div style={{ display: "grid", gap: 6 }}>
                  {stability.eigs.map((e, i) => (
                    <div key={i} style={{ display: "flex", alignItems: "center", gap: 8, padding: "6px 8px", background: C.bg, borderRadius: 5 }}>
                      <span style={{ fontSize: 11, color: C.accent3, fontWeight: 700 }}>λ{i + 1}</span>
                      <span style={{ fontSize: 11, color: e.re < 0 ? C.accent : C.accent2 }}>
                        {e.re.toFixed(4)}{Math.abs(e.im) > 1e-6 ? ` ${e.im > 0 ? "+" : "−"} ${Math.abs(e.im).toFixed(4)}i` : ""}
                      </span>
                      <span style={{ fontSize: 8, color: C.dimmer, marginLeft: "auto" }}>
                        Re {e.re < 0 ? "< 0 ✓" : "> 0 ✗"}
                      </span>
                    </div>
                  ))}
                </div>
              </Sec>

              <Sec title="Steady State" icon="⊜" color={C.dim}>
                <div style={{ display: "grid", gridTemplateColumns: "1fr 1fr", gap: 4 }}>
                  {stability.keys.map((k, i) => {
                    const parts = k.split("_");
                    const cell = cells.find(c => c.id === parts[0]);
                    const nd = template.nodes.find(n => n.id === parts[1]);
                    return (
                      <div key={k} style={{ padding: "4px 8px", background: C.bg, borderRadius: 5, fontSize: 10 }}>
                        <span style={{ color: C.dim }}>{cells.length > 1 ? `${cell?.name}:` : ""}{nd?.name} = </span>
                        <span style={{ color: C.accent }}>{(stability.ss[k] || 0).toFixed(3)} nM</span>
                      </div>
                    );
                  })}
                </div>
              </Sec>

              <Sec title="Jacobian Matrix" icon="J" color={C.dim}>
                <div style={{ overflowX: "auto" }}>
                  <table style={{ borderCollapse: "collapse", fontSize: 9, fontFamily: "inherit" }}>
                    <thead>
                      <tr>
                        <th style={{ padding: 4 }}></th>
                        {stability.keys.map(k => {
                          const nd = template.nodes.find(n => n.id === k.split("_")[1]);
                          return <th key={k} style={{ padding: "3px 6px", color: C.dim, fontWeight: 500 }}>{nd?.name}</th>;
                        })}
                      </tr>
                    </thead>
                    <tbody>
                      {stability.J.map((row, i) => (
                        <tr key={i}>
                          <td style={{ padding: "3px 6px", color: C.dim }}>{template.nodes.find(n => n.id === stability.keys[i].split("_")[1])?.name}</td>
                          {row.map((v, j) => (
                            <td key={j} style={{ padding: "3px 6px", color: Math.abs(v) < 1e-6 ? C.dimmer : v < 0 ? C.accent : C.accent2, textAlign: "right" }}>
                              {v.toFixed(3)}
                            </td>
                          ))}
                        </tr>
                      ))}
                    </tbody>
                  </table>
                </div>
              </Sec>
            </>
          )}

          <Sec title="Electronic Circuit Analog" icon="⏚" color={C.accent4}>
            <ElectronicAnalog template={template} />
            <div style={{ fontSize: 9, color: C.dimmer, textAlign: "center", marginTop: 4, lineHeight: 1.6 }}>
              🇰🇷 각 유전자/단백질 노드는 비선형 증폭기(Op-Amp)로, 분해 반응은 RC 회로의 방전으로 표현됩니다.
              피드백 루프는 전자회로의 피드백 저항과 동치예요.
            </div>
          </Sec>
        </>
      )}

      <EquationModal data={eqModal} onClose={() => setEqModal(null)} />

      <div style={{ textAlign: "center", fontSize: 7, color: C.dimmer, marginTop: 10, paddingBottom: 14 }}>
        Multi-Cell Bio-Circuit Builder · Euler ODE · Hill Kinetics · QR Eigenvalue Analysis
      </div>
    </div>
  );
}
