# CellSPICE_biocircuit

**A circuit-theoretic simulator for systems biology — from single-cell gene networks to multi-cell tissue dynamics.**

CellSPICE bridges **electronic circuit theory** and **systems biology**, letting you design, simulate, and analyze biological circuits the same way electrical engineers design chips — except your components are genes, proteins, and signaling molecules.

| Biology | ≡ Electronics |
|---|---|
| Negative feedback loop | Low-pass filter (RC) |
| Positive feedback loop | Schmitt trigger / Latch |
| Feedforward loop | Band-pass filter (RLC) |
| Repressilator | Ring oscillator |
| Hill function | Nonlinear amplifier gain |
| Protein degradation | RC discharge (τ = 1/β) |
| Gene expression | Voltage-controlled source |

Traditional systems biology tools treat biological networks as abstract graphs. CellSPICE_biocirtui treats them as **engineered circuits** — because that's what they are.

---

## Features

### Circuit Design
- Add biological components: **genes**, **proteins**, **mRNA**, **input signals**
- Wire them with **activation (→)** or **repression (⊣)** connections
- Every arrow displays its constants with units (K = 2 nM, n = 2)
- Click any node or edge to see the governing **ODE** and **Hill equation**

### Multi-Cell Tissue Mode
- Design one cell's internal circuit, then **duplicate it** as many times as needed
- Add **intercellular connections** — model quorum sensing (AHL), morphogen gradients, juxtacrine signaling
- Scale from single-cell behavior to emergent **tissue-level dynamics**

### Simulation
- Euler-method ODE solver with configurable time step
- 5 input signal types: Step, Pulse, Ramp, Sine, Square
- Interactive time-series plot with **drag-to-zoom** and **hover tooltips**
- Track all species concentrations across all cells simultaneously

### Stability Analysis
- Automatic **Jacobian matrix** computation (numerical)
- **Eigenvalue** calculation via QR iteration
- Classification: Stable Node, Stable Spiral, Unstable Node, Saddle Point, Center, etc.
- Steady-state concentration readout

### Electronic Circuit Analog
- Automatic mapping of biological network → equivalent electronic schematic
- Genes as op-amps, degradation as RC circuits, feedback as wired paths
- Side-by-side comparison of bio and electronic representations

---

## Quick Start

CellSPICE runs as a standalone React component. No backend required.

```bash
# Clone
git clone https://github.com/<your-username>/CellSPICE_biocircuit.git
cd CellSPICE_biocircuit

# Install
npm install

# Run
npm run dev
```

Or simply open `multi-cell-circuit-builder.jsx` in any React environment (Vite, Next.js, Create React App, or Claude Artifacts).

---

## Built-in Presets

| Preset | Circuit Type | Electronic Analog |
|---|---|---|
| **Negative Feedback** | Self-repressing gene | RC Low-Pass Filter |
| **Toggle Switch** | Mutual repression (A ⊣ B, B ⊣ A) | Bistable Latch |
| **Repressilator** | 3-gene cyclic repression | 3-Stage Ring Oscillator |
| **Signaling Cascade** | Ligand → Receptor → Kinase → TF | Multi-stage Amplifier |

Load any preset, modify parameters, duplicate cells, and explore.

---

## The Math

### Node dynamics (ODE)
Each non-input node follows:

```
d[X]/dt = α · ∏(repression terms) + Σ(activation terms) − β · [X]
```

### Hill activation
```
H⁺(x) = w · xⁿ / (Kⁿ + xⁿ)
```

### Hill repression
```
H⁻(x) = Kⁿ / (Kⁿ + xⁿ)
```

### Stability
Jacobian **J** is computed numerically at steady state. Eigenvalues λ of **J** determine local stability:
- All Re(λ) < 0 → **Stable** (perturbations decay)
- Any Re(λ) > 0 → **Unstable** (perturbations grow)
- Im(λ) ≠ 0 → **Oscillatory** component

---

## Use Cases

- **Synthetic biology** — Design and test gene circuits in silico before cloning
- **Systems biology courses** — Interactive teaching tool for gene regulation
- **iGEM projects** — Prototype circuit designs with quantitative predictions
- **Research** — Quick stability screening of network motifs
- **Interdisciplinary learning** — Bridge EE and biology with concrete analogies

---

## Roadmap

- [ ] Stochastic simulation (Gillespie algorithm / Chemical Langevin)
- [ ] Phase portrait visualization
- [ ] Bifurcation diagrams (parameter sweeps)
- [ ] SBML import/export
- [ ] Spatial diffusion model for morphogen gradients
- [ ] Sensitivity analysis (parameter robustness)
- [ ] Export to LaTeX equations
- [ ] Dark/light theme toggle

---

## Tech Stack

- **React** (functional components + hooks)
- **Canvas API** for high-performance plotting
- **SVG** for network diagrams and circuit schematics
- **QR iteration** for eigenvalue computation
- Zero external dependencies beyond React

---

## Contributing

Contributions are welcome! Whether it's new circuit motifs, better numerical methods, UI improvements, or documentation.

```bash
# Fork → Branch → Code → PR
git checkout -b feature/bifurcation-diagram
```

---

## License

MIT License. See [LICENSE](LICENSE) for details.

---
Code generated with assistance from Claude (Anthropic). Physics model design by Hyeonje Yang

## Future Developmet

Planned: Add frequency-response and Bode plot visualization.  
Proposed on 2026-05-02.
