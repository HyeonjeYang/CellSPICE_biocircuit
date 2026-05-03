# CellSPICE_biocircuit

CellSPICE bridges **electronic circuit theory** and **systems biology**, letting you design, simulate, and analyze biological circuits the same way electrical engineers design chips — except your components are genes, proteins, and signaling molecules.

Traditional systems biology tools treat biological networks as abstract graphs. CellSPICE treats them as **engineered circuits** — because that's what they are.

---

## Electronic-Biological Analogies

| Biology | Electronics | Equation |
|---|---|---|
| Gene (Hill activation) | Op-Amp | f(x) = w\*x^n / (K^n + x^n) |
| Mutual repression | Zener Diode | Clamps to bistable states |
| Enzyme cascade | BJT Transistor | Signal amplification with threshold |
| Irreversible reaction | Diode | One-way signal flow |
| Protein degradation | Resistor | tau = 1/beta = RC |
| Protein accumulation | Capacitor | Integrates production |
| HH Na+ channel | NMOS Transistor | Voltage-gated turn-on |
| HH K+ channel | PMOS Transistor | Delayed negative feedback |

---

## Features

**Circuit Design** — Add genes, proteins, mRNA, input signals. Wire them with activation or repression. Click any node or edge to see the governing ODE and Hill equation with all parameters.

**Multi-Cell Tissue Mode** — Design one cell's internal circuit, duplicate it, and add intercellular connections (quorum sensing, morphogens, juxtacrine signaling).

**Time-Domain Simulation** — Euler ODE solver. Step, Pulse, Ramp, Sine, Square inputs. Drag-to-zoom on time series, hover for values.

**Frequency Response** — Linearize around steady state. Bode plot (magnitude and phase). Transfer function H(s) = C(sI-A)^(-1)B with poles and DC gain.

**Stability Analysis** — Jacobian matrix, eigenvalues via QR iteration, automatic classification (Stable Node, Stable Spiral, Unstable, Saddle, Center). Steady-state readout.

**Cable Theory** — Neuron modeled as RC transmission line. Passive cable equation with configurable Rm, Cm, Ra. Optional Hodgkin-Huxley ion channels (Na+/K+). Voltage propagation visualization.

---

## Quick Start

```bash
git clone https://github.com/<your-username>/CellSPICE_biocircuit.git
cd CellSPICE_biocircuit
npm install
npm run dev
```

Open `http://localhost:5173`.

---

## Built-in Presets

| Preset | Circuit | Electronic Analog |
|---|---|---|
| Negative Feedback | Self-repressing gene | RC Low-Pass Filter |
| Toggle Switch | Mutual repression | Bistable Latch |
| Repressilator | 3-gene cyclic repression | Ring Oscillator |
| Cascade | Ligand-Receptor-Kinase-TF | Multi-stage Amplifier |

---

## The Math

**Node dynamics:**
```
d[X]/dt = alpha * prod(repression terms) + sum(activation terms) - beta * [X]
```

**Hill functions:**
```
H+(x) = w * x^n / (K^n + x^n)       (activation)
H-(x) = K^n / (K^n + x^n)            (repression)
```

**Transfer function:**
```
H(s) = C * (sI - A)^(-1) * B
```

**Cable equation:**
```
tau_m * dV/dt = lambda^2 * d2V/dx2 - (V - Vrest) + Rm * Iext
```

**Stability:** Eigenvalues of Jacobian A at steady state. Re(lambda) < 0 = stable, Re(lambda) > 0 = unstable, Im(lambda) != 0 = oscillatory.

---

## Tech Stack

- React 18, Vite
- Canvas API for plotting
- SVG for network diagrams
- QR iteration for eigenvalues
- Complex matrix inversion for Bode analysis
- No external dependencies beyond React

---

## Roadmap

- [ ] Stochastic simulation (Gillespie algorithm)
- [ ] Phase portrait visualization
- [ ] Bifurcation diagrams
- [ ] Spatial diffusion for morphogen gradients
- [ ] Sensitivity analysis

---

## License

MIT License. See [LICENSE](LICENSE) for details.

---

Code generated with assistance from [Claude] (Anthropic). Idea, review, and revision by Hyeonje Yang.
