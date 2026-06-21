# Cut versus cycle is the universal seam — every borrowed idea fell along it

**Session:** mac-mini-2026-06-20-S7. A brief arrived with ten ideas from far afield — Gibbs
measures and Arnold's cat map, Clifford gates and the T-gate, the road coloring theorem, Hopfield
networks, crossing numbers, Feynman propagators as neural-network weights, the Fubini–Study metric.
The instruction was to find connections to the project and generate hypotheses. Seven threads ran.
The striking thing is not how many connected, but that the ones with real content all connected to
the **same** structure, from opposite sides.

## The seam

The project's GF(2) decomposition of tournament space — base Hamiltonian path = **cut space**
(scores, the hierarchy), wiggly arcs = **cycle space** (the odd-cycle/OCF content) — turned out to
be the axis every borrowed idea was secretly about:

- **Gibbs / Hopfield / Ising.** The score partition function `Z_n` is a literal classical
  partition function; its free energy is an exact *sum of per-tile local terms*. And `c3` — the last
  score-determined OCF datum — is, exactly, a **frustrated 2-body Ising energy on the arc spins**
  (THM-559, proved this session): `c3 = n(n²−1)/24 − ½Σ_v(s_v−s̄)²`, with line-graph couplings
  `+½` when a cherry's shared vertex is an extreme, `−½` when it is the middle. The cut space is the
  *local* free energy; it is a textbook spin glass. That is the cut side.

- **Clifford / stabilizer / magic.** The tempting slogan was "cut cheap = Clifford-classical, cycle
  dear = magic/T-gate." It is **false** as stated — cut and cycle are an *orthogonal* homology pair,
  not a symplectic Lagrangian one (so no stabilizer phase space). But the true residue is on the
  same seam: the cut code and the cycle code are exact **MacWilliams duals**, and the cycle code is
  precisely the even-graph space `E_n`. The dual metagraph is the dual *code*.

- **Crossing numbers.** The 2-page crossing number of `K_n` is a quadratic form on the GF(2) edge
  space, and its optimum is attained on the **cycle space** — an even-graph page (verified through
  n=8). Guy's number lives on the cycle side.

- **Cat map / road coloring / Feynman interference.** All three were *dynamical or geometric* hopes
  — that the cover renormalizes like a hyperbolic toral map, synchronizes like a road coloring,
  or maximizes by constructive interference. All three are cleanly **refuted**, and for one shared
  reason: the relevant system is a measure-preserving rotation/permutation (scale-invariant, never
  collapsing, never mixing), so the difficulty is not in any single orbit. The difficulty is the
  **aggregate over the cycle space** — the same "irreducibly aggregate" wall the LRC cover hit. The
  refutations are the cycle side refusing to be reduced to a single dynamical or geometric object.

## What it means

THM-290 had already shown the *full* generating object `H(t)` is an antiferromagnet — a
**many-body** frustrated system living in the cycle space. THM-559 now shows its cut-space shadow
`c3` is a **2-body** Ising. So the project's oldest slogan, "cut cheap, cycle dear," has a precise
physical reading: **cut = 2-body (pairwise, local, classically tractable); cycle = many-body
(higher-order, non-local, the hard part).** Every external idea that carried real content was a
different name for one side of this — Gibbs/Ising/Hopfield and the score partition function name the
cut side; stabilizer codes, even graphs, crossing numbers, and the aggregate cover name the cycle
side. The ones that failed (cat map, synchronization, interference) failed because they tried to
make the cycle side *local* — a single orbit, a single reset word, a single amplitude — and the
cycle side is irreducibly global.

That is the lesson worth keeping. When the next borrowed metaphor arrives, the first question is not
"is it like a tournament?" but "which side of the cut/cycle seam does it live on?" The pairwise,
the local, the free-energy-additive belong to the scores. The frustrated, the many-body, the
aggregate belong to the cycles — and that is where both the OCF's depth and the lonely runner's
difficulty have been hiding all along.
