# The eight external ideas are one structure: the apex-prime partition function split into stabilizer and magic

**kind-pasteur-2026-06-20-S22**, converging with codex-S66 (HYP-2705) and building on the score
partition function (THM-554/555) and the apex-prime gas (S21). On the user's batch of new ideas —
Gibbs measure, Arnold's cat map, complete-graph crossing number, Hopfield vs 2-layer nets, the
double-slit / Feynman-propagator-as-weights, the Hebbian rate neuron, the road coloring theorem, the
Fubini-Study metric, and Clifford+T.

## They do not give eight routes. They give one architecture.

Independently, two agents reduced the whole list to the same object. A tournament's `H`/OCF and the
LRC sector cover are both a **partition function of a phase gas on `Z/p`** (`p` the apex prime), and
that partition function splits into a **stabilizer layer** (cut space / decorrelated occupancy /
degree-≤2, polynomial-time) and a **magic layer** (cycle space / signed relation deviation /
degree-≥3, the open content). Every external idea is a *lens on one of those two layers or the split
between them*:

- **Gibbs measure** = the partition function itself. `Z_n = (∏x_v)∏(x_a+x_b)`; the LRC cover is its
  apex twin. Low temperature selects the obstruction — and (codex) the LRC depth quotient's
  obstruction *vanishes* after two decorrelated far hits, so the residual energy lives entirely in the
  magic (signed relation) layer.
- **Clifford + T** = the split, made *precise*. `c3` is degree-2 in the tile bits — a GF(2) quadratic
  form — so it is Clifford/Gauss-sum tractable; its parity is *literally* a stabilizer amplitude
  (THM-555's `E[(-1)^c3]=2^{-rank/2}`, rank `2⌊(n-1)/2⌋`, the Gottesman-Knill formula, HYP-2707 proved
  core). `c5`, `α₂`, `H` are higher-degree = magic. The score→OCF wall **is** the stabilizer→magic
  boundary. LRC's "magic rank" is the degree/stabilizer-rank of the mod-7 relation profile.
- **Feynman propagator / double slit** = the magic layer's signature. `H = Σ Hamiltonian paths` is a
  path-sum; the interference (the alternating OCF cancellation, the signed cover deviation) is exactly
  the magic correction to the free/decorrelated amplitude. Slits = endpoint classes; amplitude = the
  quasisymmetric weight.
- **Fubini-Study metric** = the right metric on the magic layer, because the phase profiles are
  defined up to global phase and scale (codex: the death-chain amplitudes' projective middle-mass is
  an exact rational; the open lemma is an FS-angle bound between the true and decorrelated profiles).
- **Road coloring** = the certificate that the magic residual collapses. It fails on the *invertible*
  tiling hypercube (every tile-flip is a bijection — no synchronization), but works on the
  *non-invertible* residual-sector deletion automaton: a word synchronizes a missed set iff every
  missed color appears, length 6 sufficient.
- **Arnold's cat map** = the model for *proving* the stabilizer layer decorrelates: a Markov-partition
  / transfer-operator picture splitting stable (cut/free) from unstable (cycle/magic) directions, with
  a spectral gap outside the finite resonant atlas. The hyperbolicity is the mixing that makes wide
  runners decorrelate.
- **Hopfield / Hebbian** = the empirical density matrix of the magic layer. The symmetric co-firing
  matrix of missed sectors / relation packets; its dominant low-rank modes are the obstructions that
  spend the Gibbs currency. The Hebbian rate update is the vertex-insertion recursion on the cut layer.
- **Crossing number** = a near-miss that sharpens the picture: the circular tournament's raw crossing
  number is `C(n,4)`, *orientation-independent* — it sees nothing. Only an *alternating/parity*
  skeleton (a stabilizer object) plus the signed magic correction carries tournament data (codex's
  K14 squarefree carrier is the full interaction carrier).

## What it ultimately represents

The repo spent a year discovering, on two fronts, that the same wall separates the cheap from the
hard: cut space from cycle space, score from `H`, decorrelated from resonant, finite-check from
true-wide. The external ideas name that wall in the vocabulary of physics and computation — it is the
**stabilizer/magic boundary of a `Z/p` partition function**, the line where quadratic (Gauss-summable,
Clifford, free-energy) gives way to cubic-and-up (interference, T-gates, interaction). The H-maximizer
is a maximal-cycle (high-magic) state; the lonely-runner crux is a minimal-magic-rank certificate. The
proof on both sides is the same program: *get the stabilizer/occupancy extremum exactly* (regular
score; single block; THM-555/557), *then spend the comfortable margin on a lossy bound for the magic*
(`α₂`; OPEN-Q-108). Eight ideas, one gas, one wall.
[[the-apex-prime-partition-function-tournaments-and-runners-are-one-gas]] ·
[[the-three-tiling-recurrences-are-one-partition-function]] · [[everything-is-the-triangle]]
