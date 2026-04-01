# Computational Irreducibility and Tournaments

**opus-2026-03-25-S339**

## The Core Mapping

Tournament theory maps precisely onto Wolfram's framework:

| Wolfram concept | Tournament analogue |
|----------------|-------------------|
| **Ruliad** | Space of all tournaments on all n |
| **Rulial space** | The metagraph G_n (one slice of the ruliad at fixed n) |
| **Multiway system** | Tiling hypercube Q_m with waggly connections |
| **Rule application** | Tile flip (arc reversal) |
| **Causal graph** | H-gradient DAG on G_n |
| **Causal invariance** | H never decreases along any path (DAG property) |
| **Branchial space** | Complement dimension (T ↔ T^op) |
| **Branchial foliation** | G_n/Z_2 (merged metagraph) |
| **Observer** | The quotient map Q_m → G_n → G_n/Z_2 |
| **Pockets of reducibility** | Rédei, OCF, score sequences, band-limitedness |
| **Computational irreducibility** | Full isomorphism class requires enumeration |
| **Emes** | Individual arcs (no intrinsic structure, only relational) |

## Pockets of Reducibility in Tournament Theory

Wolfram's central insight: within computationally irreducible systems, there always exist **pockets of reducibility** — aspects that CAN be predicted without running the full computation. Conservation laws are the paradigmatic example.

In tournament theory, we have identified these pockets:

### Fully Reducible (predictable from simple invariants):
- **Rédei's theorem**: H(T) is always odd — regardless of T
- **Score sum**: Σ s_i = C(n,2) — always, for any tournament
- **3-cycle formula**: c₃ = C(n,3) - Σ C(s_i, 2) — determined by scores alone
- **Parity of H**: always odd — a universal conservation law
- **Complement symmetry**: H(T) = H(T^op) — pairing theorem

### Partially Reducible (predictable from richer invariants):
- **Score sequence → 97% of H variance** — the score is a powerful but incomplete predictor
- **Band-limitedness**: Walsh degree ≤ 2⌊(n-1)/2⌋ — high-frequency modes are zero
- **OCF**: H = I(Ω(T), 2) — reducible to the conflict graph
- **Fiber fraction** f(n) ~ 1/√(πn) — the average class size is predictable

### Computationally Irreducible:
- **Full isomorphism class**: requires checking all n! permutations (nauty/traces)
- **Exact H value**: no formula from polynomial invariants alone (for general T)
- **Metagraph distance**: equals minimum Hamming distance (requires search)
- **Chromatic number proof**: χ = n-1 proved computationally but no algebraic proof yet

The **boundary** between reducible and irreducible is precisely where the interesting mathematics lives. The score sequence (reducible) captures 97% but the remaining 3% (the cycle structure, the even-graph component) is where computational irreducibility begins.

## The Metagraph as a Multiway System

The tiling hypercube Q_m is a **multiway system** in Wolfram's sense:
- **States**: the 2^m tilings (binary strings of length m)
- **Rules**: tile flips (bit flips, one per tile position)
- **Multiway graph**: Q_m itself (all possible states connected by all possible rule applications)
- **Evaluation order**: which tile to flip first = which rule to apply

The **waggly layers** (d=1, d=2, ..., d=m) correspond to different "update schedules" — flipping 1 tile at a time, 2 at a time, etc.

**Causal structure**: The H-gradient on G_n is strong but NOT a strict DAG. Level edges (same H, different class) appear at n≥5: 1, 15, 136 for n=5,6,7. H-decreasing edges appear at n≥7 (962/4086 at n=7). The "arrow of time" from transitive to regular is a statistical tendency, not an absolute law. See MISTAKE-035. **CORRECTED** (opus-S1, 2026-04-01): earlier versions of this file claimed zero downhill edges — that claim was trivially true for undirected graphs and masked the real structure.

**Branchial structure**: The complement involution T → T^op creates a "branchial" dimension. Two tournaments related by complement are on "different branches" of the multiway system. The merged metagraph G_n/Z_2 collapses this branching, analogous to how an observer collapses multiway branches into a single perceived reality.

## Observer Theory and Tournament Quotients

In Wolfram's observer theory, an observer is an agent that **equivalences states** — reduces the raw complexity of the ruliad by identifying states that look the same from the observer's perspective.

In tournament theory, we have a hierarchy of observers:

1. **Labeled observer**: sees all 2^m tilings individually (no equivalencing)
2. **Isomorphism observer**: equivalences by S_n → sees G_n (A000568 classes)
3. **Complement observer**: further equivalences T ↔ T^op → sees G_n/Z_2
4. **Score observer**: only sees the score sequence → sees ~97% of structure
5. **Parity observer**: only sees H mod 2 → sees 1 bit (always 1, by Rédei)

Each observer corresponds to a different level of **computational coarsening**. The more you coarsen, the more structure you lose, but the more predictable the remaining structure becomes. This is exactly Wolfram's observation that observers create apparent simplicity from underlying complexity.

## Metamathematics of Tournament Proofs

Wolfram's metamathematics treats proofs as paths in a multiway graph of logical deductions. In our project, multiple independent proof strategies converge on the same results:

- **Rédei's parity**: proved by induction, by GF(2) cancellation, by OCF evaluation
- **Claim A (c₀ = H)**: proved by IE, by path homology, by deletion-contraction
- **Band-limitedness**: proved by polynomial degree bound, by Krawtchouk analysis

These convergent proofs are **geodesics in metamathematical space** — different paths through proof space arriving at the same destination. The fact that they exist suggests the results are "robust" features of the mathematical landscape, analogous to how conservation laws are robust features of physical spacetime.

## The Tournament Ruliad

The "tournament ruliad" is the space of ALL tournaments on ALL vertex sets:
- At each n: 2^{C(n,2)} labeled tournaments, A000568(n) isomorphism classes
- The generating functions V(s), L(x) are ways of "sampling" this ruliad
- The Fraïssé limit T_ω is the "generic point" of the tournament ruliad
- Aut(T_ω) is the "symmetry group" of the tournament ruliad

The four constants {√2, π, e, γ} that emerge from the staircase δ_{n-2} are **structural properties of the tournament ruliad** — they characterize the geometry of the space of all tournaments, just as the speed of light and Planck's constant characterize the geometry of physical spacetime.

## Key Open Problems Where Tournaments Could Help

1. **Pockets of reducibility**: Our project has systematically catalogued which tournament properties are reducible (score → H variance) and which are not (full isomorphism). This methodology could be applied to other multiway systems.

2. **Path homology for directed graphs**: GLMY path homology (active area, arXiv 2602.04140 on circulant digraphs) computes Betti numbers for tournaments. Our β₂ = 0 theorem and the seesaw conjecture β_{2k-1} × β_{2k+1} = 0 are constraints on the topology of tournament space that may have analogues in Wolfram's hypergraph rewriting.

3. **Causal invariance detection**: Our computational proof that G_n is a DAG under H (a form of causal invariance) could inform methods for detecting causal invariance in more general rewriting systems.

4. **The compression connection**: Our TC codec exploits pockets of reducibility for compression — knowing the score sequence (reducible) lets you predict most of the structure. This is a concrete engineering application of the reducibility/irreducibility dichotomy.

5. **Branchial geometry from complement symmetry**: The Z_2 complement structure of tournaments provides a simple, exactly solvable example of branchial geometry that could test conjectures about branchial space in general.
