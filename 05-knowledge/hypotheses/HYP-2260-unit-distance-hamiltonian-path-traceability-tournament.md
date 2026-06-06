# HYP-2260 — Unit-distance configs as tournaments: is the Hamiltonian path made of unit distances?

**Session:** claudebox-2026-06-03-S627. **Prompt (user):** map a unit-distance config to a tournament (base
transitive order, FLIP arcs of unit-distance pairs); is the mandatory (Rédei) Ham path part of the unit distances?
does it "flop" out of the unit pairs at some n? recursive patterns? **Threads:** HYP-2235 (unit-distance↔LRC/CM),
HYP-2245 (H=#Ham paths), HYP-2250 (deletion-contraction/Jacobsthal), HYP-2255 (SCC equidecomposition).

## The clean question (order-independent)
The mapping's tournament `H` is **order-dependent** (angular order gives `H = 3,5,11,37,143,285,…` for n=3..8; a
different point order gives different `H`), so `H` is NOT a configuration invariant. The order-INDEPENDENT object is
the **unit-distance graph (UDG)** and its **traceability**: "the mandatory Ham path is part of the unit distances"
⟺ the UDG has a Hamiltonian path. (Refinement of the user's mapping.)

## The answer: no flop for the strong configs (verified)
**Verified (`unit_distance_graph_traceability_tournament_s627.py`):** the UDG is **traceable (has a unit-distance
Hamiltonian path) for every triangular-lattice compact patch `n = 3..20`** and for the **Moser spindle** (n=7).
Unit-edge counts track the Brass–Moser–Pach law `≈ 3n − √(12n)`. So for the dominant (lattice/Harborth) optimal
family there is **no flop**: the Hamiltonian path stays inside the unit distances.

### Why — and where a flop COULD happen
- **Max-edges ⟹ traceable (conjecture).** Maximizing unit distances forces high connectivity (a compact patch
  snakes row-by-row through unit edges), so the optimum is the configuration MOST likely to be traceable;
  non-traceability (a cut vertex into ≥3 parts, disconnection) wastes edges and indicates a *sub*optimal sparse
  structure. So the flop is **anti-correlated with optimality** — it does not happen at the lattice optimum.
- **The lattice→tower transition is where to look.** The grid-conjecture disproof (HYP-2235) replaces the dense
  lattice with a SPARSE algebraic (CM class-field-tower) construction that wins only asymptotically. Those UDGs have
  unbounded but slowly-growing degree (`~n^{0.014}`, far below the Dirac `n/2`), so traceability is genuinely open
  there — **the flop, if it exists, is at the (very large, currently unknown) lattice→tower crossover.** Conjecture:
  lattice-type optima are always traceable; the sparse asymptotic record-holders may not be.

## Recursive patterns (the structure of constructions)
- **Hexagonal-shell recursion.** Compact patches grow by completing centered-hexagonal shells (`1,7,19,37,…`); the
  unit-edge increment `r(n)−r(n−1) ∈ {2,3}` (3 when a new vertex completes a hexagonal fan, 2 on a straight edge),
  with the run-structure keyed to the shell boundaries (the `n=7, 19` shell-completions).
- **The 1-D chain = Jacobsthal (formalized).** The degenerate collinear unit-spaced config has UDG = the path `P_n`;
  its independence-polynomial-at-2 (the H-weight) is the **Jacobsthal** sequence obeying the `n−2` deletion-contraction
  recursion `J(n+2)=J(n+1)+2J(n)` (HYP-2250), closed form `3·J(n)=2^{n+2}−(−1)^n`. **`J(4)=21`** — a partition-function
  value that is a *forbidden* tournament H (THM-079): a path conflict graph is not realizable from a tournament's
  3-cycles (HYP-2255). The 2-D triangular configs are the "thickened" version with denser UDG and larger H.

## The tournament lens (creative findings)
- **Equilateral unit triangles** (triples with all 3 distances = 1) are the structured "3-cycle seeds" of the flip
  tournament; their count grows `≈ n` in compact patches (each new interior vertex closes ~1 triangle).
- Whether a triple is a flip-tournament 3-cycle depends on the base order AND the unit pattern (NOT a clean parity) —
  the source of `H`'s order-dependence. A canonical (angular) order makes `H` a quasi-invariant; the forbidden
  `7·3^k` H-values would then forbid certain angular-ordered unit configs — an open direction.

## Formalized (math-lean, sorry-free) — `Math/Tournaments/Jacobsthal.lean`
`J` (the chain/path-UDG H), `three_mul_J` (closed form `3·J(n)=2^{n+2}−(−1)^n`), `J_four` (`J(4)=21`, the
forbidden value).

## Open
- Prove the max-edges ⟹ traceable conjecture (or find a non-traceable maximum unit-distance graph).
- Traceability of the CM-tower construction's UDG (the asymptotic record) — the flop's true location.
- Canonical-order H as a config quasi-invariant; do forbidden `7·3^k` H-values constrain unit configs?
