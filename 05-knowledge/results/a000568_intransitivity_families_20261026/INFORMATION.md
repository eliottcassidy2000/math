# Tournaments as information: merging the earlier work with the intransitivity observable

klein, 2026-10-26. The directive asked to merge the repo's "tournaments as information" work
(March–October 2026) with the enumeration. This note says, file by file, what each earlier
measure is in the frame of INTRANSITIVITY.md and what the new exact numbers add.

## 1. The frame: Hartley fibres after MISTAKE-231

MISTAKE-231 (audit of HYP-8875) fixed the only admissible meaning of "information" without a
probability law: `H₀(X | L) = log₂ |{X' : L(X') = L(X)}|` after fixing the ambient universe, the
equivalence and the observable `L`. Everything below uses `U_n` = tournaments on `[n]` (or the
isomorphism classes when stated), equality, and one of the observables

- `L = id`: `C(n,2)` bits, the relation itself;
- `L = C3` (intransitivity structure): fibre `N(T) = 2^{p} ∏ k_N!` (Theorem I);
- `L = scores`: fibre = number of labelled tournaments with that score function (MISTAKE-231's
  own example; not a tree quantity);
- `L = H` (Hamiltonian paths, `information_theory_deep.out`'s "ranking entropy" `log₂ H`);
- `L = (C3, scores)`: fibre `= 1` except on the balanced-prime classes (1, 2, 5, 14, 53 for
  `n = 3..7`).

The new fact is that for `L = C3` the fibre is a product over the modular tree, so the
intransitivity observable is the first one in the repo whose information loss is **structural**
rather than tabulated: it is `log₂ n!` on the transitive relation, `1` on every prime relation,
and `1.22` bits on average at `n = 10` (tending to 1).

## 2. File-by-file merge

- `07-reflections/what-we-can-and-cannot-know.md` (S307): counts `C(n,2)` bits per labelled
  tournament and `log₂ V_n ≈ C(n,2) - log₂ n!` per class; its unaudited sentence "tournament
  space has half the information content" (from `H`'s band-limit) is contradicted by its own
  estimate, as the survey noted. Exact replacement: a class carries `log₂ A000568(n)` bits and its
  intransitivity structure carries `log₂ c(n)`, which differ by `log₂(A000568/c) = 1.03` bits at
  `n = 10`; the `C3`-observer loses one bit of a class, not half.
- `05-knowledge/results/information_theory_deep.out` (ranking entropy `S = log₂ H`, the
  comparison channel, MDL of `H`): `H` is a `C3`-invariant (INTRANSITIVITY.md §3), so the ranking
  entropy is a function of the intransitivity structure alone; two relations with the same cyclic
  triples have the same number of consistent rankings. The MDL question "how many bits to specify
  `H`" is bounded by `log₂ c(n)` rather than `log₂ A000568(n)`.
- `07-reflections/information-geometry-of-tournaments.md` (S20y: Walsh–Fourier of `H` has only
  even orders because `H(T) = H(T*)`): the even-order statement is the Fourier shadow of
  `C3`-invariance restricted to the global dualization; `C3`-invariance is stronger (it also
  includes the local dualizations and the reorderings).
- THM-1420 (no `F₂`-linear invariants; every class splits into two `A_n`-classes by inversion
  parity, "the chirality that survives"): the one bit the intransitivity structure always drops is
  this chirality. For a prime `T`, `C3(T)` determines the unordered pair `{T, T*}` and nothing
  more; the `A_n`-split and the `C3`-fibre are the same bit seen from the two sides.
- THM-931 / METAGRAPH-ATLAS §8 (the 8 wild invisible pairs at `n = 9` are chirality pairs on which
  every reversal-invariant invariant ties): explained — those pairs are `C3`-classes of size 2, so
  `c3, c4, c5, |Aut|, H`, the spectral decks and every other `C3`-invariant tie by Theorem I.
  Prediction for the open `n = 10` wild hunt: wild pairs can only be `C3`-equivalent classes, and
  among primes those are exactly the 3,950,318 chirality pairs `(T, T*)` (`(7906161 - 5525)/2`).
- THM-4739 (H): first-order information = scores (the gradient part of the arc flow), second-order
  = the cyclic-triangle 2-form (curl). HYP-9351(a) asks whether the full 2-form forces the
  metagraph ports. At the level of the relation itself the answer is exact: the full 2-form as a
  set of triangles determines `T` up to the tree moves, and together with the first-order data
  (score function) up to balanced-prime dualizations. Whether this transfers to the port
  transportation problem is HYP-9351's question and is not settled here.
- HYP-3602 and `relations-not-things.md` ("order ⊕ intransitivity = cut ⊕ cycle"; 3-cycles
  generate all intransitivity, Moon): Theorem I is the quantitative form — the cycle part
  determines the relation up to the tree, and the cut part (scores) fixes the rest except at
  balanced prime nodes.
- THM-1805 (the Vandermonde is a signed sum over tournaments; the cancelling involution reverses a
  3-cycle): each cancellation step moves along one edge of `C3(T)`; intransitive tournaments are
  those with `C3 ≠ ∅`, and the survivors are the `N(T) = n!` fibre of the empty hypergraph.
- THM-559 (`c3` is a frustrated Ising energy): `c3 = |E(C3(T))|` is the crudest `C3`-invariant;
  the hypergraph is the full "frustration pattern" and Theorem I says it is almost a complete
  invariant.
- THM-462 (Kendall–Babington Smith `c3 = C(n,3) - Σ C(s_i, 2)`): the scores determine `c3` but
  not `C3(T)`; `C3(T)` determines `c3` but not the scores. The two observables are incomparable,
  and their join is almost complete (§1, last item).
- `session58_hourglass_complement_geometry_20261025/observer-profile-algebra` (an observer
  averaging over all `n!` orders is not linearly complete from `n >= 6`): the `C3`-observer is a
  different, non-linear observer; it is complete up to the finite group `∏ S_{k_N} × (Z/2)^p`.
- `session56_rmq_lca_connections_20261025/tournament-hierarchy` (species `T = SEQ(S)`,
  `S = X + P ∘ T₊`, LCA navigation of the modular tree): the tree it navigates is the tree whose
  symmetry group is the `C3`-fibre; the cycle-index data it said were missing are
  `data/prime_cycle_types_n9.tsv`.
- THM-1960 (seed census `1,1,1,0,3,15`): extended to `197, 4008, 136655, 7906161`.
- HYP-3803/3804 (flip rank, rainbow number: how many free arcs are needed to reach every class):
  an information floor `⌈log₂ A000568(n)⌉`; the `C3` analogue is `⌈log₂ c(n)⌉`, one less.

## 3. The two-line summary

A tournament is a complete antisymmetric binary relation; its intransitivity is a 3-uniform
hypergraph; the hypergraph determines the relation up to `2^{p} ∏ k_N!` choices, which is one bit
for almost every relation and `log₂ n!` bits for the transitive one. Counting the hypergraphs
(`c(n)`, `ℓ(n)`) and counting the relations (`A000568`, `2^{C(n,2)}`) are now tied by the same
prime cycle-index census through `n = 10`.
