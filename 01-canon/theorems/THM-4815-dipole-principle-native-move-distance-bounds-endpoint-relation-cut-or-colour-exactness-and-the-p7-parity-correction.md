---
id: THM-4815
title: "The dipole principle: every native-move distance lower bound in the repository is a converse-invariant functional read through K -> -K + 2A_Q (potentials, the forward-arc parity, the two-move XOR law, endpoint confinement); for score and cut functionals the change depends only on the path's two endpoints, so the bound is attained from T iff a chain of states exists in which some Hamiltonian path joins a descending pair; it fails either because a cut confines the endpoints (non-strong: endpoints = top x bottom component; strong L-local minima through order 10 all have a cut vertex, though 2-connected failures above them occur from order 7) or because every descending path is blue, which forces an anti-automorphism; exact local-minimum classifications for delta0 and for L on non-strong tournaments; every tournament has a Hamiltonian path from a maximum-score to a minimum-score vertex (classical: a corollary of Thomassen 1982; the case of a minimum vertex beating a maximum one is Lu 2018), so at odd orders the distance to the regular class is exactly R = (1/2) sum |s_v - (n-1)/2|; P_7 is exactly L = 6 moves from the transitive class (its 7 is the labelled parity)"
status: >
  PROVED: (A) the identity, its converse covariance and the instance ledger; (B1)-(B6) the endpoint principle, the
  exactness criterion, cut confinement, the delta0 and non-strong L classifications, the colour-forces-symmetry lemma,
  (B7) the max-to-min Hamiltonian path theorem (self-contained re-proof of a classical result: Thomassen 1982, case (i)
  = Lu 2018) and (B8) the regular-class distance at odd orders; (C) the forward-arc parity law and its consequences.
  FINITE-EXACT: (C) the P_7 word and labelled distances (labelled BFS over 2^21 states); (D) the class-graph census at
  orders 4-9 (order 9 in exact_n9_summary.out, completed after filing), the endpoint-relation census at orders 4-10
  (all 9,733,056 classes), random tests at orders 11-16; (B8') every regular class of orders 5-11 is exactly
  (n^2-1)/8 moves from the transitive class (regular_landing.c).
  CITED: Camion 1959; Lu, Graphs Combin. 34 (2018) 613-618 (min-max pairs: case (i) of B7 is its first theorem; abstract only); Thomassen, Proc. London Math. Soc. (3) 45 (1982) 151-168 (Theorem 5.5 of Bang-Jensen and Gutin, "Paths, trees and cycles in tournaments") (implies B7); THM-4810, THM-4811, Session74 (source/sink
  pruning lemma), THM-4729 L3, THM-4732 Lemma 1, THM-4804. OPEN: HYP-9423 (1)-(4), (7). Independently audited: PASS WITH CORRECTIONS (all proofs and census numbers reproduced
  with the auditor's own code; B7 identified as classical (Thomassen 1982); the title's failure clause and several
  slips corrected; 35 edits applied; 04-computation/experiments/dipole_exactness_20261027/audit_thm4815/REPORT.md).
session: mac-mini-2026-10-27-dipole
source: 05-knowledge/results/dipole_principle_and_prime_cycle_index_tower_20261027.md (sections 1-5)
scripts:
  - 04-computation/experiments/dipole_exactness_20261027/dipole_exactness.cpp (class-graph census; exact_n{4..8}.out, exact_n*_details.out)
  - 04-computation/experiments/dipole_exactness_20261027/endpoint_checks.c (endpoint relation; endpoint_checks_n4_9.out)
  - 04-computation/experiments/dipole_exactness_20261027/labelled_bfs.c (P_7 labelled distances; labelled_bfs_p7.out)
  - 04-computation/experiments/dipole_exactness_20261027/analyze_graph.py, median_checks.c, lu_exceptions.c, endpoint_random.c
related:
  - 01-canon/theorems/THM-4811-arc-strong-connectivity-is-the-exact-complement-distance-to-the-nonstrong-locus-black-refinement-and-the-dipole-law.md
  - 01-canon/theorems/THM-4810-h-spectrum-atlas-strong-fibres-exact-low-spectra-complement-parity-invariant-score-monotonicity-paley-budgets.md
  - 01-canon/theorems/THM-4814-order-ten-prime-cycle-index-types-order-eleven-ledger-skew-hadamard-nuclear-norm-extremality-and-black-core-distance-of-regular-tournaments.md
  - 01-canon/theorems/THM-4729-closed-and-open-tilings-hamiltonian-cycle-bipartition-of-the-metagraph-lines-and-the-forced-arrangement.md
  - 01-canon/theorems/THM-4730-strong-pure-blue-alphabet-complete-points-of-symmetry-are-k1-c3-t5.md
  - 01-canon/theorems/THM-4732-black-detours-of-closed-blue-lines-order-9-census-order-10-detours-pure-blue-layers-and-the-klein-dictionary.md
  - 01-canon/theorems/THM-4804-complement-lines-shatter-every-module-and-the-carousel-fibre-as-frobenius-orbits-of-a-normal-basis.md
  - 07-reflections/why-the-advances-work-ii-potentials-symmetry-budgets-positive-coordinates-macmini-20261027.md
  - 05-knowledge/hypotheses/HYP-9423-dipole-endpoint-conjectures-two-connected-descent-and-the-prime-tower-predictions.md
  - 05-knowledge/hypotheses/HYP-9420-strong-untouchables-1-7-21-complement-graph-components-paley-depth-quadratic-bound.md
---

# THM-4815 — The dipole principle and the cut-or-colour dichotomy

The reflection *why the advances work, II* (section 0) proposed two theses: every distance lower bound in the
repository comes from one identity (a native move negates the skew matrix and adds the path), and where such a bound
is not attained a forced symmetry explains why. This theorem makes the first precise and proves it in the form below
(A). It refutes the second as stated and proves a corrected version (B): failures come from a cut or from a symmetry,
and the endpoint relation of Hamiltonian paths decides between them. It also corrects the reflection's main example
(C).

Concurrent work. klein's THM-4818 (pushed the same morning) gives a census of the same potentials on every class of
orders 4-8, with counts that agree with (D) below. It also gives the condensation criterion for the singleton target,
which B4 restates in local form at every order. It shows that no linear converse-invariant potential is exact for the
transitive target (an LP dual to A2). It reads `L` as the `A_(n-1)` root-word metric, and finds the `P_7` class
distance 6 with the labelled count 2,520 / 2,520 (MISTAKE-682, recorded there). Not in THM-4818: the survey and
ledger of A, the endpoint criterion and the exact form of cut confinement (B1-B3; THM-4818 (C2) already uses the endpoint form of the dipole law and that every word starts in the initial and ends in the terminal component), the exact non-strong `L` classification (B5), the
colour lemma (B6), the max-to-min theorem (classical, see B7) and the regular-class distance (B7, B8), the explicit alternating half of
the parity law (C; implicit in THM-4810 Theorem P's odd toggle count and used in THM-4732 Lemma 1), and the endpoint census to order 10.

Notation. `T` is a tournament on `V`, `|V| = n`, `K` its skew matrix (`K_uv = 1` iff `u -> v`), `s` its scores. A
native move along a Hamiltonian path `Q = (q_1, ..., q_n)` of `T` from `a = q_1` to `b = q_n` is
`tau_Q(T) = T^op xor E(Q)` (THM-4810, Session74). It is black unless the positional reversal `q_i -> q_(n+1-i)` is an
anti-automorphism of `T` (blue).

## A. The identity and the bounds it produces (PROVED)

**A1.** `K(tau_Q T) = -K + 2A_Q`, `A_Q` the skew matrix of the path. Hence `tau_Q` is an involution;
`tau_(sigma Q)(sigma T) = sigma tau_Q(T)` for relabellings; `(tau_Q T)^op = tau_(Q^rev)(T^op)`; and
`tau_Q(T) = rho_(Q^rev)(T^op)`, where `rho_P` reverses the arcs of the Hamiltonian path `P` (of `T^op`). A move is a
converse followed by one path reversal.

**A2 (potential lemma).** If `Phi(T^op) = Phi(T)` and reversing one Hamiltonian path changes `Phi` by at most `c`, then
every native move changes `Phi` by at most `c`, and `d(T, Z) >= (Phi(T) - max_Z Phi)/c` for every target `Z`.

**A3 (the instances).** Reading A1 at each level of structure gives every native-move lower bound found in the
repository (survey of all THM / HYP / result packets, 2026-10-27):

| level | functional | change per move, from A1 | repository use |
|---|---|---|---|
| row sums: `d -> -d + e_a - e_b` | `L = (1/2) sum |s_(i) - i|` | at most 1 | THM-4810 Theorem L (transitive class); Session74 forest frames `beta >= L - n + 2` |
| | `delta0 = min_v min(s_v, n-1-s_v)` | at most 1 | THM-4811 (D), THM-4814 (D) for regular sources (singleton loci) |
| | `R = (1/2) sum_v |s_v - (n-1)/2|` (odd n) | at most 1 | regular class (this theorem, B8) |
| | `c3 = C(n,3) - sum C(s_v, 2)` | exactly `s(a) - s(b) - 1` | black triangle descent (black_descent_20261020c; THM-4811 (A)) |
| cut charges: `f(X) -> f(V\X) + [a in X] - [b in X]` | `lambda` | at most 1 | THM-4811 (B) (non-strong locus) |
| Hamming to a converse-closed bank | feedback number `fas`, exception counts | at most `n - 1` | THM-4810 (C), Session73/75 banks, THM-4814 (C) |
| unitarily invariant norms | `||K||_op`, `||K||_*` | `2 ||A_Q||` = `4 cos(pi/(n+1))`, `2 sum |2 cos(k pi/(n+1))|` | THM-4813 (C), Session75 section 6 (THM-4814 (C) uses the static nuclear-norm form) |
| upper triangle mod 2 (labelled) | forward arcs `F` | `F' = C(n,2) + n - 1 - F - 2 back(Q)` | THM-4810 Theorem P; (C) below |
| two moves | `K'' = K - 2A_P + 2A_Q` | toggles `E(P) xor E(Q)` | THM-4810 Lemma 1 / C2, THM-4732 Lemma 1 |
| endpoints | the endpoint relation of B3 | — | THM-4729 L3 (non-strong nodes pairwise non-adjacent), THM-4804 (modules shatter), Session74 source/sink pruning lemma |

Two remarks complete the ledger. Module shattering (THM-4804) is A1 at a module boundary: a path enters and leaves a
module `M`, so some outside vertex `v` is path-adjacent to one `m in M` and not to another `m'`, and `K'` keeps the
sign on `(v, m)` and flips it on `(v, m')`. The remaining distance bounds in the repository are not native-move bounds:
they concern other move sets (arc flips, THM-4748; circulant flips, THM-4751; triangle reversals) or are static counts
(triangle packing, trace duality, ranks) that feed the functionals above. The two colour statements (the pure-blue
alphabet THM-4730, the family E of THM-4811) are the other half of the dichotomy below. So the first thesis holds for
every native-move bound the survey found: each is A1 read through a converse-invariant functional, or A1 composed with
itself, reduced mod 2, or restricted to endpoints.

## B. Exactness: endpoints, cuts and colours (PROVED)

**B1 (endpoint principle).** If `Phi` depends only on the scores (`L`, `delta0`, `R`, `c3`) or only on the cut function
(`lambda`), then `Phi(tau_Q T) = Phi_T(a, b)` depends only on `T` and the endpoints of `Q`: by A1 the new scores are
`n - 1 - s + e_a - e_b` and the new cut function is `f(V\X) + [a in X] - [b in X]`.

**B2 (criterion).** For such an integer potential with unit tariff and a target that is exactly its zero set, the bound is
attained from `T` iff there is a chain of moves along which `Phi` drops by one each time. A state has a descending
move iff `Pi(T) ∩ D_Phi(T) ≠ ∅`, where `Pi(T) = {(a, b): some Hamiltonian path runs from a to b}` and `D_Phi(T)` is the
set of pairs with `Phi_T(a, b) = Phi(T) - 1`. It has a descending black move iff the same holds with `Pi` restricted to
black paths. A *local minimum* is a state off the target with no descending move; every local minimum is a failure of
exactness, and every failure lies above one.

**B3 (cut confinement).** If `T` is not strong, with strong components `C_1 => ... => C_k` (`k >= 2`), then
`Pi(T) = C_1 x C_k`. (Paths never go upward, so they start in `C_1` and end in `C_k`; conversely concatenate a
Hamiltonian path of `C_1` from `a`, any paths of the middle components and one of `C_k` ending at `b`; strong
components of size at least 3 are Hamiltonian, Camion.) Moreover all maximum-score vertices of `T` lie in `C_1`, all
minimum-score vertices in `C_k`, and every vertex of `C_1` scores more than every vertex of `C_k`.

**B4 (the delta0 local minima).** A tournament with `delta0 = k >= 1` has a `delta0`-descending move iff it is strong.
(If strong, start a path at a vertex of score `k` or end one at a vertex of score `n - 1 - k`. If not, by B3 the start
`a in C_1` and the end `b in C_k` cannot reach the value `k - 1`: that needs `s(a) = k` or `s(b) = n - 1 - k`, and B3
puts the scores of `C_1` above the minimum and those of `C_k` below the maximum.) So the local minima for the singleton
target are exactly the non-strong tournaments without a source or sink.

**B4' (how far the non-strong local minima are).** If `T` is not strong and `delta0(T) = k >= 1`, then every move
from `T` lands on a strong tournament (THM-4729 L3; it also follows from the cut law of A1 and B3). Choosing the endpoints `a in C_1`, `b in C_k`
away from one vertex attaining `delta0` keeps that vertex interior, so its value stays `k` and the neighbour `T_1` is
strong with `delta0(T_1) = k`. Hence `d(T, singleton locus) = k + 1` as soon as `T_1` attains its own bound. For
`k = 1` this is automatic (B4: a strong tournament with `delta0 = 1` is one move from the locus), so **every non-strong
tournament with `delta0 = 1` is exactly 2 moves from the singleton locus**, at every order. For `k >= 2` it follows
from the strong case (HYP-9421 (6); THM-4811 (D2) through order 10). The first `k = 2` case, `R5 => R5` at order 10,
has the explicit 3-move word of `r5r5_singleton.py`. This proves klein's HYP-9424 (1) for `delta0 = 1` at every order
and reduces the rest to HYP-9421 (6).

**B5 (the L local minima off the strong locus).** A non-strong, non-transitive tournament has an `L`-descending move iff
it has no source and no sink. (The sorted scores of `T` are the components' internal score sequences shifted to
disjoint blocks, so `L = sum_j L(C_j)` and the optimal matchings are blockwise. A move lowers `L` iff some optimal
matching has `s(a)` above and `s(b)` below its target. A source is never above its target `n - 1`, and a sink never
below 0. Without them, the minimum-internal-score vertex of `C_1` is above, the maximum-internal-score vertex of `C_k` is
below, and B3 joins them.) The "only if" half is Session74's source/sink pruning lemma (supplement, section 4). For strong `T` the pair (a vertex of
minimum score, a vertex of maximum score) always lies in `D_L`, so a Hamiltonian path from a minimum to a maximum
vertex gives a descent.

**B6 (colour forces symmetry).** If a state has a descending move but no descending black move (for any potential), then
every descending path `Q` is blue. Its positional reversal `sigma` is then an anti-automorphism with `sigma(a) = b`, so
`T ≅ T^op` and `s(a) + s(b) = n - 1`. This is the mechanism behind THM-4811 Theorem A (R5 and the family E) and
THM-4730 (strong pure-blue atoms `C3`, `T5`: no black move at all).

**B7 (max-to-min Hamiltonian paths).** *Every tournament on `n >= 2` vertices has a Hamiltonian path from a
vertex of maximum score to a vertex of minimum score.* Proof. Write `Delta`, `delta` for the extreme scores, `X`, `Y`
for the sets attaining them; a regular tournament is trivial. For `x in X`, `y in Y` put `T' = T - {x, y}`,
`A = N+(x) ∩ T'`, `B = N-(y) ∩ T'`; a Hamiltonian path of `T'` from `A` to `B`, with `x` prepended and `y` appended,
is the required path. If some vertex of `Y` beats some vertex of `X`, choose such a pair (case (i)); otherwise every
vertex of `X` beats every vertex of `Y` (case (ii)).

(i) *If `y -> x`.* Here `|A| = Delta` and `|B| = n - 1 - delta`. If `T'` is strong (or a single vertex), fix a
Hamiltonian cycle of `T'`. If no `u in A` had its cycle-predecessor in `B`, the predecessor map would inject `A` into
`T' \ B`, of size `delta - 1 < Delta`. So some `u in A` has its predecessor in `B`; run the cycle from `u`. If `T'` is
not strong, with components `F_1 => ... => F_q`, suppose no vertex of `F_1` is in `A`. Then every `f in F_1` beats `x`
and all of `T' \ F_1`, so `s(f) >= |T' \ F_1| + 1 > |A| = Delta`, a contradiction. Dually, if no vertex of `F_q` beats
`y`, then `s(y) >= |F_q| + 1 > s(g)` for `g in F_q`, a contradiction. So `A ∩ F_1` and `B ∩ F_q` are nonempty, and B3
joins them.

(ii) *If every vertex of `X` beats every vertex of `Y`.* Take any `x in X`, `y in Y`, so `x -> y`, `|A| = Delta - 1`
and `|B| = n - 2 - delta`. If `T'` is strong, the same count gives a path unless `Delta - 1 <= delta`, that is, unless
all scores lie in `{delta, delta + 1}`. In that case `X => Y` forces every vertex of `X` to score at least
`|Y| + (internal score)`, and summing gives `delta + 1 = |Y| + (|X| - 1)/2` and `delta = (|Y| - 1)/2`, hence `n = 2`.
If `T'` is not strong and `A ∩ F_1` is empty, the same inequality now allows equality only when `F_1 = {f}` with
`f -> x`, `y -> f` and `s(f) = Delta`. Then `f in X` is beaten by `y in Y`, which is excluded. Dually, `B ∩ F_q`
empty forces `F_q = {g}` with `g -> x` and `s(g) = delta`, so a vertex of `Y` beats a vertex of `X`, also excluded. ∎

So the failures in Lu's sense (an arc `x -> y` with no Hamiltonian path from `x` to `y`: 526 arcs in 514 classes at
order 9, never `y -> x`) are always rescued by another extreme pair, as case (ii) shows. Prior art: B7 is classical. Thomassen (Proc. London Math. Soc. (3) 45 (1982) 151-168; Theorem 5.5 of Bang-Jensen and
Gutin's survey *Paths, trees and cycles in tournaments*) characterised when a specified vertex `x` of maximum and `y` of
minimum out-degree are joined by a Hamiltonian path from `x` to `y`: always, except (a) when the terminal component of
`T - x` is one vertex `z` with `z -> x`, preceded by the component `{y}`, and `x -> y`; (b) the mirror configuration in
`T - y`; (c) for `R(1,2)` and two small exceptional tournaments. In (a) the vertex `z` also has minimum score and beats
`x`, in (b) dually, so the existential statement follows ((c) by direct check). Lu (Graphs Combin. 34 (2018) 613-618; published
abstract, full text not consulted) proves case (i) for every pair (a minimum vertex beating a maximum vertex always
gives a Hamiltonian path from the latter to the former) and the case `x -> y` up to a few exceptions. The proof above is
a self-contained re-proof; the endpoint relation agrees with the characterisation for every (max, min) pair at orders
3-9 (dipole_exactness_20261027/audit_thm4815/audit_thomassen_n3_9.out).

**B8 (the regular class at odd orders).** *For every tournament `T` of odd order `n`,
`d(T, regular class) = R(T) = (1/2) sum_v |s_v - (n-1)/2|`*, with moves of either colour. `R` is an integer potential
of unit tariff that vanishes exactly on regular tournaments (A3). A move lowers `R` iff its path starts above and ends
below the mean score, and in a non-regular tournament the maximum score is above the mean and the minimum below it. So
B7 gives a descending move at every non-regular state, and B2 gives the distance. An explicit geodesic: repeatedly move
along a Hamiltonian path from a maximum-score to a minimum-score vertex. For black moves the bound can fail (a colour
obstruction at order 7: the strong class with scores `1, 3^5, 5` and `|Aut| = 5`). At order 9 there is none: the black
distance to the regular class is `R` for every class except six black-isolated pure-blue ones.

**B8' (regular tournaments and the transitive class).** Applied to `TT_n`, B8 gives `d(TT_n, regular class) =
R(TT_n) = (n^2 - 1)/8`. Every regular tournament has `L = (n^2 - 1)/8`, and moves are involutions. So each regular class
that some `R`-descending chain from `TT_n` ends on lies at distance exactly `(n^2 - 1)/8` from the transitive class,
attaining THM-4810's lower bound. Random descending chains from `TT_n` (regular_landing.c, canonical forms) end on
every regular class at orders 5, 7 and 9 (3 of 3 at order 7, 15 of 15 at order 9) and on all 1223 of the 1,223 regular classes at order 11
(including `P_11` in that sample). So `L` is attained at every regular tournament of order 7 and 9. The
order-9 class census confirms it (all 15 regular classes at distance 10). HYP-9420 (4)'s excess is therefore 0 for
every regular tournament at these orders, not only for `P_7` and `P_11`; HYP-9423 (7) conjectures it at every odd
order.

## C. The parity law and P_7 (PROVED; FINITE-EXACT witnesses)

For a labelled tournament on `{1, ..., n}` let `F` be its number of forward arcs. Summing A1 over the upper triangle
gives `F' = C(n,2) + n - 1 - F - 2 back(Q)`, so `F' + F ≡ (n-1)(n+2)/2 (mod 2)`. For `n ≡ 1, 2 (mod 4)` the parity of
`F` is invariant (THM-4810 Theorem P). For `n ≡ 0, 3 (mod 4)` it alternates (THM-4810's toggle count is odd there; THM-4732 Lemma 1 uses this for detours), so
the labelled move graph is bipartite and
`d_lab(T, U) ≡ F(T) + F(U) (mod 2)`.

**P_7.** THM-4810 (D) found `P_7` at distance 7 from the transitive tournament `TT_id` (`i -> j` iff `i < j`). That is a
labelled distance to one fixed labelling. The Paley labelling has `F = 14` and `TT_id` has `F = 21`, so at `n = 7` the
distance is odd, at least `L + 1 = 7`. The distance to the transitive *class* is `L(P_7) = 6`. A labelled BFS over all
`2^21` states (labelled_bfs.c) finds a 6-move word, re-verified move by move, ending at a transitive tournament with
`F = 12`. Starting from `P_7 = Z_7{1,2,4}` on `0, ..., 6`, the six paths, applied in order, are

    0123456, 1203645, 0352641, 2316405, 4053162, 2163405

(output `labelled_bfs_p7.out`; the scores pass from `3^7` to `0, 1, ..., 6`).

The reflection's table entry "`P_7` no (7 > 6)", its source note, Session74's supplement ("P7 ... native distance seven
to the transitive class") and HYP-9420 (4) ("the exact excess is 1 at q = 7") read the labelled 7 as a class distance (the reflection and the body of HYP-9420 (4) are now corrected; Session74's
supplement is left as filed).
The excess of `P_7` is 0, as for `P_11` (class distance 15 by Session74's certificate; its labelled distance to `TT_id`
is at least 16 by the same parity). Recorded as MISTAKE-682.

## D. The census (FINITE-EXACT)

Exact class-level distances come from a multi-source BFS on the undirected class graph, any colour and black only, at
orders 4-8 (order 9 in exact_n9_summary.out, completed after filing). The endpoint relation `Pi(T)` comes from a subset DP at orders 4-10, all 9,733,056
classes at order 10:

| target / potential | failures at n = 4, ..., 9 | local minima at n = 4..10 | of which symmetric, `|Aut| > 1` or self-converse (n <= 9) | connectivity of the local minima |
|---|---|---|---|---|
| non-strong / `lambda` | 0, 0, 0, 0, 0, 0 | none | – | – |
| transitive / `L` | 2, 6, 33, 182, 1,431, 17,936 | 2, 6, 22, 110, 904, 13,545, 377,493 | 2, 4, 10, 32, 114, 707 | non-strong with a source or sink (2, 5, 19, 99, 855, 13,303, 376,191), or strong with a cut vertex (0, 1, 3, 11, 49, 242, 1,302) |
| singleton / `delta0` | 0, 0, 1, 3, 16, 99 | 0, 0, 1, 3, 16, 99, 915 | all through order 8; 91 of 99 at order 9 (rigid non-self-converse ones such as `S_4 => S_5'`) | all non-strong (B4), all at distance `delta0 + 1` (B4') |
| regular / `R` (odd n) | 0 (n = 5, 7, 9) | 0 (n = 5, 7, 9) | – | none (B8) |

endpoint_checks.c verifies B3, B4, B5 and B7 class by class at orders 4-10. In the black-only census, every black local
minimum that is not an unrestricted one is a colour obstruction, and these are self-converse as B6 requires (two at
order 6 for the almost-regular target defined below, one at order 7 for the regular target). At orders 4-8 the black class graph is one component plus
isolated pure-blue classes (THM-4730 with THM-4729 L4: palindromic stacks of `{K1, C3}` halves with centre in `{none, K1, C3, T5}`). So a class has no black route to a
target exactly when it is itself pure blue, or when every target class is. The transitive target is the extreme case:
its only Hamiltonian path is blue, so every other class is cut off.

The smallest rigid failures are non-strong: `K1 => S4` and `S4 => K1` at order 5 (`1111100111`, `1101111111`;
`|Aut| = 1`, not self-converse, `L = 1`, distance 2; a source or a sink, B5). The smallest strong rigid ones occur at
order 6 (10 of the 20 rigid `L`-failures there). All three strong `L`-local minima of order 6 have scores `1,1,2,3,4,4`
and `|Aut| = 1`: `111011111111101` (in gentourng order; not self-converse, `L = 1`, distance 3 to the transitive
class), its converse `101011111111110`, and the self-converse `111011111111110`. In the first, vertex 3 has the single
out-neighbour 4 and 5 the single out-neighbour 3, so a Hamiltonian path starting at 3 must end at 5, and the only one
starting at 5 is `5, 3, 4, 0, 1, 2`. As a result no path joins a minimum-score vertex (3 or 5) to a maximum-score vertex
(0 or 1). These are the only pairs in `D_L`.

At even orders the integer bound `ceil(R'/2)`, with `R' = sum_v (|2s_v - n + 1| - 1)/2` (which a move changes by at most
2), is not attained for 306 classes at order 8. Their local minima have connectivity 1 and 2; this case is not analysed
here. Random tests at orders 11-16 found no strong 2-connected tournament without an `L`-descent: uniform,
near-regular and planted-tie samples, 1,140,000 tournaments in all (endpoint_random_n11_16.out), and 595,000
near-transitive ones (transitive plus random arc reversals, where the known local minima live;
endpoint_random_neartransitive_n11_16.out).

## E. Verdict on the reflection's theses

1. *One identity.* Confirmed (A). The identity enters through converse-invariant functionals, through its square (the
   two-move XOR law), through its reduction mod 2 (parity), and through its endpoints (B1-B3).
2. *A forced symmetry explains every failure.* False as stated. Unrestricted failures are endpoint confinements: every
   local minimum of `lambda`, `delta0` and `L` through order 10 has vertex connectivity at most 1, most are not strong,
   and rigid examples exist from order 5 (strong ones from order 6); 2-connected failures that are not local minima occur
   from order 7 (1, 5, 143 classes at orders 7, 8, 9). They inherit a cut from below: the order-7 one (scores
   `2,2,2,3,4,4,4`, rigid, self-converse, `L = 3`, distance 4) has a descending region of 20 classes that never reaches
   the transitive class, and all 12 of its dead ends are cut-confined local minima. Symmetry is forced exactly when the
   obstruction is colour (B6). The corrected
   thesis is **cut or colour**. Either a cut confines the endpoints of Hamiltonian paths, or every good path is blue and
   its reversal is an anti-automorphism.
3. *Where neither can occur, the bound is exact.* Two targets have no unrestricted failure at all. For the non-strong
   locus (`lambda`, THM-4811 (B)) the Camion cycle supplies the descending endpoint pair. For the regular class at odd
   orders (`R`, B8) the descending pairs run from high to low scores, the natural direction of Hamiltonian paths, and
   B7 shows a cut can never block all of them. For the transitive class the descending pairs run from low to high
   scores, against that direction. This is why `L` fails there, and why the open case of HYP-9423 (1) asks for 2-vertex
   connectivity.

Most of the reflection's evidence (R5, the family E, the pure-blue alphabet) lies on the colour side. Its one cut-free
example, `P_7`, was not a failure.
