---
id: THM-4739
title: "Kirchhoff and Hodge laws of the tournament metagraph: intransitivity currents, second-order determination, and the consecutive-start-port ledgers. (H) The +-1 arc flow of a tournament on K_n has the Hodge decomposition gradient + curl with potential p_i = (2 s_i - (n-1))/n (the scores), gradient energy (n^2-1)/3 - 8 c3/n and curl energy (n-1)(n-2)/6 + 8 c3/n: the first-order (transitive) information is the score vector, the second-order (intransitive) information is the cyclic-triangle 2-form, whose energy is affine in c3 (the axis of THM-785). (K) Kirchhoff's law of the merged metagraph: for every class invariant with an exact line flux (c3: s(start) - s(end) - 1 per THM-785; the score variance V = n(n^2-1)/12 - 2 c3) the metagraph Laplacian applied to the potential equals the class-intrinsic current J(u) = sum over the Hamiltonian paths of the class, modulo Aut, of the flux; the total current vanishes; verified n = 4..7. (D) Determination: with only first-order port data (each closed tiling's endpoint scores and the partner's score multiset, THM-646), the class-level transportation problem has 1 solution at n = 4 but 4 at n = 5, of which one is the real metagraph; with second-order port data (the partner's curl profile, cyclic triangles per vertex) every port class is forced at n = 5 and the metagraph is determined (results at n = 6, 7 recorded in the note). (L) Consecutive-start ports of the closed endpoint-k layer, I = {2, ..., k}: the arc between b and y splits the pattern exactly in half; the half with y -> b is repaired by the inherited closed-no-chord ledger and the half with b -> y by the port-block insertion ledger Q1 = (a,b,x,U,z,J,c), Q2 = (a,x,U,z,J,b,c), Q3 = (a,x,U,J,b,c,z) (U = positions 3..n-k-2, J = n-1-k..n-3), whose directedness, path-edge cancellation and three anti-reversal witnesses are proved for all n >= k+5; exhaustive at k = 4, n = 9..12 and k = 5, n = 11, 12"
status: >
  PROVED (pre-audit): (H) elementary least squares on K_n; (K) from THM-785 (4) and THM-646 (an identity once the flux law
  is granted; its content is the list of intrinsic currents); (L) the epsilon = 1 ledger for all n >= k + 5, k >= 3 (proof
  below); the epsilon = 0 ledger is the packet's closed-no-chord ledger, proved there for k = 3 at n >= 8 (first word) and
  n >= 11 (others) and here exhaustive for k = 4, 5.
  FINITE-EXACT: electric_census.out (Hodge energies, flux law, Laplacian identities, resistance distances, n = 4..7);
  determination.out / determination_L2_n5.out / determination_L2_n6.out (CP-SAT counts of consistent metagraphs at levels
  L1-L3: n = 4: 1; n = 5: 4 / 1 / 1; n = 6: >= 169,359 (capped) / >= 200,000 (capped) / 48,000 exact); consecutive_start_ledgers_small.out / _large.out (every closed standard-blue presentation of the pattern at
  n = 9..12 for k = 4 and n = 11, 12 for k = 5: both halves repaired, no failure); invariant_collisions.py (class-level
  collisions of the invariants at n = 5, 6, 7).
  CITED: THM-785 (axis and flux laws), THM-790, THM-646 (score complement), THM-643/796 (parities), THM-4729 (closed/open
  tilings), THM-4734 (the (n-5,4) census that located the pattern), endpoint_three_repairs_20261020h (K convention, the
  closed-no-chord ledger and section 6 ledger that the port-block insertion generalises), Jiang-Lim-Yao-Ye 2011 (HodgeRank).
  NOT PROVED: that second-order port data determine the metagraph for all n (HYP-9351); the epsilon = 0 ledger for k >= 4
  at all n (the packet's permutation argument is the route); any other consecutive-port pattern.
source: mac-mini-2026-10-21-electric (owner's prompt: classify the tournament metagraph thinking in terms of intransitivity of information or electrical networks; pursue next obligations)
depends_on:
  - 01-canon/theorems/THM-785-cyclic-triangle-flow-on-blue-black-metagraph.md
  - 01-canon/theorems/THM-646-line-metagraph-score-complement-skeleton.md
  - 01-canon/theorems/THM-4734-endpoint-three-layer-closed-at-order-ten-and-the-endpoint-four-layer-at-orders-nine-and-ten.md
related:
  - 01-canon/theorems/THM-4729-closed-and-open-tilings-hamiltonian-cycle-bipartition-of-the-metagraph-lines-and-the-forced-arrangement.md
  - 05-knowledge/results/endpoint_three_repairs_20261020h/theory.md
  - 05-knowledge/hypotheses/HYP-9345-merged-metagraph-open-laws-after-thm-4729-4730.md
  - 05-knowledge/hypotheses/HYP-9351-second-order-determination-of-the-metagraph-and-the-remaining-consecutive-port-patterns.md
scripts:
  - 04-computation/experiments/electric_metagraph_20261021/electric_census.py
  - 04-computation/experiments/electric_metagraph_20261021/determination.py
  - 04-computation/experiments/electric_metagraph_20261021/determination_L2.py
  - 04-computation/experiments/electric_metagraph_20261021/invariant_collisions.py
  - 04-computation/experiments/endpoint_layers_20261020/residual_cover.py
  - 04-computation/experiments/endpoint_layers_20261020/consecutive_start_ledgers.py
outputs:
  - 04-computation/experiments/electric_metagraph_20261021/electric_census.out
  - 04-computation/experiments/electric_metagraph_20261021/determination.out
  - 04-computation/experiments/electric_metagraph_20261021/determination_L2_n5.out
  - 04-computation/experiments/electric_metagraph_20261021/determination_L2_n6.out
  - 04-computation/experiments/endpoint_layers_20261020/residual_cover.out
  - 04-computation/experiments/endpoint_layers_20261020/consecutive_start_ledgers_small.out
  - 04-computation/experiments/endpoint_layers_20261020/consecutive_start_ledgers_large.out
---

# THM-4739 -- the metagraph as a network: potentials, currents, and what information determines it

**PROVED (pre-audit) + FINITE-EXACT** as itemized in the status.

## 0. Setting

Tiling model (THM-4729 section 0): base path `n -> ... -> 1`, tilings `t in {0,1}^m`, complement
line `t -- t^c`, mirror `sigma`, merged nodes (classes with their converse), closed tilings (tile
`(n,1)` carrying `1 -> n`). Scores `s_i`, cyclic-triangle count `c3`, score variance
`V = sum (s_i - mu)^2 = n(n^2-1)/12 - 2 c3` (THM-785 (1)). THM-785 (4): along the line from a
tiling `t` with path start `P_1` (vertex `n`) and end `P_n` (vertex `1`),
`c3(t^c) - c3(t) = s(P_1) - s(P_n) - 1`. THM-646: the partner's score multiset is
`{n-1-s(v) : v interior} u {n - s(P_1), n-2-s(P_n)}`.

## 1. (H) Hodge decomposition: first- and second-order information

Write the tournament as the flow `f_{ij} = +1` if `i -> j`, `-1` otherwise, on the edges of
`K_n`. Least squares on `K_n` gives the orthogonal decomposition `f = grad p + curl part`
(there is no harmonic part on a complete graph), with

```text
p_i = (2 s_i - (n-1))/n,          E_grad = n sum_i p_i^2 = (n^2 - 1)/3 - 8 c3/n,
E_curl = C(n,2) - E_grad = (n-1)(n-2)/6 + 8 c3/n.
```

(Normal equations `n p_i - sum_j p_j = 2 s_i - (n-1)`; the energy identities follow from
`sum C(s_i,2) = C(n,3) - c3`.) So the **first-order information** of a tournament, the part
expressible as a potential, is its score vector, and the **second-order (intransitive)
information** is the curl 2-form -- the function on triangles with value `+-3` on cyclic and
`+-1` on transitive triangles -- whose energy is affine in `c3`. THM-785's "axis" is the
curl energy up to this affine change, and the transitive tournament has the minimal curl
energy `(n-1)(n-2)/6 > 0` (the `+-1` flow of a total order is not itself a gradient).
Converse negates `f` and preserves both energies; the complement line sends `f` to
`-f + 2 f_P` (`f_P` the path flow), whose gradient part is the endpoint correction of THM-646
(`p -> -p + 2 q`, `q = +-1/n` at the path ends, `0` inside).

## 2. (K) Kirchhoff's law of the metagraph

Let `M` be the merged metagraph with line multiplicities as conductances (loops ignored) and
`L_M` its Laplacian. For a class invariant `Phi` with an exact line flux `F(t) = Phi(t^c) - Phi(t)`
depending only on the tiling, summing over the lines at a node `u` gives

```text
(L_M Phi)(u) = - J_Phi(u),      J_Phi(u) := sum_{tilings t at u} F(t),
```

a **class-intrinsic current**: for `c3`, `J(u) = sum_{Hamiltonian paths P of C} (s(P_1) - s(P_n) - 1) / |Aut C|`
(twice that for a merged node), and `J_V = -2 J_c3`. Total current `sum_u J(u) = 0`. The identities
were checked exactly against the actual multigraph for `n = 4..7` (electric_census.out: zero
discrepancy at every node). The transitive node has `J = n - 2` (its single line raises `c3`
from `0` to `n-2`); `J` is a new node invariant computable without the metagraph, equal to
`[sum_v (s_v - mu)(h_v - h'_v) - H]/|Aut|` with `h_v, h'_v` the numbers of Hamiltonian paths
starting and ending at `v`. Resistance distances in `M` (unit conductance per line) are
recorded for `n <= 7` (Kirchhoff indices `3, 34.8, 157.2, 2610.6`); at `n = 7` the resistance
from the transitive node has Spearman correlation `-0.96` with `log t_u` (nodes with large
fibres are electrically close), `-0.64` with `c3`.

## 3. (D) Determination: how much information fixes the metagraph

Node-level **port data**: for each class `C`, the multiset of its closed and open tilings typed
by `(s(P_1), s(P_n), grid-symmetric?)`; a closed port of type `(a, b, sym)` must be joined to an
open port of type `(n-a, n-2-b, sym)` of a class whose invariants equal the partner invariants
computable from the tiling. With first-order partner data (the score multiset of THM-646), the
class-level transportation problem (variables: lines from each closed port class to each
compatible class; constraints: port counts, mirror symmetry of black multiplicities,
`m_black[C -> D] = m_black[C^rev -> D^rev]`) has

```text
n = 4: 1 solution (the metagraph);   n = 5: 4 solutions, one of which is the real metagraph.
```

So the scores -- first-order information -- do **not** determine the metagraph from `n = 5`
on: the colliding score sequences `(1,2,2,2,3)` (three classes) and `(1,1,2,3,3)` (two classes)
leave seven free port classes. With **second-order** partner data (score multiset plus the
curl profile, the sorted vector of cyclic-triangle counts per vertex) every port class is
forced at `n = 5` (30 forced, 0 free) and the metagraph is determined, although at the class
level two 5-tournaments still share `(scores, curl profile)`: the port structure separates
what the invariants alone do not. At `n = 6` the ladder is not yet enough (CP-SAT, class level,
56 classes):

```text
L1 (scores):                 196 of 229 port classes free;  >= 169,359 consistent metagraphs (time-capped)
L2 (+ curl profile):         153 of 355 port classes free;  >= 200,000 (enumeration cap)
L3 (+ second scores):         37 of 403 port classes free;  exactly 48,000, of which 3 class-level solutions
                                                           project to the real merged metagraph
```

so from `n = 6` on the wiring is determined neither by the transitive layer nor by the curl
profile nor by the second scores; the information the lines carry about the fibre exceeds
third-order local invariants (HYP-9351 (a)). The "intransitivity of information" reading: at
`n = 5` the wiring is invisible to the potential layer and visible to the cyclic layer; at
`n = 6` it needs more of the cyclic 2-form than its per-vertex trace.

## 4. (L) The consecutive-start-port ledgers

Closed endpoint-`k` layer (THM-4734 notation): `P = (a, b, x, 3, ..., y, c, z)`, `d(a) = n-1-k`,
`d(z) = k`, incoming ports `I` (`|I| = k-1`, `1 notin I`), `J = r_P(I)`. Take
`I = {2, ..., k}` (the first `k-1` middle positions after `b`), so `J = {n-1-k, ..., n-3}` (the
`k-1` positions before `c`; `I n J` is the centre for odd `n`), and let the middle be closed
(`c -> b`). Put `epsilon = [b -> y] = [x -> c]` (equal by blue reflection).

**Theorem.** (i) If `epsilon = 0` the packet's closed-no-chord ledger
`Q1 = (c, W, b, z, a), Q2 = (z, c, W, b, a), Q3 = (a, z, b, W, c)` (`W = (x, 3, ..., y)`) is a
three-black detour to `K_P(T)`. (ii) If `epsilon = 1` the **port-block insertion**

```text
Q1 = (a, b, x, U, z, J, c),   Q2 = (a, x, U, z, J, b, c),   Q3 = (a, x, U, J, b, c, z),
U = (3, ..., n-k-2),  J = (n-1-k, ..., n-3) in path order,
```

is a three-black detour to `K_P(T)` for every `n >= k + 5`, `k >= 3`.

*Proof of (ii).* Directedness: in `T`, `a -> b -> x -> U` are path arcs; `last(U) = n-k-2 -> z`
because every middle vertex outside `J` beats `z`; `z -> first(J)`; `J` and `last(J) = y -> c`
are path arcs. After `K_{Q1}` (non-`Q1` pairs complemented): `a -> x` (was `x -> a`, `x in I`);
`x -> U -> z -> J` retained; `y -> b` (was `b -> y`, `epsilon = 1`); `b -> c` (was `c -> b`).
After `K_{Q2}`: `a -> x`, `x -> U`, `J`, `y -> b`, `b -> c` retained; `last(U) -> first(J)` is a
path pair complemented twice (absent from `Q1` and `Q2`), hence restored; `c -> z` likewise.
Path-edge cancellation: `ab, bx` occur once (`Q1`), `x-U` and `U`-internal three times,
`last(U) z` and `z first(J)` twice, `J`-internal three times, `yc` once (`Q1`), `ax` twice,
`yb` twice, `bc` twice, `last(U) first(J)` once (`Q3`), `cz` once (`Q3`); the surviving set is
exactly `E(P)`, so the composite is the literal `K_P(T)`.
Blackness. `Q1` in `T`: the reversal of `Q1` sends `y` (position `n-2`) to `b` and `a` to `c`;
`a -> y` holds (`y notin I` as `n >= k+4`) while `b -> c` fails (`c -> b`): not an anti-automorphism.
`Q2` in `S = K_{Q1}(T)`: the reversal sends `b` to `x` and `a` to `c`; `a -> b` holds (retained
`Q1` edge) while `x -> c` fails in `S` (the pair was complemented and `epsilon = 1` gave
`x -> c` in `T`). `Q3` in `B = K_P(T)`: the reversal sends `b` to `first(U) = 3` and `a` to `z`;
`a -> b` holds (path arc) while `3 -> z` fails in `B` (`3 notin J` for `n >= k+5`, so `3 -> z` in
`T` and the non-path pair is complemented in `B`). `#`

(i) is the packet's ledger; its first word is black by the permutation argument of the
packet's section 3 for `k = 3`, and the exhaustive runs certify all of (i) for `k = 4`
(`n = 9..12`) and `k = 5` (`n = 11, 12`). In every presentation of the pattern exactly one of
the two halves applies (`residual_cover.out`: the inherited bank covers precisely the
`epsilon = 0` half, `128/256`, `1024/2048`), so the pattern `I = {2, ..., k}` is settled by the
pair. The other consecutive-port patterns (`{2, n-3, n-2}`, `{n-4, n-3, n-2}`) do not split by
a single arc (their residual halves need 3-7 templates); they are HYP-9351's next rows.

## 5. Files

`04-computation/experiments/electric_metagraph_20261021/` (`electric_census.py` / `.out`,
`determination.py` / `.out`, `determination_L2.py` / `_n5.out`, `_n6.out`,
`invariant_collisions.py`, `electric_n{4..7}.json`);
`04-computation/experiments/endpoint_layers_20261020/` (`residual_cover.py` / `.out`,
`consecutive_start_ledgers.py` / `_small.out`, `_large.out`).
