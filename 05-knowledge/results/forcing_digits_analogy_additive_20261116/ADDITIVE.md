# An additive construction: gluing without mixed odd cycles (THM-5261)

Session klein-2026-11-16-force. Scripts `scripts/glue.c`, `scripts/additive.py`, `scripts/glue_features.py`; receipts
`glue_33.out`, `glue_44.out`, `glue_55s.out`, `glue_45.out`, `glue_56s.out`, `additive.out`.

The owner asked for an additive construction "like the root graft". The repository's grafts (Session 225's single-root
graft, Session 227's set-root grafts `G(A, R; B)` with the positive cover expansion `H = Σ_l c_l(A, R) F_B(l)`) are
positive combinations of path-cover counts, not sums of counts; THM-5256 shows that one shared vertex never gives less
than the product; and THM-5312's rows `p + q − 1, p + q + 1, p + q + 3` were found by census. Here is the mechanism
behind the additive ones.

## 1. The no-mixed-odd-cycle formula (PROVED)

Let `U` be a tournament, `A`, `B` induced subtournaments with `V(A) ∪ V(B) = V(U)`, `W = V(A) ∩ V(B)`. For `X ⊆ W`
let `H_A(X)` be the weight (THM-002: `2^(number of cycles)`) of the odd-cycle packings of `A` in which the vertices of
`W` lying on cycles are exactly `X`; so `H(A) = Σ_X H_A(X)`.

**Theorem (THM-5261 (A)).** If no odd directed cycle of `U` meets both `A ∖ W` and `B ∖ W`, then

    H(U) = Σ_{X, Y ⊆ W, X ∩ Y = ∅} H_A(X) · H_B(Y).

*Proof.* Every odd cycle of `U` lies inside `A` or inside `B`; a packing of `U` is a pair (packing of `A`, packing of
`B`) whose cycles are vertex-disjoint, i.e. whose cycle-covered `W`-sets are disjoint, with the remaining `W`-vertices
singletons in both; the weight multiplies. ∎

**Corollary (additivity).** If moreover `H_A(∅) = H_B(∅) = 1` (every odd cycle of `A`, and of `B`, meets `W`) and no
nonempty `X`, `Y` with `H_A(X) H_B(Y) > 0` are disjoint — for instance every odd cycle of `A` and of `B` passes through
one common vertex `w ∈ W` — then

    H(U) = H(A) + H(B) − 1.

With `|W| = 1` the hypothesis can never hold (two cycles through `w` always splice into an odd mixed cycle, which is
THM-5256's content), so additive gluing needs `|W| ≥ 2`, and the two triangles on an arc are the first instance:
every cycle of either triangle passes through `u` and `v`, the two mixed cycles `u v a b`, `u v b a` have length 4,
and `3 + 3 − 1 = 5`. The `+1` and `+3` rows are the cases with exactly one or two disjoint pairs `(X, Y)` of weight
2 each (`2·2·1 = 4`… precisely `H(U) − H(A) − H(B) + 1 = Σ_{X,Y ≠ ∅ disjoint} H_A(X) H_B(Y)`).

## 2. Search (FINITE-EXACT)

All pairs of tournament classes glued on a shared arc with every cross orientation (`glue.c`): orders (3,3): 72
instances, 30 additive cases `H = 3 + 3 − 1 = 5` (all attaining the path-splice bound) and 12 further cases with the same
count not attaining it; (4,4): 9,216 instances, 270 + 108 cases of `p + q − 1`, 36 + 456 of `p + q + 1`, 54 + 442 of
`p + q + 3`; (4,5): 184,320 instances, 744 + 426 / 81 + 1,272 / 204 + 1,251. **Every additive case with `p, q ≥ 3`
in these ranges is two 3-cycles sharing the arc, with sources or sinks added** (`glue_features.out`: all `(3, 3) → 5`,
bound attained). The strong pairs of orders (5,5) (1,843,200 instances) and (5,6) (129,024,000 instances) have **no**
case with `H(U) − H(A) − H(B) ∈ {−3, −1, 1, 3}` at all. Pre-registered P-H is confirmed only in the weak sense (additive
gluings beyond the bare triangles exist, by padding) and refuted in spirit: no new additive mechanism appears below
order 6, exactly as §1 predicts (mixed odd paths appear as soon as a side has two vertices outside the overlap).

## 3. Reading

Additive gluing is a statement about where the odd cycles of the two pieces sit relative to the overlap, and about the
parity of the mixed cycles; the formula of §1 computes every such gluing exactly and explains THM-5312's rows:
`p + q − 1` needs all cycles of both pieces through one shared vertex, `p + q + 1` one disjoint pair of cycles, and so
on. A family with `p + q − 1` for arbitrary prime `p, q` would need pieces whose odd cycles all pass through a common
vertex and whose mixed cycles are all even; the first condition is the one-apex-like structure, the second is what
fails as soon as either side has more than one vertex outside the overlap (an odd mixed path exists), which is why
the exact additive gluings are so rare in the census.
