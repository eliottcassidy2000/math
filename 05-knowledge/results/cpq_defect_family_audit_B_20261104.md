# Independent adversarial audit B of THM-4871 statement 7 (the first-defect law for primes p < q < 2p)

**Auditor:** independent Claude session (audit B), 2026-11-04, worktree `math-wt-oddprimes-20261104`.
**Audited text:** statement 7 and its proof (**7.** in the Proofs section, facts (i)-(iii) and Steps 1-5) of
`01-canon/theorems/THM-4871-...md`, read against the Setting (sign convention, `N`-evenness, multiplier-evenness),
statements 1-4 with their proofs, THM-4721 (Theorem R, the orbit form of the sign), THM-4728 (the defect identity),
audit A (`cpq_defect_family_audit_A_20261104.md`) and the census outputs in
`04-computation/experiments/cpq_defect_20261104/`. No existing file was modified; the scripts and outputs of this
audit are in `04-computation/experiments/cpq_defect_20261104/audit_B/` (section 6).

## 0. Verdict

**Statement 7: PASS WITH CORRECTIONS.** Both halves of statement 7 are true and the proof is essentially complete:
facts (i), (ii), (iii) are correct and were confirmed computationally; Steps 1-4 (the zero-defect half, 7(a)) are
correct as written; Step 5 (the `k = q - p` half, 7(b)) has **one genuine gap** in the conjugate case (B1 below): the
sentence that excludes a conjugate `H'` with `q`-orbit `P u F` from every `N`-even graph other than `K_(q,q)` and the
empty graph "by the argument just given with `c''` in place of `g c_Q g^(-1)`" does not follow, because that argument
uses the isomorphism `X[P u F] ~= X[Q]` supplied by an automorphism, which `c''` does not supply. The claim is true and
the repair is one sentence using only fact (i) (the `p`-cycle of `H'` lies inside `Q` and fixes `q - p >= 2` points,
impossible in `AGL(1, q)`). Three minor repairs (B2-B4) and one scope remark (B5): the proof of 7(a) does not use
`q < 2p` at all and proves the zero-defect half for every prime pair `p < q` and every `k < min(p, q - p)`.

**Decisive independent confirmations.** (1) Fact (i): every nonempty non-complete circulant graph on `Z_r`,
`r = 5, 7, 11, 13, 17, 19`, is prime (brute force over all subsets for `r <= 13`, over all subsets containing `0` for
`r = 17, 19`), has `|Aut| = r |M(C)|` (so `Aut = Z_r x| M(C) <= AGL(1, r)`), has no twin pair, and `K_r` is never
multiplier-even (`audit_B_fact_i_circulants.out`). (2) Fact (ii): for `r = 3, 5, 7` and EVERY overlap pattern and
EVERY `r`-cycle on the second support (12, 720, 90,720 pairs), `<c, c'>` is 2-transitive on the union; 3000 random
pairs at `r = 11` and 1500 at `r = 13` likewise (`audit_B_fact_ii_two_cycles.out`). (3) An independent census of the
`N`-even `H`-fixed graphs (own circulant classes, own brute-force canonical forms for the coloured objects, nauty
only for the full objects) reproduces the committed rows `(3,5) k <= 3`, `(5,7) k <= 3`, `(7,11) k <= 4`,
`(11,13) k <= 2`, `(3,7) k <= 3`, with a definition-level brute force over all `2^c` objects (direct `N`-evenness,
`N` by closure) agreeing at `(3,5) k <= 2`, `(5,7) k <= 2`, `(7,11) k <= 3`, `(11,13) k <= 2`, `(3,7) k <= 2`,
`(5,11) k <= 1`, `(13,17) k = 0`; in every row with `k < q - p` every `N`-even graph is sign-even with one class, and
at `k = q - p` the only sign-odd `N`-even type is `K_(q,q)` on `Q` against `P u F` (`census_*.out`). (4) Two families
never computed before behave exactly as the proof says: `(13,17)`, `k <= 4`: defects `0, 0, 0, 0, 1`, `W = 1`
(`K_(17,17)`), all excesses `0` (statement 7 at a new pair, confirmed by the author's script too), and `(17,19)`,
`k <= 2`: `0, 0, 1` with the single witness `K_(19,19)`; `(5,11)`, `k <= 4` (`q > 2p`): all zero, which is the range
`k < min(p, q - p) = 5` that Steps 1-4 actually cover (B5). (5) The
objects of Steps 1-3 were built for every `N`-even representative of 26 rows (`M`, `D`, `F'`, the extension `g^`):
`M` and `D` are modules invariant under every nauty generator of `Aut X`, every generator is even, the sign factorises
over `M, D, F'` with the first two factors `+1`, `g^` is a colour-preserving automorphism of `gamma` and
`sgn_gamma(g^) = sgn_(X[F'])(g|F')` -- for all 24,762 representatives except the `K_(q,q)` at `k = q - p`, which is
exactly the exception the theorem names (`audit_B_steps_structure.out`). (6) The gap B1 is real: the graphs
`X_1 = (C_P = empty, C_Q a multiplier-even prime circulant, cross, F independent of colour (0, cross))` at
`k = q - p` are `N`-even, admit the `q`-cycle `c''` on `P u F` as an automorphism, are neither `K_(q,q)` nor empty,
and have `X[Q]` nonempty (18 configurations at `(5,7)`, `(7,11)`, `(11,13)`); what excludes a conjugate of `H` with
`q`-orbit `P u F` is that `Aut X[Q]` (order `q |M(C_Q)|`) contains no `p`-cycle and no non-identity element with two
fixed points, which is the repair.

## 1. Facts (i)-(iii): PASS

**(i) Prime circulants.** Connected and co-connected: `C` nonempty generates `Z_r` (`r` prime), and the complement is
the circulant of the nonempty set `Z_r^* \ C`. Not 2-transitive on vertices (a 2-transitive automorphism group forces
empty or complete), so Burnside's theorem puts `Aut` inside `AGL(1, r)`: correct. Primality via Gallai: in the
connected/co-connected case the maximal proper (strong) modules partition `V` and are permuted by `Aut X`, hence
transitively by vertex-transitivity; equal sizes dividing the prime `r` and `< r` (they are proper) force singletons,
and every proper module lies in one of them: correct. "No vertex of such a circulant is a twin of any other vertex of
`X`": inside `Q` a twin pair would be a module of size 2 of `X[Q]`; for `v` outside `Q`, `N(v) cap Q in {empty, Q}`
(uniform attachment by `H`-invariance) while `N(u) cap Q = u + C_Q` is nonempty and proper: correct (false and true
twins alike). `K_r` not multiplier-even: on `K_r` every pair is an edge, so the reversal count is the inversion
number and the reversal sign is the permutation sign; the dilation by a primitive root is an `(r-1)`-cycle, odd:
correct. Computation (`audit_B_fact_i_circulants.out`): `r = 5, 7, 11, 13` all `2, 6, 30, 62` nonempty non-complete
symmetric sets, `r = 17, 19` all `34, 58` dilation classes: zero nontrivial modules, `|Aut| = r |M(C)|` in every
case, zero twin pairs; `K_r`: reversal sign equals permutation sign for every dilation, the odd dilations are exactly
the non-squares, `K_r` not multiplier-even, the empty graph multiplier-even.

**(ii) Two overlapping cycles.** The block argument of statement 2's proof indeed uses only that `c, c'` are cycles
of prime length with distinct overlapping supports (re-derived: the per-cycle dichotomy "support inside one invariant
block / every block meets the support in at most one point" is exhaustive and the three combinations are
contradictory), so `K` is primitive on `Delta` without any parity hypothesis. `|Delta| >= r + 3`: Jordan's theorem in
Wielandt's form 13.9 (primitive, degree `m`, a prime cycle of length `r <= m - 3`) gives `A_Delta <= K`, 2-transitive
for `m >= 4`. `|Delta| = r + 1`: the stabiliser of the unique point `a` outside `O` contains `c`, transitive on
`Delta \ {a} = O`: correct. `|Delta| = r + 2`: `K_a >= <c>` is transitive on `Delta \ {a, b}`; if `K_a` fixed `b`,
then `K_a <= K_b`, equal orders by transitivity, so `K_a = K_b`, `Fix(K_a) = {a, b}` (`c` moves everything else);
`g` with `g(a) = b` satisfies `g K_a g^(-1) = K_b = K_a`, so `g` permutes `Fix(K_a)` and `g(b) = a`; then `{a, b}` is a
block (`g(a) in {a, b}` or `g(b) in {a, b}` forces `g{a, b} = {a, b}`, using `K_a = K_b`), of size `2 < r + 2`,
against primitivity: correct. 2-transitivity of `K <= Aut X` on `Delta` makes `X[Delta]` empty or complete, and
`K` fixing `V \ Delta` pointwise makes `Delta` a module: correct. Computation (`audit_B_fact_ii_two_cycles.out`):
exhaustive for `r = 3` (12 pairs, `|Delta| = 4, 5`), `r = 5` (720 pairs, `|Delta| = 6..9`), `r = 7` (90,720 pairs,
`|Delta| = 8..13`), every group 2-transitive on `Delta`; sympy orders for `r = 3, 5`: `A_m` in all cases except
20 of the 120 pairs with `|Delta| = 6 = r + 1`, where `K = PSL(2,5)` of order 60 acting 2-transitively on 6 points --
which is why the proof is right to claim only 2-transitivity (not `A_Delta`) in the case `|Delta| = r + 1`.
Sampled `r = 11` (3000 pairs, all `|Delta| = 12..21`), `r = 13` (1500 pairs), and the two block-argument cases
`|Delta| = 12, 13` at `r = 11` over all `10`- and `9`-subsets: no failure.

**(iii) Signs across invariant sets.** The orbit form (Royle et al. Lemma 4.3, re-derived: along a `<g>`-orbit of
edges of length `L` the reversal parities telescope to `[g^L` swaps the ends`]`, so an orbit contributes `1` iff it is
self-paired) shows that an edge between two `g`-invariant parts is never in a self-paired orbit, and the orbits of
edges inside a part are the orbits of the restriction: `sgn_X(g) = prod_parts sgn_(X[part])(g|part)`. Correct, and
order-free (each factor may be computed in any vertex order). Verified on every generator of every `N`-even
representative of the 26 rows of section 4 (`audit_B_steps_structure.out`, zero failures).

## 2. Steps 1-4 (statement 7(a)): PASS

**Step 1.** `|g(Q)| = q > p + k = |V \ Q|` (this is `k < q - p`) gives `g(Q) cap Q != empty`; if `g(Q) != Q`, (ii)
applied to the `q`-cycles `c_Q` and `g c_Q g^(-1)` (both in `Aut X`, distinct overlapping supports) makes
`Q u g(Q)` a module inducing the empty or complete graph; complete contradicts `C_Q` multiplier-even (fact (i): `K_q`
is not), so `Q u g(Q)` is independent, impossible when `X[Q]` is a nonempty circulant. So `X[Q]` prime circulant
`=> g(Q) = Q` for all `g`; `X[Q]` empty `=> M :=` union of the independent modules containing `Q` (the union lemma is
correct: a vertex of one outside the other is non-adjacent to the intersection, hence by uniform attachment to the
whole other module; the union is a module because an outside vertex is uniform on both and they overlap), and
`g(M) = M` by maximality (`g(M)` is an independent module meeting `M` since `g(Q)` meets `Q`). Correct. (Equivalent
description used in the computations: when `X[Q]` is empty, `M` is the false-twin class of a vertex of `Q`.)

**Step 2.** `C_P` prime circulant: `P cap M = empty` by the twin clause of (i) when `M` is independent (a vertex of
`P` in `M` would be a false twin of the vertices of `Q`), trivially when `M = Q`; `g(P) <= V \ M <= P u F` and
`|F| = k < p` give `g(P) cap P != empty`; `g(P) != P` would make `P u g(P)` empty or complete by (ii), impossible:
`g(P) = P`. `C_P` empty: `P` is a set of pairwise false twins (common neighbourhood `(Q if cross) u {f : x_f = 1}`);
if `P cap M != empty` then `P u M` is a set of pairwise false twins, hence an independent module containing `Q`, so
`P <= M` (this needs `X[Q]` empty); otherwise `D :=` union of the independent modules containing `P`, `D cap M =
empty`, and `g(D) = D` by the argument of Step 1 (`g(P) <= V \ M`, `|F| < p`, then `P u g(P)` independent `<= D`).
Correct, with one justification slip (B2): the parenthesis "(otherwise `D u M` is an independent module containing
`P` and `Q`, forcing `P <= M`)" applies only when `M` is independent; when `M = Q` is a prime circulant, `D cap Q =
empty` holds because a vertex of `Q` in the independent module `D` would be a false twin of the vertices of `P`,
against (i). The conclusion is right in both cases.

**Step 3.** `V \ (M u D) = F'` since `Q <= M` and `P <= D` in every case; `g` preserves `M`, `D`, `F'`; (iii) gives
the factorisation. `sgn_(X[M])(g|M) = +1`: `X[M]` edgeless, or `M = Q` with `g|Q in Aut X[Q] <= AGL(1, q)`, a
translation (odd order, sign `+1`) times a dilation of `M(C_Q)` (sign `+1` by multiplier-evenness; the dilation is in
`Aut X[Q]` because the translation is), using that `sgn` is a homomorphism on `Aut X[Q]`; likewise for `D = P`. The
extension `g^`: for `f in F'`, `x_f = [f ~ P] = [f ~ D]` and `y_f = [f ~ Q] = [f ~ M]` by uniform attachment to the
modules `D >= P`, `M >= Q`, both preserved; a vertex of `F cap M` is a false twin of the vertices of `Q`, so its
`F'`-neighbours are `{f in F' : y_f = 1}` (preserved) and its colour is `(cross, 0)` (`x = [u ~ P] = cross` for
`u in Q`, `y = 0` as `M` is independent); a vertex of `F cap D` (`D != M`) is a false twin of the vertices of `P`,
neighbours `{f in F' : x_f = 1}`, colour `(0, cross)`; edges between `F cap M` and `F cap D` are all present or all
absent according to `cross` (twins of a `Q`-vertex and of a `P`-vertex); none inside `F cap M`, `F cap D`; `g`
preserves the edges inside `F'`. So `g^ in Aut(gamma)` colour-preserving; colour-sign-evenness of `gamma` (the
`N`-even hypothesis via 1(c)) gives `sgn_gamma(g^) = +1`; (iii) on `F = (F cap M) u (F cap D) u F'` with `g^` the
identity on the first two parts gives `sgn_gamma(g^) = sgn_(X[F'])(g|F')`. Hence `sgn_X(g) = +1`. Correct. In the
case `D = M` (`P <= M`) the colours are `(0, 0)` and the partition is `{M, F'}`: covered. Verified on every generator
of every `N`-even representative of the 26 rows (section 4).

**Step 4.** `Q'` meets `Q` (same count as Step 1). `X[Q]` prime circulant: `Q' = Q` by (ii) (otherwise `Q u Q'`
independent), and `c''|Q`, a `q`-cycle in `Aut X[Q] <= AGL(1, q)`, lies in the unique Sylow `q`-subgroup
`<c_Q|Q>`: `c'' in <c_Q>`. `X[Q]` empty: `Q u Q' <= M` (by (ii), or `Q' = Q`), `Sym(M) <= Aut X` (independent
module) carries `Q'` to `Q`, then `Sym(Q)` (fixing `V \ Q`, hence the `p`-orbit, pointwise) conjugates `c''` into
`<c_Q>`. `P' <= V \ Q = P u F`, `|F| < p`, so `P'` meets `P`. `C_P` prime circulant: `P' = P` by (ii), `c'|P in
<c_P|P>` by Sylow in `AGL(1, p)`. `C_P` empty, `P cap M = empty`: `P u P' <= D` by (ii), `Sym(D)` fixes `Q`
pointwise (`D cap M = empty`) and carries `P'` to `P`, `Sym(P)` conjugates `c'` into `<c_P>`. `C_P` empty, `P <= M`:
`c'` preserves `M` (Step 1 applies to every automorphism), `P'` meets `M` (`|P'| = p > |F \ M|`), so `P' <= M`;
`Sym(M)` carries `(P', Q)` to `(P, Q)` and `Sym(P) x Sym(Q)` finishes. In every case `a H' a^(-1) = <c_P^i c_Q^j> = H`
by CRT (`i`, `j` units). Correct. Hence `m_X(H) = 1` for every sign-even type with `H`-symmetry (its `H`-fixed
representative is `N`-even), `E_even = 0`; with Step 3, `W = E_wit = 0`; with 3(b), `E_T = 0`; THM-4728's identity
gives `tau_T = tau_E`. Correct.

The hypothesis `q < 2p` enters 7(a) only through `k < p`, which Steps 2 and 4 (and 3(b)) use; the parenthetical
"(`k < q - p < p`)" in the statement attaches it to 3(b) alone. See B5.

## 3. Step 5 (statement 7(b)): PASS WITH ONE GAP (B1) and two minor repairs (B3, B4)

*Automorphisms.* `|V \ Q| = q`, so `g(Q) cap Q = empty` forces `g(Q) = P u F`; then `X[P u F] ~= X[Q]` via `g`, and
`g^(-1)(P)` is a module of `X[Q]` of size `p`, `1 < p < q` (`P` is a module of `X` by `H`-invariance), so by (i)
`X[Q]` is empty or complete, complete excluded by multiplier-evenness: `X[Q]` and `X[P u F]` are empty; `Q` and
`P u F = g(Q)` are modules, so `X in {K_(q,q), empty}`. Correct. For every other `N`-even `X` all automorphisms
overlap `Q`, and Steps 1-3 use only that overlap and `k < p` (`= q - p < p` iff `q < 2p`): `X` is sign-even. The
empty graph is sign-even and `K_(q,q)` is sign-odd and `N`-even (4(a)), so `W = 1`. Correct.

*Conjugates.* **B1 (gap).** Quote: "Step 4 applies to every `N`-even `X` other than `K_(q,q)` and the empty graph (a
conjugate with `Q' cap Q = empty` has `Q' = P cup F` and forces, by the argument just given with `c''` in place of
`g c_Q g^(-1)`, `X in {K_(q,q), empty}`)". With `c''` (the `q`-cycle of `H'`, support `P u F`, fixing `Q` pointwise)
in place of `g c_Q g^(-1)` the argument just given yields: `X[P u F]` is a circulant for `c''` containing the module
`P`, hence empty by (i) and the multiplier-evenness of `C_P`; `P u F` is a module (as the orbit of `c''`), `Q` is a
module. It does **not** yield `X[Q]` empty: that step used the isomorphism `X[P u F] ~= X[Q]` induced by the
automorphism `g`, and `c''` induces no such isomorphism. The configurations `X_1` (`C_P = empty`, `C_Q` a
multiplier-even prime circulant, `P u F` independent, joined to `Q` or not) are `N`-even, have `c'' in Aut X_1`, are
neither `K_(q,q)` nor empty, and have `X_1[Q]` nonempty: `audit_B_steps_structure.out` lists all 18 of them at
`(5,7)`, `(7,11)`, `(11,13)` with `|Aut X_1| = |Aut X_1[Q]| q!`. So the sentence as written does not exclude a
conjugate `H' <= Aut X_1` with `Q' = P u F`. **Repair (one sentence, using only (i)):** such an `H' = <c' c''>` has
`P' = supp(c') <= V \ Q' = Q`, so `c'|Q` is a `p`-cycle in `Aut X[Q]` fixing the `q - p >= 2` points of `Q \ P'`; if
`X[Q]` were a prime circulant, `Aut X[Q] <= AGL(1, q)` by (i), whose non-identity elements fix at most one point --
contradiction; so `X[Q]` is empty (complete is excluded by multiplier-evenness), `P u F` is independent, both are
modules, and `X in {K_(q,q), empty}`. (Equivalently: `|Aut X[Q]|` divides `q(q - 1)`, which `p` does not divide for
`q < 2p`; the fixed-point form is the one that does not need `q < 2p`.) Verified: for all 18 `X_1`, `Aut X_1[Q]`
(enumerated by closure, orders `14, 22, 26, 52`) contains no `p`-cycle and no non-identity element with more than
one fixed point. With this sentence inserted, Step 4 applies to every `N`-even `X` other than `K_(q,q)` and the empty
graph, and 7(b)'s "every `N`-even type carries a single `Aut`-class" is proved for them. (The census confirms the
conclusion: `E_even = 0` in every `k = q - p` row, i.e. the `X_1` carry one class.)

**B3 (minor, `K_(q,q)`).** Quote: "in `K_(q,q)` a `q`-cycle preserves the bipartition (an odd power of a
part-swapping automorphism swaps the parts), so `Q'` is a part, `P'` lies in the other, and `Sym(part) x Sym(part)
<= Aut X` conjugates `H'` to `H`". The parity argument is right (`c''^q = 1` cannot swap the parts, so `c''`
preserves both and its single orbit of size `q` is a whole part). But `Sym(A) x Sym(B)` preserves the parts, so it
cannot carry a conjugate with `Q' = P u F = A` to `Q = B`; the side swap `sigma in Aut X` is also needed:
`(Sym(A) x Sym(B)) x| <sigma> = Aut X` acts transitively on the pairs (`Q'` a part, `P'` a `p`-subset of the other
part). Repair: "`Aut X = (Sym(A) x Sym(B)) x| <sigma>`". Alternatively, by datum uniqueness: an `H`-fixed `K_(q,q)`
has parts that are unions of `H`-orbits (`H` of odd order fixes both parts of the canonical bipartition) of size
`q`, and the only such split of `{P, Q, {f}}` is `{Q, P u F}`; so exactly one `H`-fixed graph is isomorphic to
`K_(q,q)`, `m_X(H) = 1` (verified by brute force over all `2^c` `H`-fixed graphs at `(3,5) k = 2` and `(5,7) k = 2`:
exactly one each). Either way `E_wit = 0`.

**B4 (wording of the conclusion; statement 7 and the title).** Quote: "the least degree at which the pair `(p, q)`
violates the refined equinumerosity HYP-9340 (a) is exactly `n = 2q`". What is proved: for `H = C_pq` of orbit type
`(p, q, 1^k)` the defect is `0` for `k < q - p` and `1` at `k = q - p`, so the least `k` (least degree `p + q + k`) at
which *this family* violates HYP-9340 (a) is `k = q - p` (`n = 2q`). Other odd-order groups of order divisible by `pq`
can act on fewer than `2q` points when `q - p > 3`, e.g. `C_pq x C_3` of type `(p, q, 3)` on `p + q + 3` points
(`(7, 11)` at `n = 21 < 22`), and nothing is proved about them; for `(3, 5)` the Codex classification of all
odd-order classes through `n = 10` makes "least failure degree ten" absolute, but not for a general pair. Repair:
"the least degree at which `C_pq` of type `(p, q, 1^k)` violates HYP-9340 (a) is `2q`" (also in the title's
"(7) ... so the least degree at which `(p, q)` violates HYP-9340 (a) is `2q`").

## 4. Consistency with the census and the independent computations

Committed `.out` files versus statement 7 (columns `W, E_T, E_even, E_wit`): `(3,5)`: `k = 0, 1` all `0`, `k = 2`:
`1, 0, 0, 0`; `(5,7)`: `k = 0, 1` all `0`, `k = 2`: `1, 0, 0, 0`; `(7,11)`: `k = 0..3` all `0`, `k = 4`: `1, 0, 0, 0`;
`(11,13)`: `k = 0, 1` all `0`, `k = 2`: `1, 0, 0, 0` -- exactly 7(a) and 7(b) (`cpq_3_5_k6.out`, `cpq_5_7_k6.out`,
`cpq_7_11_k4.out`, `cpq_11_13_k2.out`). The theorem's statement 5 table and its sentence on the `(7,11) k = 4` and
`(11,13) k = 2` rows agree with the files.

Independent computations of this audit (`audit_B_statement7_census.py`; own orbit method with brute-force
canonical forms for the coloured objects, checked against Burnside's `b_k(4)` and Theorem R's `ec_k(4) = b_k(4)`;
definition-level brute force with `N` by closure, direct `N`-evenness, types by certificate, `N`-orbits as the
components under the generators):

| family | rows | orbit method | brute force (all `2^c` objects, direct `N`-evenness) |
|---|---|---|---|
| `(3,5)` | `k <= 3` | `0, 0, 1, 4` with `W = 0, 0, 1, 4`; `k = 2`: the witness is `K_(5,5)` on `Q` vs `P u F` | `k <= 2` AGREE (`|N| = 120, 120, 240`) |
| `(5,7)` | `k <= 3` | `0, 0, 1, 4`; `k = 2`: witness `K_(7,7)`, `E_T = E_even = E_wit = 0` | `k <= 2` AGREE (`|N| = 840, 840, 1680`) |
| `(7,11)` | `k <= 4` | `0, 0, 0, 0, 1`; `k = 4`: witness `K_(11,11)`, excesses `0` | `k <= 3` AGREE (`|N| = 4620, 4620, 9240, 27720`; `2^18` objects at `k = 3`) |
| `(11,13)` | `k <= 2` | `0, 0, 1`; `k = 2`: witness `K_(13,13)`, excesses `0` | `k <= 2` AGREE (`|N| = 17160, 17160, 34320`; `2^17` objects at `k = 2`) |
| `(3,7)` | `k <= 3` | `0, 0, 0, 0` with `W = 0, 0, 0, 4`, `E_T = 0, 0, 0, 8`, `E_even = 0, 0, 0, 4` | `k <= 2` AGREE |
| `(5,11)` NEW | `k <= 4` | all zero (`W = E_T = E_even = E_wit = 0`), `R = 8, 32, 128, 704, 5632` | `k <= 1` AGREE; author's `cpq_defect_census.py 5 11 4` agrees row by row |
| `(13,17)` NEW | `k <= 4` | `0, 0, 0, 0, 1`; `k = 4`: `W = 1`, the witness `K_(17,17)` (`C_P = C_Q = empty`, `cross = 1`, all four fixed points of colour `(0, 1)`, 289 edges), excesses `0`; `R = 192, 768, 3072, 16896, 135168` | `k = 0` AGREE (`|N| = 42432`); author's script agrees row by row |
| `(17,19)` NEW | `k <= 2` | `0, 0, 1`; `k = 2`: `W = 1`, the witness `K_(19,19)` (361 edges, `|Aut| = 2 (19!)^2`), excesses `0`; `R = 960, 3840, 15360` (`t_17 = 16`, `t_19 = 30`) | -- |

All rows satisfy Theorem R and THM-4728's identity. The structural checks (`audit_B_steps_structure.out`) over the
24,762 `N`-even representatives of `(3,5) k <= 2`, `(5,7) k <= 2`, `(7,11) k <= 4`, `(11,13) k <= 2`, `(3,7) k <= 2`,
`(5,11) k <= 4`, `(13,17) k <= 2` find every object of Steps 1-3 as the proof describes it (case counts of
`M = Q` / `M` independent and `D = P` / `D = M` / `D` independent are printed), with exactly 12 "problems", all for
the `K_(q,q)` at `k = q - p` of the four prime pairs (generators move `M` and `D`, one generator is odd): the theorem's
own exception. The only representatives with a nauty generator carrying `Q` off itself are these `K_(q,q)`.

## 5. Scope and citations (item 10)

**B5 (the proof gives more than the statement; the Scope bullet is misleading).** Steps 1-4 use exactly two
inequalities: `p + k < q` (Step 1 and the first line of Step 4) and `k < p` (Steps 2 and 4, and 3(b)). Neither needs
`q < 2p`. Hence **7(a) holds for every pair of primes `p < q` and every `k < min(p, q - p)`**: every `N`-even
`H`-fixed graph is sign-even with one class and `tau_T(H) = tau_E(H)`. For `q < 2p` this is `k < q - p` (the stated
range); for `q > 2p` it is `k < p`, which the census confirms in every existing row (`(3,7)`, `(3,11)`, `(3,13)`
`k <= 2`; `(7,13) k <= 3`) and in the new family `(5,11) k <= 4` computed here (two implementations). The Scope bullet
"Statement 7's Step 2 uses `|F| < p` and Step 1 uses `|V \ Q| < q`, i.e. `k < min(p, q - p) = q - p`; the proof does
not extend to `q - p >= p`, where `g(Q)` can miss `Q` without `X` being `K_(q,q)` or empty" conflates the two halves:
the zero-defect half extends verbatim to `k < min(p, q - p)` for `q > 2p`; only the identification of the first
defect (Step 5, which needs `k = q - p < p`) does not, and there the first defect at `k = p <= q - p` is indeed a net
value (statement 4(b), `(3,7) k = 4`). Suggested repair: add "(a') for all primes `p < q` and `k < min(p, q - p)`
the conclusions of (a) hold (same proof)" and reword the bullet; this upgrades the cited `q > 2p` census rows from
FINITE-EXACT to PROVED and makes `(5,11) k <= 4` a proved row. Nothing in Steps 1-5 silently needs `q - p < p`
beyond `k < p` at `k = q - p` (B1's repair in its fixed-point form also needs only `q - p >= 2`).

**Citations.** Burnside's prime-degree theorem (transitive of prime degree: 2-transitive or inside `AGL(1, r)`) is
used correctly in (i) and in the Sylow steps of Step 4 (the unique Sylow `q`-subgroup of a subgroup of `AGL(1, q)`
containing a `q`-cycle is the translation group). Jordan (Wielandt 13.9) is used with its exact hypothesis in (ii)
(`r <= |Delta| - 3`, primitivity established first). Gallai's modular decomposition is used in its classical form
in (i). The orbit form of the sign (Royle et al. Lemma 4.3) is correctly stated in 4(a)'s proof and used in (iii).
Neither Feit-Thompson nor Jones is needed for statement 7 (fact (ii) avoids the `|Delta| = r + 2` affine case by the
stabiliser argument). The `(iii)`-free parts of the Setting (homomorphism and order-independence of the sign) are as in
THM-4721 and audit A.

**Audit A corrections (quick verification).** In THM-4871: title, statement 3(c) and the proof of 3 now read
`E_T(p) = 4 t_p^2 t_q` and `8 t_p^2 t_q` `N`-orbits (the only remaining "`4 t_p t_q`" is the status line describing
the correction); "62 rows, thirteen families" in the title, statement 5 and its closing sentence; Jones = Bull. Aust.
Math. Soc. 89 (2014) 159-165 in the related list; the status describes the brute-force agreements as audit A
corrected them. HYP-9463 (1) says "PROVED for primes p < q < 2p (THM-4871 statement 7 ...)" in the title and the
status, with the unique-witness description restricted to the proved range. **Residuals** (outside statement 7, for
the author): HYP-9463's title still says in statement (3) "now with 64 further rows" (its status and body say 62),
and the note `cpq_defect_family_mirror_folds_20261104.md` line 10 still types "FINITE-EXACT (64 census rows)".

## 6. Scripts and outputs (`04-computation/experiments/cpq_defect_20261104/audit_B/`)

| file | content |
|---|---|
| `auditb_common.py` | own helpers: symmetric/skew sets and dilation classes, multiplier stabiliser and multiplier-evenness, brute-force canonical forms and colour-sign-evenness of 4-coloured graphs/tournaments on `k <= 4` points (min over `S_k`), Burnside's `b_k(4)`, object builders, nauty wrappers |
| `audit_B_fact_i_circulants.py/.out` | fact (i): modules by brute force, `|Aut| = r |M(C)|`, twins, `K_r` signs, `r = 3..19` |
| `audit_B_fact_ii_two_cycles.py/.out` | fact (ii): 2-transitivity of `<c, c'>` for every overlap and every second cycle, `r = 3, 5, 7` exhaustive, `r = 11, 13` sampled; sympy orders for `r = 3, 5` |
| `audit_B_statement7_census.py` | independent orbit census and definition-level brute force (`--brute kb`) |
| `census_3_5_k3.out`, `census_5_7_k3.out`, `census_7_11_k4.out`, `census_11_13_k2.out`, `census_3_7_k3.out` | the committed families (brute force at `k <= 2, 2, 3, 2, 2`) |
| `census_5_11_k4.out`, `census_13_17_k4.out`, `census_17_19_k2.out` | the new families (`q > 2p` zero range through `k = 4`; two new `p < q < 2p` pairs with their first defect `1` at `k = q - p` and the witness `K_(q,q)`) |
| `author_script_cpq_5_11_k4.out`, `author_script_cpq_13_17_k4.out` | the author's `cpq_defect_census.py` on the new families (cross-check) |
| `audit_B_steps_structure.py/.out` | Steps 1-3 objects on every `N`-even representative of 26 rows; the Step 5 gap configurations `X_1`; datum uniqueness of `K_(q,q)` |

(Display note: nauty returns `|Aut|` as mantissa x `10^exponent`; the census outputs print the mantissa for groups
beyond `10^10` -- e.g. `7755157608` stands for `2 (13!)^2 = 7.755 x 10^19` -- the structure output prints the product.)

## 7. Summary of required corrections (for the author)

1. **Step 5, conjugate case (B1):** insert the fixed-point sentence: a conjugate `H' = <c' c''>` with `Q' = P u F`
   has `P' <= Q`, so `c'|Q in Aut X[Q]` is a `p`-cycle fixing `q - p >= 2` points, impossible in `AGL(1, q)` if
   `X[Q]` were a prime circulant (fact (i)); hence `X[Q]` is empty, and the rest of the sentence follows.
2. **Step 2 (B2):** "`D cap M = empty`" is justified by the union lemma only when `M` is independent; when `M = Q`
   is a prime circulant cite (i) (no vertex of `Q` is a twin of a vertex of `P`).
3. **Step 5, `K_(q,q)` (B3):** "`Sym(part) x Sym(part)`" must include the side swap (or argue by the uniqueness of the
   `H`-fixed `K_(q,q)`).
4. **Statement 7 / title (B4):** "the least degree at which the pair `(p, q)` violates HYP-9340 (a) is `2q`" ->
   "... at which `C_pq` of type `(p, q, 1^k)` violates ...".
5. **Scope (B5, optional but recommended):** record that Steps 1-4 prove (a) for all primes `p < q` and
   `k < min(p, q - p)`; reword the Scope bullet; the `q > 2p` zero rows (and `(5,11) k <= 4`) become PROVED.
6. **Residual "64"** in HYP-9463's title statement (3) and in the note's truth-typing line.
