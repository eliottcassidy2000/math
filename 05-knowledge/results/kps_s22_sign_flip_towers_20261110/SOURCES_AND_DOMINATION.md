# Source interfaces, the domination-defect inequality, and an Euler-function criterion

kind-pasteur-2026-11-10-S22, continuation of 2026-11-11. Owner directive: *keep going; ponder how the universal
domination-defect inequality reframes the constructive task of a Chan-Pak proof through a uniform supply of reusable
source interfaces, including their initial cases; relate the criterion "RH iff n/phi(n) < e^gamma log log n +
e^gamma (4 + gamma - log 4 pi)/sqrt(log n) for all n >= p_120569#" to areas the repo has touched; two papers attached
(Sulyma, arXiv:2611.12043; Chen-Mogielnicki-Ono, arXiv:2607.12873).*

Checks: `scripts/sources_check.py` (S1-S8, all pass, no input files). Canon: THM-5130 (H). Chan-Pak remains OPEN.

## 0. Summary

| statement | status |
|---|---|
| Every interface (pair `(B, O)`) splits as a domination skeleton plus a defect: `d(G_t) = rho_dom + 2^t delta_dom`, `delta_dom = sum over nonempty S in O of 2^(number of out-neighbours of S outside O)`, and `(r - h - rho_dom + 1) + 2^t (M - delta_dom) >= 0` | PROVED (lemma here + the inequality of Session261) |
| `M = delta_dom` exactly when `B` is transitive | FINITE-EXACT (3,978 pairs of order <= 6; 32,767 labelled pairs of order <= 5) |
| The zero-defect interfaces are the one-apex ones: `s(TT_n, O) = (1 - w(O), sum 2^(n-o), 1, [n in O])` | PROVED (S7) |
| One-vertex gadget on a seed of slope `2^k`: source `(D, U) = (1 + 2^(a+k), 2^(a-p)(1+h) + 2^(p+k) 3^(a-p) - [p=1]E)`; three of the five 9-grid sources of Session261 are of this form | PROVED (corollary of THM-5130 (A)) |
| A one-vertex source of the `d`-grid at level `s` has slope `d 2^v` with `s - floor(log2 d) <= v`; its excess `U - H` is one of `floor(log2 d) + 1` numbers `(d - 2^i) 2^(s-i) - 1` | PROVED (from `U_B(o) >= H(B)`) |
| 9-grid: sources at levels 7, 8, 9, 10, 11, so `1 + 2^s + 9 * 2^(r+s)` is a path count for all `r >= 1`, `s <= 11` | FINITE-EXACT (levels <= 6: Session261) |
| 43-grid complete: `1 + 2^s + 43 * 2^(r+s)` is a path count for all `r, s >= 1` | PROVED from one carrier + four sources + one head (after the grids 3, 21, 35 of other sessions) |
| Abundancy index `a(X) = Pi_2(X)/H(X)`: multiplicative over ordered joins, a subset sum with all terms at most 1, `n + 1 <= a(X) <= 2^n` with equality on the right only for transitive `X`; the closing law is `N = H(X) (r + M a(X))` | PROVED (lower bound KNOWN, Session197) |
| A uniform source family for the 9-grid; a seed for every odd `c`; Chan-Pak | OPEN |

## 1. Dictionary and credit

Sessions 258-261 and the multiplicand packet write `(H, W, U)` for my `(h, h - r, M)` and `(D, U)` for `(r, M)`; a
*source* or *carrier* is a selected base used through `H = D + 2^t U`; a *port* is an unselected vertex where a module
is substituted; an *interface* is the list of statistics retained. Their stated target (multiplicand REVIEW.md):
"a uniform family of physical source interfaces and their initial cases. Composition through arbitrary modules must
retain the output's Q2 as well as H."

Already there when this note was written (first or concurrent; not claimed here):

- unit sets lie in one strong core with transitive padding (Session250, Session258); `U >= |O| H` (Session258);
  `|Aut(B, O)|` divides the offset (odd_c_seed_symmetry);
- the twin identity `W = 2 H(T/xy)`, the bracket `3^c - 2^(c+1) + 1` for a selected pair inside a transitive module
  of one five-vertex carrier, unit iff `c` is 0 or 1, and the slopes `15, 25, 39, 65` times `3^b` (Session260
  coverage/PROOF.md, 2026-11-11); the lines `1 + (b 3^a - 1) 2^s`, `b` = 30, 50, 78 (endpoint_coverage_compiler);
- the carriers with slope families `16 + 24 * 2^t`, `4 + 42 * 2^t`, `32 + 35 * 2^t` and their arbitrary-module form
  `H(M) + (a + b 2^t) Q_2(M)` (multiplicand OFFSET_TRANSPORT.md); the complete 3-, 21- and 35-grids;
- order-eleven roots for the 33 lines of HYP-9658 (B) (odd_c_seed_symmetry, first);
- the domination inequality `H(T) >= d(G_z)` with equality iff `T - z` transitive (Session260), its quantitative form
  `d(G) - 1 <= 2 P_even <= H(T) - H(T - z)` with nonnegative fibres (Session261 EXTENSION4.md) and a second proof
  (multiplicand exchange/PROOF.md).

Specific to this session (not found elsewhere by the sweep of 2026-11-11): the transport law for the whole state
through a new closer with re-selection inside the inserted tournament (THM-5130 (A)); the sign flip for an arbitrary
base and the equation `2^delta - 3^(delta-1) = 1` (pre-registration pushed 2026-11-10, before Session260's carrier
instance); the tower recurrence; the anti-seeds; and the items of this note.

## 2. An interface is a domination skeleton plus a defect

For a pair `(B, O)` with `O` nonempty, root the closed tournament `T_t = (B => TT_t => v, v -> O)` at its closer.
The graph of Sessions 260-261 is `G_t`: left side `O`, right side `(B - O) + tail`, `o ~ x` iff `o -> x`.

**Lemma 1 (domination ray).** `d(G_t) = rho_dom + 2^t delta_dom` for `t >= 1`, where

    delta_dom = sum over nonempty S in O of 2^(number of vertices outside O beaten by some vertex of S),
    rho_dom   = d(G_0) - iota - delta_dom + 1,

`G_0` the graph without the tail and `iota = 1` if every vertex of `O` beats a vertex outside `O`, else 0.

*Proof.* Every tail vertex is adjacent to all of `O`. Split a dominating set `S` by its tail part. If it contains
the whole tail, `O` is dominated and the vertices of `B - O` not in `S` must be beaten by `S ∩ O`:
`D_1 = sum over all S_O in O of 2^|N(S_O)|` sets. If the tail part is proper and nonempty (`2^t - 2` ways), the
missing tail vertices force `S_O` nonempty: `D_1 - 1`. If it is empty: the dominating sets of `G_0` with `S_O`
nonempty, `d(G_0) - iota`. Sum. ∎

**Corollary 2.** By the inequality of Session261 (`H(T) - H(T - z) >= d(G_z) - 1`) applied to `T_t` and `z = v`,

    (r - h - rho_dom + 1) + 2^t (M - delta_dom) >= 0   for all t >= 1,   hence   M >= delta_dom.

For one selected vertex this is `U_B(o) >= 2^(outdeg of o)`. FINITE-EXACT: `M = delta_dom` holds exactly for
transitive `B` (S5).

Reading. The skeleton `(rho_dom, delta_dom)` sees only which outside vertices each selected vertex beats. Everything
that depends on the paths of `B` sits in the *defect state* `(r - h - rho_dom + 1, M - delta_dom)`, which is
nonnegative as a ray. In a tower the skeleton slope of a level is `2 alpha 2^m`, the same for every base; the
recursion of THM-5130 (C) is entirely a recursion for defects.

## 3. Initial cases

The zero-defect interfaces are the one-apex ones. In state form (S7), for `O` in `TT_n`:

    r = 1 - w(O),   M = sum over o in O of 2^(n-o),   E = [n in O],   w(O) = sum over pairs of O of max(2, 2^(distance-1)).

So the initial supply is indexed by the positive integers `m` (binary digits = `O`): slope `m 2^j`, offset
`1 - w(m)`, and the one-apex spectrum is `{1 - w(m) + 2^t m}`. Their offsets are 1 (at most one selected vertex),
-1 (one pair at distance 1 or 2), -3, -5, -7, -11, -13, ... (`w` = 10 does not occur for at most four selected
vertices); the sign flip of THM-5130 (B) sends each to its negative with an explicit slope. This is the precise sense in which seeds have initial cases: the equality cases of the domination
inequality are the start states of the towers, and every other interface pays a defect of at least `H(B) - 1`.

## 4. Sources from seeds, and two grids

**Proposition 3 (one-vertex gadget).** For a pair with state `(r, M, h, E)` and the gadget `(TT_a, {p})`:

    r+ = h+ = r + 2^a M,   M+ = 2^(a-p) (r + h) + 2^p 3^(a-p) M - [p = 1] E,   E+ = [p = a] (r - h) + 2^(p-1) M.

(THM-5130 (A) with `w = w3 = 0`.) For a seed of slope `2^k` this is a source with Fermat-type offset
`D = 1 + 2^(a+k)`. The sources `(5, 9)`, `(17, 18)`, `(33, 36)` of Session261's 9-grid are the cases
(point, `a = 2, p = 1`), (point, `a = 4, p = 4`), (seed `h = 5, M = 16` of order 4, `a = p = 1`) (S2).

**Proposition 4 (excess types).** If `H(B) = 1 + 2^s` and `U_B(o) = d 2^v` (a one-vertex source of the `d`-grid at
level `s`), then `d 2^v >= 2^s + 1`, so `v >= s - floor(log2 d)`, and the excess
`U_B(o) - H(B) = sum over out-neighbours x of o of Q_B(x, o)` is `(d - 2^i) 2^(s-i) - 1` with `i = s - v`.
(From THM-5130 (F1).) For `d = 9` the excess is `2^(s-3) - 1`, `5 * 2^(s-2) - 1`, `7 * 2^(s-1) - 1` or `2^(s+3) - 1`.

So a uniform supply of sources for the 9-grid is a family, one for each `s`, of tournaments with `2^s + 1` paths and
a vertex with one of four prescribed excess counts. The least-order one-vertex sources have the first type at
`s` = 3, 4, 5, 8, 9, 11 (excess 0, 1, 3, 31, 63, 255), the second at `s` = 2, 7 (excess 4, 159) and the third at
`s` = 6 (excess 223).

**9-grid (FINITE-EXACT).** Census of every selected set over every tournament of order at most 9 (`src10.c`,
`runs/src9grid.txt`): sources exist at every level `2 <= s <= 11`. New levels (S3):

| s | base (order, code) | O | offset | slope | H(B) |
|---|---|---|---|---|---|
| 7 | 8, 11026436 | {2} | 129 | 288 | 129 |
| 8 | 8, 6562209 | {2} | 257 | 288 | 257 |
| 9 | 9, 442533377 | {0} | 513 | 576 | 513 |
| 10 | 9, 6793626886 | {1, 2} | 1025 | 4608 | 1565 |
| 11 | 9, 194216770 | {0} | 2049 | 2304 | 2049 |

So `1 + 2^s + 9 * 2^(r+s)` is a path count for all `r >= 1` and `s <= 11`. Level 10 has no one-vertex source of
order at most 9; it needs a two-vertex selected set. No uniform family: among the unit roots of order at most 9 no
port has a slope family `2^i + 9 * 2^j * 2^a` (the constants found are `d` = 3, 21, 35 and 43), and the gadget route
through seeds with `H(B) = 2^e - 1` stops at `e = 4` at these orders (pre-registered G1 failed, G3 confirmed).

**43-grid (complete).** The unit root of order 9 with code 45174751748, root 2, has a port (vertex 8, in-degree 1,
not selected) with slope family `16 + 172 * 2^a`; a module at such a port is visited at most twice, so `H`, `W`, `U`
are affine in `2^a` and the unit offset at `a = 1, 2` holds for every `a`. This gives all `s >= 5`, `r >= 1`. Sources
at levels 2, 3, 4 give all `r >= 1`, a source at level 1 gives `r >= 2`, and `N(43, 1, 1) = 175` is a path count of
order 7 (S8). Hence `1 + 2^s + 43 * 2^(r+s)` is a path count for all `r, s >= 1`.

`runs/grid_levels_singletons.txt` lists, for every odd `d <= 129`, the levels with a one-vertex source of order at
most 9.

## 5. The constructive task in interface form

1. **Grid normal form (KNOWN, seed_recursion_births).** Every odd `N` other than `1 + 2^s` is `1 + 2^s + d 2^(r+s)`
   for one triple `(d, r, s)`. Chan-Pak says: every grid point is a path count except `(1, 1, 1)` and `(1, 2, 2)`.
2. **What one grid needs.** One carrier (a unit pair with a port of in-degree 1 whose slope family is
   `2^i + d 2^j 2^a`) gives all `s > i`; then `i` sources (one per level, each all `r`), finitely many seeds or heads
   for small `r`. The grids 3, 21, 35 (other sessions) and 43 (here) are finished exactly this way. So the supply
   problem is: *for every odd `d`, one carrier and finitely many sources.*
3. **Where the domination inequality bites.**
   - It grades interfaces. Grade 0 (zero defect) is one-apex, and its slopes are sums of powers of two,
     `delta_dom`. Whatever makes a slope non-dyadic is defect. A source for an odd `d` that is not a sum of few
     powers of two must take `d 2^v - delta_dom` from the defect.
   - It bounds slopes from below (`M >= delta_dom`, `U >= H`), which is what leaves only `floor(log2 d) + 1` excess
     types per level (Proposition 4) and what forbids seeds of slope 7 or anti-seeds of slope 8 (THM-5130, the
     remark on 15).
   - It gives a budget for heads: for every vertex `z` of a realisation, `N >= d(G_z) + H(T - z) - 1`, and along any
     order of deleting vertices `N - 1 >= sum of (d_i - 1)`. A small target leaves a finite list of skeletons.
4. **Composition.** The requirement that a composition retain `Q_2` as well as `H` is met in one direction by the
   gadget law: the four numbers `(r, M, h, E)` are closed under insertion of a new closer with re-selection, for
   every inserted tournament `X`; no further statistic of the base is needed. What is not closed is the use of a
   closed tournament as a module elsewhere: that needs its `Pi_2`, which depends on three-path data of the base.
5. **A proposal.** Split the search for a carrier of a given `d` into a skeleton problem and a defect problem: choose
   the out-neighbourhood pattern of the selected set so that `delta_dom` is the dyadic part of the wanted family, and
   look for bases whose defect state is the remaining part. The towers show that defect states of a fixed shape
   propagate linearly; the open step is a base whose defect is prescribed. Test for a later session: for
   `d = 9`, enumerate the pairs of order 10 whose defect slope family is `9 * 2^j 2^a + 2^i - delta_dom`.

## 6. The Euler-function criterion

The statement pasted by the owner is Nicolas's theorem (2012); `n/phi(n)` is extremal at primorials, and Robin's
criterion is its analogue for `sigma(n)/n`. Typed bridges to the repo:

- **[PROVED] An abundancy index for tournaments.** Put `a(X) = Pi_2(X)/H(X)`. Then `a(X => Y) = a(X) a(Y)`;
  `a(X) = sum over vertex sets S of H(X[S]) H(X - S)/H(X)` with every term at most 1 (two disjoint paths of a
  tournament merge, keeping both orders, into one path; the merge is injective); `n + 1 <= a(X) <= 2^n`, with
  `a = 2^n` only for transitive `X` (S6; the lower bound is Session197's). The closing law is
  `H(B => X => v) = H(X) (r + M a(X))`, and the slope of a padded seed is `M_core` times the product of `a(C)` over
  the components after the core (THM-5130 (F2)): an Euler product whose local factor at a single vertex is 2.
  Dictionary: strong components ↔ primes; `a(X)` ↔ `sigma(n)/n` (a divisor sum with terms at most 1) and `n/phi(n)`
  (an Euler product of local factors); transitive tournaments ↔ primorials (the extremal chain); rays ↔ the values
  of the abundancy law on the extremal chain.
- **[CITED, in canon]** THM-878 already uses the unconditional Rosser-Schoenfeld bound
  `phi(q) > q/(e^gamma log log q + 3/log log q)`; the pasted inequality is its sharp form under RH, with the
  second-order constant `e^gamma (2 + beta)`, `beta = sum over zeros of 1/(rho (1 - rho))`.
- **[CITED, in canon]** THM-4666: for Paley tournaments the ranking excess is at least
  `(e^gamma/pi + o(1)) p^(3/2) log log p` for suitable primes, through `L(1, chi_-p)`. Same constant, same mechanism:
  `L(1, chi) = prod (1 - chi(q)/q)^(-1)` is extremal when all small primes are residues, exactly as `n/phi(n)` is
  extremal when all small primes divide `n`; the opposite extreme (all small primes non-residues) is the Heegner
  side (163; opus-S90, THM-5125).
- **[ANALOGY with exact halves]** `n/phi(n) = sum over squarefree d | n of 1/phi(d)` is a free gas of the primes of
  `n`; `H = I(Omega, 2)` is the hard-core gas of odd cycles (THM-510, THM-517); `d(G_z)` is a sum of
  `2^(components)` over admissible sets (Session261). The domination-defect inequality, like Nicolas's, is one-sided,
  exact on an extremal chain (one-apex ↔ primorials), and silent about the size of the defect, which is where the
  arithmetic is.
- **[STRUCTURE]** The threshold `p_120569#` plays the role of our ray thresholds: a uniform statement beyond it,
  finitely many initial cases checked apart. Chan-Pak has the same architecture, with two initial cases that fail.
- **[NUMEROLOGY]** T234: `prod over p <= 19 of (p+1)/(p-1) = 21`, which is `(n/phi(n)) (sigma(n)/n)` at `n = 19#`:
  the product of the Nicolas and Robin ratios at that primorial; the partial products are integers only at 2, 3, 5,
  7 and 19. No mechanism.
- **[ANALOGY]** Session260: a seed of odd slope has `H = 2M + 1` (quasiperfect shape), an anti-seed `2M - 1`
  (almost perfect shape).

## 7. The two papers

- **Chen-Mogielnicki-Ono (record compositions, noncommutative symmetric functions).** [ANALOGY] Their count is a
  product over blocks of `binomial(2 s_j - 1, 2 l_j - 1) E_(2 l_j - 1)`, blocks cut at records. The constants of a
  transitive gadget are counts of the same kind: `u3(j) = 2^(j-1) 3^(a-j) - 2^(a-j+1) + [j = 1]` counts ordered
  triples of blocks whose marked block has least element `j`, and a tower slope is a product over levels of such
  factors. [QUESTION] The state `(r, M, h, E)` is a pointed, ordered refinement of the Redei-Berge function at
  `(1)`, `(1, 1)`; the natural home of ordered path covers with marked ends is noncommutative. Is there a
  noncommutative lift of the odd-cycle formula whose coefficients in products of power sums of the first kind count
  packings with a record structure, as their Theorem 1.3 does for alternating permutations?
- **Sulyma (RO(G)-graded norms, de Rham-Witt).** [ANALOGY] A ray `r + 2^t M` is a section over the Rees algebra of
  the 2-adic filtration, and multi-module families with bases 2 and 3 are multi-Rees (two filtrations: two-block
  and three-block covers). The doubling `K -> K[TT_2 at y]` behaves like a norm for the group of order 2: the offset
  of the twins is `H(K_2) - 2 H(K)`, norm minus transfer, and the two unit offsets `+1`, `-1` are exchanged by the
  pair gadget as the two real characters are by the sign representation. [QUESTION] `Z_2^* = {+1, -1} x (1 + 4 Z_2)`:
  the two kinds of line are the two cosets and `s` is the depth in the unit filtration. Is the sign flip the
  shadow of a genuine grading of the family of towers by the representation ring of the group of order 2?

## 8. Pre-registered checks E4 (PREREGISTRATION.md)

G1 FAILED (the pairs with `(h, M) = (15, 64)` have `E = 6`, slope 138; no level-7 source from the gadget on these
seeds); G2 CONFIRMED (levels 7 and 8, and also 9 and 11, by one-vertex sources of order at most 9); G3 CONFIRMED
(no pair of order at most 9 with `H(B) = 31` and slope 64 or 128); G4 CONFIRMED.
