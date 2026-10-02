# Edge multiset dimension of hypercubes III: Open Problem 3, one estimate from d = 10

**Provenance:** Claude thread session (project thread, branch `claude/project-thread-96889h`), 2026-10-02;
temporary identity claude-96889h, no `.machine-id`. Independent audit owed.
**Source problem:** J. Allikvere, *The edge multiset dimension of hypercubes*, arXiv:2608.09983v1, Section 10,
Open Problem 3, as recorded in
[THM-4534](../../01-canon/theorems/THM-4534-edge-multiset-dimension-grows-like-exp-theta-cube-root-d.md) and the lane
note [`procgen_edim2_20261001_growth_and_uniform_bounds.md`](procgen_edim2_20261001_growth_and_uniform_bounds.md)
(section 4). The paper proves that a uniformly random vertex set of Q_d resolves with positive probability for every
d >= 11. For 11 <= d <= 50 it does this with a separate computer certificate of its union bound U_d at each d; for
d >= 51 it uses an analytic tail estimate. Open Problem 3 asks for one analytic estimate valid from d = 11. THM-4534
lists the strict form as OPEN and proves a closed-form estimate for every d >= 17 (its Theorem D).
**Builds on:** the forest lemma (the paper's Lemma 11), the pair types and cells, and the antipodal halving of
[`edge_multiset_dimension_growth_20261002.md`](edge_multiset_dimension_growth_20261002.md) (sections 1.1, 1.2, 1.4,
1.6), which are also in THM-4534 and
[THM-4525](../../01-canon/theorems/THM-4525-edge-multiset-dimension-of-q6-is-15.md) (L2).
**Runner:** `04-computation/experiments/edge_multiset_dimension_op3_20261002_run.py` (pure Python, exact integer and
rational arithmetic, 9 s); stdout in `05-knowledge/results/edge_multiset_dimension_op3_20261002.out`: 573,465
checks, ALL CHECKS PASSED.
**Labels used:** PROVED, FINITE-EXACT.

## Status header

| # | Claim | Label |
|---|-------|-------|
| 1 | **Open Problem 3 is settled in its strict form, and from d = 10.** For every d >= 10, a uniformly random subset of V(Q_d) is edge-multiset resolving with probability at least 0.32. The proof is one estimate: the recursion V_(d+1) <= V_d/2 + A_(d+1) for the paper's own union bound (with antipodal halving), started from the single exact value V_10 = 0.66023. No per-d certificate is used for d >= 11. | PROVED (computer-assisted in two finite places: the exact rational V_10, and a rational check of the ratio inequality for 10 <= d <= 40; analytic beyond d = 40) |
| 2 | **Transfer lemma (Theorem 1).** For every pair type t' of Q_(d+1) other than the antipodal one there is a pair type t of Q_d, and an explicit near-central cell size N*(t'), with B(t') <= B(t) beta(N*(t')). Here B is the potential-forest bound of a pair type. | PROVED by hand; also machine-checked sub-cell by sub-cell for every type, 3 <= d <= 40 (41,116 transferred forest edges) |
| 3 | **Ratio bound (Lemma 4).** rho_d <= 1/2 for every d >= 10 (0.4959 at d = 10, 0.0875 at d = 16, 0.00006 at d = 40). | PROVED (rational upper bounds for 10 <= d <= 40; a closed-form bound below 0.001 for d >= 41) |
| 4 | With antipodal halving, the paper's union bound is already below 1 at d = 10: V_10 = 0.66023 (1.30745 without halving, the value the lane note quotes). For 6 <= d <= 12, the potential forests reach the maximum-weight-forest value exactly. | FINITE-EXACT |
| 5 | V_d -> 0, and V_d <= 0.6767 for every d >= 10. | PROVED |

## 1. Setting

All of this is at density q = 1/2. S is a uniformly random subset of V(Q_d).

* **Pair types.** Unordered pairs {e, f} of distinct edges fall into the types P(h) (parallel; L = d - 1 remaining
  coordinates, h of them separate the base points, nu = L - h; 1 <= h <= L) and X(h) (crossing; L = d - 2,
  nu = L - h, 0 <= h <= L). Write a type as (kind, h, nu). P(h) contains d 2^(d-2) C(d-1, h) pairs, X(h) contains
  d(d-1) 2^(d-1) C(d-2, h) pairs (growth note 1.1; the counts add up to C(E, 2)).
* **Sub-cells and the level graph.** For P(h) a vertex s has coordinates (c, j): c disagreements on the nu common
  coordinates, j on the h separating ones. The sub-cell (c, j) has 2 C(nu, c) C(h, j) vertices and sends them from
  level c + j (the value d(e, s)) to level c + h - j (the value d(f, s)). For X(h) the sub-cell (x, y, c, j), with
  two bits x, y, has C(nu, c) C(h, j) vertices and goes from level x + c + j to level y + c + h - j. The levels run
  from 0 to T, where T = h + nu for P and T = h + nu + 1 for X. The *weight* N_ab of a pair of levels a != b is the
  total size of the sub-cells joining a and b, in either direction. Diagonal sub-cells (equal levels) play no role.
  In closed form, for P(h): N_ab = 4 C(nu, (a+b-h)/2) C(h, (a-b+h)/2) when a + b - h is even and both binomials are
  defined, and 0 otherwise.
* **Forest lemma (the paper's Lemma 11; growth note 1.2 at q = 1/2).** For every forest F on the levels,
  P(H_e = H_f) <= prod over ab in F of beta(N_ab), where beta(N) = C(N, floor(N/2))/2^N. The reason: the flow through
  a forest edge is X - Y with X ~ Bin(M_1, 1/2), Y ~ Bin(M_2, 1/2), M_1 + M_2 = N_ab, and the largest atom of X - Y
  is beta(N_ab).
* **beta.** It is non-increasing, since beta(2m-1) = beta(2m) >= beta(2m+1). Also beta(N) <= sqrt(2/(pi N))
  (C(2m, m) <= 4^m/sqrt(pi m)).
* **Antipodal halving (growth note 1.4; THM-4525 L2).** H_(e-bar) is H_e reversed. So {e, f} and {e-bar, f-bar}
  collide together, and the involution fixes exactly the antipodal pairs, which form the type P(d-1) (nu = 0). In
  the union bound every type gets weight w = 1/2, except P(d-1), which gets weight 1.

## 2. Potential forests and the bound V_d

**Centrality.** On the levels 0..T of a type, b is *more central* than a, written b < a, when |2b - T| < |2a - T|,
or when |2b - T| = |2a - T| and b > a. (Of two mirror levels, the upper one counts as more central.) This is a strict
total order.

**Potential forest.** Every level a that has a neighbour b < a with N_ab > 0 *chooses* one such neighbour of largest
weight; write N(a) for that weight. Levels with no such neighbour are *roots*.

**Lemma 1 (PROVED).** The chosen edges form a forest. *Proof.* Orient each chosen edge away from its chooser. Two
levels cannot choose the same edge, since each would be more central than the other. Every level has out-degree at
most 1, and centrality strictly increases along an oriented edge. A cycle of chosen edges, with as many edges as
levels, would have to be a directed cycle, which is impossible. QED

The forest is the same idea as the toward-the-centre forest of the growth note (Lemma 3 there), now with the
heaviest edge at every level and with ties broken. Define
* B(t) = prod over the choosers a of t of beta(N(a)),
* V_d = sum over the types t of Q_d of w(t) |t| B(t), the union bound with halving,
* A_d = w |t| B(t) for the antipodal type t = P(d-1): A_d = d 2^(d-2) prod over 0 <= j < (d-1)/2 of beta(4 C(d-1, j)).
  (Its level graph is the matching {j, d-1-j}, each lower level chooses its mirror, so B is the exact collision
  probability.)

By the forest lemma and the halving, **P(S is not resolving) <= V_d.**

**FINITE-EXACT (runner S2).** V_6, ..., V_12 = 41.80, 29.21, 13.00, 3.843, **0.66023**, 0.08596, 0.005561. These
equal the maximum-weight-forest (Kruskal) values with halving, exactly, for every 6 <= d <= 12. Without halving
the Kruskal values are 1.30745 (d = 10), 0.15618 (d = 11) and 0.01075 (d = 12); these are the values the lane note
quotes for the paper's bound. The weights and counts agree with the independent `cells` and `type_count` of
`procgen_edim2_20261001_lib.py` for 6 <= d <= 12, and with brute force over representative edge pairs for
3 <= d <= 8 (runner S1).

## 3. The transfer lemma

Fix d >= 3 and a type t' = (kind, h', nu') of Q_(d+1) other than the antipodal type P(d). Its *predecessor* t is a
type of Q_d:
* kind P: always the **nu-step** t = (P, h', nu' - 1) (here nu' >= 1);
* kind X: the **h-step** t = (X, h' - 1, nu') if h' >= nu' and h' >= 1; otherwise the nu-step t = (X, h', nu' - 1).

A nu-step keeps h and adds one common coordinate. An h-step keeps nu and adds one separating coordinate. Either
way the top level rises from T to T + 1.

**Level map.** For a level a of t with 2a != T put phi(a) = a if 2a < T (a *lower* level) and phi(a) = a + 1 if
2a > T (an *upper* level). For a chooser a with chosen neighbour b put
* psi(b) = b (a lower) or b + 1 (a upper) for a nu-step;
* psi(b) = b + 1 (a lower) or b (a upper) for an h-step.
The most central level of t is never a chooser, so phi and psi are defined wherever they are used.

**The extra chooser.** Define the level a* of t' and its target b*:

| case | condition | a* | b* | N*(t') = N'_(a* b*) |
|---|---|---|---|---|
| X, or P with h odd | T even | T/2 | T/2 + 1 | see below |
| X, or P with h odd | T odd | (T+3)/2 | (T+1)/2 | see below |
| P, h even (B1) | T even | T/2 + 2 | T/2 | see below |
| P, h even (B2) | T odd | (T-1)/2 | (T+3)/2 | see below |

with, in terms of t' = (kind, h', nu'):
* kind P: **N*(t') = 4 C(nu', floor(nu'/2)) C(h', ceil(h'/2) - 1)**;
* kind X, h' odd: **N*(t') = 2 C(nu'+1, floor((nu'+1)/2)) C(h', (h'-1)/2)**;
* kind X, h' even: **N*(t') = 2 C(nu', floor(nu'/2)) C(h'+1, h'/2)**.

**Theorem 1 (transfer lemma, PROVED).** B(t') <= B(t) beta(N*(t')).

*Proof.* Write N for the weights of t and N' for those of t'.

*Step 1: the image of a chosen edge is at least as heavy.* Let a be a chooser of t with chosen neighbour b. Map
each sub-cell of t that joins a and b (in either direction) to a sub-cell of t' as follows. "From a" means the
sub-cell goes from level a to level b; "from b" means the reverse.

| step, side of a | sub-cell from a | sub-cell from b |
|---|---|---|
| nu-step, a lower | unchanged | unchanged |
| nu-step, a upper | c -> c + 1 | c -> c + 1 |
| h-step, a lower | unchanged | j -> j + 1 |
| h-step, a upper | j -> j + 1 | unchanged |

(The bits x, y of a crossing sub-cell are kept.) A direct substitution into the level formulas (c + j, c + h - j),
resp. (x + c + j, y + c + h - j), shows that the image goes from phi(a) to psi(b), resp. from psi(b) to phi(a). For
example, in an h-step with a lower, a sub-cell from b has b = c + j and a = c + h - j. Its image (c, j + 1) in t' goes
from c + j + 1 = b + 1 = psi(b) to c + (h + 1) - (j + 1) = a = phi(a). The image is never smaller, by Pascal's rule:
its size has C(nu + 1, c) or C(nu + 1, c + 1) in place of C(nu, c), or C(h + 1, j) or C(h + 1, j + 1) in place of
C(h, j). The map is injective, because within one direction it is a translation, and the two directions land on
sub-cells that start at different levels (phi(a) != psi(b) by Step 2). Hence N'(phi(a) psi(b)) >= N_ab.

*Step 2: psi(b) is more central than phi(a) in t'.* Put delta = 2a - T and gamma = 2b - T, so that b < a means
|gamma| < |delta|, or |gamma| = |delta| with gamma > 0 > delta. In t' the top is T + 1.
* a lower (delta < 0): 2 phi(a) - (T+1) = delta - 1, of absolute value |delta| + 1. For a nu-step,
  2 psi(b) - (T+1) = gamma - 1. If |gamma| < |delta|, then |gamma - 1| < |delta| + 1. If gamma = -delta > 0, then
  |gamma - 1| = |delta| - 1. For an h-step, 2 psi(b) - (T+1) = gamma + 1. If |gamma| < |delta|, then
  |gamma + 1| < |delta| + 1. If gamma = -delta > 0, then gamma + 1 = |delta| + 1: a tie, but psi(b) is the upper level
  of the two, so it is more central.
* a upper (delta > 0): a tie |gamma| = delta is impossible, because a would win it. So |gamma| < delta, and
  2 phi(a) - (T+1) = delta + 1, while 2 psi(b) - (T+1) = gamma + 1 (nu-step) or gamma - 1 (h-step). Both have
  absolute value at most |gamma| + 1 < delta + 1.

*Step 3: the extra chooser.* phi is injective, and a* is not phi(a) for any chooser a of t:
* T even, X or P with h odd: a* = T/2 is not in the image of phi at all.
* T odd, X or P with h odd: a* = phi((T+1)/2), and (T+1)/2 is the most central level of t, a root.
* B1 (P, h even, T even): a* = phi(T/2 + 1). Since h is even, every edge of t joins levels of the same parity.
  The only level more central than T/2 + 1 is T/2, of the other parity, so T/2 + 1 is a root.
* B2 (P, h even, T odd): a* = phi((T-1)/2). The only level more central than (T-1)/2 is (T+1)/2, of the other
  parity, so (T-1)/2 is a root.
b* is more central than a* in t' (top T + 1): in the four rows the pairs (|2a* - T - 1|, |2b* - T - 1|) are (1, 1)
with b* upper, (2, 0), (3, 1) and (2, 2) with b* upper. Finally, N'(a* b*) = N*(t'):
* P: substitute into N'_ab = 4 C(nu', c) C(h', j) with c = (a+b-h')/2, j = (a-b+h')/2, h' = h and nu' = nu + 1.
  - Row A1 (h odd, T even, so nu' is even): c = (T+1-h)/2 = nu'/2 and j = (h-1)/2.
  - Row A2 (h odd, T odd, nu' odd): c = (T+2-h)/2 = (nu'+1)/2 and j = (h+1)/2.
  - Row B1 (h even, T even, nu' odd): c = (T+2-h)/2 = (nu'+1)/2 and j = h/2 + 1.
  - Row B2 (h even, T odd, nu' even): c = (T+1-h)/2 = nu'/2 and j = h/2 - 1.
  In every row C(nu', c) = C(nu', floor(nu'/2)) and C(h', j) = C(h', ceil(h'/2) - 1), by the symmetry C(n, k) = C(n, n-k).
* X: write {a*, b*} = {m, m + 1}, the two most central adjacent levels of t' (m = (T'-1)/2 for T' = T + 1 odd and
  m = T'/2 for T' even). For h' odd only the sub-cells with x = y join m and m + 1: j = (h' - 1)/2 in one direction
  and j = (h' + 1)/2 in the other, with x = 0, 1. Pascal's rule merges the two values of x, and the total is
  2 C(nu' + 1, m - (h'-1)/2) C(h', (h'-1)/2). For h' even only x != y contributes: j = h'/2 and h'/2 - 1 in one
  direction, h'/2 + 1 and h'/2 in the other. The total is C(nu', m - h'/2) (2 C(h', h'/2) + C(h', h'/2 - 1) +
  C(h', h'/2 + 1)) = 2 C(nu', m - h'/2) C(h' + 1, h'/2). In both cases the index of the nu'-binomial is central
  (floor((nu'+1)/2), resp. floor(nu'/2)), which gives the table.

*Step 4: conclusion.* Every phi(a) has the neighbour psi(b), which is more central (Step 2), so phi(a) is a chooser
of t' with N'(phi(a)) >= N'(phi(a) psi(b)) >= N(a) (Step 1). The level a* is a further chooser, with
N'(a*) >= N*(t') (Step 3). The other levels of t' contribute factors at most 1, and beta is non-increasing, so
B(t') <= prod over choosers a of t of beta(N(a)), times beta(N*(t')), which equals B(t) beta(N*(t')). QED

(Runner S3: for every type and 3 <= d <= 40, every one of the 41,116 chosen edges is transferred sub-cell by
sub-cell, and the maps, sizes, targets, extra choosers and N*(t') are checked as stated. For d <= 13 the
inequality B(t') <= B(t) beta(N*(t')) is also checked in exact arithmetic.)

*Where the gain comes from.* Moving to Q_(d+1) multiplies the number of pairs of a type by about 2(d+1)/nu' (or
2(d+1)/h'). The old cells only grow, but the gain that pays for this is one new forest edge at the centre, with
N*(t') of order 2^d/d. The weight of the antipodal type halves under the nu-step that leaves it, and the
crossing types use whichever step has the larger denominator.

## 4. The recursion

**Corollary 2 (PROVED).** V_(d+1) <= rho_d V_d + A_(d+1), where
rho_d = max over types t of Q_d of the sum, over the types t' of Q_(d+1) whose predecessor is t, of
r(t') = (w(t') |t'|)/(w(t) |t|) beta(N*(t')).

*Proof.* Split off the antipodal type of Q_(d+1). By Theorem 1,
w(t')|t'| B(t') <= r(t') w(t)|t| B(t) for every other t'. Sum, grouping the t' by their predecessor. QED

The count ratios are |t'|/|t| = 2(d+1)/nu' for a nu-step and 2(d+1)/h' for an h-step (for both kinds), and the
weight ratio is 1, except 1/2 when t is antipodal. A parallel type has one successor. A crossing type (h, nu) has
at most two: (h, nu+1) when nu >= h, and (h+1, nu) when h + 1 >= nu.

**Lemma 3 (binomial bounds).** C(m, floor(m/2)) >= 2^m/(m+1) for m >= 0 (it is the largest of m + 1 terms summing to
2^m), and C(h, ceil(h/2) - 1) >= 2^h/(2(h+1)) for h >= 1 (for odd h it is central; for even h it is
C(h, h/2) (h/2)/(h/2 + 1) >= C(h, h/2)/2).

**Lemma 4 (ratio bound).** rho_d <= 1/2 for every d >= 10.
* *10 <= d <= 40 (FINITE-EXACT, runner S4).* Every r(t') is bounded above by a rational number: beta(N) itself
  for N <= 4000, else an integer-square-root upper bound for sqrt(2/(3.14159 N)). The loads are summed exactly.
  The maxima are

  | d | 10 | 11 | 12 | 13 | 14 | 15 | 16 | 24 | 32 | 40 |
  |---|---|---|---|---|---|---|---|---|---|---|
  | rho_d <= | 0.49587 | 0.35085 | 0.26127 | 0.19845 | 0.15036 | 0.11648 | 0.08751 | 0.00877 | 0.00077 | 0.00006 |

  For d = 10..12 the largest load sits on a central crossing type ((X, 4, 4), (X, 4, 5), (X, 5, 5)). From d = 13 on
  it sits on the parallel types next to the antipodal one (nu = 0 or 1). For comparison, rho_8 <= 0.923 and
  rho_9 <= 0.654, while rho_6, rho_7 > 1.
* *d >= 41 (PROVED).* Parallel: h' + nu' = d, and N*(t') >= 4 (2^(nu')/(nu'+1)) (2^(h')/(2(h'+1)))
  = 2^(d+1)/((nu'+1)(h'+1)) by Lemma 3. With beta(N) <= sqrt(2/(pi N)) and (nu'+1)(h'+1)/nu'^2 <= 2(d+1),
  r(t') <= (2(d+1)/nu') sqrt((nu'+1)(h'+1)/(pi 2^d)) <= 2 (d+1)^(3/2) sqrt(2/pi) 2^(-d/2).
  Crossing: h' + nu' = d - 1, the denominator is max(h', nu') >= (d-1)/2, (nu'+1)(h'+1) <= ((d+1)/2)^2 and
  N*(t') >= 2 C(nu', floor(nu'/2)) C(h', floor(h'/2)) >= 2^d/((nu'+1)(h'+1)). So each r(t') is at most
  (2(d+1)^2/(d-1)) sqrt(2/pi) 2^(-d/2), and a load (two successors at most) is at most twice that. Both bounds
  decrease for d >= 41, where their larger value is 2.9e-4. (The runner also checks Lemma 3 for m <= 300, and that
  this crude bound dominates the exact one for 41 <= d <= 64.) QED

## 5. Antipodal terms

A_d = d 2^(d-2) prod over 0 <= j < (d-1)/2 of beta(4 C(d-1, j)).
* *11 <= d <= 60 (FINITE-EXACT, runner S5).* Rational upper bounds as in Lemma 4: A_11 <= 0.01575,
  A_12 <= 3.75e-4, A_13 <= 3.84e-4, A_14 <= 4.0e-6, and the sum over 11 <= d <= 60 is at most 0.016518.
* *d >= 61 (PROVED).* beta(4) = 3/8, and beta(4 C(d-1, j)) <= (2 pi (d-1))^(-1/2) for 1 <= j <= d-2. There are at
  least (d-3)/2 factors with j >= 1, so log2 A_d <= log2 d + d - 2 + log2(3/8) - ((d-3)/4) log2(2 pi (d-1)).
  For d >= 61, log2(2 pi (d-1)) >= 8.55, so log2 A_d + d/2 <= log2 d - 0.637 d + 3.0. This is negative at d = 61
  and decreasing beyond, so A_d <= 2^(-d/2), and the tail sum is at most 2.2e-9.

## 6. Open Problem 3

**Theorem 5 (PROVED; Open Problem 3).** For every d >= 10,
P(a uniformly random subset of V(Q_d) is not edge-multiset resolving) <= V_d <= V_10 + sum over k >= 11 of A_k
<= 0.66023 + 0.01652 < 0.6767.
Moreover V_(d+1) <= V_d/2 + A_(d+1), so V_d -> 0.

*Proof.* Section 2 gives the first inequality. Corollary 2 and Lemma 4 give V_(d+1) <= V_d/2 + A_(d+1) <= V_d + A_(d+1)
for d >= 10. Induct from V_10 (runner S2: an exact fraction below 0.6603) and use section 5. QED

So existence of resolving sets for all d >= 10 follows from one estimate: the recursion, plus one base value. The
paper's union bound needs an evaluation at a single d only, and that d can be 10 instead of 11 thanks to the
halving. Together with explicit resolving sets for 6 <= d <= 9 (THM-4525 and the first lane note
[`procgen_edim_20261001_edge_multiset_dimension.md`](procgen_edim_20261001_edge_multiset_dimension.md)), this gives
edim_m(Q_d) < infinity for every d >= 6.

**How this answers the question as posed.** THM-4534 and its lane note record OP3 as "one analytic estimate valid from
d = 11", and the lane note suggested exactly this route (section 4.3 there): a monotonicity lemma V_(d+1) <= V_d for
the potential-forest bound would reduce the problem to one value. Theorem 1 and Lemma 4 prove a quantitative form
of that lemma (with the antipodal term separated), for d >= 10. If one insists on no computer evaluation at all, two
finite inputs remain: the base value V_10 (an exact rational, one d) and the ratio check for 10 <= d <= 40. The
paper's own tail estimate also rests on a single evaluation, at d = 51.

## 7. Relation to other work

* THM-4534 (audited canon; procgen_edim2 lane): Theorem D proves U_d < 1 in closed form for every d >= 17. Its
  Fourier-Hoelder multi-forest lemma certifies U_10 <= 0.2649 and U_11 <= 0.0448, computed at each d. Its OPEN item
  "OP3 in the strict form" is what Theorem 5 settles. Proposition K there (weighted spanning trees of the parallel
  level graph) is not used here.
* Growth note ([`edge_multiset_dimension_growth_20261002.md`](edge_multiset_dimension_growth_20261002.md), section 6)
  first recorded OP3 as OPEN and asked for "a lifting lemma from Q_d to Q_(d+1)". Theorem 1 is such a lemma, for the
  union bound rather than for resolving sets; that section now points here. Its own sparse-regime estimate (Theorem A
  there) is effective from d = 41.

## 8. Not done

* A base below d = 10. V_9 = 3.843 > 1, although rho_9 <= 0.654. The true expected number of colliding pairs at
  d = 9 is about 0.8 (EMPIRICAL, lane note 4.3), so a first-moment argument at d = 9 would need much more than
  forests.
* The Fourier-Hoelder lemma of THM-4534 was not combined with the transfer. Doing so would need edge-disjoint extra
  edges, one for each forest.
* An independent audit of this note.

## 9. Reproduction

```
python3 04-computation/experiments/edge_multiset_dimension_op3_20261002_run.py [DMAX_TRANSFER=40]
```
Sections: S1 cells and counts against brute force; S2 the union bound, the Kruskal comparison and V_10; S3 the
transfer lemma, sub-cell by sub-cell; S4 the ratio bound; S5 the antipodal terms and the conclusion. 573,465
checks, about 9 seconds, ALL CHECKS PASSED. Timing goes to stderr; stdout is deterministic.
