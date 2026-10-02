# Edge multiset dimension of hypercubes III: Open Problem 3, a recursion from one base value

**Provenance:** Claude thread session (project thread, branch `claude/project-thread-96889h`), 2026-10-02;
temporary identity claude-96889h, no `.machine-id`.
**Audit and revision:** a blind independent audit (2026-10-02) of the first version (commit ada1e7bc) reproduced
every number and found the mathematics correct (transfer lemma, recursion, ratio bound, antipodal bounds, V_10). It
also found the headline overstated. That version claimed "one estimate, no per-d certificate for d >= 11", but its
induction used an exact ratio check at each 10 <= d <= 40 and numerical A_d for 11 <= d <= 60, and several of its
printed upper bounds were rounded down. This version follows the audit's fix. Theorem 5 now starts from the single
exact value V_11 and uses closed forms, valid for every d >= 3, at every step. Every printed upper bound is rounded
up, and the runner asserts each constant quoted here exactly. The auditor checked the new route of Theorem 5 in its
own interval arithmetic; this revised text has not had a second audit.
**Source problem:** J. Allikvere, *The edge multiset dimension of hypercubes*, arXiv:2608.09983v1, Section 10,
Open Problem 3, as recorded in
[THM-4534](../../01-canon/theorems/THM-4534-edge-multiset-dimension-grows-like-exp-theta-cube-root-d.md) and the lane
note [`procgen_edim2_20261001_growth_and_uniform_bounds.md`](procgen_edim2_20261001_growth_and_uniform_bounds.md)
(section 4). The paper itself could not be fetched here (the container's network policy blocks arxiv.org). The paper
proves that a uniformly random vertex set of Q_d resolves with positive probability for every d >= 11. For
11 <= d <= 50 it uses a separate computer certificate of its union bound U_d at each d. For d >= 51 it uses an
analytic tail estimate, whose constant rests on one evaluation (B_51; lane note section 4). Open Problem 3 asks for
one analytic estimate valid from d = 11. THM-4534 lists this strict form as OPEN. Its Theorem D proves a closed-form
estimate for every d >= 17, whose only numerical inputs are interval evaluations at d = 17, 18 and 19.
**Builds on:** the forest lemma (the paper's Lemma 11; growth note 1.2 at q = 1/2) and the pair types and cells
(growth note 1.1) of [`edge_multiset_dimension_growth_20261002.md`](edge_multiset_dimension_growth_20261002.md); the
antipodal halving (growth note 1.4, which rests on Lemma L2 of
[THM-4525](../../01-canon/theorems/THM-4525-edge-multiset-dimension-of-q6-is-15.md)).
**Runner:** `04-computation/experiments/edge_multiset_dimension_op3_20261002_run.py` (pure Python, exact integer and
rational arithmetic, about 9 s); stdout in `05-knowledge/results/edge_multiset_dimension_op3_20261002.out`:
573,986 checks, ALL CHECKS PASSED.
**Labels used:** PROVED, FINITE-EXACT, OPEN.

## Status header

| # | Claim | Label |
|---|-------|-------|
| 1 | **Open Problem 3: one recursion with closed-form coefficients, from one base value.** For every d >= 11, a uniformly random subset of V(Q_d) is edge-multiset resolving with positive probability. The proof is one estimate for the paper's union bound with antipodal halving, V_(d+1) <= R(d) V_d + g(d+1), whose coefficients R and g are explicit elementary functions, valid for every d >= 3. It starts from the single exact value V_11 < 0.0859648. It gives V_12, ..., V_16 <= 0.4606, 0.7406, 0.8026, 0.6490, 0.4039 and V_d <= 1/2 for every d >= 16 (Theorem 5). The case d = 10 follows from a second exact value, V_10 < 0.660228. | PROVED, computer-assisted at one point: the exact rational V_11 is FINITE-EXACT (runner S2); the rest is closed-form, re-checked in exact arithmetic by the runner (S4-S6). The case d = 10 uses V_10 in the same way |
| 2 | **V_d <= V_10 < 0.660228 for every d >= 10**, so a uniformly random subset of V(Q_d) is resolving with probability more than 0.3397 for every d >= 10 (Corollary 6). Also V_d -> 0. | PROVED, computer-assisted. Its finite inputs are FINITE-EXACT: V_10, the ratio bounds rho_10, ..., rho_15, and A_11, A_12 (runner S2, S4, S5). V_d -> 0 needs only Theorem 5 |
| 3 | **Transfer lemma (Theorem 1).** For every pair type t' of Q_(d+1) other than the antipodal one there is a pair type t of Q_d, and an explicit near-central cell size N*(t'), with B(t') <= B(t) beta(N*(t')). Here B is the potential-forest bound of a pair type. | PROVED by hand; also machine-checked sub-cell by sub-cell for every type, 3 <= d <= 40 (41,116 transferred forest edges) |
| 4 | **Closed forms (Lemmas 4 and 5).** rho_d <= R(d) = max(2 (d+1)^(3/2), 4 (d+1)^2/(d-1)) sqrt(2/pi) 2^(-d/2) and A_d <= g(d) = d 2^(d-2) (3/8) (2 pi (d-1))^(-(d-3)/4), for every d >= 3. R decreases from d = 7 on and is below 1/2 from d = 16 on; g decreases from d = 6 on. | PROVED (the runner checks, as a sanity check, that they dominate the exact values for 3 <= d <= 64, resp. 3 <= d <= 60) |
| 5 | Exact values: V_10 < 0.660228 and V_11 < 0.0859648 (without halving the bounds are at most 1.30746 and 0.156177); rho_d <= 1/2 for every 10 <= d <= 40, with rho_10 <= 0.49587. For 6 <= d <= 12 the potential forests reach the maximum-weight-forest value exactly. | FINITE-EXACT |
| 6 | Open Problem 3 with no evaluation at all (no exact base value). | OPEN (section 8) |

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
* **beta.** It is non-increasing, since beta(2m-1) = beta(2m) >= beta(2m+1). Also beta(N) <= sqrt(2/(pi N)) for
  N >= 1 (C(2m, m) <= 4^m/sqrt(pi m), and beta(2m-1) = beta(2m)).
* **Antipodal halving (growth note 1.4; THM-4525 L2).** H_(e-bar) is H_e reversed. So {e, f} and {e-bar, f-bar}
  collide together, and the involution fixes exactly the antipodal pairs, which form the type P(d-1) (nu = 0). In
  the union bound every type gets weight w = 1/2, except P(d-1), which gets weight 1.

## 2. Potential forests and the bound V_d

**Centrality.** On the levels 0..T of a type, b is *more central* than a, written b ≺ a, when |2b - T| < |2a - T|,
or when |2b - T| = |2a - T| and b > a. (Of two mirror levels, the upper one counts as more central.) This is a strict
total order.

**Potential forest.** Every level a that has a neighbour b ≺ a with N_ab > 0 *chooses* one such neighbour of largest
weight; write N(a) for that weight. Levels with no such neighbour are *roots*.

**Lemma 1 (PROVED).** The chosen edges form a forest. *Proof.* Orient each chosen edge away from its chooser. Two
levels cannot choose the same edge, since each would be more central than the other. Every level has out-degree at
most 1, and centrality strictly increases along an oriented edge. A cycle of chosen edges, with as many edges as
levels, would have to be a directed cycle, which is impossible. QED

The forest is the same idea as the toward-the-centre forest of the growth note (section 1.6 there), now with the
heaviest edge at every level and with ties broken. Define
* B(t) = prod over the choosers a of t of beta(N(a)),
* V_d = sum over the types t of Q_d of w(t) |t| B(t), the union bound with halving,
* A_d = w |t| B(t) for the antipodal type t = P(d-1): A_d = d 2^(d-2) prod over 0 <= j < (d-1)/2 of beta(4 C(d-1, j)).
  (Its level graph is the matching {j, d-1-j}, each lower level chooses its mirror, so B is the exact collision
  probability.)

By the forest lemma and the halving, **P(S is not resolving) <= V_d.**

**FINITE-EXACT (runner S2).** V_6, ..., V_12 <= 41.7956, 29.2119, 13.0044, 3.84260, **0.660228**, **0.0859648**,
0.00556133 (exact fractions, rounded up). These equal the maximum-weight-forest (Kruskal) values with halving,
exactly, for every 6 <= d <= 12. Without halving the Kruskal values are at most 1.30746 (d = 10), 0.156177 (d = 11)
and 0.0107481 (d = 12). These are the values the lane note quotes for the paper's bound (1.307, 0.156, 0.0107); its
section 4.3 quotes 0.1564 at d = 11 for its own version of the potential forests. The weights and counts were also
compared with the independent `cells` and `type_count` of `procgen_edim2_20261001_lib.py` for 6 <= d <= 12, during
development and again by the blind audit; the deposited runner checks them against brute force over representative
edge pairs for 3 <= d <= 8 (S1).

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

*Step 2: psi(b) is more central than phi(a) in t'.* Put delta = 2a - T and gamma = 2b - T, so that b ≺ a means
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

## 4. The recursion and the ratio bound

**Corollary 2 (PROVED).** For every d >= 3, V_(d+1) <= rho_d V_d + A_(d+1), where
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

**Lemma 4 (ratio bound, PROVED).** For every d >= 3,
rho_d <= R(d) := max(2 (d+1)^(3/2), 4 (d+1)^2/(d-1)) sqrt(2/pi) 2^(-d/2).
For d >= 7 the first term is the larger one, and R(d+1) < R(d). In particular R(15) <= 0.56419 and
R(16) <= 0.43693, so R(d) < 1/2 for every d >= 16.

*Proof.* Parallel t' = (P, h', nu'): h' + nu' = d with h', nu' >= 1, and t' is the only successor of its
predecessor. By Lemma 3, N*(t') >= 4 (2^(nu')/(nu'+1)) (2^(h')/(2(h'+1))) = 2^(d+1)/((nu'+1)(h'+1)). With
beta(N) <= sqrt(2/(pi N)), a weight ratio at most 1, and (nu'+1)(h'+1)/nu'^2 <= 2d <= 2(d+1) (as (nu'+1)/nu' <= 2
and (h'+1)/nu' <= d),
r(t') <= (2(d+1)/nu') sqrt((nu'+1)(h'+1)/(pi 2^d)) <= 2 (d+1)^(3/2) sqrt(2/pi) 2^(-d/2).
Crossing t' = (X, h', nu'): h' + nu' = d - 1, the count ratio has the denominator max(h', nu') >= (d-1)/2, the weight
ratio is 1, (nu'+1)(h'+1) <= ((d+1)/2)^2, and N*(t') >= 2 C(nu', floor(nu'/2)) C(h', floor(h'/2))
>= 2^d/((nu'+1)(h'+1)). So each r(t') is at most (2(d+1)^2/(d-1)) sqrt(2/pi) 2^(-d/2), and a load (two successors
at most) is at most twice that. For the last sentence: comparing squares, the first term is the larger one exactly
when (d-1)^2 >= 4(d+1), which holds for d >= 7; and the first term decreases because ((d+2)/(d+1))^3 <= (9/8)^3 < 2.
QED (Runner S4 checks Lemma 3 for m <= 300, the two inequalities of the last step for 7 <= d < 400, and, as a sanity
check, that R(d), rounded up with pi > 3.14159, dominates the exact bound below for every 3 <= d <= 64.)

The closed form at d = 10, ..., 16 is R(d) <= 1.8194, 1.4659, 1.1688, 0.92357, 0.72427, 0.56419, 0.43693.

**Exact ratio bounds (FINITE-EXACT, runner S4; used only in Corollary 6).** Every r(t') is bounded above by a
rational number: beta(N) itself for N <= 4000, else an integer-square-root upper bound for sqrt(2/(3.14159 N)). The
loads are summed exactly, and the table rounds up:

| d | 10 | 11 | 12 | 13 | 14 | 15 | 16 | 24 | 32 | 40 |
|---|---|---|---|---|---|---|---|---|---|---|
| rho_d <= | 0.49587 | 0.35085 | 0.26127 | 0.19845 | 0.15036 | 0.11649 | 0.087512 | 0.0087701 | 0.00077196 | 0.000063118 |
| largest load on | (X,4,4) | (X,4,5) | (X,5,5) | (P,12,0) | (P,12,1) | (P,14,0) | (P,14,1) | (P,22,1) | (P,30,1) | (P,38,1) |

Every value for 10 <= d <= 40 is at most 1/2. For d = 10..12 the largest load sits on a central crossing type. In the
rows from d = 13 on it sits on the antipodal type (P, d-1, 0) for odd d and on its neighbour (P, d-2, 1) for even d.
For comparison, rho_6, ..., rho_9 <= 1.7311, 1.2326, 0.92321, 0.65417.

## 5. Antipodal terms

**Lemma 5 (PROVED).** For every d >= 3, A_d <= g(d) := d 2^(d-2) (3/8) (2 pi (d-1))^(-(d-3)/4), and g(d+1) < g(d)
for every d >= 6.

*Proof.* A_d = d 2^(d-2) prod over 0 <= j < (d-1)/2 of beta(4 C(d-1, j)). The factor j = 0 is beta(4) = 3/8. For
1 <= j < (d-1)/2 (so j <= d - 2), C(d-1, j) >= d - 1, hence beta(4 C(d-1, j)) <= (2 pi (d-1))^(-1/2) < 1. There are
at least (d-3)/2 such factors. For the monotonicity,
g(d+1)/g(d) = (2(d+1)/d) (2 pi d)^(-1/4) ((d-1)/d)^((d-3)/4) <= (2(d+1)/d) (2 pi d)^(-1/4),
and (2(d+1)/d)^4 < 2 pi d for d >= 6 (the left side decreases, the right side increases, and at d = 6 they are
2401/81 < 29.7 and 12 pi > 37.6). QED (Runner S5 checks the last inequality for 6 <= d < 400 and, as a sanity check, that g(d)
dominates the exact rational upper bound for A_d for every 3 <= d <= 60.)

The closed form at d = 11, ..., 17 is g(d) <= 0.53498, 0.33457, 0.20226, 0.11863, 0.067701, 0.037687, 0.020507.
The exact rational upper bounds (FINITE-EXACT, runner S5; A_11 and A_12 are used only in Corollary 6) are
A_11 <= 0.015753, A_12 <= 0.00037463, A_13 <= 0.00038347 and A_14 <= 4.0273e-6.

## 6. Open Problem 3

**Theorem 5 (PROVED, computer-assisted at one point; Open Problem 3).** For every d >= 3,
V_(d+1) <= R(d) V_d + g(d+1).
Started from the exact value V_11 < 0.0859648 (runner S2), this gives V_12, ..., V_16 <= 0.4606, 0.7406, 0.8026,
0.6490, 0.4039, and V_d <= 1/2 for every d >= 16. Hence, for every d >= 11,
P(a uniformly random subset of V(Q_d) is not edge-multiset resolving) <= V_d < 1,
and V_d -> 0. With the exact V_10 < 0.660228 the same holds at d = 10.

*Proof.* Corollary 2 with Lemmas 4 and 5 gives the recursion for every d >= 3, and its right side increases with V_d.
The five steps from V_11 are evaluated in exact rational arithmetic, with R and g rounded up (runner S6). For
d >= 16, if V_d <= 1/2 then V_(d+1) <= R(16)/2 + g(17) <= 0.2390 < 1/2, because R decreases from d = 7 on and g from
d = 6 on. So V_d <= 1/2 for every d >= 16, by induction, and then V_(d+1) <= R(d)/2 + g(d+1) -> 0. Section 2 gives
the probability bound. QED

**Corollary 6 (PROVED, computer-assisted).** V_d <= V_10 < 0.660228 for every d >= 10. So a uniformly random subset
of V(Q_d) is resolving with probability more than 0.3397 for every d >= 10.

*Proof.* Induct on d. The step from d to d + 1 needs rho_d <= 1/2 and A_(d+1) <= V_10/2. The first holds by the exact
values for 10 <= d <= 15 (section 4) and by R(d) <= R(16) < 0.437 for d >= 16 (Lemma 4). For the second,
V_10 > 0.6602 (an exact fraction), so V_10/2 > 0.33; A_11 <= 0.015753 and A_12 <= 0.00037463 (section 5); and
A_k <= g(k) <= g(13) <= 0.20226 for every k >= 13 (Lemma 5). Then V_(d+1) <= V_10/2 + V_10/2 = V_10. QED

**How this answers the question as posed.** THM-4534 and its lane note record Open Problem 3 as "one analytic
estimate valid from d = 11". The lane note suggested this route (section 4.3 there): "a monotonicity lemma
V_(d+1) <= V_d for the potential-forest bound would reduce Open Problem 3 to the single value V_11". Theorem 1 and
Lemmas 4 and 5 make it rigorous in a quantitative form, with the antipodal term separated. The result is not
V_(d+1) <= V_d itself, but V_(d+1) <= R(d) V_d + g(d+1) with closed-form coefficients. The estimate is one recursion,
valid at every d >= 3, and its one numerical input is the exact value V_11 (with antipodal halving). This is the
same shape as the paper's tail (one closed form, one evaluation at d = 51) and THM-4534's Theorem D (one closed form,
evaluations at d = 17, 18, 19), and it starts at d = 11, where the paper's computations start. Read this way, it
answers Open Problem 3. If the problem is read as forbidding any evaluation of the union bound, it stays open
(section 8).

Explicit resolving sets for 6 <= d <= 9 (THM-4525 and the first lane note
[`procgen_edim_20261001_edge_multiset_dimension.md`](procgen_edim_20261001_edge_multiset_dimension.md)), with
Theorem 5 from d = 10 on, give edim_m(Q_d) < infinity for every d >= 6. THM-4534 already proves this, with explicit
sets up to d = 16 and its Theorem D beyond.

## 7. Relation to other work

* THM-4534 (audited canon; procgen_edim2 lane): Theorem D proves U_d < 1 in closed form for every d >= 17. Its only
  numerical inputs are interval evaluations at d = 17, 18 and 19 (bounds 0.6175, 0.01186 and 0.00297), and from
  d = 20 on it is elementary. Its Fourier-Hoelder multi-forest lemma certifies U_10 <= 0.2649 and U_11 <= 0.0448
  without halving, so positive probability at d = 10 was already known there (lane-level numerics, as THM-4534
  labels them). What is new here is the halved potential-forest bound V_10 < 1, the transfer lemma, and the recursion
  that reduces every d >= 11 to the single value V_11. THM-4534 lists "OP3 in the strict form" as OPEN; Theorem 5
  answers it in the sense of section 6. Proposition K there (weighted spanning trees of the parallel level graph) is
  not used here.
* Growth note ([`edge_multiset_dimension_growth_20261002.md`](edge_multiset_dimension_growth_20261002.md), section 6):
  an earlier version recorded Open Problem 3 as OPEN and asked for "a lifting lemma from Q_d to Q_(d+1)". Theorem 1
  is such a lemma, for the union bound rather than for resolving sets; that section now points here. Its own
  sparse-regime estimate (Theorem A there) is effective from d = 41.

## 8. Not done

* **Open Problem 3 with no evaluation at all (OPEN).** The closed forms alone cannot start below d = 16, since
  R(15) > 1/2. From V_10 alone they give only V_11 <= 1.74. The audit also found that from V_10, with the exact rho_10
  and A_11, they fail at d = 13 and 14 (bounds 1.18 and 1.21). A proof with no exact base value would need closed forms for
  rho_d, or for V_11 itself, that are much sharper at d <= 15.
* A base below d = 10. V_9 <= 3.84260 > 1, although rho_9 <= 0.65417. The lane note's Monte Carlo estimate (section
  4.3 there; not a proof) of the expected number of colliding pairs at d = 9 is about 0.8, so a first-moment
  argument at d = 9 would need much more than forests.
* The Fourier-Hoelder lemma of THM-4534 was not combined with the transfer. Doing so would need edge-disjoint extra
  edges, one for each forest.
* A second independent audit of this revised version.

## 9. Reproduction

```
python3 04-computation/experiments/edge_multiset_dimension_op3_20261002_run.py [DMAX_TRANSFER=40]
```
Sections: S1 cells and counts against brute force; S2 the union bound, the Kruskal comparison, V_10 and V_11; S3 the
transfer lemma, sub-cell by sub-cell; S4 the exact ratio bounds and the closed form R(d); S5 the exact antipodal
terms and the closed form g(d); S6 Theorem 5 and Corollary 6. 573,986 checks, about 9 seconds, ALL CHECKS PASSED.
Timing goes to stderr; stdout is deterministic.
