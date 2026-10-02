# Edge multiset dimension of hypercubes II: ln edim_m(Q_d) = Theta(d^(1/3))

**Provenance:** Claude thread session (project thread, branch `claude/project-thread-96889h`), 2026-10-01/02;
temporary identity, no `.machine-id`. This is the "edim2" lane named in the wave-24 letter of
collatz-procgen-20260922 (Allikvere's Open Problems 2-4), continued here; no edim2 result of that session was on
`origin/main` at 05791bbd. It continues the first note,
[`procgen_edim_20261001_edge_multiset_dimension.md`](procgen_edim_20261001_edge_multiset_dimension.md).
**Source problem:** J. Allikvere, *The edge multiset dimension of hypercubes*, arXiv:2608.09983v1,
Section 10: Open Problems 2, 3 and 4, and the bounds for Q_7. This session could not fetch the paper
(the container's network policy blocks arxiv.org downloads), so the problems are used as the first lane
note and [THM-4525](../../01-canon/theorems/THM-4525-edge-multiset-dimension-of-q6-is-15.md) record them:
OP2 asks for the growth of edim_m(Q_d), OP3 for one analytic estimate valid from d = 11, OP4 compares
density-1/2 landmark sets with sparse ones.
**Builds on:** [THM-4525](../../01-canon/theorems/THM-4525-edge-multiset-dimension-of-q6-is-15.md)
(edim_m(Q_6) = 15; Lemmas L1-L5; explicit sets for Q_7..Q_12) and
[HYP-9169](../hypotheses/HYP-9169-edge-multiset-dimension-of-hypercubes-grows-like-exp-cube-root.md)
(the growth conjecture; proved here and, independently, in THM-4534, which upgraded it).
**Audit:** a blind independent audit (2026-10-02) found no false claim and no substantive gap, and re-certified
every row of the table with its own code; its label, credit and wording fixes, and one monotonicity line for
d > 10^5, are applied here.
**Runner:** `04-computation/experiments/edge_multiset_dimension_growth_20261002_run.py`; stdout in
`05-knowledge/results/edge_multiset_dimension_growth_20261002.out`, ending with ALL CHECKS PASSED.
**Parallel work.** The mac-mini lane note
[`procgen_edim2_20261001_growth_and_uniform_bounds.md`](procgen_edim2_20261001_growth_and_uniform_bounds.md)
landed on `origin/main` (12d7b60c, 00:49:53 UTC) while this note was being finished, and was then audited and
promoted as [THM-4534](../../01-canon/theorems/THM-4534-edge-multiset-dimension-grows-like-exp-theta-cube-root-d.md).
The two were written independently: this note's first commit (fc801d67, 00:38:50 UTC) already had Theorems A-D and
both constants. It reaches the same constants kappa_+ = 2.0526 and kappa_- = 0.8146, with a different forest
construction (a star lemma with path and zigzag forests, against Lemma H and the toward-the-centre forest here). For
Open Problem 4 it states the upper half of the density dichotomy (its item 10); the lower half, P -> 0 below
kappa_-, is not stated there, although it follows at once from its Theorem B and Markov's inequality. It goes
further on two points. (i) Open Problem 3: a closed-form estimate at density 1/2 for every d >= 17 (its Theorem D),
against d >= 41 here in the sparse regime (Theorem A; section 6). The companion note
[`edge_multiset_dimension_op3_20261002.md`](edge_multiset_dimension_op3_20261002.md) now settles Open Problem 3
from d = 10. (ii) The sparse table: its Lemma FH gives M_11 = 361, M_12 = 370, M_13 = 457 and M_14 = 490, against
492, 543, 622 and 724 here. Its deposited output is a `--quick` run, which certifies only d = 11..14; its larger
values (M_16 = 636, M_64 = 31808) come from a full run that was not deposited, and THM-4534 lists them as lane-only.
It also has a result with no counterpart here: its Theorem C shows that no choice of forests makes the forest-lemma
union bound work below kappa_+. New here: the lower half of the dichotomy (Theorem D(2)), Lemma H and the
toward-the-centre forest (item 5), the explicit bound edim_m(Q_d) < 2 exp(5 d^(1/3)) for every d >= 6 (Theorem C),
the certified table for every 11 <= d <= 64 with its certificates deposited, and the Q_7 search (13 <= edim_m(Q_7),
against 11 there). So items 1 and 4 (upper half) below are an independent second derivation of results stated
there.

## Status header

| # | Claim | Label |
|---|-------|-------|
| 1 | **Open Problem 2 (the growth of edim_m(Q_d), as THM-4525 and the first note record it): ln edim_m(Q_d) = Theta(d^(1/3)).** More precisely, kappa_- <= liminf ln edim_m(Q_d)/d^(1/3) <= limsup <= kappa_+, where kappa_+ = (3 sqrt2 ln 2)^(2/3) = 2.0526 and kappa_- = ((3 sqrt2/4) ln 2)^(2/3) = kappa_+/4^(2/3) = 0.8146. This proves HYP-9169, independently of THM-4534 (same constants). The paper's own wording of Open Problem 2 could not be read here; if it asks for the constant, that part is open (section 9). | PROVED (Theorems A and B); blind independent audit 2026-10-02 found no false claim |
| 2 | **Explicit form, every d >= 6:** edim_m(Q_d) < 2 exp(5 d^(1/3)). | PROVED, computer-assisted (Theorem C). Its finite inputs are FINITE-EXACT: the explicit sets of THM-4525 for 6 <= d <= 12, the certified table for 11 <= d <= 64 (check S5), and Theorem A in interval arithmetic for 64 <= d <= 10^5 (check S6); by hand beyond |
| 3 | **Certified sparse table.** For every 11 <= d <= 64, a random set of rational density q_d is resolving and has at most M_d points with positive probability, so edim_m(Q_d) <= M_d. Examples: Q_11 <= 492, Q_16 <= 978, Q_32 <= 6242, Q_64 <= 57738 (full table in section 4). ln M_d / d^(1/3) lies in [2.731, 2.788]. | FINITE-EXACT (outward-rounded interval arithmetic). Supersedes item 7 of the first note, which was double precision only. THM-4534 is sharper where its deposited output certifies (d = 11..14, Lemma FH) |
| 4 | **Open Problem 4 (density 1/2 versus sparse).** For random sets S_q: P(S_q resolves) -> 1 when ln(q 2^d) >= (kappa_+ + eps) d^(1/3), uniformly in q <= 1/2; P(S_q resolves) -> 0 when ln(q 2^d) <= (kappa_- - eps) d^(1/3). Minimum resolving sets have density 2^(-d) e^(Theta(d^(1/3))), so density-1/2 sets are larger by a factor 2^(d - O(d^(1/3))). | PROVED (Theorem D). THM-4534 (lane note item 10) states the upper half |
| 5 | Lemma H (a "valid cell" lemma for the hypergeometric law) and the toward-the-centre forest. For every pair type with h, nu >= 1, every level with \|t\| > 1 has a flow edge toward the middle whose cell carries at least 1/(3(min(h,nu)+1)) of the level. | PROVED; exact check for all L <= 60 |
| 6 | Open Problem 3 (one analytic estimate valid from d = 11). Here Theorem A is effective from d = 41 (lambda chosen per d) and from d = 64 with lambda = 5 d^(1/3). THM-4534 reaches d = 17 at density 1/2. | Settled from d = 10 in the companion note [`edge_multiset_dimension_op3_20261002.md`](edge_multiset_dimension_op3_20261002.md) (PROVED there, computer-assisted at its base; audit owed). Here: PROVED for d >= 41 (check S6) |
| 7 | **Q_7: no resolving set of size <= 12, so 13 <= edim_m(Q_7) <= 19** (was 8..19; 11..19 in procgen_edim2_20261001, k <= 10). Exhaustive search up to Aut(Q_7), validated end to end on Q_6 (section 7; full write-up in [the Q_7 note](edge_multiset_dimension_q7_20261002.md)). | FINITE-EXACT (computer-assisted); a blind independent audit (2026-10-02) of the Q_7 note found no false claim and no gap |

## 0. Definitions and notation

* Q_d has vertex set {0,1}^d and edges e = {u, u + e_i}. For a vertex s, d(e, s) = min(d(u,s), d(u+e_i, s)).
  For a vertex set S, the histogram of e is H_e(r) = #{s in S : d(e,s) = r}, r = 0..d-1.
  S is (edge-multiset) resolving when the d 2^(d-1) histograms are pairwise distinct; edim_m(Q_d) is
  the least |S|. edim_m(Q_d) is infinite exactly for 2 <= d <= 5 (the paper's Theorem 6, as cited in the first
  note): it is finite for every d >= 6, and trivially for d = 1, which has a single edge.
* Projection lemma (paper, Lemma 3): for e = {u, u + e_i}, d(e, s) = d_H(u', s'), where ' deletes
  coordinate i.
* b_n(k) = C(n,k)/2^n. D(p||q) = p ln(p/q) + (1-p) ln((1-p)/(1-q)) is the binary divergence, so
  D(p||1/2) = p ln(2p) + (1-p) ln(2(1-p)). Logarithms are natural throughout.
* A random landmark set S_q contains each vertex independently with probability q. Write m = q 2^d and
  lambda = ln m.
* E = d 2^(d-1) is the number of edges.

## 1. Tools

### 1.1 Pair types and cells (PROVED; brute-force checked)
Aut(Q_d) acts on unordered pairs of distinct edges. Two edges are *parallel* (same direction i) or
*crossing* (directions i != k). Put L = d - 1 (parallel) or L = d - 2 (crossing). Let h be the
number of the remaining L coordinates in which the two edges' base points differ, and nu = L - h.

* **Parallel, 1 <= h <= d-1.** There are N_P(h) = d 2^(d-2) C(d-1,h) such pairs. A vertex s
  determines two numbers: c, its disagreements with u' on the nu common coordinates, and j, its
  disagreements with u' on the h differing coordinates. Then (d(e,s), d(f,s)) = (c + j, c + h - j).
  The cell (c, j) has M(c,j) = 2 C(nu,c) C(h,j) = 2^d b_nu(c) b_h(j) vertices.
* **Crossing, 0 <= h <= d-2.** There are N_X(h) = d(d-1) 2^(d-1) C(d-2,h) such pairs. Besides c and j
  there are two bits: x (coordinate k against e) and y (coordinate i against f). Then
  (d(e,s), d(f,s)) = (x + c + j, y + c + h - j), and the sub-cell (x, y, c, j) has
  C(nu,c) C(h,j) = 2^d b_nu(c) b_h(j)/4 vertices.
* The counts add up to C(E, 2). The runner re-derives every cell matrix and count by brute force over
  all pairs of edges for d = 3, 4, 5 (check S1).

### 1.2 The forest lemma at any density (PROVED; the paper's Lemma 11 at general q)
Fix edges e != f and let X_ab = |S_q intersect {s : d(e,s) = a, d(f,s) = b}|. These are independent
binomials Bin(M_ab, q). H_e = H_f holds exactly when the "flow" D_ab = X_ab - X_ba has zero net
outflow at every level r: the sum over b != r of D_rb is 0 (diagonal cells cancel).

**Lemma 1.** For every forest F on the level set {0, ..., d-1},
P(H_e = H_f) <= prod over {a,b} in F of max_k P(D_ab = k).

*Proof.* Condition on all D_ab with {a,b} not in F; they involve cells disjoint from those of F. List
the edges of F as g_1, ..., g_t so that g_j is a leaf edge of the forest {g_j, ..., g_t}, with leaf
level r_j. Among the forest variables, the conservation equation at r_j involves only
D_(g_1), ..., D_(g_j). So the equation fixes D_(g_j) as a function of the conditioned variables and of
D_(g_1), ..., D_(g_(j-1)). The D_(g_j) are independent, so the probability that all t equations hold
is at most the product of the largest atoms. QED

### 1.3 Atom bounds (PROVED)
**Lemma 2.** Let D = X - Y with X ~ Bin(M_1, q), Y ~ Bin(M_2, q) independent, N = M_1 + M_2 and
x = N q (1 - q). Then for every integer k,
P(D = k) <= e^(-x) I_0(x) = (1/pi) int_0^pi e^(-x(1 - cos t)) dt.
Moreover:
1. e^(-x) I_0(x) <= sqrt(pi/(8x)) for all x > 0;
2. e^(-x) I_0(x) <= (1 + 1/(8x) + 9/(128 x^2) + C_4/x^3)/sqrt(2 pi x) + e^(-x)/2 for all x > 0, with
   C_4 = (75/1024) 2^(7/2) = 0.8286.

*Proof.* Each of the N Bernoulli steps contributes a factor of modulus
|1 - q + q e^(+-it)| = (1 - 2q(1-q)(1 - cos t))^(1/2) to the characteristic function, so
P(D = k) <= (1/2pi) int_(-pi)^(pi) (1 - 2q(1-q)(1-cos t))^(N/2) dt <= (1/2pi) int e^(-x(1-cos t)) dt.
(1) Since 1 - cos t = 2 sin^2(t/2) >= 2t^2/pi^2 on [-pi, pi], the integral is at most
(1/2pi) int_R e^(-2x t^2/pi^2) dt = sqrt(pi/(8x)).
(2) Substitute u = 1 - cos t: the integral is (1/pi) int_0^2 e^(-xu) du / sqrt(u(2-u)). On [1,2], bound
e^(-xu) by e^(-x); the remaining integral is pi/2, giving e^(-x)/2. On [0,1], write
(2-u)^(-1/2) = 2^(-1/2) (1-y)^(-1/2) with y = u/2 <= 1/2. Taylor's theorem with the Lagrange
remainder gives (1-y)^(-1/2) <= 1 + y/2 + 3y^2/8 + (5/16) 2^(7/2) y^3 on [0, 1/2]. Extend the
integral to [0, infinity) and use int_0^inf e^(-xu) u^(k - 1/2) du = Gamma(k + 1/2) x^(-k-1/2). QED

The runner also encloses e^(-x) I_0(x) directly, from the series sum (x^2/4)^k/(k!)^2 with a geometric
tail bound. It uses the series for x < 60 and bound (2) for x >= 60. Bound (2) exceeds the true value by
a relative 2.8e-5 at x = 30 and 7.5e-7 at x = 100 (check S2).

### 1.4 Antipodal halving (PROVED)
For the antipodal edge e-bar, d(e-bar, s) = d - 1 - d(e, s), so H_(e-bar) is H_e reversed (L2 of
THM-4525). Hence the collision events of the pairs {e,f} and {e-bar, f-bar} coincide. The map
{e,f} -> {e-bar, f-bar} is an involution on pairs. It preserves every pair type, and its only fixed
pairs are the antipodal pairs {e, e-bar}, which form the type P(d-1). So a union bound may count every
type except P(d-1) with weight 1/2.

### 1.5 Lemma H: a valid cell toward the middle (PROVED)
Fix h, nu >= 1, L = h + nu, and a level a with t := a - L/2 > 1. Let
p(k) = C(h,k) C(nu, a-k)/C(L,a). This is the hypergeometric law of the h-coordinate count k given the
level a, supported on S_a = [max(0, a - nu), min(h, a)]. Call k *valid* if h/2 < k < h/2 + t.

**Lemma H.** Some valid k in S_a has p(k) >= 1/(3 |S_a|) >= 1/(3 (min(h, nu) + 1)).

*Proof.*
1. *Unimodality.* For k, k+1 in S_a, p(k+1)/p(k) = (h-k)(a-k)/((k+1)(nu-a+k+1)), and
   p(k+1) >= p(k) iff k + 1 <= rho := (a+1)(h+1)/(L+2). So the maximizers of p are floor(rho) and,
   when rho is an integer, also rho - 1. One checks rho - 1 < h, rho - 1 < a and
   (a+1)(h+1) - (a-nu)(L+2) = nu(L-a+1) + (L-a) + 1 > 0, so floor(rho) lies in S_a.
2. *Location.* With mu = ah/L, rho - mu = (a nu + h(L-a) + L)/(L(L+2)) lies in (0, 1). So every
   maximizer k* satisfies mu - 1 < k* < mu + 1. Here mu = h/2 + th/L with 0 < th/L < t.
3. *Case: some maximizer k* is valid.* Then p(k*) >= 1/|S_a|.
4. *Case k* <= h/2.* Then k* = floor(h/2), because k* > mu - 1 > h/2 - 1. The point k* + 1 is valid
   (it is <= h/2 + 1 < h/2 + t) and lies in S_a. Its ratio is
   p(k*+1)/p(k*) = [(h-k*)/(k*+1)] [(a-k*)/(nu-a+k*+1)] >= [h/(h+2)] [1] >= 1/3,
   since h - k* >= h/2, k* + 1 <= h/2 + 1, a - k* >= nu/2 + t and 1 <= nu - a + k* + 1 <= nu/2 - t + 1.
5. *Case k* >= h/2 + t.* Then k* = ceil(h/2 + t) (as k* < mu + 1 < h/2 + t + 1), and c* = a - k* <= nu/2.
   The point k* - 1 is valid and lies in S_a (its c-value is c* + 1 <= floor(nu/2) + 1 <= nu). Its ratio is
   p(k*-1)/p(k*) = [(nu-c*)/(c*+1)] [k*/(h-k*+1)] >= [nu/(nu+2)] [1] >= 1/3.
In every case some valid k has p(k) >= (max p)/3 >= 1/(3|S_a|), and |S_a| <= min(h, nu) + 1. QED

The runner verifies the inequality exactly, in rational arithmetic, for every 2 <= L <= 60, every h and
every level with |t| > 1: 69310 cases (check S3). The smallest ratio of the best valid p to the bound
is 2.45.

### 1.6 The toward-the-centre forest (PROVED)
In the flow graph of a pair type with h, nu >= 1, the cell (c, j) joins the levels c + j and
c + h - j. Take a level a with t = a - L/2 > 1 and a valid k. The two cells (a - k, k) and
(a - k, h - k) join a to b = a - 2(k - h/2), in the two orientations, and
|b - L/2| = |t - 2(k - h/2)| < |t|. Levels with t < -1 are handled by the reflection a -> L - a, k -> h - k.

**Lemma 3.** Choose for every level with |t| > 1 one such edge. The chosen edges form a forest.

*Proof.* Each level chooses one edge, toward a strictly smaller |t|. Two different levels cannot
choose the same edge, since each would then be strictly more central than the other. Orient every chosen
edge away from its chooser. Then every level has out-degree at most 1, and |t| strictly decreases along
oriented edges. A cycle in the chosen edges would have to be oriented cyclically, which is impossible. QED

### 1.7 Binomial estimates (standard)
* b_n(k) >= e^(-n D(k/n || 1/2))/(n+1): the term k of Bin(n, k/n) is its largest, hence >= 1/(n+1).
* b_n(k) <= e^(-n D(k/n || 1/2)) <= e^(-2(k - n/2)^2/n).
* n D(1/2 + t/n || 1/2) <= phi_n(t) := 2t^2/n + R_n(t) with R_n(t) = (4t^4/(3n^3))/(1 - 4t^2/n^2), for
  |t| < n/2. This follows from D(1/2+e || 1/2) = sum over j >= 1 of (2e)^(2j)/(2j(2j-1)).

### 1.8 A lattice-sum inequality (PROVED)
**Lemma 4.** Let f be decreasing on [t_1, tau], and let t_1 < t_1 + 1 < ... be the lattice points.
Then the sum of max(0, f(t_1 + i)) over the lattice points t_1 + i <= tau is at least int_(t_1)^(tau) f.

*Proof.* If f(t_1) < 0, then f < 0 on [t_1, tau] and the integral is negative. Otherwise let t_0 be the
supremum of the points in [t_1, tau] where f >= 0. Then int_(t_1)^(tau) f <= int_(t_1)^(t_0) f. For each
lattice point t_1 + i <= t_0, monotonicity gives f(t_1 + i) >= the integral of f over
[t_1 + i, min(t_1 + i + 1, t_0)], and these intervals cover [t_1, t_0]. QED

## 2. Theorem A: the upper bound

For L in {d-1, d-2} and kappa in {2, 1/2}, write
Lambda(L, kappa) = lambda - ln(3 pi (L+1)^2/(8 kappa (1-q))).
With tau = sqrt(L Lambda/2), set
g(L, Lambda) = (2/3) Lambda tau - 3 Lambda - (tau/15 + 2/3) Lambda^2/(L - 2 Lambda) and
g'(L, Lambda) = (1/3) Lambda tau - Lambda - (tau/15 + 2/3) Lambda^2/(2(L - 2 Lambda)).

**Theorem A (explicit union bound).** Let d >= 6, lambda > 0, m = e^lambda and q = m 2^(-d) <= 1/4.
Put Lambda_P = lambda - ln(pi d/(16(1-q))) and Lambda_X = lambda - ln(pi(d-1)/(4(1-q))). Assume that every
Lambda used below satisfies 2 Lambda < L and tau >= 3. Then
P(S_q is not resolving) <= U := (E^2/4) e^(-g_*) + d 2^(d-2) e^(-g'(d-1, Lambda_P)) + d(d-1) 2^(d-2) e^(-g'(d-2, Lambda_X)),
where g_* is the least of the four values g(L, Lambda(L, kappa)) and g(d-2, Lambda_X).
If U + e^(-m/3) < 1, then edim_m(Q_d) < 2m.

*Proof.*
1. *Generic types* (parallel with 1 <= h <= d-2, crossing with 1 <= h <= d-3). Use the forest of
   Lemma 3 (Lemma H supplies the valid k). For crossing types use only the sub-cells with x = y = 0;
   they have the same (c, j) geometry with L = d - 2, and further cells only enlarge N. The edge
   {a, b} chosen at level a uses the two cells (a - k, k) and (a - k, h - k). Each has
   2 C(nu, a-k) C(h, k) vertices (parallel) or C(nu, a-k) C(h, k) vertices (crossing, x = y = 0).
   Since C(nu, a-k) C(h, k) = 2^L b_L(a) p(k), this gives N = 2^(d+1) b_L(a) p(k) (parallel) and
   N >= 2^(d-1) b_L(a) p(k) (crossing). Hence
   x_e = N q (1-q) >= kappa (1-q) m b_L(a) p(k), with kappa = 2 (parallel) or 1/2 (crossing).
   By Lemma H and 1.7, x_e >= kappa (1-q) m e^(-phi_L(t))/(3 (L+1)^2). Every atom is at most 1, so
   Lemmas 1 and 2(1) give
   -ln P(H_e = H_f) >= (1/2) sum over levels with |t| > 1 of max(0, Lambda - phi_L(t)), with Lambda = Lambda(L, kappa).
2. *Lattice sum.* The function f(t) = Lambda - phi_L(t) is even and decreasing in |t| on [0, L/2), and
   tau < L/2. On each side, the first lattice point t_1 with t_1 > 1 satisfies t_1 <= 2. By Lemma 4 and f <= Lambda,
   that side contributes at least int_(t_1)^(tau) f >= int_0^tau f - 2 Lambda. Now
   int_0^tau (Lambda - 2t^2/L) dt = Lambda tau - 2 tau^3/(3L) = (2/3) Lambda tau, because tau^2 = L Lambda/2. Also
   R_L(t) <= (4t^4/(3L^3)) L/(L - 2 Lambda) for t <= tau, so int_0^tau R_L <= tau Lambda^2/(15(L - 2 Lambda)).
   Adding both sides and halving,
   -ln P(H_e = H_f) >= (2/3) Lambda tau - 2 Lambda - tau Lambda^2/(15(L - 2 Lambda)) >= g(L, Lambda).
3. *Crossing, h = 0.* The flow graph is the path {c, c+1}, c = 0..d-2, and the edge {c, c+1} has
   N = 2 C(d-2, c), so x = (1-q) m b_(d-2)(c)/2. Give every level with |s| > 1 the path edge toward the centre
   (d-1)/2.
   A level at offset s from (d-1)/2 then uses b_(d-2)(c) with c at distance |s| - 1/2 from (d-2)/2, and
   phi_(d-2)(|s| - 1/2) <= phi_(d-2)(s), so (1/2) ln(8x/pi) >= (1/2)(Lambda_X - phi_(d-2)(s)). The
   computation of step 2 with L = d - 2 then gives -ln P >= g(d-2, Lambda_X).
4. *Antipodal types.* For P(d-1) (nu = 0, L = d-1) the cells (0, j) give the matching {j, L-j}, j < L/2,
   with x = 2(1-q) m b_L(j) (when d is odd the middle level L/2 is left unmatched). It is a forest, and step 2 applied to one side (the matching pairs t with -t)
   gives -ln P >= (1/3) Lambda_P tau - Lambda_P - tau Lambda_P^2/(30(L - 2 Lambda_P)) >= g'(d-1, Lambda_P).
   For X(d-2) the x = y = 0 sub-cells give the matching {j, L - j} with x >= (1-q) m b_L(j)/2, L = d-2,
   and -ln P >= g'(d-2, Lambda_X) in the same way.
5. *Union bound with halving (1.4).* The generic types and X(0) contain at most C(E,2) <= E^2/2 pairs,
   counted with weight 1/2. P(d-1) has d 2^(d-2) pairs, weight 1. X(d-2) has d(d-1) 2^(d-1) pairs, weight 1/2.
6. *Size.* |S_q| ~ Bin(2^d, q), and P(|S_q| >= 2m) <= e^(-m/3) (Chernoff). So with probability at least
   1 - U - e^(-m/3) > 0, the set S_q is resolving and has fewer than 2m points. QED

The runner checks numerically (check S4) that the actual toward-the-centre forest sums dominate
g(L, Lambda) for d in {20, 30, 40, 64, 100, 150}, lambda = C d^(1/3) with C in {3, 4, 5, 6}, and every
generic type: 3042 cases, smallest slack 19.9 nats. The same check for the special types X(0), P(d-1) and
X(d-2) against g(d-2, Lambda_X), g'(d-1, Lambda_P) and g'(d-2, Lambda_X) covers 63 cases, smallest slack 7.3 nats.

**Corollary A (asymptotics).** For every eps > 0 there is d_0(eps) such that
edim_m(Q_d) < 2 exp((kappa_+ + eps) d^(1/3)) for all d >= d_0(eps). Hence
limsup ln edim_m(Q_d)/d^(1/3) <= kappa_+ = (3 sqrt2 ln 2)^(2/3) = 2.0526.

*Proof.* Take lambda = (kappa_+ + eps) d^(1/3). Every Lambda equals lambda - O(ln d) and
(2/3) Lambda tau = (sqrt2/3) L^(1/2) Lambda^(3/2). The other terms of g are O(d^(1/3)), so
g_* >= (sqrt2/3) d^(1/2) lambda^(3/2) - O(d^(2/3) ln d) = (1 + eps/kappa_+)^(3/2) (2 ln 2) d - o(d),
because (sqrt2/3) kappa_+^(3/2) = 2 ln 2. Since ln(E^2/4) = 2d ln 2 + O(ln d), the first term of U tends
to 0. Likewise g' >= (1 + eps/kappa_+)^(3/2) (ln 2) d - o(d), against ln(d^2 2^d) = d ln 2 + O(ln d).
The side conditions hold for large d, since q -> 0 and Lambda = Theta(d^(1/3)) = o(L). QED

*Where kappa_+ comes from.* A typical pair has about 2 tau = sqrt(2 L lambda) central levels. Each pays
about (1/2) ln(m b_L) in the forest product, so the sum is about (sqrt2/3) L^(1/2) lambda^(3/2), and it
must beat ln(number of pairs) = 2 d ln 2.

## 3. Theorem B: the lower bound

**Theorem B.** For every eps > 0 and all large d, every resolving set S of Q_d has
ln|S| >= (kappa_- - eps) d^(1/3), with kappa_- = ((3 sqrt2/4) ln 2)^(2/3) = 0.8146.
This improves the constant c = (ln 2/sqrt 2)^(2/3) = 0.6216 of L5 in THM-4525 (written 0.6215 there).
The inequality used is L5's own; only its evaluation is sharper.

*Proof.* Let m = |S|, lambda = ln m, L = d - 1, mu_r = m b_L(r) and t = r - L/2. In nats, L5 says
(d-1) ln 2 + ln d <= sum_(r=0)^(L-1) G(mu_r), where G(mu) = (mu+1) ln(mu+1) - mu ln mu. Adding the nonnegative
term r = L only weakens it, and the proof uses that weaker form, with the sum running to L.
The right side increases with m, so it suffices to show that it is below (d-1) ln 2 when
lambda = (kappa_- - eps) d^(1/3) and d is large.
1. G(mu) = ln(1+mu) + mu ln(1 + 1/mu) <= ln(1 + mu) + 1. For mu <= 1, G(mu) <= mu(1 + ln 2 + ln(1/mu)).
2. mu_r <= m e^(-2t^2/L), by 1.7.
3. Let tau_0 = sqrt(L lambda/2). For |t| <= tau_0, m e^(-2t^2/L) >= 1, so
   G(mu_r) <= lambda - 2t^2/L + 1 + ln 2. For an even function that decreases in |t|, the lattice sum is at
   most the integral plus twice its maximum. So these levels contribute at most
   (4/3) lambda tau_0 + 2 lambda + (2 tau_0 + 2)(1 + ln 2).
4. For |t| > tau_0, mu_r <= e^(-y) with y = 2(t^2 - tau_0^2)/L >= 4 tau_0 (|t| - tau_0)/L. The bound
   mu(1 + ln 2 + ln(1/mu)) increases on (0, 1], and e^(-y)(1 + ln 2 + y) decreases in y. So these
   levels contribute at most 2(1 + ln 2) + (2 + ln 2) L/(2 tau_0).
5. With lambda = (kappa_- - eps) d^(1/3) we have tau_0 = O(d^(2/3)) and L/tau_0 = O(d^(1/3)). Altogether
   sum G(mu_r) <= (2 sqrt2/3) L^(1/2) lambda^(3/2) + O(d^(2/3)) = (1 - eps/kappa_-)^(3/2) d ln 2 + O(d^(2/3)),
   because (2 sqrt2/3) kappa_-^(3/2) = ln 2. This is below (d-1) ln 2 for large d. QED

Check S9 solves L5 numerically for d up to 10^6 (exact G; the binomials come from `gammaln` in double precision,
so S9 is a sanity check, not part of the proof). The least lambda that L5 allows,
divided by d^(1/3), is 1.11, 1.07, 0.99, 0.92, 0.876 at d = 10^2, ..., 10^6, so it approaches kappa_- = 0.8146
slowly. The excess is of the size of one lower-order term: the factor sqrt(2/(pi L)) in b_L shifts lambda by
about (1/2) ln(pi d/2), which is 0.071 d^(1/3) at d = 10^6, against an excess of 0.061 d^(1/3).

## 4. Theorem C: an explicit bound for every d, and the certified table

**Theorem C.** edim_m(Q_d) < 2 exp(5 d^(1/3)) for every d >= 6.

*Proof.*
* 6 <= d <= 12: the explicit resolving sets of THM-4525 have sizes 15, 19, 26, 38, 48, 65, 76. All are
  below 2 e^(5 d^(1/3)), which exceeds 16000 at d = 6.
* 11 <= d <= 64: the certified table below gives M_d <= e^(2.79 d^(1/3)).
* 64 <= d <= 100000: Theorem A with lambda = 5 d^(1/3). Every quantity is enclosed in outward-rounded
  interval arithmetic (mpmath.iv), and U + e^(-m/3) < 1 at every integer d in the range (check S6;
  the largest value is 0.661, at d = 64).
* d > 100000: by hand. Here q <= e^(-0.6 d). Every Lambda is at least lambda - 2 ln d - 0.87
  (the worst constant is ln(3 pi/4) = 0.857), and this is at least 4.4 d^(1/3): the function
  0.6 d^(1/3) - 2 ln d - 0.87 is positive at 10^5 and increasing for d > 1000. Also Lambda <= 5 d^(1/3),
  tau <= 1.6 d^(2/3) and L - 2 Lambda >= 0.99 d. Then
  g_* >= (sqrt2/3) (d-2)^(1/2) (4.4 d^(1/3))^(3/2) - 15 d^(1/3) - 3 d^(1/3) >= 4.2 d and
  g' >= 2.1 d, so U <= e^(2d ln 2 + 2 ln d - 4.2 d) + 2 e^(d ln 2 + 2 ln d - 2.1 d) < e^(-1.4 d).
  Each step holds for every d > 10^5, by monotonicity. The third term of g is at most
  (1.6 d^(2/3)/15 + 2/3) 25 d^(2/3)/(0.99 d) = (2.694 + 16.84 d^(-2/3)) d^(1/3) <= 3 d^(1/3). Since
  (d-2)^(1/2) d^(1/2) >= d - 2, g_* >= 4.3508 (d-2) - 18 d^(1/3), and 0.1508 d - 8.71 - 18 d^(1/3) is positive at
  10^5 and increasing for d > 251; so g_* >= 4.2 d. Likewise g' >= 2.1754 (d-2) - 6.5 d^(1/3) >= 2.1 d (the
  difference is positive at 10^5 and increasing for d > 155). Then U e^(1.4 d) <= d^2 e^(-1.41 d) + 2 d^2 e^(-0.0068 d),
  and both terms decrease for d > 2/0.0068 and are below e^(-600) at 10^5. Finally m = e^(5 d^(1/3)) > e^230, so
  e^(-m/3) is negligible and U + e^(-m/3) < 1. Check S7 evaluates the estimates at d = 10^5, 10^6, 10^8 and
  10^12 and checks each monotonicity step at 10^5. QED

**Certified table (FINITE-EXACT).** For each d, the union bound uses, for every pair type, a
maximum-weight spanning forest of the flow graph (Kruskal on the exact integer cell counts). The other
ingredients are the atom bounds of Lemma 2 (series enclosure for x < 60, bound (2) for x >= 60), the
antipodal halving, a rational q, and the Chernoff tail
P(Bin(2^d, q) >= M + 1) <= exp(-2^d D((M+1)/2^d || q)). All arithmetic is outward-rounded interval
arithmetic, and the table lists rigorous upper endpoints. Then edim_m(Q_d) <= M_d whenever U + tail < 1.

| d | q_d | M_d | U <= | tail <= | ln M_d / d^(1/3) |
|---|---|---:|---:|---:|---:|
| 11 | 222861/1024000 | 492 | 0.95349086 | 0.04382090 | 2.7871 |
| 12 | 494559/4096000 | 543 | 0.93396591 | 0.06497403 | 2.7505 |
| 13 | 5681/81920 | 622 | 0.93430088 | 0.06270534 | 2.7359 |
| 14 | 10403/256000 | 724 | 0.93068529 | 0.06928310 | 2.7321 |
| 15 | 155287/6553600 | 843 | 0.94643315 | 0.05342962 | 2.7317 |
| 16 | 90673/6553600 | 978 | 0.94144631 | 0.05802286 | 2.7325 |
| 17 | 1052691/131072000 | 1129 | 0.93884429 | 0.06113297 | 2.7337 |
| 18 | 1223837/262144000 | 1306 | 0.93725654 | 0.06223008 | 2.7377 |
| 19 | 1400923/524288000 | 1490 | 0.94151823 | 0.05822074 | 2.7382 |
| 20 | 1614059/1048576000 | 1708 | 0.93527653 | 0.06435201 | 2.7421 |
| 21 | 459133/524288000 | 1936 | 0.93286757 | 0.06709204 | 2.7432 |
| 22 | 2082059/4194304000 | 2193 | 0.94770116 | 0.05191233 | 2.7455 |
| 23 | 1176267/4194304000 | 2464 | 0.92644587 | 0.07084907 | 2.7461 |
| 24 | 331207/2097152000 | 2774 | 0.94593670 | 0.05395847 | 2.7485 |
| 25 | 1481233/16777216000 | 3095 | 0.94809234 | 0.05153052 | 2.7488 |
| 26 | 1656007/33554432000 | 3454 | 0.95230840 | 0.04768395 | 2.7501 |
| 27 | 1845761/67108864000 | 3833 | 0.93214052 | 0.06621092 | 2.7505 |
| 28 | 1/65536 | 4254 | 0.95187451 | 0.04750444 | 2.7516 |
| 29 | 22599/2684354560 | 4688 | 0.95607469 | 0.04379412 | 2.7513 |
| 30 | 5004947/1073741824000 | 5179 | 0.95137631 | 0.04849663 | 2.7524 |
| 31 | 5497383/2147483648000 | 5681 | 0.95288147 | 0.04661436 | 2.7519 |
| 32 | 6046957/4294967296000 | 6242 | 0.95632316 | 0.04310674 | 2.7526 |
| 33 | 52939/68719476736 | 6823 | 0.95852123 | 0.04105533 | 2.7523 |
| 34 | 7246457/17179869184000 | 7463 | 0.96036725 | 0.03943140 | 2.7527 |
| 35 | 7901081/34359738368000 | 8127 | 0.96033344 | 0.03964454 | 2.7523 |
| 36 | 8621929/68719476736000 | 8857 | 0.95866492 | 0.04065031 | 2.7526 |
| 37 | 9358381/137438953472000 | 9604 | 0.95997012 | 0.03990039 | 2.7519 |
| 38 | 203499/5497558138880 | 10431 | 0.95993254 | 0.03995793 | 2.7521 |
| 39 | 11017147/549755813888000 | 11283 | 0.95923865 | 0.04051745 | 2.7515 |
| 40 | 5969909/549755813888000 | 12216 | 0.95869994 | 0.04105807 | 2.7516 |
| 41 | 2576067/439804651110400 | 13172 | 0.96269144 | 0.03688155 | 2.7509 |
| 42 | 6956827/2199023255552000 | 14221 | 0.96591537 | 0.03365056 | 2.7510 |
| 43 | 14973433/8796093022208000 | 15294 | 0.96724905 | 0.03243543 | 2.7502 |
| 44 | 197/214748364800 | 16469 | 0.96555499 | 0.03381236 | 2.7502 |
| 45 | 17331657/35184372088832000 | 17674 | 0.96541431 | 0.03409503 | 2.7495 |
| 46 | 18622661/70368744177664000 | 18978 | 0.96595190 | 0.03378695 | 2.7493 |
| 47 | 4987297/35184372088832000 | 20317 | 0.96620906 | 0.03376261 | 2.7486 |
| 48 | 1069591/14073748835532800 | 21777 | 0.96820855 | 0.03127281 | 2.7485 |
| 49 | 5720211/140737488355328000 | 23268 | 0.96214135 | 0.03785760 | 2.7477 |
| 50 | 24485831/1125899906842624000 | 24884 | 0.96013215 | 0.03931705 | 2.7475 |
| 51 | 26096877/2251799813685248000 | 26536 | 0.97471244 | 0.02495601 | 2.7468 |
| 52 | 27884361/4503599627370496000 | 28320 | 0.96636642 | 0.03334185 | 2.7465 |
| 53 | 7429391/2251799813685248000 | 30150 | 0.95681340 | 0.04304046 | 2.7458 |
| 54 | 3988127/2251799813685248000 | 32262 | 0.86332502 | 0.13521689 | 2.7466 |
| 55 | 6754181/7205759403792793600 | 34186 | 0.92144304 | 0.07785459 | 2.7451 |
| 56 | 35812997/72057594037927936000 | 36308 | 0.96712996 | 0.03274437 | 2.7444 |
| 57 | 9586113/36028797018963968000 | 38726 | 0.85065821 | 0.14927675 | 2.7450 |
| 58 | 2017107/14411518807585587200 | 40911 | 0.98161153 | 0.01820588 | 2.7433 |
| 59 | 42805259/576460752303423488000 | 43347 | 0.96710224 | 0.03250573 | 2.7426 |
| 60 | 45414131/1152921504606846976000 | 45966 | 0.96465642 | 0.03502393 | 2.7423 |
| 61 | 48263011/2305843009213693952000 | 48740 | 0.90467709 | 0.09449701 | 2.7421 |
| 62 | 639839/57646075230342348800 | 51658 | 0.88476896 | 0.11435897 | 2.7420 |
| 63 | 13446593/2305843009213693952000 | 54415 | 0.97452018 | 0.02545243 | 2.7404 |
| 64 | 28397833/9223372036854775808000 | 57738 | 0.99958142 | 0.00041343 | 2.7409 |

Each M_d is chosen as small as the certificate allows, so some margins 1 - (U + tail) are tiny: 1.05e-6 at d = 49
and 5.2e-6 at d = 64. They are rigorous (outward-rounded upper endpoints; the blind audit re-certified all 54 rows
with its own code and a tighter enclosure), and nothing else depends on them: Theorem C needs only
M_d <= e^(2.79 d^(1/3)), and with a larger M the tail, and hence the margin, improves at once.

The certified bounds are smaller than the first note's double-precision values (Q_11 <= 511,
Q_16 <= 1056, Q_32 <= 6638). The union bound here also uses the antipodal halving of section 1.4, which
the first note did not.

## 5. Theorem D: Open Problem 4, density 1/2 versus sparse

**Theorem D.** Let S_q be a random set of density q in Q_d and m = q 2^d.
1. If ln m >= (kappa_+ + eps) d^(1/3), then P(S_q resolving) >= 1 - o(1), uniformly in q <= 1/2.
   By complementation the same holds for q >= 1/2 when ln((1-q) 2^d) >= (kappa_+ + eps) d^(1/3).
2. If ln m <= (kappa_- - eps) d^(1/3), then P(S_q resolving) -> 0.

*Proof.* (1) For fixed forests, every factor of the union bound decreases in x = N q(1-q), and q(1-q)
increases on (0, 1/2]. So for q_0 <= q <= 1/2, with q_0 = e^((kappa_+ + eps) d^(1/3)) 2^(-d), the bound at
q is at most the bound at q_0, which tends to 0 by the proof of Corollary A. Complementation:
H^(V - S)_e(r) = 2 C(d-1, r) - H^S_e(r), so S resolves iff its complement does (paper, Prop. 2).
(2) By Theorem B, every resolving set has at least M := e^((kappa_- - eps/2) d^(1/3)) points for large d.
By Markov's inequality, P(|S_q| >= M) <= m/M <= e^(-(eps/2) d^(1/3)). QED

**Quantitative comparison.** A density-1/2 set has about 2^(d-1) points, while the best sizes are
exp(Theta(d^(1/3))). For every d >= 6 the ratio 2^(d-1)/edim_m(Q_d) is at least 2^(d-1)/(2 e^(5 d^(1/3))).

| d | density 1/2 (about 2^(d-1)) | paper's certificate | certified sparse M_d | explicit set |
|---|---:|---:|---:|---:|
| 7 | 64 | 63 | - | 19 |
| 10 | 512 | 492 | - | 48 |
| 11 | 1024 | existence | 492 | 65 |
| 12 | 2048 | existence | 543 | 76 |
| 16 | 32768 | existence | 978 | - |
| 32 | 2.1e9 | existence | 6242 | - |
| 64 | 9.2e18 | existence | 57738 | - |

## 6. Open Problem 3: settled in the companion note
The paper proves existence for d >= 11 at q = 1/2, by exact computation for 11 <= d <= 50 and an
elementary estimate beyond. Theorem A is one analytic estimate. With lambda = 5 d^(1/3) it becomes effective at
d = 64 (check S6, first part). With lambda chosen for each d it is effective from d = 41: lambda = (d + 37)/4 works
for every 41 <= d <= 63 (check S6, second part; the largest value of U + e^(-m/3) there is 0.054, at d = 41). At
d = 40 no lambda on a grid of step 0.01 in [15, 22] works (the best value is 1.37, at lambda = 19.1). Its losses are the factor 1/(3(L+1)) of Lemma H, the factor 1/(L+1) of the method-of-types bound,
and the crude atom constant sqrt(pi/8). The parallel note behind THM-4534 goes further at density 1/2: a closed-form
estimate for every d >= 17 (its Theorem D) and, with its Fourier-Hoelder lemma, a certified U_10 <= 0.2649.

An earlier version of this section asked for a different idea, such as a lifting lemma from Q_d to Q_(d+1). The
companion note [`edge_multiset_dimension_op3_20261002.md`](edge_multiset_dimension_op3_20261002.md) supplies one.
It does not lift landmark sets
(the naive lift S x {0} fails: after two lifts the set is fixed by swapping the two new coordinates, so L1 excludes
it). It lifts the union bound itself: each pair type of Q_(d+1) inherits a potential forest from a pair type of Q_d
plus one extra near-central edge, so the paper's union bound at density 1/2 with antipodal halving satisfies
V_(d+1) <= V_d/2 + A_(d+1). Started from the exact value V_10 = 0.66023, this gives V_d <= 0.6767 for every d >= 10,
which settles Open Problem 3 from d = 10 (PROVED there, computer-assisted at its base; audit owed).

## 7. Q_7
**Theorem 7.1 (FINITE-EXACT, computer-assisted; blind-audited 2026-10-02).** Q_7 has no edge-multiset resolving set of any
size k <= 12. Hence 13 <= edim_m(Q_7) <= 19; the previous bounds were 8 <= edim_m(Q_7) <= 19 (THM-4525).

Full write-up: [`edge_multiset_dimension_q7_20261002.md`](edge_multiset_dimension_q7_20261002.md). Code and run
records: `04-computation/edge_multiset_dimension_q7_20261002/`.
* **Method.** Every k-set is Aut(Q_7)-equivalent to A x {0} u B x {1}, where the last coordinate has the least
  imbalance a - b >= 0 and A is one of the Aut(Q_6)-orbit representatives of a-subsets (Lemma 1 there, proved).
  For every (k, a), every representative A and every admissible B is tested exactly. Nothing is pruned on
  partial collisions, because resolvability is not monotone. A weight lemma (Lemma 2 there) shows that the pairs
  with a >= 11 hold no resolving set (in the k <= 12 run it removed (12, 12); the a = 11 cases were searched anyway
  and have no leaves). The representatives come from orderly generation and match the Burnside counts 1, 1, 6, 16,
  103, 497, 3253, 19735, 120843, 681474, 3561696, 16938566 (a = 0..11).
* **Size.** 303,583,126,680 leaves for k <= 12, in 73 CPU minutes. For every (k, a) the leaf count equals an
  independent dynamic-programming count of the domain.
* **Validation.** The same code reproduces the Q_6 result of THM-4525 end to end: nothing for k <= 14, and for
  k = 15 exactly the 229 deposited orbits, with both search engines. A verification mode re-tested 154 million
  Q_7 leaves from scratch with no discrepancy. The C and pure-Python checkers agree on 1,636 sets.
* **Annealing** (EMPIRICAL) found a second resolving 19-set of Q_7, {2, 4, 21, 22, 32, 38, 42, 44, 47, 54, 70,
  79, 84, 90, 110, 114, 116, 120, 126}, inequivalent to the known one. Sizes 17 and 18 were not tried.
* **k = 13** (1.51e12 leaves, about 6 CPU hours) was still running when this note was committed and is not
  claimed. k = 14 would take about 55 CPU hours by the same method.
* Runner section S10 recomputes the two 19-sets and the orbit counts from the definitions and audits the stored
  records; with `--q7` it rebuilds the search and re-runs Q_6 (k <= 13) and Q_7 (k <= 10).

## 8. Where the truth may lie (heuristic, not claimed)
* Random sets cannot do much better than kappa_+. A typical pair of edges collides with probability
  about the product of its central atoms, so the expected number of colliding pairs is about U, and
  U -> infinity for ln m <= (kappa_+ - eps) d^(1/3). A second-moment argument would turn this into a
  threshold at kappa_+ for random sets; it was not attempted. (procgen_edim2_20261001 Theorem C proves the
  related statement that no choice of forests makes the union bound work below kappa_+.)
* The entropy bound treats each level count as if it could take about mu_r values; a random set only
  spreads it over about sqrt(mu_r) values. And the union bound pays for pairs (2 d ln 2) rather than
  edges (d ln 2). These two factors of 2 are exactly the ratio kappa_+/kappa_- = 4^(2/3). Whether a
  structured landmark set can beat random sets is the remaining question of Open Problem 2.

## 9. Open
* The constant: does lim ln edim_m(Q_d)/d^(1/3) exist, and where in [0.8146, 2.0526] is it?
* Open Problem 3 is settled from d = 10 in the companion note. What remains open there is a version without the
  finite computations at its base (the exact V_10 and the ratio check for 10 <= d <= 40).
* Q_7: 13 <= edim_m(Q_7) <= 19. Is a resolving set of size 13 to 18 possible? (k = 13 and a first annealing run
  at size 18 were running at commit time.)

## 10. Reproduction
`python3 -u 04-computation/experiments/edge_multiset_dimension_growth_20261002_run.py` prints its results to stdout
(deposited as `05-knowledge/results/edge_multiset_dimension_growth_20261002.out`) and section timings to stderr,
and ends with ALL CHECKS PASSED (470144 checks). It takes about 4 to 5 minutes on one core (231 s on the
original machine, 295 s in the audit re-run), mostly S5 and S6.
Needs numpy, scipy and mpmath. With `--q7`, S10 also builds the Q_7 search (gcc) and re-runs the fast subset in a
temporary directory: 161 s more with 2 of 4 shared cores, 455 checks in S10.
