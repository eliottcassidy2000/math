# Busch's staircase, Fermat-block cuts, and where 7 and 21 live

mac-mini-2026-11-13-cuts · **Chan–Pak OPEN.** Owner prompt (2026-11-13): show that an arbitrary prescribed odd
target admits one of the cuts by thinking of it as something other than an arbitrary odd number (it is never 7 or
21); possibly prove that a further impossible value would be isolated or part of an infinite family; relate Busch's
bound to Neményi's theorem and the Fermat primes.

Status of what follows: **PROVED** (§3 hub-port identity, §4 transfer at the new roots, the port sets Λ(1), Λ(3));
**FINITE-EXACT** (§1 staircase through order 9, §2 the 398,131 horizon, §3 coverage); **CONJECTURE** (Staircase-Run,
Shallow-Cut). Scripts: [`scripts/`](scripts/README_SCRIPTS.md).

Inherited: the gluing law N = H(A) + 2^j L(A,v) is Session 324's
[ordered Busch overlap](../session324_busch_inverse_and_gap_structure_20261113/BUSCH_INVERSE.md) §3, and its closure
countermodels are in [GAP_STRUCTURE](../session324_busch_inverse_and_gap_structure_20261113/GAP_STRUCTURE.md). Also
used: the strong inventory below 625 in THM-4744, the one-centre family of THM-4876 / HYP-9469, the 2-adic lines of
THM-5126 / HYP-9658, and the threshold-95 ternary rule in
[decomposition16](../decomposition16_odd_spectrum_20261110/induction/PROOF.md). These notes were mined while
deleting ten stale worktrees.

## 0. The picture in one paragraph

Busch's floor b(n) = 3, 5, 9, 15, 25, 45, 75, 125, … consists of products of the first two Fermat primes,
F0 = 3 and F1 = 5. Every strong realizer of N has order at most n(N), the largest n with b(n) ≤ N. So each odd
target lives on one *step* [b(n), b(n+1)) of a staircase, and only strong tournaments of order ≤ n can serve it
directly. Through order 9 every step is filled, with exactly two gaps: **7 in the step [F1, F0²) = [5, 9)** and
**21 in the step [F0·F1, F1²) = [15, 25)**. These are the two steps where the strong spectrum is still smaller than
the step: order 4 has the single count F1, and order 6 has 13 counts. From order 7 on, each order's strong spectrum
contains a long unbroken run of consecutive odd counts, and the runs overlap from one order to the next. On the
constructive side, a target is better thought of as a *2-adic lift* N = h + 2^j·L of a smaller count h than as an
arbitrary odd number. Here L is a "port" of some realizer of h, read off at one root vertex. Every odd target from
69 to 398,129 has such a lift with a source on at most 9 vertices.

## 1. The Busch staircase (FINITE-EXACT through order 9)

Busch (EJC 13, 2006, Theorem 1) gives b(n) = 3β^(n−3), β^(n−1), 9β^(n−5) for n ≡ 0, 1, 2 (mod 3), where β = 5^(1/3).
Thus b(n) ∈ {3·5^k, 5^k, 9·5^k}. The floor is attained by zero-defect Fermat-block chains (Moon's interval chains,
THM-4744): gluing R_j at a root of defect 0 multiplies H by 2^j + 1 and adds j + 1 vertices. The most efficient
block is R_2, with count F1 = 5 for 3 new vertices. This is Busch's base, since
5^(1/3) < 3^(1/2) = 9^(1/4) < 17^(1/5) < 33^(1/6) < …

`staircase.py` (all tournaments of orders 3–9 from `gentourng`; Sp(n) = strong counts at order exactly n):

| n | step [b(n), b(n+1)) | odd targets | missed at order exactly n | missed by orders ≤ n | longest odd run in Sp(n) | gaps of Sp(n) below its run |
|---|---|---:|---|---|---|---|
| 3 | [3, 5) | 1 | – | – | [3, 3] | – |
| 4 | [5, 9) | 2 | 7 | **7** | [5, 5] | – |
| 5 | [9, 15) | 3 | – | – | [9, 15] | – |
| 6 | [15, 25) | 5 | 21 | **21** | [23, 33] | 21 |
| 7 | [25, 45) | 10 | – | – | [65, 105] | 63 |
| 8 | [45, 75) | 15 | 47, 55, 67 | – | [69, 609] | 47, 55, 67 |
| 9 | [75, 125) | 25 | 77, 79, 83 | – | [85, 2881] | 77, 79, 83 |

The runs chain: 105 > 69 and 609 > 85, and run(9) reaches 2881, past b(15) = 1875. The gaps below each run sit in
a thin boundary layer just above b(n), and the previous order covers them.

**Staircase-Run Conjecture.** For every n ≥ 8, Sp(n) contains every odd number in [ℓ(n), r(n)], with
ℓ(n) − b(n) = o(b(n)) and r(n) ≥ b(n+2). Together with THM-4744's inventory, it implies the strong form of Chan–Pak:
every odd N ∉ {1, 7, 21} is the count of a strong tournament. The route to a proof is to induct on runs at fixed
order: this is the decomposition16 threshold-95 rule with order bookkeeping. Suppose every count h in run(n) extends
by one vertex to 3h+2 and 3h+4, and 3h is supplied. Then order n+1 contains every odd number in
[3ℓ(n), 3r(n)+4]. The boundary layer [b(n+1), 3ℓ(n)) lies inside run(n) whenever r(n) ≥ 3ℓ(n), which holds at
n = 8 and 9 (609 ≥ 207 and 2881 ≥ 255) and is preserved by the step. So the conjecture reduces to the ternary rule on
runs.

**Fermat remark (elementary, PROVED).** 7 is the least prime dividing no 2^j + 1: ord_7(2) = 3 is odd, so −1 is
not a power of 2 mod 7. Hence 7 divides no Fermat-block count 2^j + 1, no Busch floor value, and no product of
block counts. It also misses every ray h(1 + 2^s) of a source whose ports are 2-power multiples of its count; §2
shows that sources of count 1 and 3 have only such ports. The order-4 step contains exactly one non-Fermat odd
number, 7 = Φ3(2). The step [15, 25) loses 21 = Φ3(4) = 3·7. Primes with ord_p(2) odd (7, 23, 31, 47, 71, 73, 79,
89, …; density 7/24 by Hasse) are the natural places for Fermat-built constructions to thin out. 47 is the extra
miss of the span-2 ribbons (THM-4880).

**Neményi.** Session 324 reads the request as P. F. Neményi's inverse problem for plane stress trajectories: the
trajectories fix the stress field up to five real constants (Neményi–Sáenz, *Acta Mech.* 208, 2009). In the
staircase picture Busch's theorem plays the same role for H. The observation N fixes the strong order to at most
n(N), so the strong inverse fibre is finite. The rate of the staircase is F1^(1/3); the coincidence of the "5" in
the five-constant theorem with the Fermat base 5 is an analogy, not a theorem. See Session 324 for the exact
five-coordinate padding fibre.

## 2. An arbitrary target as a 2-adic lift: Fermat-block cuts

Glue the near-transitive block R_j (order j+2, count 2^j+1) at a root v of any tournament A, with A − v beating the
block. The count is **N = H(A) + 2^j·L(A,v)** (Session 324). Here L(A,v) counts ordered pairs (P,Q) of disjoint
paths covering A − v, either possibly empty, with Q empty or v → first(Q). Call N a *cut* of the source (A, v).
So N ≡ h (mod 2^j) with h = H(A) ≤ N/3: the low binary digits of the target are those of a smaller count, and
(N − h)/2^j must be one of h's ports. Write Λ(h) = {L(A,v) : H(A) = h} for the port set.

**FINITE-EXACT (`cover_cuts.py`).** Sources of order ≤ 9 give 1,273,591 distinct pairs (h, L). Every odd N with
69 ≤ N ≤ 398,129 is a cut of one of them. The odd non-cuts below 69 are exactly 1, 7, 11, 13, 21, 47, 67, and of
these all but 7 and 21 are counts of small tournaments. The first miss above 67 is **398,131**. Beyond it the misses
grow as counting predicts: 1,903,018 odd N ≤ 10^7 are not order-≤9 cuts. So a target is one Fermat cut away from a
tournament on at most nine vertices far beyond any scale those sources were designed for.

**Port sets of the smallest counts (PROVED).** Λ(1) = {2^t} and Λ(3) = {3·2^t}. Realizers of 1 are transitive, and
a vertex at position i of TT_k has port 2^(k−i). Realizers of 3 are C3 plus padding. With b right-padding vertices,
each C3 vertex has port 3·2^b. Left-padding roots see q(C3)·2^t = 12·2^t, and right-padding roots see 3·2^t. The
products obey L(A1⇒A2, v∈A2) = H(A1)·L(A2,v) and L(A1⇒A2, v∈A1) = L(A1,v)·q(A2), where q counts ordered
2-path covers, empty paths allowed. Realizers of 5 are R_2 plus padding, since R_2 is the only strong 4-tournament.
Its core ports are 5, 6, 7, 9 times 2^b, and its padding ports are 5·2^t or at least q(R_2) ≥ 25. Hence
Λ(5) ∩ [5, 16] = {5, 6, 7, 9, 10, 12, 14}.
Reading the exceptions: 7 − 1 = 2·3 needs an odd port 3 at h = 1; 7 − 3 = 4 and 7 − 5 = 2 are below the ports.
21 − 1 = 4·5 needs an odd port at h = 1, 21 − 3 = 2·9 needs port 9 at h = 3, and 21 − 5 = 2·8 or 16·1 needs
port 8 or 16 at h = 5. All fail, and h = 7 does not exist. Session 324 §5 reads the same exclusion in defect
coordinates.

**Shallow-Cut Conjecture.** Every odd N ≥ 69 is a Fermat-block cut H(A) + 2^j L(A,v). A strengthened form is that
the order of A can be taken O(log N). It is not a reformulation of Chan–Pak: it asserts that one specific gadget
shape suffices.

## 3. The target as a binary word with a cut position (hub words)

Let A_x be TT_m plus a hub u with u → i iff x_i = 1. Then H = 1 + 2W(x) with W the weighted inversion count of
THM-4876 (density 1/2, HYP-9469).

**PROVED (hub port = binary value).** L(A_x, u) = 1 + X(x), where X(x) is the binary number x_1x_2…x_m. Every
ordered 2-partition of the spine is one cover; Q = ∅ gives 1, and min Q = i with x_i = 1 gives 2^(m−i). Gluing R_j
at u returns a hub tournament with word x·1^j·0, hence the identity

    H(x 1^j 0) = H(x) + 2^j (1 + X(x)),   equivalently   W(x 1^j 0) = W(x) + 2^(j−1)(1 + X(x)).

So gluing at the hub never leaves the half-dense hub family. All new values come from cutting at a *spine* vertex
k. There L(x,k) is computed by a 12-state transfer over the spine (u in P with slots Pa/Pb/Q, or u in Q with
slots P/Qa/Qb). `hubL.py` matches brute force on all 9,218 spine roots with m ≤ 10.

**FINITE-EXACT (`hubcut`, m ≤ 20).** The hub-cut family {1 + 2W(x) + 2^j L(x,k)} ∪ {1 + 2W(x)} misses only 1,843 of
the 100,000 odd N ≤ 2·10^5 (1.8%). Hub words alone miss 50%. The first misses are
7, 13, 21, 43, 67, 83, 231, 233, 313, 317, 327, 373, 379, 533, 631, 667, … Composite misses are mostly filled
by products (231 = 3·77, 327 = 3·109, 533 = 13·41). The first prime misses above 100 are 233, 313, 317, 373, 379.
A second cut, or a second hub, is the natural completion to try. This is an explicit two-parameter numeration:
target = (binary word, cut position, 2-adic shift).

## 4. Transfer at the new roots (PROVED; `transfer_check.py`, 300 random checks)

For a rooted (A, v), track H, the defect D = L − H, h = H(A − v), E = end_v(A) and S = start_v(A). Glue R_1
(new vertices r1 → r2, v → r1, r2 → v, A − v beats both). The new count is 3H + 2D at every root, and:

| new root | D′ | h′ | E′ | S′ |
|---|---|---|---|---|
| v | H − h + 2D | h | h | S |
| r1 | H − E | 2H + D − E | H + D | S |
| r2 | 0 | H | H | S |

S is invariant at all three roots, and the system closes linearly on (H, D, h, E, S). This is the j = 1 instance of
the changed-port compilers (endpoint_context_recursion has a four-coordinate version). Restricted to these three
roots, there are three choices per step against growth at least 3, with equality only along zero-defect roots, so
this sublanguage alone should not be expected to be complete. Codex proved a nine-rule analogue sparse. A complete
scheme must also cut at old roots, whose ports need finer statistics.

## 5. The dichotomy question: what is and is not provable

- Not provable from closure alone. Session 324's countermodels delete any divisor-downset of 7m from the odd
  integers and stay multiplicatively closed, so ordinal products permit several finite extra holes. Fermat-block
  rays do not repair this abstractly: a ray is a property of a concrete rooted realizer, and an abstract set can
  avoid sparse rays.
- Staircase form (exact). An extra exception m is a value missed by every strong tournament of order ≤ n(m), and
  m is not a product of attained values. THM-4875 (D) and THM-4744 put any such m above 2^32, so n(m) ≥ 42.
- Heuristic dichotomy. A step of index n has about 0.36·b(n) odd targets but about 2^(n²/2)/n! strong tournaments of
  order n, whose counts are already observed to fill long runs. A sporadic, random-like gap therefore has
  super-exponentially small expectation beyond the first few steps, and those steps are exhausted (7 and 21 are
  their only gaps). So a further exception could not be sporadic. It would need an arithmetic cause acting on a
  whole family: a congruence class, a multiplicative ray 7p^j (Session 324 shows this is exactly what an infinite
  set of 7p-holes must look like), a 2-adic line 1 + c·2^s (THM-5126 / HYP-9658), or the multiples of a prime with
  ord_p(2) odd. In that sense a third impossible value would be part of an infinite family, or would not exist.
  Every candidate family checked so far is complete beyond 7 and 21.
- Conditional theorem. The Staircase-Run Conjecture implies that the exceptional set is exactly {7, 21}.

## 6. Next moves

1. Prove run chaining at order n+1 from a run at order n, using ternary extensions on a positive proportion of the
   run (§1).
2. Close the hub-cut numeration with one more cut or a second hub, and attack its completeness as a digit problem
   in (x, k, j) (§3).
3. Port-wedge induction. Use the §4 transfer, together with old-root ports, to show that the set of pairs (h, L) at
   order m is dense enough in a cone that its rays cover [69, X(m)] with X(m) → ∞ (§2).
