# Generation beyond order 11, and the spanning-tree twin of the 7/21 proof

THM-4985 and THM-4986 (with HYP-9565). Session mac-mini-2026-11-09-genesis.

Code and outputs are in
[`04-computation/experiments/genesis_20261109/`](../../../04-computation/experiments/genesis_20261109/).
Pre-registration: [PREREGISTRATION.md](../../../04-computation/experiments/genesis_20261109/PREREGISTRATION.md), commit
fb31d2c38a, pushed before `sa13.c` and `tau_share.c` ran.

## 0. Definitions and inherited facts

- **H(T)** is the number of Hamiltonian paths. **S(n)** is the set of values H(T) over n-tournaments.
- S(n - 1) is contained in S(n): pad with a sink (THM-1370).
- klein's exact censuses (`05-knowledge/results/path_spectrum_two_angles_metagraph_20261106/data`) give the following.

  | order | contents of S(n) |
  | --- | --- |
  | 10 | every odd h <= 15,745 except 7, 21 and 208 values, all at least 14,685 |
  | 11 | every odd h <= 95,095 except 7, 21 and 3,635 values, all at least 80,453 |
  | <= 9 | read from `census_n.txt` |
- A **first creator** (FC) of order n is a class with H not in S(n - 1) (THM-4922, THM-4942).
- **lambda(T)** is arc-strong connectivity: the fewest arc reversals that make T non-strong.
- A **module** is a vertex set M, with 2 <= |M| <= n - 1, such that every outside vertex beats all of M or loses to
  all of M. T is **prime** if it has no module, and T is a substitution iff it is not prime.
- **THM-4942's five-operation grammar** (kind-pasteur S17, finite-exact and exhaustive through order 11). Every class
  is one of the following:
  - **FC**: a first creator;
  - **FCV**: a first creator plus a vertex, i.e. some T - v is a first creator;
  - **SUM**: non-strong;
  - **SUB**: non-prime;
  - **SUM2**: lambda <= 2. Being within r reversals of a non-strong tournament is the same as lambda <= r.
- An **escaper** is a tournament in none of the five classes. Equivalently, an escaper is:
  - strong and prime;
  - lambda >= 3;
  - H(T) in S(n - 1);
  - H(T - v) in S(n - 2) for every v.

  Escapers are exactly the lambda >= 3 part of S17's REST category (THM-4942 PROOF §3.1: strong, prime, not a first
  creator, no first-creator deletion). HYP-9532 (F) asked whether escapers exist at any order.
- **Defects of T**:
  - **module defect**: the fewest reversals creating a module. It equals the minimum over M of
    sum over x not in M of min(|N+(x) & M|, |N-(x) & M|). Reversing arcs inside M or outside M never helps.
  - **twin distance**: the module defect restricted to |M| = 2, i.e. the fewest distinguishers of a pair.
  - **decomposability defect**: dec(T) = min(lambda(T), module defect), the fewest reversals reaching a
    decomposable (non-strong or non-prime) tournament. For n >= 3 the minimum is redundant (audit G1, m1). A
    non-strong tournament is non-prime, because the side of size >= 2 of a dominating cut is a module. So the module
    defect never exceeds lambda, and dec equals the module defect.

## 1. Order 12: the grammar fails (FINITE-EXACT)

### 1.1 Extensions of the seven 3-arc-strong order-11 non-creators (`ext12.py`)

S17 listed the seven classes of 3-arc-strong 11-tournaments with H in S(10) (`lambda3_classes_order11.txt`): H = 11,475,
15,213 (2 classes), 15,273 (2) and 15,561 (2). All seven are substitutions.

Extending each by a 12th vertex in all 2^11 ways gives 14,336 tournaments. These were classified with exact H for T and
every deletion (`tprop.c`), the exact spectra S(10) and S(11), lambda by cut enumeration, and primality by module
enumeration.

- Exactly 36 of the extensions are escapers.
- They form **12 isomorphism classes** (`canon.c`, 128-bit individualization-refinement), each found 3 times.
- The classes are 6 converse pairs, with no self-converse class.
- Each converse pair has one H value: 48,997, 49,215, 49,935, 50,221, 51,451 or 52,303.
- All 12 are rigid (|Aut| = 1), have lambda = 3 and twin distance 1.

An independent pure-Python check (`verify_escapers.py`) reconfirms all 12 classes. It shares no code with the C tools:
- a separate subset DP for H;
- S(10) and S(11) read from the census T-lines, not from the hole lists;
- lambda by enumerating all 4,094 nonempty proper vertex subsets;
- primality by module closure of every pair.

Since THM-4942 (3) is exhaustive through order 11, **order 12 is the first order at which the five-operation grammar
fails**. This answers HYP-9532 (F). It also gives k*(12) >= 3, where k*(n) is the largest lambda of a non-creator
(HYP-9532 (D)).

### 1.2 Annealing (`sa12.c`; EXPLORATORY search, exact verification of every hit)

The annealer works over 12-tournaments with minimum semidegree >= 3. At n = 12 <= 4·3 + 1 this is equivalent to
lambda >= 3 (THM-4811 (B)).
- Its energy is zero exactly when H(T) is in S(11) and every H(T - v) is in S(10).
- Moves are single arc reversals that keep minimum semidegree >= 3.
- Run A (seed 11) has no further constraint. Run B (seed 22) requires every vertex pair to have at least 2
  distinguishers.
- Both runs use 1.5M steps.

Every hit is re-checked by `analyze_escapers.py`: canonical forms, re-classification with `grammar.py`, lambda, and
module defect by enumeration of the 4,082 vertex sets M with 2 <= |M| <= 11.

Run A recorded 240 hits and run B recorded 7 (each hit ends a restart). After canonical forms:
- **181 escaper classes**, none among the 12 classes of 1.1, for **193 classes in all**
  (`escapers12_all_classes.txt`).
- H runs from 32,897 to 58,017: 150 distinct values, 43 of which carry two classes.
- 17 classes are self-converse, and 99 classes had their converse class found too.
- The (lambda, module defect, twin distance) triples:

  | (lambda, module defect, twin distance) | classes |
  | --- | --- |
  | (3, 1, 1) | 190 |
  | (3, 2, 2) | 1 (H = 57,643) |
  | (3, 3, 3) | 2 (H = 48,537 and 48,849) |

  The two classes with module defect 3 are self-converse with |Aut| = 3 (`canon.c`; audit G1 agrees).
  Each has scores 3^6 8^6: six vertices beat 8 others and
  send only 3 arcs back to the six vertices that beat 3 others. So each is an ordinal sum of two 6-tournaments
  re-connected by exactly 3 arcs.
- **Creator depth** is the least number of vertex deletions reaching a first creator; it is >= 2 for escapers by
  definition (`depth.py`).
  - Depth is 2 for 192 classes: they contain a first-creator 10-subtournament.
  - Depth is 3 for one class: H = 33,021, scores 3^6 8^6, with 54 first-creator 9-subtournaments and no order-10 one.

**How the escapers are built from first creators** (exact per class, over the 193 classes found;
`cuts.py`, `fcpairs.py`, `twinquot.py`, `subfc.py`):

| description | classes |
| --- | --- |
| **re-connected sum of two first creators**: a cut with exactly 3 back arcs whose two sides are both first creators | 106 |
| **substitution of first creators plus one vertex**: some T - v is a strong substitution whose quotient over the maximal modules and every module are first creators, like S17's parent C3[R5, R5, 1] (3-cycle, two copies of R5, one vertex) | 47 |
| both | 13 |
| neither (all of creator depth 2, i.e. an order-10 first creator plus two vertices) | 53 |

- In all 106 re-connected sums the sides are two order-6 first creators, with H in {33, 37, 41, 43, 45} (the largest
  values new at order 6). The most frequent pair is (43, 43), in 15 classes.
- These 106 are exactly the classes with scores 3^6 8^6 (audit G1, m3). Their halves have minimum out-degree >= 2, or
  minimum in-degree >= 2. Every 6-tournament with minimum out-degree >= 2 has H in {33, 37, 41, 43, 45}, except the
  source-plus-R5 tournament with H = 15. So "both halves are first creators", and the value set, follow almost
  automatically from the score sequence. The census describes these escapers; it does not explain them.
- The two defect-3 classes are re-connected sums of two first creators of 45, with the 3 back arcs forming a matching.
  The depth-3 class is a re-connected sum of two first creators of 33.
- "An order-11 first creator with one vertex doubled into a transitive pair and one arc reversed" cannot describe an
  escaper (0 classes, as it must be). Its quotient is a vertex deletion of T, and no vertex deletion of an escaper is a
  first creator (audit G1, m2).

Every escaper found has lambda = 3 exactly. In general a non-creator of order n has lambda <= k*(n). So the grammar
with "ordinal sum re-connected by at most k*(n) arcs" is trivially complete at order n, and HYP-9532 (F)'s question
"how many re-connecting arcs?" is the question of how k*(n) grows (section 4).

### 1.3 What the order-12 escapers say about repairs

- The 12 classes of 1.1, and most annealing classes, are **one arc reversal from a substitution**: some pair has a
  single distinguisher. Adding "a substitution re-connected by one arc" to the grammar covers them.
- The two annealing classes with **decomposability defect 3** are not within two reversals of any non-strong or
  non-prime tournament. So no grammar of the form {FC, FCV, everything within 2 reversals of a decomposable tournament}
  covers order 12. They are within three reversals: lambda = 3.
- **Every escaper found has creator depth 2 or 3, and lambda = 3.** So each found escaper is covered by one of two
  repairs:
  - "a first creator plus at most two vertices" covers 192 of the 193 classes;
  - "an ordinal sum re-connected by at most three arcs" covers all 193.

  Neither repair survives to all orders under (S-inf) (section 4): the bands need unboundedly many added vertices and
  re-connecting arcs.
- Every escaper has a deletion with lambda <= 2. This is automatic: delete an in-neighbour of a vertex of in-degree 3,
  or an out-neighbour of a vertex of out-degree 3. At n = 12, lambda = 3 means minimum semidegree 3 (THM-4811 (B)).
  It carries no information.

## 2. Band tournaments (PROVED)

**Definition.** band_s(n) has vertices 0, ..., n - 1. For i < j the arc is j -> i when j - i <= s (a *back arc*), and
i -> j otherwise. band_s(2s + 1) is the rotational tournament. The reflection i -> n - 1 - i is an isomorphism from
band_s(n) to its converse.

**Credit (audit G1, M2).** band_s(n) is THM-4880's ribbon tournament R(a^n) with the constant letter a = {1, ..., s}
(mac-mini-2026-11-04-ribbon). THM-4880 §F already gives H(band_2(n)) (letter 12, with its recurrence) and H(band_3(n))
for n <= 10 (letter 123: 1, 1, 1, 1, 9, 41, 175, 629, 2121, 6889, continued here by 22,143, 72,073, 237,817). The proof of
Lemma B3 is THM-4880 Theorem A's window argument made explicit. New here: Lemmas B1 and B2, the escaper certificates,
and the (S-inf) theorem.

**Distinguishers.** Fix u < w and x not in {u, w}. In band_s(n), x -> u iff either x > u and x - u <= s, or x < u and
u - x >= s + 1. Hence:

| position of x | x distinguishes u and w exactly when |
| --- | --- |
| x > w | x - w <= s < x - u |
| x < u | u - x <= s < w - x |
| u < x < w | both distances x - u and w - x are <= s, or both are >= s + 1 |

**Lemma W (two-step walks).** Let gcd(a, b) = 1 and n >= a + b. The digraph on {0, ..., n - 1} with arcs i -> i + a
and i -> i - b (whenever in range) is strongly connected.

*Proof.* For n = a + b every vertex has exactly one out-arc. The out-map is i -> i + a (mod a + b), which is one cycle
because gcd(a, a + b) = 1. So the digraph is a Hamiltonian cycle. For n > a + b, the digraph restricted to
{0, ..., n - 2} is the digraph for n - 1, strongly connected by induction. Vertex n - 1 has an out-arc to n - 1 - b and
an in-arc from n - 1 - a, both inside {0, ..., n - 2}. ∎

`band_lemmas.py` checks Lemma W for all coprime a, b <= 25 and n <= 80.

**Lemma B1 (primality).** band_s(n) is prime for every n >= 2s + 3. band_s(2s + 2) is not prime: {0, 2s + 1} is a
module.

*Proof.* Claims 1 and 2 below give the result for n >= 2s + 3: by Claim 2 every module contains a consecutive pair, and
by Claim 1 it is then all of V.

*Claim 1.* If n >= 2s + 2, a module M containing two consecutive vertices u, u + 1 is all of V.

- First, an interval [p, q] in M with q - p >= s forces M = V:
  - if q < n - 1, then x = q + 1 distinguishes q - s and q (x - q = 1 <= s < s + 1), so q + 1 is in M;
  - symmetrically, if p > 0, then x = p - 1 distinguishes p and p + s, so p - 1 is in M.
- Now D(u, u + 1) = {u + s + 1, u - s}, restricted to V. Two cases:
  - **u + s + 1 <= n - 1.** Then u + s + 1 is in M. The pair (u + 1, u + s + 1) has span s, so every vertex strictly
    between them distinguishes it. Hence the interval [u, u + s + 1] lies in M.
  - **u + s + 1 > n - 1.** Then u >= n - s - 1 >= s + 1, so u - s >= 1 is in M. Every vertex strictly between u - s and
    u distinguishes that pair. Hence [u - s, u + 1] lies in M.

  In both cases the first bullet applies.

*Claim 2.* If n >= 2s + 3, every module M contains two consecutive vertices. Take u < w in M with d = w - u minimal, and
suppose d >= 2.
- **2 <= d <= 2s.** Some x with u < x < w has x - u <= s and w - x <= s: the interval
  [max(u + 1, w - s), min(w - 1, u + s)] is nonempty when d <= 2s. That x is in M and is closer to u, a contradiction.
- **d >= 2s + 2.** x = u + s + 1 has both distances >= s + 1, so x is in M and closer, a contradiction.
- **d = 2s + 1.** If w < n - 1, then w + 1 distinguishes u and w (1 <= s < 2s + 2), so w, w + 1 are in M. If u > 0, then
  u - 1, u are in M likewise. Otherwise u = 0 and w = n - 1 = 2s + 1, which means n = 2s + 2.

*Non-primality at n = 2s + 2.* For {0, 2s + 1}, every x strictly between has exactly one of the two distances <= s, so
no x distinguishes the pair. ∎

**Lemma B2 (arc-strong connectivity).** lambda(band_s(n)) = s for every n >= 2s + 1.

*Proof.* The in-degree of i is min(s, n - 1 - i) + max(0, i - s), and its out-degree is min(s, i) + max(0, n - 1 - i - s).
Both are >= s when n >= 2s + 1. Vertex 0 has in-degree exactly s, so lambda <= s.

For the lower bound there are two ranges:
- **n <= 4s + 1.** Minimum semidegree s gives lambda >= s by THM-4811 (B) (re-proved in Session230 §5).
- **n >= 3s.** Let sigma be a bijection from {1, ..., s} to {s + 1, ..., 2s} with gcd(r, sigma(r)) = 1. Such a bijection
  exists by the coprime-mapping theorem. This case, {1..s} -> {s + 1..2s}, is due to Daykin and Baines (Mathematika
  10 (1963) 132-136); the general theorem is Pomerance and Selfridge (*Proof of D. J. Newman's coprime mapping
  conjecture*, Mathematika 27 (1980) 69-83). Both CITED.
  - Put D_r = {i -> i - r} ∪ {i -> i + sigma(r)} for r = 1, ..., s.
  - These are back arcs of span r <= s and forward arcs of span sigma(r) >= s + 1, so each D_r lies in band_s(n), and
    different r use different spans, so the D_r are arc-disjoint.
  - Each D_r is spanning and strongly connected by Lemma W, since n >= 3s >= r + sigma(r).
  - Every cut is crossed by each D_r in each direction, so lambda >= s.

The two ranges cover every n >= 2s + 1, because 3s <= 4s + 1.

Without the citation, sigma(r) = rs + 1 works for n >= s^2 + s + 1. The (S-inf) theorem of section 4 needs only large
n, so it uses only this elementary range. ∎

`band_lemmas.py` checks lambda = s by max-flow, primality exactly for n >= 2s + 3, and that each D_r with
sigma(r) = rs + 1 is strongly connected, for s <= 6 and n <= 30.

**Lemma B3 (bounded back-span).** Suppose some vertex order of T puts every back arc (an arc from a later vertex to an
earlier one) at span <= s. Then H(T) <= (s + 3)^(2n).

*Proof.* Let P be a Hamiltonian path. Process the vertices in the order 0, 1, ..., n - 1.

- After vertex v - 1, the vertices <= v - 1 form maximal runs of consecutive vertices of P.
- The head h of any run, other than the run containing the first vertex of P, has a predecessor >= v in P. The arc
  from that predecessor to h is a back arc, so h >= v - s.
- So there are at most s + 1 runs.
- Vertex v either has no predecessor, or its predecessor is the tail of one of these runs, or its predecessor is a
  vertex > v. That is at most s + 3 options. The same holds for its successor.
- The sequence of choices determines P.

Hence there are at most (s + 3)^(2n) paths. ∎

## 3. Certified band escapers (FINITE-EXACT, witness-certified)

Membership h in S(m) for m >= 12 is certified by an explicit m-tournament with exactly h paths.
- `wit.c` finds such tournaments with the attachment form of a random base (HYP-2879 / THM-4097 / THM-4114, restated
  in THM-4942 (1)), and re-checks every hit by a direct DP.
- `verify_bands.py` recomputes every witness's H in pure Python.

| band | H | certificate for H | deletion values | certificate for deletions | lambda | prime |
| --- | --- | --- | --- | --- | --- | --- |
| band_3(13) | 237,817 | in S(12), witness | 20,485 .. 72,073 (7 values) | all <= 80,451, so in S(11) | 3 | yes |
| band_4(16) | 64,053,917 | in S(15), witness | 4,330,581 .. 15,978,923 (8 values) | in S(14), witnesses | 4 | yes |

Both are escapers. band_4(16) needs four re-connecting arcs, so the grammar fails at order 16 even with "ordinal sum
re-connected by at most 3 arcs". Both bands have twin distance 1: the pair {0, 1} is distinguished only by s + 1. So
both are one reversal from a substitution.

Witness files: `wit_band3_13_H.out`, `wit_band4_16_dels.out`, `wit_band4_16_H.out`.

The deletion witness for 15,927,743 = H(band_4(15)) also shows that band_4(15) is a non-creator. Its lambda is 4 (Lemma B2), so k*(15) >= 4. band_4(15) is not an escaper: its end deletion band_4(14) has H = 4,013,709, which exceeds 3,719,831, the largest known count of order 13 (THM-329; the Z_13 circulant gives 3,711,175). So band_4(14) may well be a first creator, but this is not proved.

## 4. No bounded grammar of this kind works at all orders (PROVED, conditional)

**Hypothesis (S-inf).** S(m) contains every odd h <= g(m) except 7 and 21, for some function g with g(m)^(1/m) -> inf.

HYP-9481 (a) (n(h) <= m(h) + 1) implies (S-inf) with g(m) = max H(m - 1), which is >= (m - 1)!/2^(m - 2) (Szele).
Chan-Pak's conjecture alone does not imply (S-inf), because it does not control the order.

**Theorem.** Assume (S-inf). Then for all integers s >= 3 and k >= 0, band_s(n) is, for every sufficiently large n:
- not a first creator;
- not a first creator plus at most k vertices;
- strong, prime, and of arc-strong connectivity exactly s.

So it escapes the grammar {FC; FC plus at most k vertices; SUM; SUB; ordinal sum re-connected by at most s - 1 arcs}.

Consequently no grammar of this form with fixed (k, r) generates all tournaments, and k*(n) -> inf (HYP-9532 (D)).

*Proof.*
- WLOG g is nondecreasing, because S(m) is contained in S(m + 1).
- Pick n with (s + 3)^(2n) < g(n - k - 1).
- Let J be a set of at most k vertices and U = band_s(n) - J.
- By the insertion lemma (THM-1370 (2)) and Lemma B3, H(U) <= H(band_s(n)) <= (s + 3)^(2n) < g(|U| - 1).
- H(U) is not 7 or 21 (THM-1370), so H(U) is in S(|U| - 1) and U is not a first creator.
- Strength, primality and lambda = s are Lemmas B1 and B2. ∎

The bands do not settle the question with near-substitutions allowed, because they have twin distance 1. Whether
dec(T) is bounded on escapers is HYP-9565 (A). It is 3 at order 12 (section 1.2).

## 5. First creators are the typical tournaments

**FINITE-EXACT (`fc_share.py`, klein's census).** The share of first-creator classes at orders 3..11 is:

| order | 3 | 4 | 5 | 6 | 7 | 8 | 9 | 10 | 11 |
| --- | --- | --- | --- | --- | --- | --- | --- | --- | --- |
| share | 0.5000 | 0.2500 | 0.5000 | 0.5893 | 0.7193 | 0.7932 | 0.8545 | 0.8884 | **0.9200** |

Orders 7-10 agree with HYP-9532 (C). Order 11 (831,487,872 of 903,753,248 classes) is new here.

**Mechanism (FINITE-EXACT).** Let R_H(n) = max H(n - 1)/E[H(n)], with E[H(n)] = n!/2^(n - 1) (Szele).

| n | 5 | 6 | 7 | 8 | 9 | 10 | 11 | 12 |
| --- | --- | --- | --- | --- | --- | --- | --- | --- |
| R_H(n) | 0.667 | 0.667 | 0.571 | 0.600 | 0.466 | 0.474 | 0.404 | 0.407 |

For n = 13 and 14 the maxima are not known: OEIS A038375 lists a(1..11) only (audit G1, M1).
- THM-329 (opus-2026-05-27-S6) found tournaments with 531,205 paths (order 12) and 3,719,831 paths (order 13) by
  search; neither is proved maximal.
- Mapping a Hamiltonian path to its last vertex v and a Hamiltonian path of T - v is injective. So
  max H(m) <= m·max H(m - 1), and max H(12) <= 1,141,140.
- Hence 0.349 <= R_H(13) <= 0.751 (PROVED), and R_H(14) >= 0.3495; no upper bound below 1 is known for R_H(14).

So for every n <= 13 the largest value of order n - 1 lies below the mean of order n, and every n-tournament with H
above max H(n - 1) is a first creator. For larger n this is conditional on (M_eps) below.

**Theorem (PROVED, conditional).** Assume:

> **(M_eps)** max H(m) <= (1 - eps)·((m + 1)/2)·m!/2^(m - 1) for all large m.

Then a uniformly random labelled n-tournament is a first creator with probability at least
1 - (W(n)/n! - 1)/eps^2 = 1 - O(1/n).

*Proof.*
- (M_eps) at m = n - 1 says max H(n - 1) <= (1 - eps)E[H(n)].
- Chebyshev's inequality and THM-589 (both legacy files: `THM-589-metagraph-H-variance-closed-form-nud-succession.md`
  and `THM-589-cvH-even-run-closed-form-and-2-over-n-rate.md`) give Var H/(E H)^2 = W(n)/n! - 1 ~ 2/n. ∎

**Corollary (conditional on (M_eps)).** Every tournament T is an induced subtournament of a first creator.

*Proof.* In a random N-tournament, the floor(N/|T|) disjoint |T|-sets each induce a copy of T independently, each with
probability >= 2^(-C(|T|, 2)). So a copy exists with probability -> 1. A first creator also occurs with probability
-> 1, so for large N some first creator contains T. ∎

(M_eps) follows from max H(m) = o(m · m!/2^(m - 1)); it is weaker, allowing c·m·m!/2^(m - 1) with c < 1/2. The known upper bounds do not reach it:
- Alon (1990) proved max H(m) <= c·m^(3/2)·m!/2^(m - 1).
- Friedgut and Kahn improved the exponent to 3/2 - xi for a small xi.
- The lower bound is (e - o(1))·m!/2^(m - 1) (Adler, Alon and Ross).

So (M_eps) is open, and its truth would follow from max H(m) = O(m!/2^(m - 1)).

## 6. The spanning-tree twin of the 7/21 proof (THM-4986)

**Lemma T1 (block floor; PROVED here, elementary).** Every 2-connected simple graph B on m >= 3 vertices has
tau(B) >= m, with equality iff B is the cycle C_m.

*Proof.*
- Write B as a cycle C plus open ears P_1, ..., P_k (Whitney).
- Adding an ear P with l >= 1 edges between distinct vertices x, y of a connected graph G gives
  tau(G ∪ P) = l·tau(G) + tau_2(G; x|y). Here tau_2 counts spanning forests of G with two trees separating x from y.
  - A spanning tree either omits exactly one edge of P and restricts to a spanning tree of G (l choices each),
  - or contains all of P and restricts to such a forest.
- tau_2 >= 1: delete an edge on the x-y path of a spanning tree.
- So tau(G ∪ P) >= l·tau(G) + 1 > tau(G) + (l - 1), and the ear adds l - 1 vertices.
- Starting from tau(C) = |C|, induction gives the claim, with strict inequality once there is an ear. ∎

**Theorem T2 (KNOWN, folklore; the setting of Sedláček 1970; proof transcribed).** The spanning-tree counts of connected simple graphs are exactly
{1} ∪ {N >= 3}.

The proof has the four steps of THM-1370:
- **Multiplicativity.** tau is the product of tau over the blocks.
- **Trivial blocks.** Bridges have tau = 1.
- **Floor.** By Lemma T1, a 2-connected block on m >= 3 vertices has tau >= m >= 3, so tau = 2 needs a block of order
  <= 2. Such a block is a bridge, with tau = 1. 2 is prime, so no product rescues it.
- **Completeness.** tau(C_N) = N for every N >= 3.

**Lemma T3 (crevasse lemma; PROVED, elementary).** Suppose:
- every structure decomposes into irreducible blocks;
- the count is the product of the block counts and equals 1 exactly on trivial blocks;
- every finite family of blocks occurs as the block family of some structure (ordinal sum; gluing along a path);
- the floor F(m) = min{c(B) : B irreducible nontrivial, |B| = m} is nondecreasing.

Then:
- (i) the set of counts is the multiplicative monoid generated by Irr = the set of irreducible counts (this needs only
  multiplicativity and realizability);
- (ii) since F is nondecreasing and unbounded, deciding whether v is in Irr needs only the finitely many orders
  m <= B(v) = max{m : F(m) <= v}.

(Audit G1, m11: the other hypotheses are not needed for (i), and (ii) is where the floor does its work.)

Both omission proofs are instances:

| world | blocks | floor F | the gaps |
| --- | --- | --- | --- |
| tournaments | strong components | Busch's f | 7 in (f(4), f(5)) = (5, 9), 7 prime; 21 = 3·7 and 21 not in Irr(<= 6) (finite check, THM-1370), where B(21) = 6 since f(7) = 25 |
| simple graphs | 2-connected blocks | F(m) = m | 2 in (tau(K_2), tau(C_3)) = (1, 3); 2 prime |

**The dictionary (ANALOGY with exact anchors on both sides).** See the README table.

**The contrast (PROVED / FINITE-EXACT).** Let R_tau(m) = max over (m - 1)-vertex graphs of tau, divided by
E[tau(G(m, 1/2))]. The maximum is tau(K_(m - 1)) = (m - 1)^(m - 3) (Cayley, and adding edges never lowers tau). The
mean is m^(m - 2)/2^(m - 1), since each spanning tree of K_m survives with probability 2^(1 - m). So

R_tau(m) = 2^(m - 1)·(m - 1)^(m - 3)/m^(m - 2) >= 2^(m - 1)/(e·m) -> inf,

whereas R_H(n) < 1 at every computed n and trends down (0.667 at n = 5, about 0.40 at n = 11 and 12; section 5).

**FINITE-EXACT tau-spectra (`tau_share.c`, every labelled graph on m <= 8 vertices, exact Bareiss determinants).**

| m | connected labelled graphs | values abs(S_tau(m)) | max tau | new values | labelled share new |
| --- | --- | --- | --- | --- | --- |
| 4 | 38 | 5 | 16 | 3 | 0.2632 |
| 5 | 728 | 16 | 125 | 11 | 0.3475 |
| 6 | 26,704 | 65 | 1,296 | 49 | 0.4190 |
| 7 | 1,866,256 | 386 | 16,807 | 321 | 0.5150 |
| 8 | 251,548,592 | 3,700 | 262,144 | 3,314 | 0.5719 |

abs(S_tau(m)) reproduces OEIS A182290 (1, 1, 2, 5, 16, 65, 386, 3700, ...) and the connected counts reproduce A001187, which validates `tau_share.c`. The spectrum is sparse: 3,700 values below 262,144. Tournament spectra are nearly solid below their maximum (1,520 of the 1,679 odd values up to 3,357 at order 9).

**Pre-registered Q3 (REFUTED).** The labelled share of graphs whose count is new rises from 0.4190 (m = 6) to 0.5150 and 0.5719; I had predicted a decrease. So new counts are common in both worlds at small orders, for opposite reasons:
- tournaments overshoot an almost solid old spectrum;
- graphs fill a very sparse one.

The first orders of alpha(N) for N <= 120, from `ALPHA=1 ./tau_alpha 8` (`tau_alpha_le8.out`; m <= 7 in
`tau_alpha_le7.out`), agree with OEIS A182258:
- alpha(N) = N for N = 3..7;
- 10, 13 and 22 do not occur below 9 vertices;
- 21 already occurs at 5 vertices.

The fixed points of alpha, 3, 4, 5, 6, 7, 10, 13, 22, are Nebeský's theorem (Časopis pro pěstování matematiky 98
(1973) 95-97); Azarija and Škrekovski completed the case 22.

**First creators across the dictionary (KNOWN, cited).**
- Sedláček (1970) introduced alpha(N), the least number of vertices of a graph with exactly N spanning trees. It is the
  graph analogue of n(h).
- Azarija and Škrekovski (Math. Bohem. 138 (2013) 121-131; arXiv:1210.6335) proved that alpha(N) <= (N + 4)/3 and
  beta(N) <= (N + 7)/3 hold together (beta = the least number of edges) iff N is not in
  {3, 4, 5, 6, 7, 9, 10, 13, 18, 22}. This set lies inside Euler's idoneal numbers, via the representation
  N = ab + bc + ca of theta graphs. They conjectured alpha(N) = o(log N).
- Stong (Australas. J. Combin. 82 (2022) 182-196) proved alpha(N) = O((log N)^(3/2)/log log N).

On the tournament side, n(h) >= m(h), and HYP-9481 (a) predicts n(h) <= m(h) + 1, with m(h) ~ log h/log log h.
