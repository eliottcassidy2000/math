# Exact strong/decomposable pairs: census, R-factors, one-apex form, a prime pair

kind-pasteur-2026-11-13-S26, THM-5151. **Chan–Pak remains open.** This note studies the "exact pairs" behind
Session 310's neutral refill. In such a pair a strong tournament and a decomposable one have the same full slot
polynomial. The note also asks what such pairs can and cannot supply to a repeatable return.

## 0. Conventions and inherited facts

- Q_T(t) = sum_j F_j(T) binom(t, j). F_j counts ordered covers of V(T) by j nonempty vertex-disjoint directed
  paths, so Q_T(t) counts path covers in t labelled slots with empty slots allowed. H(T) = F_1 = Q_T(1).
- Rows are outgoing-neighbour bit masks with labels from 0.
- TT_m is the transitive tournament. X => Y is the ordinal sum (every vertex of X beats every vertex of Y).
- C(T) counts directed Hamiltonian cycles.

The following are inherited and used without re-proof:

- (I1) Q is multiplicative over ordinal sums (THM-3166).
- (I2) Q_T(t) = sum over packings P of vertex-disjoint directed odd cycles of 2^|P| t^(u(P)+|P|), where u(P) is
  the number of uncovered vertices. This is the all-ones case of Grinberg–Stanley Theorem 1.39, re-derived in
  Session 310, POSITIVE_CLOCK.md.
- (I3) The component budget (Session 310, POSITIVE_CLOCK.md, last section). ord_0(Q_T), the order of vanishing
  at t = 0, is the number of singleton and odd strong components plus twice the number of even strong components.
- (I4) Equal Q permits replacement of a uniformly attached module in every exterior context (Session 310,
  NEUTRAL_TAILS.md).
- (I5) The frame Π has count H = 60 P_M(3) P_N(2) + H(M) H(N) for modules M (left) and N (right). Moving a
  trailing singleton of N to the front of M is one arc flip (Session 309, PHYSICAL_REFILL.md; Session 310,
  NEUTRAL_TAILS.md). Here P = Q.

## 1. Shape of a strong/decomposable collision (corollary of (I3))

**Proposition 1.** Let X be strong of order n >= 3, and let Z be non-strong with Q_Z = Q_X. Then:

- n is even;
- Z = Z_1 => Z_2 has exactly two strong components;
- each Z_i is a singleton or a strong tournament of odd order;
- at most one Z_i is a singleton.

*Proof.* By (I3), ord_0(Q_X) is 1 for odd n and 2 for even n. Z has at least two components, each contributing at
least 1, so ord_0(Q_Z) >= 2. Equality therefore forces n even, exactly two components, and contribution 1 from
each, which means a singleton or an odd strong component. Two singletons would make n = 2, which is impossible for
a strong X. QED.

This is immediate from Session 310's budget. It is recorded because it fixes the search space of section 3: a
detachable partner Y -> TT1 needs Y odd strong, and a block partner Y_1 => Y_2 needs both blocks odd strong.

## 2. R-factors: where neutral sites can come from

For a strong tournament S put R_S(t) = Q_S(t) / t^(ord_0(S)), and R_TT1 = 1.

- R_S is a polynomial with R_S(0) > 0.
- For odd |S| >= 3 by (I2), R_S(0) = 2C(S).
- For even |S| >= 4, R_S(0) = 2 c_(|S|-1)(S) + 4 d(S), where c_(|S|-1)(S) counts cycles of length |S| - 1 and
  d(S) counts unordered spanning pairs of disjoint odd cycles.

Write the strong components of T as S_1, ..., S_k. By (I1) and (I3),

    Q_T(t) = t^beta(T) * prod_i R_(S_i)(t),    beta(T) = s(T) + o(T) + 2 e(T),

where s, o, e count the singleton, the odd nontrivial and the even strong components.

**Proposition 2 (site bookkeeping).**

(a) Two tournaments have the same Q if and only if they have the same beta and the same product of R-factors.

(b) If Q_T = Q_T', then s(T') - s(T) = (o(T) - o(T')) + 2 (e(T) - e(T')).

(c) cap(Q) := max{ s(T) : Q_T = Q } equals beta - mu(R). Here mu(R) is the least value of o + 2e over all ways of
writing R = prod R_(S_j) as a product of R-factors of nontrivial strong tournaments, with mu(1) = 0.

Because Q is commutative over ordinal sums (I1), the singleton components of a realization can be placed last.
So cap(Q_N) is the largest number of trailing transfer sites that a Q-neutral rebase of a right module N can
expose.

*Proof.* (a) Both directions follow from the display: beta is the t-adic valuation (every R-factor is nonzero at
0), and division by t^beta leaves the R-product. (b) beta(T) = beta(T'). (c) Maximise s = beta - (o + 2e). QED.

**Reading.** A neutral rebase creates sites only by recombining R-factors:

- Split an even strong X into an odd strong Y and a singleton. This needs R_X = R_Y and gives +1 site. Session 310's
  B -> A -> TT1 is of this kind.
- Merge, meaning any recombination that lowers o + 2e. For example R_Z = R_(Y1) R_(Y2) for odd strong Z gives
  Y_1 => Y_2 -> Z -> TT1, +1 site. Three odd cores merging into one even component (R_X = R_(Y1) R_(Y2) R_(Y3))
  also gives +1 site.
- Recombinations through even components, such as R_X = R_(Y1) R_(Y2) (the H = 45 type below), change no site
  count by themselves.

For N = B^r, beta = 2r and the R-product is R_A^r, so cap <= 2r - 1. This is Session 310's ceiling. Their A^r -> TT_r
realizes r. Exceeding r requires a realization with o + 2e < r, that is, a merge.

**Proposition 8 (site budget; unconditional).** If Q/t is the slot polynomial of some tournament, then

    cap(Q/t) = cap(Q) - 1.

Consequently, consider any sequence of Π transfers (a trailing singleton of the right module N moves to the
left module) interleaved with Q-neutral rebases of N. The number of transfers is at most cap(Q_(N_0)). A neutral
rebase leaves cap(Q_N) unchanged, and a transfer lowers it by exactly one.

*Proof.* If T' realizes Q/t with s' = cap(Q/t) singleton components, then T' -> TT1 realizes Q with s' + 1. So
cap(Q) >= cap(Q/t) + 1, and in particular cap(Q) >= 1.

Conversely, take T realizing Q with s = cap(Q) >= 1. Reordering the strong components keeps Q by (I1), so we may
place a singleton component last and write T = T'' -> TT1. Then T'' realizes Q/t with s - 1 singletons, so
cap(Q/t) >= cap(Q) - 1.

A transfer replaces N = N' -> TT1 by N', so Q_N becomes Q_N / t. A neutral rebase keeps Q_N. Every transfer needs
a singleton in the current realization, and s <= cap. QED.

So a clock fed only by neutral refills runs at most cap(Q_(N_0)) steps. Through order 9 this budget is
s(N_0) + d(N_0) (Proposition 3' below), and HYP-9678 (C) conjectures it in general. For B^r this would be exactly
r, the number Session 310 realizes. Unconditionally, an unbounded return must raise the budget by non-neutral
moves.

**Proposition 3 (no small odd merge; FINITE-EXACT).** Take odd strong Y_1, Y_2 with |Y_1| + |Y_2| <= 10. Then no
tournament Z has Q_(Y_1 => Y_2) = Q_(Z -> TT1). This covers all pairs of orders (3,3), (3,5), (3,7), (5,5).

A candidate Z would have ord_0(Q_Z) = 2 - 1 = 1. By (I3) it would then be a single strong component of odd order
|Y_1| + |Y_2| - 1 <= 9. Comparing coefficients, it would also have triangle count c_3(Y_1) + c_3(Y_2) (from
t^(n-2)) and C(Z) = 2 C(Y_1) C(Y_2) (from t^2).

*Check.* All strong classes of orders 3, 5, 7, 9 (1, 6, 353, 178,133 classes) were taken from the complete
class lists (classgen, first_creators_genealogy_20261107). Every product Q_(Y1) Q_(Y2) was compared, on t = 1..n,
with t Q_Z for every strong Z of order n - 1. There are no matches (`scripts/refill2.py`).

So a right module that has spent its sites, consisting of odd cores of small total order, has no hidden
capacity. The order-11 case is pre-registered search S2 (section 7).

**Proposition 3' (capacity formula; FINITE-EXACT through order 9).** For a tournament T let s(T) count its
singleton components. Let d(T) count its even strong components that have a detachable partner, that is, are
Q-equal to Y -> TT1 for some odd strong Y. Then for every tournament T of order <= 9,

    cap(Q_T) = s(T) + d(T).

Equivalently, mu(R) is additive over strong components at these orders: an odd strong core hides no site, an even
component hides at most one, and it hides one exactly when it is detachable.

*Check.* `scripts/capcheck.py` runs over all 198,946 classes of orders 3-9.

- It computes the strong components of each class.
- It identifies the detachable even components by profile. Through order 8 these are the order-6 and order-8
  classes of Proposition 4.
- In every one of the 52,432 profile groups, s + d is the same for all members and equals the largest s.

The general statement is HYP-9678 (C).

## 3. Census of exact strong/decomposable pairs

**Precursor.** Session 261 constructors/MODULE_REBASE.md §3 ran full cover-profile collision censuses at orders 5
and 6 (56 classes and 24 profiles at order 6), without separating strong from non-strong.

**Proposition 4 (FINITE-EXACT, all orders <= 9).** Exactly 7 isomorphism classes of strong tournaments of order
<= 9 share their full profile with a non-strong tournament. All 7 have C = 1 (uniquely Hamiltonian).

| order | H | strong classes | non-strong partners | trade type |
|---:|---:|---:|---|---|
| 6 | 15 | 2 (converse pair) | Y5 -> TT1, TT1 -> Y5 (Y5 = classgen code 152, H = 15, C = 3) | none: identical odd-cycle census {3: 4, 5: 3}, no disjoint pair |
| 8 | 45 | 3 | C3 => Y5', Y5' => C3 (Y5' = code 739, H = 15, C = 2; a different profile from code 152) | block pair, no singleton site |
| 8 | 95 | 2 = Session 310's B and its converse | A -> TT1, TT1 -> A, with A any of 5 strong order-7 classes of one profile | six triangle pairs <-> twelve pentagons (Session 310) |

There are none at orders 7 and 9, by Proposition 1. Method: fprof.c computes C and all F_j for each of the
1+1+2+4+12+56+456+6880+191,536 classes of orders 1-9. Classes are then grouped by profile, and the groups that
mix strong and non-strong classes are reported (`scripts/census.py`).

**Proposition 5 (FINITE-EXACT, uniquely Hamiltonian slice to order 12).**

(a) Generation. Extending each C = 1 strong class of order m by one vertex in all 2^m - 2 strong ways keeps
exactly the extensions with C = 1. Canonical forms come from the unique Hamiltonian cycle (the lexicographically
least pair word over its rotations). The procedure gives 1, 3, 8, 21, 55, 144, 377, 987, 2584 classes at orders
4..12. These equal F_(2n-6), so the generation is complete: it is enough that every generated class is C = 1
and strong, and that the true count is F_(2n-6).

The count is KNOWN: Garey 1972 (also Moon 1982), re-proved in Session 197's C1_COMPLETE / FIBONACCI_CLASS and in
this machine's THM-4868 (G) (S10, which enumerated all F_(2n-6) classes through order 18). The structure is also
KNOWN: Shelor 1977, opus-S74's THM-4867 (an in-degree-1 vertex and a one-arc cut) and THM-4868 (E) (rigidity). Only
the slot profiles of these classes are new here.

(b) Order 10. Exactly 6 of the 377 classes collide with a non-strong tournament:

- H = 153: 2 classes, partners C3 => Y7 (4 strong order-7 classes of one profile);
- H = 219: 2 classes, both one-apex (O(001100001) and its converse O(011110011)), partners C3 => Y7 (6 classes);
- H = 317: 2 classes, partners Y9 -> TT1 (2 strong order-9 classes of one profile). **317 is prime.**

This is complete within C = 1, since all two-block partners of Proposition 1 have block orders <= 9.

(c) Order 12. None of the 2584 classes collides with a two-block sum of type (3,9) or (5,7). Type (1,11) is
pre-registered search S1 (section 7).

**Proposition 5(d) (outside C = 1 at order 10; FINITE-EXACT in its range; pre-registered S3).** Search every
strong order-10 tournament with at most 5 backward arcs relative to 0<1<...<9 (1,385,979 candidates). Compare each
against all two-block sums of Proposition 1: TT1 + Y9, C3 => Y7 and Y5 => Y5' over all odd strong blocks (48,908
profiles). The matches have four profiles:

- H = 153, 219, 317, all with C = 1;
- **H = 333, with C = 5.**

The H = 333 profile equals Q(C3 => Y7) for Y7 = classgen order-7 code 43104 (H = 111, C = 11). One member is
X = (1002,1020,1017,1008,481,704,896,768,544,16). Among the 120 labelled H = 333 matches there are exactly 2
isomorphism classes, a converse pair. The author's vertex-deletion fingerprints separate 2; the audit's canonical
forms show there are no more. So the uniquely Hamiltonian property of Proposition 4 is a small-order fact and fails at order 10. Within
this range, every detachable match (partner Y9 -> TT1) is the H = 317 profile, which has C = 1.

**Remark (C = 1 is not forced by the construction at orders <= 9).** Proposition 4 searched every class. The C = 1
property of all 7 strong members is therefore an observation, not a restriction of the search. Proposition 5(d)
shows it does not persist: HYP-9678 (A) is refuted and replaced by the detachable form (A'). A complete order-10
census would need the 9.7 million order-10 classes and was not attempted.

## 4. One-apex form, and what Session 310's B is

A one-apex tournament O(w) has an apex x and a transitive body 1 -> 2 -> ... -> m (i -> j for i < j). So body
vertex 1 beats the rest of the body and is its top-ranked vertex. The word w in {0,1}^m has w_i = 1 iff body vertex i
beats x.

Session 191 (one_apex/PROOF.md) uses the complementary letters: there x -> i iff b_i = 1. So O(w) is Session 191's
T(w̄), and the same word names different tournaments in the two packets.

**Proposition 6.**

    Q_(O(w))(t) = t^(m+1) + 2 * sum over a < b with w_a = 0, w_b = 1 of sum over l >= 1 of
                  binom(b - a - 1, 2l - 2) t^(m + 1 - 2l).

O(w) is strong iff w_1 = 0 and w_m = 1, and then C = 1.

*Proof.* The body is acyclic, so every directed cycle passes through x. It is x -> a, then an increasing chain
from a to b in the body, then b -> x, with w_a = 0 and w_b = 1. Its length is one plus the chain size, so it is
odd iff the chain has 2l vertices; there are binom(b - a - 1, 2l - 2) such chains. All cycles share x, so the only
packings are the empty one and single cycles. Apply (I2). A Hamiltonian cycle needs the chain to be the whole
body: this exists iff w_1 = 0 and w_m = 1, and it is unique. A Hamiltonian cycle makes O(w) strong. Conversely, a
strong tournament has one by Camion's theorem, which forces w_1 = 0 and w_m = 1. QED.

**Credit.** "Odd cycles of a one-apex tournament are increasing chains through the apex" is KNOWN:

- feedback_word_hamilton_spectrum_20261020/intercyclic_rays.md §4 (Oct 20);
- THM-5134 (G), "every odd cycle of a one-apex tournament passes through its apex";
- THM-5070 (2), this machine's S18.

The t-graded form is (I2) applied to it, and the H formula is the t = 1 case. Write O(m,k) = O(0^m 1^k): the apex
beats the first, top-ranked, m body vertices and loses to the last k. Its value H = 2 + (2^m - 1)(2^k - 1) is
Session 191's U(o,z), formula (5).

**Identification (new; verifier step 3).** Session 310's eight-vertex B (rows 2,80,3,7,173,143,189,15) is
isomorphic to O(0011111) = O(2,5). Its apex beats the two top-ranked vertices of TT7 and loses to the other five.
Its converse is O(0000011) = O(5,2). Among all 128 words of length 7, B matches only O(0011111) (audit).

- Proposition 6 explains Session 310's cycle table by counting chains: 10 = 2·5 triangles, 30 pentagons and 7
  heptagons. Q_B = t^8 + 20 t^6 + 60 t^4 + 14 t^2 and H = 2 + 3·31 = 95.
- The order-6 pair member is O(01001).
- The order-6 member's converse is O(01101). Both order-10 members at H = 219 are one-apex: O(001100001) and its
  converse O(011110011).
- The members at H = 45, 153, 317 are not one-apex (no vertex deletion is transitive).
- No partner block other than C3 is one-apex. This was checked for the order-5, order-7 and order-9 partner
  classes listed above.
- Strong O(m,k) are C = 1, so they lie in the complete slice. Most of them have no partner: at order 6 none does
  (H = 17, 23); at order 8 only O(2,5) and O(5,2) do (H = 95; O(1,6), O(3,4) have H = 65, 107); at order 10 none
  does (H = 257, 383, 443, 467).

**Proposition 6' (Hamiltonian cycles of partners).** Compare the t^2 coefficients through (I2).

- If an even strong X of order n has Q_X = Q_(Y -> TT1), then C(Y) = c_(n-1)(X) + 2 d(X).
- If instead Q_X = Q(Y1 => Y2) with both blocks odd strong, then c_(n-1)(X) + 2 d(X) = 2 C(Y1) C(Y2).

For a strong one-apex X = O(w) with n >= 6, d(X) = 0. A cycle on n - 1 vertices omits exactly one body vertex.
Omitting an interior one always leaves an admissible chain. Omitting body vertex 1 (the body's source) is
admissible iff w_2 = 0, and omitting body vertex n - 1 (its sink) iff w_(n-2) = 1. Hence

    c_(n-1)(O(w)) = (n - 3) + [w_2 = 0] + [w_(n-2) = 1] >= 3.

So the detachable partner of a one-apex member has at least n - 3 Hamiltonian cycles. It is never uniquely
Hamiltonian, and in particular never one-apex. Checks against the census:

- the order-6 member O(01001) gives C(Y5) = 3;
- B = O(0011111) gives C(A) = 7, matching Session 310's seven heptagons;
- the H = 219 member O(001100001) gives C(Y7) = 8/2 = 4.

The H = 317 member is not one-apex, and its partner Y9 has C = 8.

## 5. A prime exact pair, and a Π run from a prime refill to 7p

Take X = (254,508,120,112,224,192,128,268,637,255) and Y = (378,380,505,496,224,448,384,259,16).

- X is strong of order 10 with C = 1. Y is strong of order 9.
- H(X) = H(Y) = 317, which is prime.
- F(X) = (317, 17926, 284076, 2092344, 8620200, 21511440, 33304320, 31288320, 16329600, 3628800).
- Q_X(t) = t Q_Y(t) = Q_(Y -> TT1)(t), checked at t = 0..11.
- Q_Y(2) = 9280.

By (I4), X^r may be replaced by Y^r -> TT_r in every uniform context. In Π with N = X^r and any nonempty M, the r
transfers of Session 310's argument apply verbatim. That argument uses only the trailing singletons, not the
core. They give

    H_k = 60 P_M(3) · 9280^r · 2^(r-k) · 3^k + H(M) · 317^r,   k = 0..r,   H_(k+1) = (3 H_k - H(M) 317^r) / 2.

Unlike Session 310's B^r clock, whose values are all divisible by 5, the 5-adic and 3-adic residues here are those
of 317^r H(M).

**Examples** (deterministic Miller–Rabin; checked in verify_s26.py):

- M = TT3, r = 2: 558,047,332,489 (prime) -> 837,070,948,489 -> 1,255,606,372,489 = 7 × 179,372,338,927 (prime).
  The neutral refill sits at a prime, exposes two sites, and the two literal transfers end at a 7p target.
- M = TT1: r = 2, k = 2 gives the prime 139,511,908,489; r = 3, k = 3 gives the prime 3,884,008,766,575,013.

**Proposition 7 (where a clock can be divisible by a prime; elementary).** Let H_k = I + c 2^(r-k) 3^k and let l
be a prime not dividing 6c. Then l | H_k iff (3/2)^k ≡ -I / (c 2^r) (mod l). So the k with l | H_k form either
the empty set or one residue class modulo ord_l(3/2).

For l = 7, 3/2 ≡ 5 is a primitive root, so ord_7(3/2) = 6:

- if 7 does not divide c·I, exactly one of any 6 consecutive k has 7 | H_k;
- if 7 | I but 7 does not divide c, no term is divisible by 7.

Session 309's clock has I = 499, c = 1,303,200 and r = 6, so c 2^r = 83,404,800. It hits 7 at k = 2:
187,661,299 = 7 × 26,808,757. A prime l with -I/(c 2^r) outside the subgroup generated by 3/2 never divides any
term. For that clock these are exactly 13, 19, 29, 43, 47, 53 below 60 (audit); the next are 61, 67, 71, 97,
101.

So a long run cannot avoid 7 unless 7 | c·I. That is good for 7p targets, which recur with period 6 in k, and bad
for a run of consecutive primes.

**Credit.** The mod-7 mechanism is already in Session 309's MIGRATING_CLOCK.md §2, on the same day and for the same
clock family ("Modulo 7 the ratio ... has order 6 ... 7 | H_n iff n ≡ 0 mod 6"). That 3/2 is a primitive root mod 7
is also in THM-4800. Proposition 7 is the elementary general form. The history sweep's comparison with
Rabinowitsch's prime-producing polynomials is an analogy only, with no class-number content.

## 6. What this does and does not give toward a return rule

**Supplied.**

- The complete exact-pair stock through order 9.
- The C = 1 stock through order 10: detachable partners at H = 15, 95, 317, and block partners at H = 45, 153, 219.
- Outside C = 1, in S3's range at order 10, the block pair H = 333 with C = 5.
- A prime-valued detachable pair whose Π clock is not stuck at a fixed prime divisor.
- The bookkeeping of Proposition 2: every new site in a neutral rebase is an R-factor recombination.

**Not supplied.**

- An infinite family of primitive pairs.
- A second refill after a run. After the r transfers the right module is Y^r, a sum of odd cores. By Proposition
  2 its hidden capacity is zero unless merges exist (any recombination lowering o + 2e). Proposition 3 excludes
two-core odd merges with |Y1| + |Y2| <= 10. Proposition 3' excludes every merge when the module has total order at
most 9.
- A target selector. By Proposition 8, an unbounded return needs non-neutral moves. If, as HYP-9678 (C)
  conjectures, merges never occur, each refill needs a fresh even component with a partner, brought in by a
  non-neutral move (for example Session 309's migration and recharge) before the neutral rebase that splits it.
  HYP-9678 records this and the conjectures behind it.

## 7. Pre-registered searches (PREREGISTRATION.md; results filled in after the runs)

- S1: order-12 detachable partners among strong order-11 tournaments with at most 5 backward arcs. RESULT: no
  match, so the prediction FAILED. 1,291 targets; 3,847,591 candidates, 1,243,347 strong, 312,197 profiled.
- S2: order-11 odd merges, targets Q(Y5 => Y7)/t, same restriction. RESULT: no match, prediction CONFIRMED. 655
  targets; 189,829 profiled.
- Not pre-registered: biclique backward patterns L × U with L below U and |L|, |U| <= 5 (7,469 order-11
  tournaments) contain no S1 partner either.
