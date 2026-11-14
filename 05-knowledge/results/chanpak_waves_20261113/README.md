# Chan–Pak waves: prediction, test, verdict (session mac-mini-2026-11-13-cuts, sixth prompt)

**Chan–Pak OPEN.** Each wave states predictions, runs the cheapest decisive test, and records the verdict. Data:
complete strong spectra Sp(n) of orders ≤ 9 (`gentourng` + `hl`, see the busch_staircase note). Notation: in base
b, I(Ω(T), b) = 1 + b·m with m = c1 + b·c2 + b²·c3 + …, where c_k counts sets of k pairwise vertex-disjoint odd
cycles (THM-5070). Gap_b is the set of unattained m.

## Wave 1: pairwise-meeting odd cycles

- **Prediction P11.** Every attained H is attained by a tournament whose odd cycles pairwise meet (Ω complete).
- **Verdict: REFUTED.** Below 125, the odd counts with no Ω-complete realizer are exactly:

  7, 21, 29, 35, 41, 43, 55, 57, 59, 67, 71, 75, 77, 83, 89, 93, 101, 105, 109, 111, 113, 119 (22 values).

  The list is FINITE-EXACT. Ω is complete iff the tournament has no two vertex-disjoint 3-cycles, since any
  non-transitive set contains a triangle. Such a realizer has a single nontrivial strong component, and Busch's floor
  b(10) = 125 caps its order at 9.

- **Theorem A (PROVED, elementary).** For m < b, m ∈ Gap_b ⟺ 2m + 1 is not the count of an Ω-complete tournament.
  *Proof.* If c2 ≥ 1 then m ≥ b. So m < b forces c2 = c3 = … = 0, i.e. Ω ≅ K_m, and then H = 1 + 2m. ∎
- **Corollary.** For every base b ≥ 60, Gap_b contains the 22 m-values of the list above:
  {3, 10, 14, 17, 20, 21, 27, 28, 29, 33, 35, 37, 38, 41, 44, 46, 50, 52, 54, 55, 56, 59}.
  This is consistent with THM-5070's tables for b ≤ 8, including m = 14 (N = 29) in bases 4, 6, 7, m = 27 (55) in
  bases 4, 5, and m = 21 (43) in base 6. The number of base-b gaps grows with b; HYP-9530 (B) asks only for
  finiteness at each b.
- **Prediction P13 (OPEN).** The Ω-complete spectrum misses a positive proportion of odd numbers (about 35% below
  125), so |Gap_b| ≥ c·b. The route is a structure theorem for tournaments with no two disjoint triangles. They
  include the one-centre (hub) family, but not only it: for example the regular 5-tournament, H = 15, has no hub.

## Wave 2: the universal gaps

- **Theorem B (PROVED).** m = 3 ∈ Gap_b for every b ≥ 2.
  *Proof.* m = 3 forces c1 = 3 and c2 = 0; in base 2 the alternative c1 = 1, c2 = 1 is impossible. So there are three
  pairwise-meeting odd cycles and no others. They lie in one strong component C, and every other component is
  trivial. But C has 1 odd cycle if |C| = 3 and 2 if |C| = 4 (the unique strong 4-tournament R_2). If |C| ≥ 5 it has
  at least |C| − 2 ≥ 3 triangles (Moon) plus a 5-cycle (pancyclicity), so at least 4 odd cycles. ∎
- **Theorem C (PROVED; uses 21 ∉ S).** m = 10 ∈ Gap_b for every b ≥ 9.
  *Proof.* For b ≥ 11, Theorem A applies and N = 21 is not attained at all. For b = 9 or 10, the decompositions
  10 = c1 + b·c2 force c2 = 0, since a disjoint pair needs c1 ≥ 2. So Ω ≅ K_10 again. ∎
  - With THM-4940 / THM-5070 for b ≤ 8: m = 10 is a gap for every b except **4 and 8**. In base 8 the exception is
    C3 ⇒ C3, with m = 2 + b = 10. In base 4 it is a profile c1 = 6, c2 = 1.
  - *Correction.* The brick_moves note listed base 3 as an exception; base 3 does miss N = 31 (m = 10).

## Wave 3: the brick law

- **Prediction P1** (from the brick_moves note). A count shared by two-component compositions (a, b) and (c, d) of n,
  with a ≤ b and c ≤ d, is an S-unit re-association with S = primes(Sp(a)) ∪ primes(Sp(c)).
  - **Verdict: CONFIRMED at n = 13.** (4,9)/(5,8): 373 shared counts; (4,9)/(6,7): 249; (5,8)/(6,7): 320. 0
    violations. Together with n = 11, 12 that is 1,481 shared counts with 0 violations.
- **Stronger prediction:** the fixed set {3, 5, 11, 13} suffices for all n.
  - **Verdict: REFUTED.** At n = 13, 425 = 5·85 = 17·25 moves a 17, a prime of Sp(6). At n = 14 there are 360 and
    363 failures. The movable bricks grow with the smallest core size, exactly as P1 says.
- **Prediction P7:** 273 = Φ3(16) is a Fermat-block cut 17 + 2^j·L.
  - **Verdict: REFUTED for sources of order ≤ 9.** No realizer of 17 up to order 9 has a power-of-two port. The
    ports of R_k are 2^(k+1) − 1, then a descending interior sequence, then 2^k + 2 and 2^k + 1; for example R_4 has
    31, 69, 44, 28, 18, 17. Powers of two occur only as q(R_2) = 32 (a left-padding root) and the port 16 of R_3. 273
    is attained by other means.

## Wave 4: structure of pairwise-meeting tournaments (FINITE-EXACT, orders 5–9)

Strong tournaments with no two vertex-disjoint triangles (Ω complete):

| order n | total | with a hub (one vertex on every triangle) | hub-free | path counts of hub-free ones |
|---|---:|---:|---:|---|
| 5 | 6 | 3 | 3 | 13, 15 |
| 6 | 11 | 7 | 4 | 25, 33 |
| 7 | 20 | 15 | 5 | 49, 69 |
| 8 | 37 | 31 | 6 | 97, 141 |
| 9 | 70 | 63 | 7 | 193, 285 |

**Exact fits.** The total is 2^(n−3) + n − 3, of which 2^(n−3) − 1 have a hub and **n − 2 are hub-free**. Every
hub-free one has a triangle transversal of size 2, and its count is **1 + 3·2^(n−3)** or **9·2^(n−4) − 3**. The first
sequence is THM-5126's line 1 + 3·2^s for s ≥ 2; the line's only missing point is 7.

- **Prediction P16.** The pattern holds for every n ≥ 5. Then the pairwise-meeting spectrum is the hub (one-centre)
  spectrum W together with {1 + 3·2^s} and {9·2^s − 3}. By Theorem A, for large b the base-b gaps below b are the
  m with 2m + 1 outside that set: about half of all m < b, if HYP-9469's one-centre density 1/2 holds. Base 2 has
  only two gaps; large bases have about b/2 below b.
  - *Test:* order 10, where the prediction is 8 hub-free classes with counts 385 and 573. **CONFIRMED EXACTLY**
    (`scripts/omega10.c` over all 9,733,056 order-10 tournaments): 135 = 2^7 + 7 strong tournaments with no two
    disjoint triangles, 127 = 2^7 − 1 with a hub, and 8 hub-free. All 8 have τ = 2, and their counts are 385 and 573.
    So the counts for orders 5–10 are FINITE-EXACT and match the formulas.
  - The τ ≤ 2 part should be a directed-triangle packing/covering (Erdős–Pósa type) fact; this has not been checked
    against the literature.

## Wave 5: the Two-Hub Conjecture (headline)

- **Prediction.** Two hubs suffice for Chan–Pak, while one hub misses about half.
- **Verdict: CONFIRMED far beyond expectation.** Order ≤ 9 tournaments with a triangle transversal of size ≤ 2
  attain every odd N < 700 except 7 and 21. The exact counter `scripts/twohub.c` (below) then shows that every odd
  N ≤ **22,156,805** except 7 and 21 is the path count of a two-hub tournament on at most 18 vertices (2,327,673 at
  16 vertices). The first miss
  at word length ≤ 14 is 2,327,675, a length-cap effect. Length ≤ 12 already covers every odd N ≤ 200,000 except 7
  and 21.

**The counter.** The spine is TT_m. Hubs u1 and u2 have words x, y ∈ {0,1}^m (x_i = [u1 → i]) and ε = [u1 → u2]. A
Hamiltonian path has the form S0 h1 S1 h2 S2, with increasing spine segments S0, S1, S2 and junction arcs at the
hubs. So H is a **rational series in the four-letter alphabet (x_i, y_i)**: a 36-state transfer (2 hub orders × 3 × 3
× 2 segment states), with H = u_ε·M_w·v0. It agrees with brute force on 300 random configurations (0 mismatches).

**Two-Hub Conjecture (TWO-HUB CHAN–PAK).** Every odd N ∉ {7, 21} is H(T) for a tournament T with two vertices
whose deletion leaves a transitive tournament. This is a strengthening of Chan–Pak and a reduction of it to one
explicit 4-letter numeration.

**Proof route.** By OCF, H = 1 + 2(c1 + c2 + c12) + 4b. Here c1, c2 count odd cycles through u1 only or u2 only,
c12 those through both, and b the disjoint pairs. When every u1-cycle is disjoint from every u2-cycle, b = c1·c2 and

    H = (1 + 2c1)(1 + 2c2) + 2·c12,

a product of two one-hub counts (each from the half-dense W) plus an even cross-cycle correction. The aim is a
sub-automaton in which c12 moves in steps of 1 over an interval while c1 and c2 stay fixed; combined with W·W this
gives a digit argument for every large odd N.

**Predictions.**
- **P17. CONFIRMED.** The horizon grows about tenfold per two extra spine vertices. At word length ≤ 16 (18
  vertices), every odd N ≤ **22,156,805** except 7 and 21 is a two-hub count; the first miss is 22,156,807. The run
  took 4 minutes.
- **P18.** A sub-automaton exists with H′ = 3H + d for d ∈ {0, 2, 4}, closed under itself: the clone-clock ternary
  digits realized inside two hubs. *Test:* search 1–2 letter suffixes with affine action on H for a fixed prefix
  class.
- **P19.** No odd N other than 7 and 21 needs more than two hubs. *Test:* the two-hub counter at word length 16
  combined with the 2^32 verification of THM-4875 (D).

## Wave 6: the digit structure of two-hub counts

**Horizon law (FINITE-EXACT).** X(L) is the first odd miss above 21 using words of length ≤ L. The ratio X(L)/3^L
rises toward balanced-ternary scale:

| L | 6 | 7 | 8 | 9 | 10 | 11 | 12 | 13 | 14 | 16 |
|---|---|---|---|---|---|---|---|---|---|---|
| X(L)/3^L | 0.41 | 0.33 | 0.40 | 0.38 | 0.41 | 0.44 | 0.43 | 0.49 | 0.487 | 0.515 |

L-digit balanced ternary covers (3^L − 1)/2. At L = 14 the first miss is 2,327,675 = 5²·7·47·283, against
3^14/2 = 2,391,484; the last covered value 2,327,673 = 3·601·1291. The length cap binds exactly where an L-digit
ternary register runs out.

**Mechanism (PROVED, elementary): two-hub counts are sums of 3-smooth frames.** In a path S0 h1 S1 h2 S2, fix the
extreme elements: a = max S0, b = min S1, c = max S1 and d = min S2. Each of these must carry the right hub label.
Every other spine vertex i may then join any segment whose range contains it, so it contributes a factor
r(i) = [i < a] + [b < i < c] + [i > d] ∈ {0, 1, 2, 3}. Hence

    H = Σ over hub orders and admissible frames (a, b, c, d) of 2^α·3^β,

where α and β count the vertices covered by exactly two or exactly three ranges (frames with a zero factor drop out).
The ternary digits are the triply covered vertices.

**Nested-frame lemma (PROVED, elementary).** Take the hub order (u1, u2) and an admissible frame with positions
d < b < c < a. Then x_a = 0, x_b = 1, y_c = 0 and y_d = 1. The coverage multiplicity is:
- 1 on [0, d) and on (a, m);
- 2 on (d, b) and on (c, a);
- 3 on (b, c).

Every vertex is covered, so the frame contributes exactly

    2^((b − d − 1) + (a − c − 1)) · 3^(c − b − 1).

Example: if u1 beats only b0 and u2 is beaten only by c0, then b and c are forced. Summing the free positions a and d
gives geometric series, so the counts contain 3^(c0 − b0 − 1)·(2^k − 1)(2^l − 1) plus the contributions of the other
frame orders. THM-4876's one-centre family (2^i − 1)(2^l − 1) + 2 is the β = 0 shadow of this. The ternary register
β is set by the gap between u1's chosen out-neighbour and u2's chosen in-neighbour. Many admissible frames together
give the dense numeration.

**P18 in its local form (append one vertex, get H ↦ 3H + d): REFUTED.** Over all 43,690 configurations with words of
length ≤ 7, appending one spine vertex gives all three branches 3H, 3H+2, 3H+4 for only 38 of them, and none for
39,276. Prepending gives identical statistics, as the converse symmetry predicts. So the ternary horizon is a global
frame effect, not a vertex-by-vertex digit recursion.

**Revised proof route.** Show that frame sums realize every odd N in [69, c·3^L] for some c > 0. Since each frame is a
3-smooth term 2^α·3^β, a natural lemma would realize sums of 3-smooth numbers with a divisibility-antichain pattern.
That compares with the classical d-completeness of {2^a·3^b}: Erdős's question whether every positive integer is a
sum of 3-smooth numbers none dividing another (attribution and exact statement still to be checked against the
literature). Two hub orders and the four frame positions supply the signed or overlapping choices.

**FRACTRAN / PRIMEGAME reading.** Conway's PRIMEGAME is 14 fractions over the 10 primes 2–29. Four of them are
registers (2, 3, 5, 7), and six are state primes (11, 13, 17, 19, 23, 29). Minsky's two-counter universality encodes
a configuration as 2^a·3^b, exactly the shape of one two-hub frame. Read this way, the two-hub tournament is a
two-register machine: the 36 transfer states play the role of the state primes, and each frame is one register
configuration 2^α·3^β. That PRIMEGAME's program length 14 equals the word length 14 whose horizon first exceeds 2×10^6
is a coincidence. The structural match is "few registers × finite control".

**Next predictions.**
- **P20.** Single-frame-dominated words realize every number 2^α·3^β·(small cofactor) + 1. *Test:* enumerate words
  whose admissible frames are at most 3, and list their counts.
- **P21.** X(L)/3^L converges to a constant c ∈ (0.5, 0.6). *Test:* L = 17, 18 with a C run on a 64-bit coverage
  bitmap.
- **P22.** Every odd N ≤ X(L) has a realizer of length ≤ L that uses both hub orders. *Test:* count realizers by
  order.

## Wave 7: the absolute-target problem, trijections, and 69/34/45

**Absolute-target decoder (`scripts/decode2hub.c`).** Given N, it builds a two-hub word letter by letter. For every
remaining length r it precomputes upper and lower bounds on what each of the 36 automaton states can still
contribute, and it prunes any prefix whose range misses N. The search is exact: no false pruning, and witnesses are
verified at the leaf.
- It solved 20/20 random odd targets near 10^3 and 20/20 near 10^6, with explicit words. Example: N = 1,359,795 is
  the word 10213223200100, with letters 2x + y and e = 0.
- The search cost grows roughly like N (about 10^6 nodes at 10^6), because the per-state bounds are loose. A run at
  10^9 was stopped.
- Simulated annealing on |H − N| solved **0/20** even at 10^6: a single bit flip moves H by a large factor, so the
  landscape is rugged.

Status: the absolute-target problem is constructive at these scales. The next step is tight per-prefix ranges, using
the frame structure or a low-rank split of the 36-dim bilinear form H = c(suffix)·v(prefix) for meet-in-the-middle.
This is the step that goes beyond coverage tables and witness counts.

**Trijection link (CITED: [tournament_graph_trijections](../tournament_graph_trijections_20261020/README.md)).**
A constraint graph with m edges and c components charts a 2^−(m−n+c) fraction of switching classes. It is one-to-one
exactly for spanning trees, that is n − 1 edges and connected.
- A two-hub Hamiltonian path is such a tree.
- Deleting the hubs leaves a 3-component linear forest on the spine, the three runs, with cyclomatic deficiency
  c − 1 = 2.
- The triply covered spine vertices, the ternary digits, are exactly those free to join any of the three
  components.
So the frame count is a count of 3-component linear forests with prescribed endpoint labels.

**69, 34, 45 (connections found).**
- 69 is where the order-8 strong run [69, 609] begins. It is also the first value from which every odd N is a single
  Fermat-block cut of a ≤ 9-vertex source.
- 34 = m(69) is the cycle-content coordinate of 69, and it equals A059735(6), the converse-merged tournament classes
  on 6 vertices (trees5_tournaments4 note).
- 45 = H_max(6) = b(8) = 9·5, a Busch floor and a product of the bricks.
- 69 = 45 + 24 = 34 + 35 (m + r). 69 = 3·23, and 23 is the next Fermat-blind prime after 7, with P_23 = Cay(Z_23, ⟨2⟩).
- Two-hub small horizons: words of length 3, 4, 5 first miss 17, 35 = 5·7, and 61 / 63 = 7·9. The 7-lane is the
  last thing small two-hub words reach.

**The 3-smooth sum lemma: status.**
- **Needed (Frame Completeness).** Every odd N ≥ 23, other than 21, is an admissible frame sum
  Σ 2^α·3^β of some two-hub word. By the frame identity this is equivalent to the Two-Hub Conjecture.
- **Proved here.** The frame identity and the nested-frame closed form (wave 6).
- **Not proved.** The completeness itself; an Erdős-type sum-of-3-smooth-numbers representation would also need the
  frames to be simultaneously admissible.
- **The decoder shows the gap concretely.** Every tested target has a witness, but the witness is found by search,
  not by a digit rule.

## Wave 8: microcosm and macrocosm, decoder bounds, and the summand graph

**Microcosm 1 + 2^a·3^b (a ≥ 3) and macrocosm 1 + 36·2^k·3^l (`scripts/blocks.c`, FINITE-EXACT).** Restrict to
low-block words: x has two constant blocks and y at most three, length ≤ 26.
- They realize only **94 of 419** microcosm targets up to 3×10^11. They miss 1 + 8·3^b for every b ≥ 5.
- They realize **46 of 349** macrocosm targets.
- They do realize 1297 = 1 + 6^4 = 1 + 36·2²·3², at length 9, the owner's wheel example.

So single-frame words do not carry the microcosm. Its members need several admissible frames, and the 3-smooth
targets are not the frame-simple ones. (All of them below 22,156,805 are two-hub counts by wave 5.) Note that 36 is
also the number of states of the counter: macrocosm values are 1 + (#states)·2^k·3^l.

**Decoder bounds (diagnosis).** The per-state lower bound collapses to 0 at almost every state, because some suffix
always kills every path: S1 is left unclosed, or a junction label is wrong. So only the upper bound prunes, which is
why the cost is about N. Tightening needs a target-aware lower bound: the least count over suffixes that keep the
frame admissible. A cleaner route is interval arithmetic on frames (wave 6): a prefix fixes which frame endpoints
are already committed, and the remaining multiplicities are 1–3 per vertex. That gives two-sided bounds of the form
Σ committed 2^α·3^β·[3^r or the least ratio].

**Summand graph and three strands (ANALOGY with an exact core).**
- The summand graph has arcs X → Z and Y → Z whenever X + Y = Z. Its defining pattern, two arrows into a common
  target generated by the pair, is the additive shadow of transitivity: x → y → z forces x → z.
- Built one node at a time, it splits into 3 strands by how each integer is pulled into existence. In the two-hub
  path S0 h1 S1 h2 S2, the spine is likewise split into 3 increasing strands, and a vertex covered by all three
  ranges is a free ternary digit.
- The Zeckendorf three-diagonal decomposition of triangular numbers is the third lens: T_n split along three
  Fibonacci diagonals.
- The common mechanism: a quantity decomposes along exactly three monotone strands, and the freedom lives where the
  strands overlap. On the tournament side this is exact: H = Σ 2^α·3^β with α and β the doubly and triply
  overlapped vertices.

**Next.**
1. A frame-interval decoder with two-sided bounds, aiming at near-linear search in log N.
2. Find the 3-frame words that realize 1 + 8·3^b.
3. Recompute the owner's wheel rotation 1297 → 357 inside the switching class. The trijection chart predicts that
   rotation preserves the size of the switching-class fibre while moving the frame multiset.

## Next waves (queued predictions)

- **P13.** Classify tournaments with no two disjoint triangles; prove the Ω-complete spectrum has density < 1.
- **P14.** P1 at n = 14 with the full prime set: (5,9)/(6,8), (5,9)/(7,7), (6,8)/(7,7).
- **P15.** The exponent pattern of Gap_b. For b ≥ 60, every m < b with 2m + 1 in the 22-list is a gap. Predict which
  of them survive for m ≥ b, extending THM-5070's 30,000 bound to b = 9–12.
