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

## Wave 9 (overnight): horizons, Berggren, logarithmic length, and the Erdős decoder step

- **P21, horizon ratio (FINITE-EXACT).** At L = 17 the first miss is 73,505,585, so X(17)/3^17 = **0.569**. The
  ratios at L = 14, 16, 17 are 0.487, 0.515, 0.569, still inside the predicted band (0.5, 0.6). L = 18 is running.
- **P23, Berggren depth.** For the 395 primitive hypotenuses c ≤ 2500, minimal two-hub length minus Berggren depth
  ranges over [−25, 4]: deep outer-branch triples are cheap, and middle-branch ones cost at most 4 extra. So the
  prediction holds with constant 4, but the depth is not the right scale.
- **P24, logarithmic length (FINITE-EXACT to 2.2×10^7).** Every odd N ∉ {7, 21} is a two-hub count with word length
  ≤ ⌈log₃(2N)⌉ + 1.
  - Max excess is exactly 1 on the hypotenuse sample.
  - In general it follows from the horizons: X(L − 1) ≥ 0.33·3^(L−1) > 3^L/9 for L ≥ 7. So N > X(L − 1) forces
    log₃(2N) > L − 2.
- **Erdős decoder step.** The 3-smooth decoder lemma (census note) writes an odd n as 2m + 3^k, choosing the parity
  of k so that m is odd.
  - Over all 10,922 configurations with words of length ≤ 6, a one-vertex extension (append or prepend, any labels,
    either hub arc) realizes **2H + 3^k for some k in 47.8%** of configurations.
  - It realizes 2H + 1 in 16.2% and 2H − 1 in 16.2%.
  - For comparison, all three ternary branches occur in only 0.09%.
  So the doubling step of the decoder is locally available about half the time, unlike the ternary step.
- **Prediction P25 (decoder closure, the next proof target).** For every odd n ≥ 23 other than 21, there are k and an
  odd m = (n − 3^k)/2 such that some two-hub realizer of m has a one-vertex extension with count n.
  - With about (4/3)^L realizers per target, this is plausible. Proved, it gives the Two-Hub Conjecture by
    induction on n.
  - *Test:* for n ≤ 10^4, keep all realizers of each m and check the extensions.
  - **REFUTED in the one-vertex form** (`scripts/decclose.c`). With every realizer of length ≤ 10 (these cover
    all m ≤ 10^4), **5,777** of the odd n in [23, 20,000] are not reached. The first are 43, 151, 153, 175 and 251.
    Local availability, about half of all configurations, does not give per-target closure. The surviving form to
    test is closure under bounded multi-vertex extensions, for example two or three vertices, or two Erdős steps
    composed.
  - **Also REFUTED with two-vertex extensions** (`scripts/decclose2.c`: append two, prepend two, or one at each end;
    realizers of length ≤ 9). **3,912** of the odd n in [23, 14,000] are still unreached, with the same first cases.
    For example 43 = 2·17 + 9 is the only odd-m decomposition (2·21 + 1 would need 21), and no realizer of 17
    extends to 43. Conclusion: the 3-smooth antichain decoder does not lift to a local tournament move. Frame
    admissibility is non-local, so a proof must choose whole words, as the branch-and-bound decoder does, rather
    than extend realizers of smaller counts.

## Wave 10: the owner's wheel rotation, reproduced

The owner reported that rotating one two-hub switching class takes the count from 1,297 to 357 without changing its
size. In the word model, a rotation is a cyclic shift of both hub words along the spine: the first spine vertex moves
to the end. `scripts/wheel.c` checks this.
- Among the 724 two-hub words with H = 1297 = 1 + 6^4 = 1 + 36·2²·3² at lengths 8–10, **24** have a cyclic rotation
  with H = 357 = 3·7·17 at the same length.
- Example: L = 8, e = 1, x = 10000110, y = 10100000.

The rotation moves frames, not counts: it changes which (a, b, c, d) are admissible while preserving the size. So
the wheel is an exact frame-shuffling symmetry of the word space, and a natural generator for moving between counts
at fixed size. **P26:** the cyclic group's orbits on words of length L hit each count class in a pattern governed by
the frame multiset; for example, rotation orbits of microcosm counts 1 + 2^a·3^b contain non-smooth counts.

## Wave 11: the six lens and the Collatz ghosts

**The r-tree is the Collatz/ghost tree (PROVED, elementary).** In r = (H + 1)/2 coordinates, the three ternary
children of an attained count H are:

| H-branch | r-branch | map | N mod 6 | realized by |
|---|---|---|---|---|
| 3H | r ↦ 3r − 1 | the ghost: Collatz on the negatives, since −n ↦ 3(−n) + 1 = −(3n − 1) | 3 | C3 product (always available) |
| 3H + 2 | r ↦ 3r | tripling | 5 | one-vertex branch |
| 3H + 4 | r ↦ 3r + 1 | Collatz | 1 | one-vertex branch |

Every r ≥ 2 has the unique balanced-ternary parent round(r/3). So Chan–Pak says this tree, pruned to the branches
that tournaments realize, reaches every r ∉ {4, 11}. The two exceptions are exactly:
- **r = 4 = 3·1 + 1**, the Collatz step from 1 that closes the trivial cycle 1 → 4 → 2 → 1. A single vertex cannot
  take the +1 branch, since 7 is impossible.
- **r = 11 = 3·4 − 1**, the ghost step from the dead node 4.

So "3n + 1 on the positives and 3n − 1 on the negatives, glued at 0" is exactly the r-tree. The ghosts are the
branches each side does not use.

**The prime lens (6k ± 1).** Primes p ≥ 5 are ≡ ±1 mod 6, which are exactly the classes of the two one-vertex
branches. The prime lane of Chan–Pak (THM-4094: Chan–Pak ⟺ every prime p ≠ 7 and every 7p is attained) therefore
lives entirely in the non-ghost branches r ↦ 3r and r ↦ 3r + 1. The ghost branch (the C3 product) only produces
multiples of 3, which are never prime targets. The 36 = 6² states of the two-hub counter are 2 hub orders × 3 × 3 × 2
segment states. The two hub orders are the two sides glued at the hub arc e, which matters only when the hubs are
adjacent in the path (S1 empty): that is the "0" of the gluing.

**Decoder at scale (negative).** Starting at the length law L* = ⌈log₃ 2N⌉ + 1 with randomized child order
(`scripts/decode2hub_fast.c`) still costs about N nodes: 1.4×10^6 at N ≈ 10^6, and no hit at 10^9 within 1.5×10^8
nodes. Bound pruning is intrinsically linear in N. The next idea is meet-in-the-middle on H = c(suffix)·v(prefix),
bucketing prefixes by their dominant coordinates.

**Better lenses: summary and predictions.**
- **P27, prime-branch rule.** For every attained h ≥ 23 other than the dead nodes, the one-vertex branch whose class
  is ≡ ±1 mod 6 is available. This is the threshold-95 rule restricted to the prime classes; it suffices for the
  prime lane. *Test:* from the decomposition16 data.
- **P28, ghost bound.** Of the three r-branches, only the ghost branch (C3 product) is unconditionally available. A
  proof needs the Collatz branch 3r + 1 and the tripling branch 3r at every r ≥ 12. *Test:* realizer bank counts of
  3r and 3r + 1 from two-hub words, which are FINITE-EXACT to 2.2×10^7.

## Wave 12: the glue at 0 as an adjuster (negative), and the route ledger

**Hub-arc flip** (`scripts/jstat.c`, all 1,398,100 word pairs of length ≤ 10). J = H(e=1) − H(e=0) counts the paths
with adjacent hubs, signed by direction; this is the "0" where the two orders are glued.
- |J| = 2 in only **1.0%** of word pairs, |J| = 0 in 0.7%, and |J| ≤ 16 in 8.3%.
- The histogram over small even values is flat.

So the glue is not a fine ±2 adjuster: it scales with the count. A covering induction cannot rely on it except
probabilistically, through some realizer of N − 2 with J = 2.

**Route ledger for full Chan–Pak coverage (status after this session).**

| route | status |
|---|---|
| Busch staircase runs (ternary rule on runs) | open; reduces to the threshold-95 rule |
| Fermat-block cuts from ≤ 9-vertex sources | complete to 398,129, then counting-limited |
| two-hub numeration (Two-Hub Conjecture) | complete to 22,156,805; length law ≤ ⌈log₃ 2N⌉ + 1 |
| frame identity (2^α·3^β sums) | proved |
| Erdős 3-smooth antichain decoder | lemma proved; lift to one- or two-vertex extensions refuted |
| local ternary digits (append one vertex) | refuted (38 of 43,690) |
| exact absolute-target decoder | works; cost about N (no hit at 10^9) |
| hub-arc glue as ±2 adjuster | refuted (1%) |
| Collatz/ghost r-tree lens | reframing; exceptions are r = 4 (Collatz image of 1) and r = 11 (ghost image of 4) |
| pairwise-meeting structure | formulas exact for orders 5–10 |

**Rank of the prefix states (FINITE-EXACT).** The prefix state vectors v_p span exactly **27 = 3³ of the 36**
coordinates, at prefix lengths 4, 6 and 7. Singular values decay quickly: s2/s1 ≈ 0.46, s5/s1 ≈ 0.08 and
s10/s1 ≈ 0.04 at length 7. So H = c·v is far from rank 1, but approximately rank ≈ 5.

**Most promising open route.** Meet-in-the-middle on the bilinear form H = c(suffix)·v(prefix). Compute the rank of
the prefix state vectors {v_p}. Exact divisor matching is ruled out (rank 27). The design that remains is coarse-to-fine
meet-in-the-middle: bucket prefixes by their projections on the top ~5 singular directions; for each suffix,
compute the target cell and check candidates exactly. That would give an exact absolute-target solver for N ~ 10^15
and large-scale evidence for the Two-Hub Conjecture.

## Wave 13: the coarse-to-fine solver, built and measured (negative, with the obstruction identified)

`scripts/beam2hub.c` is a beam search over two-hub words. At each spine position it keeps the B prefixes whose
sampled random completions best straddle N (fraction above versus below, plus log distance). It prunes with the exact
per-state upper bound and checks exactly at the last layer, trying lengths L* to L* + 3.
- With B = 400 and 8 samples it solved **1 of 10** random targets near 10^6, worse than the exact branch-and-bound
  (10 of 10, about 1.4×10^6 nodes).

**Why exact large-N decoding resists generic search.**
- An exact hit must land on N among the completions of a prefix. For the last R letters there are about 2·4^R
  completion values spread over a range of order N, so the hit density per prefix is about 4^R/N. A beam of width B
  then succeeds with probability about B·4^R/N. With 3^L ≈ 2N, the total work is about N/(4/3)^L ≈ N^0.74 at best.
  That is the "many witnesses" regime and still super-polylogarithmic.
- Meet-in-the-middle would give about N^0.63, from 2^L prefixes against 2^L suffixes, but only with a fast exact
  hyperplane query "find p with c·v_p = N". Since v_p spans 27 dimensions (wave 12), no such query is available.
- **Conclusion.** Exact absolute targets at 10^12 and beyond need a *structural* decoder that chooses frames, not
  letters: for example, realize the Erdős 3-smooth antichain of N (census note) frame by frame, with a correction
  term handled by bounded search. The B&B decoder remains the exact tool up to about 10^7, and coverage tables reach
  2.2×10^7.

**Prediction P29.** A frame-level decoder that places one nested frame 2^α·3^β per Erdős-decoder step, plus O(1)
correction letters, solves every odd N ≥ 23 other than 21 with word length O(log N) in polynomial time. This is
the constructive form of the Two-Hub Conjecture, and it would settle Chan–Pak.

## Wave 14: where X^0.98 comes from, the Dirichlet-abscissa lens, and canonical representatives

**Source.** Session 337, [TERRACE_DENSITY](../chanpak_target_first_terrace_density_20261114/TERRACE_DENSITY.md).
Strong tournaments with the fewest directed triangles (Moon's n − 2) form the "first triangle terrace". Their count
is an interval-packing sum over a gap word, and the spectrum is exactly the multiplicative closure of the connected
values. The note proves that Σ H^−s converges at s = 49/50 via R(s)·Q(s) < 1, where
- R(s) = Σ_r (1 + 2^r)^−s runs over Fermat-block counts, and
- Q(s) = Σ_q E_q^−s runs over a kernel with series (9 + 8z + 5z²)/(9(1 − z − z² − z³)).

So the number of distinct counts up to X is O(X^0.98), even though every odd residue mod 2^k is hit.

**Where 0.98 comes from (computed here).** The true critical exponent is the root of R(s)·Q(s) = 1:
**s\* = 0.963577…** At 0.98, R·Q = 0.949 < 1; 0.98 = 49/50 is a convenient rational above s\*, verified there with
integer arithmetic. So the sharp form of the method gives O(X^(0.9636 + ε)).

**Connections inside it.**
- 9E_q = 17, 31, 57, 105, 193, 355, 653, 1201 is a **tribonacci** sequence: each term is the sum of the previous
  three. Its growth rate is the tribonacci constant τ = 1.839286…, the root of the denominator 1 − z − z² − z³.
- 17 = F2 is a Fermat prime. 31 = 2^5 − 1 is a Mersenne prime. **193 = 1 + 3·2^6** is exactly the hub-free
  pairwise-meeting count at order 9 (wave 4).
- R(s) is a **Fermat zeta function** over the counts 1 + 2^r of the Fermat blocks R_r.
- The ratio E_q ≥ (9/5)·E_(q−1) used in the proof is Busch's largest step ratio, b(n+1)/b(n) = 9/5.
- The bit probability 2/3 is w/(1 + w) at w = 2, the C3 block.
- The bound has the shape of a renewal equation: words in the family factor as alternating Fermat blocks and
  tribonacci corridors, Σ_k R^k·Q^(k−1).

**Dirichlet-abscissa lens.** For a construction family F, let σ(F) = inf{s : Σ_(T∈F) H(T)^−s < ∞}. Any family with
σ(F) < 1 realizes O(X^(σ+ε)) counts, so it can never prove Chan–Pak, however well it covers residues.
- The minimum-triangle terrace has σ = 0.9636.
- Two-hub words have σ = log 4 / log 3 = **1.2619**, from 4^L words against counts of size about 3^L. That excess
  (4/3)^L is exactly why their coverage is complete up to 2.2×10^7.

So residue completeness mod 2^k is not a valid local-to-global principle for Chan–Pak. The deciding invariant is σ
together with the distribution of counts.

**Canonical representative (proposal).** For each odd N ∉ {7, 21}, select the **two-hub shortlex normal form**: the
two-hub word of minimal length, then lexicographically least (x, y, e), with the tournament class it determines.
- It exists for every N ≤ 22,156,805, with length ≤ ⌈log₃ 2N⌉ + 1 (wave 9).
- The selection is natural because it lives in the σ > 1 family with an exact 36-state counter.
- The proof route this suggests is **heredity**: show that the normal form of N is obtained from the normal form of
  a smaller count by a bounded frame move (a frame-level decoder). This failed locally for arbitrary realizers
  (wave 9), but normal forms are a much smaller, canonical set.
- *Test (P30):* compute the normal forms for N ≤ 10^5 and check whether removing the last spine vertex of a normal
  form yields a normal form, or a bounded edit of one.
- **P30 REFUTED for shortlex** (`scripts/normalform.c`). Among the 14,722 normal forms with N ≤ 29,523, deleting the
  last spine vertex leaves a normal form in only **32.6%** of cases (first failures 11, 13, 15, 19, 23, 25). Deleting
  the first vertex does so in only **25.5%** (first failures 5, 9, 13, 15, 17). Shortlex minimality is not
  hereditary in either direction. The canonical selection should instead come from the construction: for example
  the word produced by a greedy balanced-ternary or Erdős frame rule. That choice is hereditary by design, and the
  remaining task is to prove that it never gets stuck except at r = 4 and r = 11.

**Frame-level decoder status.** Designed, not built as a provably complete algorithm. Each Erdős step
n = 2m + 3^k would place one nested frame, contributing 2^α·3^β by the nested-frame lemma. Wave 8 showed that
admissible frames come in interacting groups, since low-block words realize only 94 of 419 microcosm targets. So a
correct frame decoder must track the cross-frames each placement creates. The exact branch-and-bound decoder remains
the constructive tool to about 10^7.

## Wave 15: the greedy canonical decoder, the corridor-13 rank-one coincidence, MRB and Leech

**Greedy canonical decoder (`scripts/greedy2hub.c`).** This is the exact branch-and-bound with a fixed rule: at each
spine position, try the child whose upper bound is the smallest value still ≥ N. The first word found, at the
minimal length that has a solution, is the **greedy canonical representative** of N. The rule is deterministic,
and it is hereditary in the sense that the canonical word's prefixes are exactly the greedy choices along its path.
- It solved 20/20 random targets near 10^3 (mean 623 nodes, against 1,727 for the plain decoder) and 20/20 near
  10^6 (mean 3.6×10^5 nodes, against 1.05×10^6). That is about 3× faster, though still super-polynomial in log N.
- ~~What remains for a proof: show that the greedy search never needs to backtrack past a bounded depth.~~
  **Refuted in wave 16:** the final word leaves the first greedy descent within its first 1–3 letters.

**Why only corridors of length 0 or 13 force the terrace count to factor (FINITE-EXACT, q < 30).** In the terrace
interval packing (Session 337), a zero-corridor of length q contributes the 2×2 kernel [[A_q, B_q], [B_q, D_q]]. Its
series over 1 − z − z² − z³ are (1+z²), (1+z) and (1+z+z²). The count factors across the corridor exactly when the
kernel has rank one, i.e. det_q = A_q·D_q − B_q² = 0.

det_q = 0, −2, 3, −1, −4, 8, −5, −7, 20, −18, −9, **47, −56, 0**, 103, −159, …

The zeros are exactly **q = 0 and q = 13**. At q = 13 the kernel is a perfect square:
- A_13 = 2209 = 47², B_13 = 2632 = 47·56, D_13 = 3136 = 56²;
- the two determinants just before it are 47 and −56, and det_14 = 103 = 47 + 56.

**P31:** det_q ≠ 0 for every q ≠ 0, 13. **PROVED in wave 16:** det_q = −T_(−q−3), so the corridor determinant is
the tribonacci sequence run backward.
Connections:
- 13 is in the owner's set {4, 2, 5, 3, 7, 13} (THM-4838).
- 47 is the extra miss of the span-2 ribbons (THM-4880). It is also a prime where 2 has odd order (23), with
  Paley tournament P_47.
- 56 = C(8, 3) = 7·8.

**MRB (prior work, merged).** The repository's MRB product hierarchy (mrb_seven_cutoffs; Session 144 and 162 MRB
notes) shows that the signed product with weights (−1)^n·λ·(n^(1/n) − 1)^α survives iff α > 1/2. That is a sharp
convergence threshold of the same type as the terrace's Dirichlet abscissa s\* = 0.9636 (wave 14). In both, a
sparse-versus-dense dichotomy is decided by where a weighted series crosses 1. A natural test is the MRB-weighted
terrace series Σ (−1)^q·E_q^(−s) and its abscissa.

**Leech (pointer).** The repository's Leech threads are THM-4760 (Fricke twinings of the Leech VOA), THM-484
(involution modulus 24) and THM-4752 (tournament scores in II_{25,1} and the A24 deep hole). A testable link: the
terrace kernels at rank one (q = 13) give a 2-dimensional rank-one quadratic form (47, 56). Its norm
47² + 56² = 5345 = 5·1069 can be compared with the norms of II_{25,1} vectors attached to tournament scores
(THM-4752). This is ANALOGY until such a map is exhibited. (Wave 16: 5345 = T_(−30), forced by the doubling identity.)

## Wave 16: the corridor theorem, the backward tribonacci zero, Singer 13, and the decoder's limits

Script: [`scripts/corridor_tribonacci.py`](scripts/corridor_tribonacci.py) (output in `corridor_tribonacci.out`).
T denotes the tribonacci numbers T_0 = T_1 = 0, T_2 = 1, extended to all of Z by T_(n−3) = T_n − T_(n−1) − T_(n−2).
Put S_n = T_(−n), the backward sequence: 0, 1, −1, 0, 2, −3, 1, 4, −8, 5, 7, −20, 18, 9, −47, 56, **0**, −103, …

**Theorem C (corridor factorization, PROVED).** The corridor kernel has entries A_q = T_(q+3) − T_(q+1),
B_q = T_(q+3) − T_q and D_q = T_(q+3), and

det K_q = T_(q+3)·T_(q−3) − T_q² = **−T_(−q−3)**.

For q ≥ 0 the kernel is therefore rank one exactly when q ∈ {0, 13}.

*Proof.*
1. The entries follow from the generating functions: D − B = z²/(1 − z − z² − z³) and D − A = z/(1 − z − z² − z³).
2. The recurrence gives 2T_q − T_(q+1) = T_(q−3). This turns A·D − B² into T_(q+3)T_(q−3) − T_q².
3. Write T_n = Σ c_k α_k^n. In T_(q−3)T_(q+3) − T_q² the diagonal terms (i = j) cancel. The product of the three roots
   is 1, so α_iα_j = 1/α_k, and what remains is a combination of α_k^(−q), i.e. a backward tribonacci sequence. Two
   such sequences agree once they agree at q = 0, 1, 2.
4. The two-sided sequence vanishes exactly at n ∈ {−16, −3, 0, 1}. This is CITED: García–Gómez–Luca,
   [arXiv:2211.00248](https://arxiv.org/abs/2211.00248), give the exact zero-multiplicity, and Mane,
   [arXiv:2507.11596](https://arxiv.org/abs/2507.11596), names T_(−16) = 0 as the zero not explained by
   identically vanishing polynomials. An independent certificate (Skolem's method at p = 83):
   - the tribonacci matrix has order 287 = 7·41 modulo 83;
   - T ≢ 0 mod 83 outside the four zero classes;
   - in each zero class, v_83(T_(z+287) − T_z) = 1.

   So, by Strassmann (p ≥ 5), the 83-adic interpolation k ↦ T_(z+287k) has exactly one zero in each class. ∎

So "the terrace decoder works backward" is literal: the corridor determinant *is* the backward tribonacci sequence.
Corridor 13 is its sporadic zero T_(−16) = 0, and corridor 0 is the zero T_(−3) = 0.

**The backward zero forces every known tribonacci square (PROVED identities, checked |q| < 400).**

| identity | at a zero of S gives |
|---|---|
| T_(q−1)T_(q+1) − T_q² = −T_(1−q) | the progressions (1, 2, 4) and (3136, 5768, 10609) |
| T_(q−3)T_(q+3) − T_q² = −T_(−q−3) | the corridor progression (81, 504, 3136) |
| T_(q−4)T_(q+4) − T_q² = −4T_(−q−2) | the progression (81, 927, 10609) |
| S_(q−1)S_(q+1) − S_q² = −T_(q+1) | T_16 = S_15² = 56², T_18 = S_17² = 103², T_3 = 1, T_5 = 4 |
| S_(q−3)S_(q+3) − S_q² = −T_(q−3) | T_10 = S_13² = 9² |
| S_(2n) = S_(n−1)² + S_n² + 2S_nS_(n+1) | S_30 = 47² + 56² = 5345 |

The square tribonacci numbers in [OEIS A128911](https://oeis.org/A128911) are 0, 1, 4, 81, 3136 and 10609. None other
exists below index 47,000, but completeness is open. Every one of them is forced this way, either by S_3 = 0 or by
S_16 = 0. The only three-term geometric progressions T_a, T_b, T_c with 3 ≤ a < b < c < 300 are (3,4,5),
(10,13,16), (10,14,18) and (16,17,18), each an identity above evaluated at a zero. The corridor progression
(3^4, 2^3·3^2·7, 2^6·7^2) is {2,3,7}-smooth, and the {2,3,7}-smooth T_n below index 3000 are only 1, 2, 4, 7, 24,
81, 504 and 3136.

**P32:** every square in the two-sided tribonacci sequence is forced by a zero through such an identity. This would
make A128911 complete.

**One sporadic zero in the whole k-bonacci family (FINITE-EXACT, k = 3..15, indices ≥ −4000).** Since
(x − 1)(x^k − x^(k−1) − … − 1) = x^(k+1) − 2x^k + 1, each two-sided k-bonacci sequence has a structured family of
negative zeros: −3 − j(k+1) − i, for 0 ≤ j ≤ k−3 and 0 ≤ i ≤ k−3−j. These number (k−1)(k−2)/2. Every negative zero
found is structured except one, tribonacci's T_(−16). So the corridor-13 factorization is the single sporadic event
of the family.

**Why 13: a Singer cycle of PG(2,3) (ELEMENTARY).**
- x³ − x² − x − 1 is irreducible over F_3, and its roots have norm αβγ = 1. They therefore lie in the norm-one
  subgroup of F_27^×, which has prime order (27 − 1)/2 = 13.
- So the tribonacci period mod 3 is 13 = 1 + 3 + 9 = |PG(2,3)| = Φ3(3), and the matrix acts on the projective plane
  of order 3 as a Singer cycle.
- The zeros −3 and −16 lie in one class mod 13: corridor 13 is exactly one 3-adic period.
- Mod 5 the period is 31 = Φ3(5) = |PG(2,5)|, the base-3 miss of wave 2.
- ANALOGY: the Chan–Pak exceptions are 7 = Φ3(2) = |PG(2,2)| (Fano, P_7) and 21 = Φ3(4) = |PG(2,4)|. So the
  corridor's 13 = Φ3(3) is the middle term of 7, 13, 21.
- In the ghost language of wave 11, the backward ("ghost") side of the recurrence has a sporadic zero that controls
  the forward side: the corridor factorization and all three large squares.

**{4,2,5,3,7,13} and the Leech lattice.** Already in the repository: {4,2,5,3,7,13} = {d + 1 : d | 12}, with primes
vSC(12) (THM-4838). Today's joint_metagraph_dyadic note shows these are the levels whose Hauptmodul eta exponent
24/(N − 1) is even. New exact facts (ELEMENTARY; checked for 2k ≤ 60):
- θ_Leech = E_12 − (65520/691)Δ with 65520 = 24·2730 = 2^4·(2^12 − 1) = 2^4·3^2·5·7·13. The kissing number is
  196560 = 2^4·3^3·5·7·13, so the owner's set is read straight off it: 2^4, then 3, 5, 7, 13.
- The numerator of the Eisenstein constant 4k/|B_2k| equals 2^a(2^(2k) − 1) exactly for 2k ∈ {2, 4, 6, 12}, giving
  24, 240, 504 and 65520.
- These are the weights where every prime of 2^(2k) − 1 has (p − 1) | 2k. Their primes 2k + 1 = 3, 5, 7, 13 are
  Φ2(2) = F0, Φ4(2) = F1, Φ3(2) and Φ12(2): the two Busch bricks, the first Chan–Pak exception, and the corridor.
- NUMEROLOGY: the middle of the corridor progression, T_13 = 504 = 2^3(2^6 − 1), is also the E_6 constant.

**MRB (papers given to earlier sessions, merged).** The PDFs are not in this checkout. Their summaries are in
[fermat_self_duality §8](../fermat_self_duality_fewcycle_coincidences_20261104.md) and session31's contour audit.
- *Greedy residue expansions with an inverse-state map* (MRBaid2, MRBgreedy3, Burns' "A Convergent Greedy Residue
  with Geometry"). There the initialisation question stays open. The two-hub analogue is measured below. Greedy
  order is not an autonomous map, and everything rests on initialisation, as in MRBaid2. Unlike MRB, the state space
  is finite (36 states), so the search is exact.
- *Mirror pair* (NUMEROLOGY): the only rational solution of x^(1/x) = y^(1/y) with x ≠ y is (2, 4), with
  2^4 = 4^2 = 16, which gives r_2 = 4 in the mirror-base constant B. The sporadic zero sits at index −16.
- *Threshold 1/2* (ANALOGY): the MRB product survives iff α > 1/2 (mrb_seven_cutoffs). |det K_q| = |T_(−q−3)| grows
  like τ^(q/2), because the complex roots have modulus τ^(−1/2). In both, 1/2 is a square-root cancellation
  exponent.

**The greedy decoder is a search order, not a decoder (FINITE-EXACT; `scripts/greedy_tail.c`,
`scripts/greedy_logcenter.c`, 1000 consecutive odd targets per window).**

| child order | window | first-descent success | mean nodes at final length | correct greedy head |
|---|---|---:|---:|---|
| tightest upper bound ≥ N | [23, 2023) | 0.2% | 287 | usually 0–2 letters |
| tightest upper bound ≥ N | [10^4, 10^4 + 2000) | 0% | 3,886 | 0–2 of 9–10 |
| tightest upper bound ≥ N | [10^6, 10^6 + 2000) | 0% | 260,860 | 0–2 of 13–14 |
| N nearest the log-centre of [lo, hi] | [23, 2023) | 2.1% | 335 | 1–4 |
| N nearest the log-centre of [lo, hi] | [10^4, 10^4 + 2000) | 0% | 2,887 | 2–3 |

The max and min bounds over the 36 states are too coarse to steer the search. A real decoder needs finer
invariants of the attainable set, such as its residues modulo 3^j from each state.

**P33:** steering by the attainable residues mod 9 from each state (a finite automaton) raises first-descent success
above 50% at 10^4.

## Wave 17: the residue-steered decoder, horizon 2.26×10^8, and four local proof routes tested

**Horizon at word length 18 (FINITE-EXACT, `twohub_bits`).** Every odd N < **226,079,703** other than 7 and 21 is a
two-hub count with a word of length ≤ 18. This is ten times wave 5's bound. The ratio X(L)/3^L keeps rising, which
fits P21: 0.487 (L = 14), 0.515 (16), 0.569 (17), 0.584 (18).

**The residue-steered decoder (`scripts/resdec2hub.c`; P33 REFUTED, but the residues pay as a filter).**
- *Oracle study* (`scripts/oracle_steering.c`). All words of length 9 are enumerated and grouped by their first 2 or
  3 letters. For every odd N in the horizon we ask whether the group a statistic ranks first contains a solution:

  | statistic of the subtree | first 2 letters | first 3 letters |
  |---|---:|---:|
  | random feasible prefix | 91.5% | 73.9% |
  | exact residue density mod 9 | 83.5% | 58.2% |
  | exact residue density mod 27 | 79.3% | 62.7% |
  | log-centre of [min, max] | 98.9% | 84.1% |
  | normal density at N (exact mean and variance) | **99.1%** | **88.3%** |

  Residues are worse than random as a steering signal. A subtree whose values sit in few residue classes is a rigid
  subtree, and rigid subtrees rarely contain N. The magnitude statistics carry the signal.
- *Built decoder.* Children are ordered by the normal density of N. The exact mean is v·μ_r and the exact second
  moment is vᵀQ_r v, with Q_r = ¼ Σ_ℓ M_ℓ Q_(r−1) M_ℓᵀ. In the last 6 letters, every prefix is checked against the
  deduplicated table of suffix responses modulo 216 = 8·27, an exact filter with no false pruning.

  | targets | interval bounds only | + residue filter | residue-steered (moment + filter) |
  |---|---|---|---|
  | 100 near 10^6 | 757,891 nodes, 26.7 s | 288,527, 21.4 s | **99,376, 10.9 s** |
  | 30 near 10^7 | 8,433,406, 89.0 s | 3,234,428, 74.4 s | **898,797, 35.5 s** |
  | 10 near 10^8 | – | 47,389,003, 419.6 s | **14,689,904, 226.3 s** |

  The decoder is 1.9–2.5× faster, with 3–8× fewer nodes. First-descent success is still **0%** from 10^4 upward
  (10.7% below 623), so P33 is refuted: the decoder is still a search.
- *Why residues cannot steer from the top.* The two-hub series has minimal linear dimension **17**: Hankel rank 17
  over Q, mod 3 and mod 5, with prefix span 28 and suffix span 20 (`scripts/minimal_rank.py`). The prefix automaton
  has 129 states mod 2 and 458,754 mod 4, but more than 2·10^7 mod 3, 5, 7, 8, 9, 16 or 27
  (`scripts/residue_automaton.c`). Exact residue densities are therefore only available near the leaves.

**Four local proof routes, tested.**
1. *Per-value ternary fertility.* Call n fertile if some realizer w has single-letter insertions, at any position,
   realizing all three of 3n, 3n+2 and 3n+4. If every n ≥ 23 were fertile, induction would prove the conjecture
   (balanced ternary in r = (n+1)/2). **FALSE:** only 323 of 1231 values are fertile with realizers of length ≤ 9
   (`scripts/fertility.c`). For example, no realizer of 35 has an insertion child 107 or 109. Binary fertility
   (2n ± 1) fails too: 256 of 1231. This strengthens P18's refutation from configurations to values.
2. *Affine numeration by a common eigenvector.* Each letter matrix has spectrum 0⁸ 1¹⁶ 2¹⁰ 3². A left vector
   tripled by every letter would make H affine in one register, x ↦ 3x + digit. **None exists:** the common left
   eigenspaces have dimension 0 for λ = 2 and 3, and 4 for λ = 1, spanned by the goodness symmetries
   (`scripts/common_eigen.py`).
3. *Smaller families.* Measured misses among odd numbers:
   - one hub: 45.6% of odd N ≤ 2000;
   - x arbitrary, y with a single switch: 24.9% of odd N ≤ 10^5 (m ≤ 16);
   - three-letter alphabets: {00, 01, 10} covers every odd N in [23, 102,083) at L ≤ 13, but its relative horizon
     X/3^L drifts down (0.153, 0.086, 0.052, 0.064, 0.070, 0.064, 0.065, 0.037 for L = 8..15). By converse symmetry
     {01, 10, 11} behaves identically, while {00, 01, 11} and {00, 10, 11} already miss 175.

   The fourth letter is what keeps the horizon at a fixed fraction of 3^L.
4. *Sparse digits.* Decoded realizers of N ≈ 10^5–10^6 have about 1,000–2,000 frames
   (`scripts/frames.py`, checked against brute force on 150 random words). Their weights lie on the staircase
   α + β = m − 4, and the multiplicities grow with α. So a realizer is a dense frame sum, not a digit expansion.

   One lens remains: frames nested inside a fixed outer pair (a, d) weigh 2^(A−g)·3^(g−1) by their gap g = c − b, a
   rational base-3/2 digit system in the sense of Akiyama–Frougny–Sakarovitch, with digit = number of (b, c) pairs at
   gap g. Whether two-hub words realize enough digit vectors is open (**P34**).

**Route ledger after wave 17.** Closed:
- single-vertex ternary and binary recursions (P18; per-value);
- decoder closure (P25), shortlex heredity (P30), bounded backtracking (wave 16), residue steering (P33);
- annealing, beam search, the hub-arc ±2 adjuster;
- uniform affine numeration in one register.

Open:
- the horizon law X(L) ≥ c·3^L (P21; data 0.58 at L = 18);
- the rational-base frame lens (P34);
- port-wedge induction (busch_staircase §6.3);
- the Shallow-Cut Conjecture.

## Wave 18: the interval-hub family and the bi-quadratic law

**Bi-quadratic frame law (PROVED, elementary).** For two-hub words, H(x, y, e) is a multilinear polynomial of
degree ≤ 2 in the x-bits and degree ≤ 2 in the y-bits, so total degree ≤ 4.

*Proof.* In the frame identity each admissible frame has a weight ∏ r(i) that does not depend on the labels, times
an admissibility indicator. For the order (u1, u2) the indicator is (1 − x_a)·x_b·(1 − y_c)·y_d; for (u2, u1) it is
(1 − y_a)·y_b·(1 − x_c)·x_d. Empty runs drop literals. ∎

All 1024 words with m = 5 attain exactly the degrees (2, 2, 4) (`scripts/ihub_anf.py`). Example, with one hub, an
empty interval and e = 0, where X(x) is the binary value with x_0 most significant:

H = 1 + 2X(x) − Σ_(i<j) 2^max(1, j−i−1) x_i x_j.

For instance x = 11 gives 1 + 6 − 2 = 5, checked by listing the 5 paths.

**Interval-Hub Conjecture (NEW; FINITE-EXACT to 26,610,535).** Every odd N other than 7 and 21 is the count of
TT_m + u1 + u2, where u1's out-set on the spine is arbitrary and u2's out-set is a single interval [p, p+q)
(`scripts/intervalhub.c`).

| m ≤ | 12 | 14 | 16 | 18 | 19 | 20 | 21 | 22 | 23 |
|---|---|---|---|---|---|---|---|---|---|
| first odd miss above 21 | 17,789 | 84,271 | 345,143 | 1,467,007 | 2,064,419 | 4,007,635 | 8,070,555 | 15,785,407 | 26,610,535 |
| that miss ÷ 2^m | 4.34 | 5.14 | 5.27 | 5.60 | 3.94 | 3.82 | 3.85 | 3.76 | 3.17 |

The horizon grows like 3–5·2^m: **binary** scale, against the ternary 3^L/2 of unrestricted two-hub words. So m ≈
log2 N − 2 spine vertices suffice. The family is far smaller (2^m·m² words, not 4^m) and still complete through
2.66·10^7. The last ratio (3.17 at m = 23) may signal a slow decline; P35 needs it only to stay above 2.

Robustness (`scripts/intervalhub_mult.c`): with m ≤ 18 and four or more spare spine vertices, every odd N has at least about
50 realizers in each dyadic window up to 2^15. The minimum falls to 1 only at the horizon edge.

Which part of the family matters (`scripts/intervalhub_qcap.c`):
- *Intervals of length ≤ 1* (u2 a near-sink): the coverage **stops at 7,477** for every m ≤ 20.
- *Lengths ≤ 2 and ≤ 3*: the horizons reach 1,174,699 and 2,028,565 at m ≤ 20.
- *An interval touching an end of the spine*: misses 25%.
- *No interval* (one hub): misses 45.6%.

The interior interval is what completes the spectrum.

**Structure for a proof (PROVED / FINITE-EXACT).**
- For fixed (I, e), H is a quadratic pseudo-Boolean function of x (`scripts/ihub_coef.py`).
- Before the interval the prefix vectors span only **6** dimensions (`scripts/ihub_phase_rank.py`), so that part of
  x acts through six one-hub statistics.
- Example, m = 8 and I = [3, 5):
  - the quadratic coefficients among the first three bits are 23 times the one-hub ones (−46 = 23·(−2));
  - their linear coefficients are 37·2^(5−i).

  The interval and the suffix set the scales, and the prefix supplies a scaled one-hub numeral.

**No local binary recursion (FALSE).** Single-vertex insertions inside the family realize only 4,917 of the 42,113
odd N′ ≤ 84,271 from a realizer of their odd parent n with N′ = 2n ± 1 (`scripts/intervalhub_inherit.c`). The
doubling is global, as everywhere else. HYP-9469's increment law (one added vertex either keeps h or raises it to at
least (3h+1)/2) explains why local steps are coarse.

**P35.** The Interval-Hub Conjecture holds with m ≤ ⌈log2 N⌉ − 1.

**One interval alone gives only short runs (FINITE-EXACT, `scripts/intervalhub_runs.c`).** For each single
(p, q, e), the values over all x ∈ {0,1}^m contain at most 15, 15, 21, 21 and 23 consecutive odd numbers at
m = 8, 10, 12, 14 and 16. Example: [14103, 14147] at m = 16, p = 13, q = 2. Completeness is a property of the
*union* over intervals, so a proof must move the interval as a digit, not fix it.

**P36 (proof route).** Write each H(·; I, e) as A + B·(one-hub numeral of the prefix) plus quadratic carries, with
explicit scale factors (B, C) from the 6-dimensional prefix phase. Then prove that the union over intervals of the
images of these scaled numerals covers every odd N ≥ 23, by a digit decoder with carries. This is a binary-scale,
low-dimensional numeration problem: the smallest complete family found so far.

## Wave 19: a residue-steered source decoder for root-change words, and what the final equality needs

**Prior art first (KNOWN; re-derived this wave before the repository search caught up).** The family built here as
"caterpillar tournaments" is the repository's **O/Z/J root-change grammar** (THM-4904; Sessions 197–205; THM-4888,
THM-5148, HYP-9494, HYP-9677). Its letters are:
- O or Z: add a spine vertex that the current apex beats (O), or that beats the apex (Z);
- J: add a new apex that beats the old one.

These are this wave's moves a, b, c. The following are all KNOWN:
- the root triple evolves linearly (join law THM-4866 (B); block matrices THM-5148 (A));
- the one-apex count is g(0) = 1, g(2v+1) = g(v), g(2v) = g(v) + v + (v mod 2), with the binary pair energy of
  THM-5093;
- JZ = OJ (Session 205, THM-4888 (D) 1) is the identity called "factor migration" in this session's scratch work;
- the slope (port) set of a fixed count is a finite union of geometric progressions (THM-5129 (A));
- the grammar misses exactly 29 odd numbers, the largest 2647, verified through 10^9 (HYP-9494 (D)). The depth-first
  search [`scripts/caterpillar.c`](scripts/caterpillar.c) reproduces this through 2·10^5.

**Corrections to earlier waves.**
- Wave 18's one-hub closed form is THM-4876 (6).
- Wave 17's Hankel rank 17 of the two-hub series reproduces THM-4878 (G(2,0), Session 191).
- Wave 18's phrase "the smallest complete family found so far" is withdrawn. The root-change grammar is smaller: a
  4-dimensional state and three letters.

**NEW: the residue-steered source decoder ([`scripts/catdec.c`](scripts/catdec.c)).** This is the owner's request
made literal. The reader handles its factors only by handing them to the source; the source then chooses its digits
by residues; and the run ends in one exact equality.
1. *Peel the readers.* Reader blocks J O^i Z^j are removed from the end of the word. They act on the source's root
   triple (H, E, λ) by integer matrices, so the target reads N = A·g(U) + C·E(U) + B·(1 + U). Here U is the
   one-apex source word read as a binary number (O = 1).
2. *Put the factors into the source.* The equation A·g(U) + B·U = T is solved from the lowest bit of U upward by
   g's recursion. A last letter O sends (B, T) to (2B, T − B). A last letter Z sends it to (2B + A, T − A·[the
   preceding letter is O]). The reader's factor B is carried inside the source problem, and the parity of the next letter is
   the residue requirement passed down.
3. *Final absolute equality.* The recursion closes only at U = 0 with T = A exactly. Once the window
   [(T − A)/(A + B), (T − A)/B] holds at most 8 integers, these are checked directly (g(U) ≤ U + 1). Every output
   word is re-run through the root-state recurrence.

| targets | decoded | time per target | solver nodes per target |
|---|---|---|---|
| every odd N ≤ 3001 | all but the grammar's 29 misses, 1877 and 1897 (depth 7: decoded once 7 readers are allowed), and 2377 (needs a reader outside J O^i Z^j) | – | – |
| 100 random, near 10^6 | 100 | 1 ms | 1.6·10^4 |
| 100 random, near 10^9 | 100 | 31 ms | 3.5·10^5 |
| 50 random, near 10^12 | 50 | 1.6 s | 1.4·10^7 |
| 12 random, near 10^15 | 12 | 13 s | 9.8·10^7 |

The cost grows like √N. A reader chain with factors (A, B) covers values at density about 1/(A + B), and its
source search costs about √(NA)/B. The random samples at 10^12 and 10^15 lie beyond the 10^9 census, so they are
new (sampled, not exhaustive) evidence for HYP-9494 (D).

**NEW (FINITE-EXACT): monotone readers suffice.** Take a one-apex source of any length followed by readers
J O^i Z^j with i, j ≤ 3 (16 block types). This family attains every odd N in [2649, 10^6]
([`scripts/readeralpha.c`](scripts/readeralpha.c)). It has 30 misses below 2649: the grammar's 29 and 2377.
Smaller block sets cover less:
- i, j ≤ 2: 4969 is also missed;
- i, j ≤ 1: 71 misses below 20,000, the largest 19,781;
- a single block type: thousands of misses.

The extension to 10^7 was still running at this commit (P37 below).

**NEW: the multiplicity regime ([`scripts/repcount.c`](scripts/repcount.c)).** Allow at most 4 readers with
i, j ≤ 3 and count representations (chain, source word) of each target over 40 consecutive odd targets per scale:

| scale | 3·10^3 | 10^4 | 3·10^4 | 10^5 |
|---|---:|---:|---:|---:|
| mean | 157 | 186 | 174 | 197 |
| minimum | 13 | 12 | 3 | 41 |

The mean is flat. A fixed finite set of readers is a critical regime: each chain contributes a fixed density
1/(A + B), so the mean cannot grow. Growth must come from chains whose length grows with log N. HYP-9494 (D)'s depth
table points the same way: 1877 and 1897 need 7 root changes.

**NEW: inside the grammar, runs of root changes are the wave-16 corridors
([`scripts/corridor_factor.py`](scripts/corridor_factor.py)).** A run of k root changes acts on (H, E, λ) by J^k,
where J = [[1, −1, 1], [1, −1, 0], [1, 0, 1]] has characteristic polynomial x³ − x² − x − 1.
- *Determinant (PROVED).* The (H, λ)-block of J^k has determinant −T(−k−2), where T is the two-sided tribonacci
  sequence with T(0) = T(1) = 0, T(2) = 1. This is checked for k ≤ 300; both sides satisfy linear recurrences of
  order at most 12, which makes the check a proof. It is wave 16's corridor determinant −T(−q−3) at q = k − 1.
- *Rank one.* So the block is rank one exactly at k = 1 and k = 14, by the cited zero set {−16, −3, 0, 1}. At
  k = 14 the block is (47, 56)ᵀ(47, 56), built from the tribonacci squares 2209 = 47² and 3136 = 56².
- *Consequence.* After an O, where E = 0, a run of exactly 14 root changes gives H = 47·(47H₀ + 56λ₀) and
  λ = 56·(47H₀ + 56λ₀). This was checked on 1000 random words. The corridor reads the whole source through a single
  linear form.
- *Forced factors.* Runs of other lengths entered after an O force smaller common factors of H (k ≤ 60):
  - 3 at k = 3, 29, 42, 55;
  - 5 at k = 9, 40;
  - 47 at k = 14, 60;
  - 159 at k = 16;
  - 499 at k = 54.

  A decoder aiming at a prime target must avoid these run lengths after an O.

**The final absolute equality: what remains (ANALYSIS).** With the pair energy Ψ (THM-5093's Q), the one-apex
count is g(U) = 1 + U − u₀ − Ψ(U). So one chain's equation is

N = (A + B)(1 + U) − A·Ψ(U) + C − (A + C)·u₀.

Split U = 2^m V + W by Session 193's concatenation law. This gives

N = V·[(A + B)·2^m − A·rev_m(W)] − A·Ψ(V) + c(W).

Every target is therefore a *near-product* of the high source bits V with a quantity read from the reversed low
bits W. The source chooses its digits from both ends at once, and the absolute equality becomes a fixed-point
condition for V inside a window of O(A) integers. Consequences:
- the event "chain R represents N" has probability about 1/(A_R + B_R), so no single chain is complete;
- a proof must control the union over chains, and the critical-regime data above say it must use chains whose length
  grows like log N;
- by THM-5129 (A), no proof can fix a source count and free its port: the (count, port) pairs of sources are always
  graph-like.

This matches the repository's reading of Chan–Pak as cofiniteness of an explicit matrix-semigroup orbit (HYP-9640,
HYP-9677). There thin-orbit methods (Bourgain–Kontorovich) give density one at best
([tilted ensemble §7](../tilted_ensemble_local_limit_20261110/README.md#7-typed-connections-and-named-principles)).
The decoder settles each instance of the final equality cheaply; it does not prove the equality in general.

**P37 (pre-registered).** The sub-language "one-apex source, then readers J O^i Z^j with i, j ≤ 3" attains every
odd N in [2649, 10^7].

## Wave 20: the two-root-change conjecture, HYP-9677 to 2^24, and a run-cap trichotomy

**The two-root-change family F_≤2 (HYP-9494 (D)).**
- *Census.* [`scripts/f2miss.c`](scripts/f2miss.c) and the bitset version [`scripts/f2census.c`](scripts/f2census.c)
  enumerate root-change words with at most two root changes. They reproduce HYP-9494 (D) exactly: 140 misses below
  10^8, the largest 3,780,967, none in (5·10^7, 10^8].
  The extension to 4·10^9 (P38) was still running at this commit.
- *Decoder.* [`scripts/f2dec.c`](scripts/f2dec.c) allows arbitrary reader blocks J w and at most two root changes. It
  enumerates the readers by total length and solves the source with the readers' factors carried inside.

  | random odd targets | decoded with at most two root changes | time per target |
  |---|---|---|
  | 10 near 10^9 | 10 | 0.04 s |
  | 50 near 10^12 | 50 | 1.4 s |
  | 20 near 10^15 | 20 | 30 s |

  About half of random targets are already one-apex values (no root change), so random samples mainly test the
  decoder. They support HYP-9494's conjecture only where they land.

**Bounded runs (HYP-9677).** An exact pruning makes the census cheap. A word J Z… reaches the same root state as
O J… (JZ = OJ, Session 205), and the second form leaves more run budget. So a Z directly after J may be skipped
whenever the run before that J was below the cap. With this, [`scripts/runcap2.c`](scripts/runcap2.c) runs about 45×
faster at 2^19 with an identical result.

*HYP-9677 extended (FINITE-EXACT).* With every run of spine letters of length at most 5 (the 63 matrices), every odd
N in (10,391, 2^24 = 16,777,216] is attained. The misses are exactly the 49 known ones; the earlier record was 2^19.

*A trichotomy in the run cap* ([`scripts/runtail.c`](scripts/runtail.c)). Misses per dyadic window:

| max run | blocks | misses | behaviour of the per-window misses |
|---|---:|---:|---|
| 1 | 3 | 1,940,377 below 2^22 | the missed share rises to 93.5% of the odd numbers of [2^21, 2^22): density tends to 0 |
| 2 | 7 | 350,882 below 2^22 | the missed share falls slowly (27.5%, 23.8%, 20.9%, 19.7%, 19.4%, 17.5%, 14.4% for 2^15 to 2^21) while the count grows 1.65–2× per doubling |
| 3 | 15 | 1,867 below 2^26, largest 26,538,691 | peak 285 in [2^15, 2^16), then 277, 188, 136, 82, 56, 36, 12, 1, 1, 0 (2^16 to 2^25) |
| 4 | 31 | 127 below 2^25, largest 130,247 | none in (130,247, 2^25] (P39, CONFIRMED) |
| 5 | 63 | 49 below 2^24, largest 10,391 | none in (10,391, 2^24] |

So the run cap shows the same trichotomy that HYP-9480 (D) found for the number of root changes:
- runs ≤ 1: the density tends to 0;
- runs ≤ 2: the density apparently tends to 1, with infinitely many misses expected;
- runs ≥ 3: apparently finitely many misses.

The middle regime is not visible in word counts. Words per odd number (trailing O's excluded) grow at every cap: from
2^14 to 2^22 they rise from 38 to 212 for runs ≤ 2, from 104 to 1091 for runs ≤ 3, and from 169 to 2391 for runs ≤ 4
([`scripts/rundim.c`](scripts/rundim.c)). Yet runs ≤ 2 still leaves about 14% of the odd numbers of [2^21, 2^22)
unhit. Collisions dominate, so word counts, and the Poisson heuristics built on them, are a poor proxy for coverage
in these structured families.

**Fewer matrices ([`scripts/runcap_noZJ.c`](scripts/runcap_noZJ.c), [`scripts/runwhite.c`](scripts/runwhite.c)).**
Forbid Z directly after J; by JZ = OJ, such a word can be rewritten only if the run before the J has room.
- *Cap 4.* This leaves 16 blocks: w empty or beginning with O, |w| ≤ 4, after a first run of any word of length ≤ 4.
  Through 2^24 they have 147 misses, with the same largest miss, 130,247, as all 31 cap-4 blocks.
- *Greedy removal (searched at 2^21).* OZOO and OOZO can be removed together, leaving 14 blocks that still have the
  largest miss 130,247 (checked to 2^24: 170 misses, the largest 130,247). Removing OOZZ as well moves the largest miss to 145,147. More generally, removing any one
  of the six remaining length-4 blocks (13 blocks) keeps the census apparently cofinite through 2^22, with the
  largest miss between 145,147 (OOZZ removed) and 418,429 (OZZO or OOOZ removed).
- *Greedy to 11 blocks ([`scripts/greedy_blocks.py`](scripts/greedy_blocks.py)).* Continuing the greedy removal,
  OOZZ, OOOO and OZZO go, leaving 11 blocks: J followed by ε, O, OO, OZ, OOO, OOZ, OZO, OZZ, OOOZ, OZOZ or OZZZ,
  after a first run of length at most 4. Through 2^26 = 67,108,864 they have 408 misses, the largest 418,429, and
  none in [2^19, 2^26]. Every further single removal leaves misses in [2^19, 2^21].
- *The singular blocks are load-bearing ([`scripts/runterm.c`](scripts/runterm.c)).* Allow J only after a complete
  block. With all 11 blocks as complete blocks the census is unchanged (408 misses, largest 418,429, below 2^21).
  Keeping only the invertible ones (J alone, det 1, or J w with w beginning with O and ending with Z, det 2^|w|)
  misses 34% of the odd numbers below 2^21 (7 blocks), or 24% with OOZZ added. So the target cannot be moved into a
  semigroup of invertible matrices this way. This is one more reason thin-group methods do not transfer directly.
- *Cap 3.* The 8 Z-free blocks sit near the critical line: per-window misses 579, 566, 517, 477, 452, 378 for 2^16
  to 2^21.

So the repository's 63-matrix conjecture (HYP-9677) appears to hold with 14 matrices after the first factor beyond
130,247, and with 11 beyond 418,429. In the owner's language this is the reader's factor put into the source: the
Z-free normal form moves every leading Z of a reader into the preceding run as an O (JZ = OJ), and at cap 4 this
costs nothing in coverage.

*Heuristic for the threshold.* A block J w with |w| = n multiplies a typical state's count by about c·2^n. The
pressure sum Σ_w (c·2^|w|)^(−s) at s = 1 is (L + 1)/c, so it crosses 1 near L + 1 ≈ c. With c ≈ 3–4 this puts the
critical cap at L ≈ 2–3. This is the Hausdorff-dimension condition of Zaremba-type problems; it is a heuristic, not
a theorem. Two facts block a direct transfer of Bourgain–Kontorovich: most block matrices are singular (J w is invertible only when w
is empty, with det 1, or begins with O and ends with Z, with det 2^|w|; checked for |w| ≤ 7), so the semigroup is not
inside a group, and density one would still fall short of
cofiniteness.

**Corridors force primes along tribonacci periods (FINITE-EXACT, [`scripts/corridor_primes.py`](scripts/corridor_primes.py)).**
After an O (E = 0), a run of k root changes gives H = A_k·H₀ + B_k·λ₀. Here (A_k, B_k) are the H- and λ-entries of
the first row of J^k. So a prime p divides every such count exactly when that row is a multiple of (0, 1, 0) mod p.
Scanning one full period of the row for each prime p < 400:
- the forcing lengths, when they exist, are a single residue class modulo the period, which is the tribonacci period
  of p. Examples: 3 at k ≡ 3 (mod 13); 5 at k ≡ 9 (mod 31); 11 at k ≡ 72 (mod 110); 47 at k ≡ 14 (mod 46); 53 at
  k ≡ 16 (mod 52);
- 31 primes below 400 can be forced (3, 5, 11, 23, 31, 43, 47, 53, 59, 71, 89, 101, ...), and 47 never are (2, 7,
  13, 17, 19, 29, 37, 41, ...).

The modulus 13 for p = 3 is wave 16's Singer cycle of PG(2, 3). If J mod p has projective order p² + p + 1 (below
60: p = 3, 5, 23, 59), it permutes all points of PG(2, p) in one cycle, so the orbit of [e₁] meets [e₂] and p is
forced. Irreducibility of x³ − x² − x − 1 alone is not enough: at p = 31 and p = 37 the projective order is a third
of p² + p + 1, and 31 is forced while 37 is not.

**Pre-registered.**
- **P38.** F_≤2 misses no odd N in (3,780,967, 4·10^9].
- **P39 (CONFIRMED).** Runs ≤ 4 (31 matrices) miss no odd N in (130,247, 2^25]: 127 misses below 2^25, the largest
  130,247. This is HYP-9677 with 31 matrices in place of 63.
- **P40 (REFUTED, by one value).** Runs ≤ 3 were predicted to miss no odd N in (13,009,019, 2^26]. The census finds
  exactly one such miss, 26,538,691, and none in (26,538,691, 2^26].
- **P40′.** Runs ≤ 3 miss no odd N in (26,538,691, 2^27].

## Wave 21: combining strategies — what covers what, and where the 11-matrix orbit is weakest

Owner prompt: work toward universal placement, combine strategies so that each covers the others' gaps, and think about
the cofiniteness of the 11-matrix orbit. **Chan–Pak OPEN.** The short answer of this wave: the strategies in hand fail
on the same targets (primes and numbers coprime to 3 and 5), so combining them improves thresholds and constants but
not the regime. Unit-step windows cannot be stacked or widened by single back arcs.

**One formula behind both placement programs (KNOWN, THM-4866 (B)).** Join a rooted source (A, a) to a rooted reader
(B, b) by one back arc b → a, with A beating B otherwise. Then

N = H_A·H_B − E_A(a)·S_B(b) + λ_A(a)·ρ_B(b),

where E is the number of Hamiltonian paths ending at the root, S the number starting at it, and λ, ρ the left and right
ports (pairs of paths covering the tournament, split at the root).
- The root-change decoder of waves 19–20 is this formula iterated, with one-apex blocks as readers.
- The other session's free-address windows are this formula with a one-apex wing as source. For example,
  908u + 9g(u) + {752, …, 758} uses four readers with core count 9 whose slopes coincide at 908 after reindexing the
  source address.
- "Universal placement" asks that this pairing, over realizable source and reader profiles, reach every large odd N.

**Three obstructions to building exact windows by single back arcs (PROVED, elementary).**
- *(a) Ports dominate counts.* For every rooted tournament, ρ_B(b) ≥ H_B and λ_A(a) ≥ H_A. Split each Hamiltonian
  path at the root: this gives distinct pairs.
- *(b) Width of a fixed-profile window.* Fix the source and vary the reader over tournaments with the same
  (H_B, ρ_B). Then N moves only through −E_A·S_B. Consecutive odd counts need E_A = 1. Since S_B lies in [0, H_B] and
  Rédei fixes its parity, the window has at most ⌊H_B/2⌋ + 1 values. The other session's windows use a different,
  tilted mechanism: readers with different profiles whose coefficients coincide after a reindexing, or cancel on a
  special line of sources. Bound (b) does not cover those.
- *(c) No stacking.* In a chain P_1 ⇒ … ⇒ P_k with one back arc into each part from the next, vary the start
  statistic of one part with everything else fixed. The count changes by E of the part before it times the product of
  the counts of all parts outside that junction. So a unit step at a junction forces every other part to be
  transitive, and two independent unit-step windows cannot coexist.
- *(d) Counting corollary.* One-apex sources with E_A = 1 (the word ends in Z) and one fixed-(H_B, ρ_B) bank cover at
  most a fraction 1/2 + 1/H_B of the odd numbers. By (a), the base H_B·g(U) + ρ_B·(1 + U) exceeds ρ_B·U ≥ H_B·U, so
  at most X/(2H_B) sources fall below X, and each carries at most ⌊H_B/2⌋ + 1 values.

So windows can only add a bounded number of counts on top of a base that is already dense. They cannot turn a
one-parameter base into a covering. The two ways around (b)–(c) are the other session's tilted readers and interfaces
with several back arcs (seam ranks 13, 62, … in THM-4878).

**Products do not cover the gaps of the root-change families (FINITE-EXACT).** Multiplicativity realizes every
composite that has a factorization into realized counts, so Chan–Pak reduces to primes and numbers like 7p. Those are
exactly where the root-change families are weakest: local factors favour the residue class 0 modulo small primes
(THM-5128 (A)), so composites are over-covered and primes under-covered.

| family | primes among misses | primes among odd numbers | window |
|---|---:|---:|---|
| one root change (F_≤1) | 18.8% | 11.3% | [2^25, 2^26) |
| runs ≤ 2 | 29.0% | 12.8% | [2^22, 2^23) |
| runs ≤ 3, Z-free blocks | 49.8% | 12.8% | [2^22, 2^23) |
| invertible blocks only (8) | 17.0% | 13.4% | [2^21, 2^22) |

([`scripts/missclass.c`](scripts/missclass.c), [`scripts/classwhite.c`](scripts/classwhite.c),
[`scripts/classterm.c`](scripts/classterm.c).) The closer a family is to complete, the more its last misses are primes.

**Where the 11-matrix orbit is weakest ([`scripts/rep11.c`](scripts/rep11.c), [`scripts/primemiss.c`](scripts/primemiss.c)).**
Count representations of each odd N in [X/2, X] (strict block words of the 11-block set):

| X | mean | primes: mean | multiples of 3: mean | class 1 or 2 mod 3: mean | minimum (where) |
|---|---:|---:|---:|---:|---|
| 2^18 | 321 | 76 | 660 | 150 | 0 (misses remain below 418,429) |
| 2^20 | 545 | 109 | 1,159 | 238 | 2 (at a prime) |
| 2^22 | 925 | 152 | 2,026 | 375 | 5 (at the prime 2,115,427) |
| 2^24 | 1,571 | 216 | 3,529 | 592 | 18 (at 8,431,181 = 11 × 766,471; the weakest prime has 23) |

The class-0 boost is about 6.0 modulo 3, 3.7 modulo 5 and 1.6 modulo 7 at 2^24. The mean grows like X^0.38 and the
prime mean like X^0.25. The minimum sits at a rough number coprime to 3, 5 and 7: a prime up to 2^22, and 11 × 766,471
at 2^24. With products added, the orbit's misses are its 208 prime misses (the
largest 386,641) and 36 numbers whose only factorizations use 7 or 21. Its 164 other misses are products of two
counts. So products lower the orbit's threshold only from 418,429 to 386,641.

**Approximate placement (FINITE-EXACT, [`scripts/gaps.c`](scripts/gaps.c), [`scripts/gmult.c`](scripts/gmult.c)).**
The one-apex values alone hit 0.498 of the odd numbers, with maximal gaps 30, 42, 46, 52, 50 on the windows starting
at 2^16, 2^19, 2^23, 2^24, 2^25 (about 2 log₂ N). There are 2 ln 2 = 1.386 words per odd number, the largest
multiplicity grows slowly (6 near 2^8, 18 near 2^24), and the Cauchy–Schwarz density bound settles at 0.406. Every odd
target is therefore within O(log N) of a one-apex count. Exact placement from there is blocked by (b)–(c). A
positive-density theorem from there is the repository's open linear-energy statement (THM-4886, THM-4887).

**No local fine corrections (PROVED for added vertices, FINITE-EXACT for flips).** Approximate placement would become
exact if some local move could shift a one-apex count by small exact amounts. None of the tested moves can:
- *adding vertices:* one new vertex either keeps h or raises it to at least (3h+1)/2 (HYP-9469 (1), proved);
- *attaching readers:* obstructions (a)–(d) above;
- *flipping one spine arc of a one-apex tournament* ([`scripts/flipstat.c`](scripts/flipstat.c),
  [`scripts/flips.c`](scripts/flips.c)): the smallest nonzero change averages 20% of H at spine length 10 and 18% at
  12. Flips changing H by at most 20 are 5.8% of all flips at length 10 and 1.4% at length 12. Over 40 random words,
  no word realizes every even change in [−20, 20].

So universal placement cannot be local: it has to come from the global variety of a dense family.

**A singular-series check for the 11-matrix orbit ([`scripts/singser.c`](scripts/singser.c)).** On [2^21, 2^22], take
for each small prime p the ratio between the mean multiplicity over N not divisible by p and the overall mean:

| p | 3 | 5 | 7 | 11 | 13 | 17 | 19 | 23 | 29 | 31 | 37 | 41 | 43 |
|---|---|---|---|---|---|---|---|---|---|---|---|---|---|
| ratio | 0.405 | 0.672 | 0.930 | 0.958 | 0.973 | 0.898 | 0.958 | 0.995 | 0.940 | 0.945 | 0.988 | 0.993 | 0.991 |

Their product over p ≤ 43 is 0.175, and the observed ratio for primes is 0.164, so the prime deficit is consistent with
an Euler product of local factors. The factors are not monotone in p: there are strong class-0 boosts at 17, 29 and
31 and almost none at 23, 41 and 43. Numbers coprime to every p ≤ 43 have ratio 0.219, above the product, so the
local factors are positively correlated. Any local limit theorem for this orbit would carry this irregular singular
series.

**The two-dimensional shadow ([`scripts/pairs2d.c`](scripts/pairs2d.c), [`scripts/shadowcap.c`](scripts/shadowcap.c)).**
Could the full root-change grammar (unbounded runs) admit a Euclid-like descent, as continued fractions with unbounded
partial quotients do? That would need its reachable (count, port) pairs to fill a positive share of a cone. They do
not. In the cone X/2 ≤ H ≤ X, H ≤ λ ≤ 4H, the reached share is 3.3%, 1.6%, 0.77%, 0.37% at X = 2^8, 2^10, 2^12,
2^14, so the number of reached pairs grows like X^1.47. For bounded runs, the same pair exponent (fit over 2^10 to
2^16) and the share of odd counts in [2^15, 2^16) that are hit are:

| run cap | 1 | 2 | 3 | 4 | 5 | unbounded |
|---|---|---|---|---|---|---|
| pair exponent | about 0.9 | about 1.15 | about 1.26 | about 1.34 | about 1.39 | 1.47 |
| odd counts hit in [2^15, 2^16) | 13% | 72% | 98.3% | 99.98% | 100% | 100% |

The regimes of wave 20 line up with this "shadow dimension" as in Zaremba-type problems, but the threshold is not
1. Cap 2 already has exponent above 1 and still misses infinitely often, and local factors move the effective
cofiniteness threshold to roughly 1.25. Since the reachable set is a thin, about 1.5-dimensional subset of the plane,
cofiniteness is a fixed-direction projection statement for a thin set. Marstrand-type theorems speak only of almost
every direction, so a proof would need non-concentration along vertical lines.

**What the combinations say about cofiniteness.** All routes end at the same place: a local limit theorem for a
supercritical family on its weakest residue classes, which are the primes. The 11-matrix orbit is supercritical even
there (prime multiplicity growing like X^0.25), which is why it looks cofinite. Two facts keep thin-orbit methods
away: the singular blocks are load-bearing (wave 20), and even with expansion those methods prove density one, not
cofiniteness. A proof needs a new ingredient, either an exact numeration or a prime-specific mechanism. None of the
strategies combined here supplies one.

**A combination that does change something: a route to positive density.** The other session's free-address
window is proved for every u ≥ 1: 908u + 9g(u) + {752, 754, 756, 758} are Hamiltonian-path counts
(LINEAR_WINDOWS.md). So every distinct value of F(u) = 908u + 9g(u) gives a distinct count. If F is boundedly
many-to-one, the Hamiltonian-path spectrum has positive lower density. Positive density is open: the repository's
best proved lower bound is X^(3/4) up to logarithms (THM-5141).
- *Data (FINITE-EXACT, [`scripts/fwin.c`](scripts/fwin.c), [`scripts/finj.c`](scripts/finj.c),
  [`scripts/fcoll.py`](scripts/fcoll.py)).* This one family realizes 0.873% of the odd numbers up to 10^8, 10^9 and
  10^10, which is 8/917 to three digits. For u ≤ 2·10^7, F takes 19,977,121 distinct values: 0.23% of arguments sit
  in repeated values, and the largest multiplicity is 3. Colliding pairs share a long high prefix and differ only in
  their low bits, for example u = 3232 = 110010100000 and u′ = 3259 = 110010111011.
- *Structure (PROVED, elementary; checked on 5,000 random cases, [`scripts/phirec.py`](scripts/phirec.py)).* Write
  Ψ for the binary pair energy, so g(u) = 1 + u − u₀ − Ψ(u). Then F(u) = 917u − 9Ψ(u) − 9u₀ + 9. Splitting
  u = 2^m V + W gives F(u) = 917·2^m V − 9Ψ(V) + 9 + Φ_V(W), where
  Φ_V(W) = 917W − 9Ψ(W) − 9V·rev_m(W) − 9·v₀·w_(m−1) − 9·w₀.
  Peeling the top bit b of W gives Φ_V(b·2^(m−1) + W) = Φ_(2V+b)(W) + b·(917·2^(m−1) − 9V − 9v₀). Hence a
  collision F(u) = F(u′) whose highest differing bit is m − 1 is a "cross" coincidence
  Φ_(2V)(W) = Φ_(2V+1)(W′) + (917·2^(m−1) − 9V − 9v₀) between (m−1)-bit words.
- *Two partial lemmas (PROVED).*
  - For fixed (W, W′), the cross equation is linear in V with coefficient 2(rev W − rev W′) − 1, which is odd and
    nonzero. So at most two values of V (one per parity) give a cross collision.
  - Collisions also need 9 | (u′ − u). So at a level with 4^m ≤ U there are at most 2·4^(m−1)/9 + 2^m colliding
    pairs, and all such levels together contribute at most about U/13.
- *The gap.* At the high levels (4^m > U, high part V smaller than 2^m) only crude bounds hold. Injectivity of Φ_V'
  on the bottom half of the bits gives multiplicity at most about 2^(m/2)·(103/V)^(1/2), hence about U^(3/2) pairs.
  The data (Φ_V at most 2-to-one for all sampled V ≤ 103·2^m and m ≤ 14, [`scripts/phimult.c`](scripts/phimult.c))
  point to O(U log U).
- *What would follow.* Collision pairs O(U^(5/4 − ε)) would already beat X^(3/4). O(U log U) gives at least c·X/log X
  counts below X. O(U) gives positive lower density. This is the one place in this wave where two strategies combine
  into something neither has alone: the other session's exactness and this wave's distinctness analysis.

**Pre-registered.**
- **P41 (registered in the session draft before the 2^24 run finished; two of three parts CONFIRMED).** For the
  11-block orbit on [2^23, 2^24]:
  - mean prime multiplicity at least 200: yes, 216;
  - minimum multiplicity at least 7: yes, 18;
  - minimum attained at a prime: no. It is attained at 11 × 766,471, a rough composite; the weakest prime has 23.
- **P42 (REFUTED).** The maximal gap between one-apex values on [2^26, 2^27) was predicted to lie between 48 and 58.
  It is 44. The maximal gaps on the windows starting at 2^23, 2^24, 2^25, 2^26 are 46, 52, 50, 44: they fluctuate
  and show no clear growth, so "about 2 log₂ N" overstates the trend. Bounded gaps are not excluded by the data.

## Wave 22: collisions of 908u + 9g(u) — Poisson statistics, a spacer factorization, and the 3/4 barrier

This wave continues wave 21's route to positive density. **Chan–Pak OPEN; positive density OPEN.** The values of
F(u) = 908u + 9g(u) collide as a Poisson process would, at every scale tested: up to 2^30, 0.11% of the arguments
share a value and no value is hit four times. A proof is blocked by the wide-window case that also stops THM-5141 at
exponent 3/4.

**Verdicts.**
- **P37 CONFIRMED.** The sub-language "one-apex source, then readers J O^i Z^j with i, j ≤ 3" attains every odd N in
  [2649, 10^7]. Its 30 misses are all at most 2,647 ([`scripts/readeralpha.c`](scripts/readeralpha.c), 2 h 22 min).
- **P38 CONFIRMED.** F_≤2 misses no odd N in (3,780,967, 4·10^9]. It has 140 misses below 4·10^9
  ([`scripts/f2census.c`](scripts/f2census.c), 1 h 49 min of CPU).
- **P43 and P44 CONFIRMED** (below).

**F is always odd, and its collisions look Poisson (FINITE-EXACT).**
- g(2v) = g(v) + v + v₀, and v + v₀ is even. Also g(2v+1) = g(v). So g is always odd, and so is F(u).
- Take φ(u), the binary radical inverse, to be uniform. Then the values of F near X have density
  ρ = ln(917/908)/9 = 1/912.5 per integer, or 2ρ per odd integer. In that model the number of pairs u < u′ ≤ U with
  |F(u) − F(u′)| ≤ L (L even) is U(L + 1)/912.5.
- Measured pairs divided by U(L + 1) ([`scripts/fnear.c`](scripts/fnear.c)); the model value is 0.001096:

  | L | U = 2^18 | 2^20 | 2^22 | 2^24 |
  |---:|---:|---:|---:|---:|
  | 0 | 0.00123 | 0.00110 | 0.00112 | 0.00115 |
  | 2 | 0.00114 | 0.00108 | 0.00113 | 0.00112 |
  | 8 | 0.00106 | 0.00104 | 0.00107 | 0.00109 |
  | 98 | 0.00109 | 0.00108 | 0.00109 | 0.00109 |
  | 9,998 | 0.00108 | 0.00109 | 0.00109 | 0.00109 |

- Exact collisions at larger U ([`scripts/fmult9.c`](scripts/fmult9.c)). A collision forces 9 | u′ − u, so each
  residue class of u modulo 9 is sorted separately.

  | U | colliding pairs / U | values hit twice | values hit three times | Poisson model for triples |
  |---|---:|---:|---:|---:|
  | 2^26 | 0.001122 | 75,150 | 45 | 54 |
  | 2^28 | 0.001120 | 299,925 | 225 | 215 |
  | 2^30 | 0.001104 | 1,182,425 | 910 | 860 |

  No value is hit four times up to 2^30 (model: 0.47 expected), and 99.78% of the arguments u ≤ 2^30 have a value
  of their own. The largest departure from the model in the first table is 12% (L = 0 at U = 2^18); for L ≥ 8
  every entry is within about 5%.

**Exact structure of a colliding pair (PROVED, elementary).**
- *The collision equation is THM-5141 (C) transported to F.* Split a colliding pair u < u′ into:
  - the common top V, above the highest differing bit;
  - the common bottom L, the j bits below the lowest differing bit;
  - the windows M, M′ of n₁ bits between them.

  Put A = M′ − M, B = rev(M′) − rev(M) and ρ = rev_j(L). Then 9 | A and

  Aρ + BV = 917·2^j·(A/9) − (Ψ(M′) − Ψ(M)) − v₀ − σε,

  where σ = m′₀ − m₀, ε = l_(j−1) for j ≥ 1 and ε = 1 for j = 0. This is a lattice line in the outer data (ρ, V).
  The weight 917·2^j − 9ρ plays the role of THM-5141's suffix weight c, and V that of its prefix P. Checked on all
  19,275 colliding pairs below 2^24 ([`scripts/fcoll_windows.py`](scripts/fcoll_windows.py), with the pairs listed by
  [`scripts/fcollide.c`](scripts/fcollide.c)).
- *Monotone residue classes.* For u = 2^k a + r with r < 2^k and a′ > a,

  F(2^k a′ + r) − F(2^k a + r) = (a′ − a)(917·2^k − 9 rev_k r) − 9(Ψ(a′) − Ψ(a)) − 9(a′₀ − a₀)r_(k−1) ≥ 908·2^k − 9a′,

  because Ψ(a′) ≤ a′. So F is strictly increasing on every residue class modulo 2^k inside [0, 100·4^k]. A
  colliding pair therefore satisfies 100·4^j < u′, where j is its lowest differing bit: j is below half the length
  of u′. Data below 2^24: j ≤ 6 (the bound allows 8). The counts for j = 0, 1, 2, 3 are 9,593, 4,885, 2,380 and
  1,325, halving each time, and most pairs have V/2^j between 2^4 and 2^8.
- *The spacer factorization.* Let q be the number of bits just below the highest differing bit where u has 1 and u′
  has 0. Then u = 2^m V + (2^q − 1)2^p + c and u′ = 2^m V + 2^(q+p) + c′, with m = q + p + 1 and p-bit words c, c′
  whose top bits are not (1, 0). With D = rev_p(c′) − rev_p(c) − 1 and 9t = u′ − u = 2^p + c′ − c,

  (2V + 1)·(2^(q+1)·D + 3) = 1834t + 2(Ψ(c) − Ψ(c′)) − 2 rev_p(c) + 2(c₀ − c′₀) − 2v₀ + κ,

  where κ = 2c_(p−1) − 1 if q ≥ 1 and κ = 1 − 2c′_(p−1) if q = 0. Checked on all 19,275 pairs
  ([`scripts/spacer_identity.py`](scripts/spacer_identity.py)).
  - *D ≠ 0.* If D = 0, then c = 0X and c′ = 1X, so u′ − u = 3·2^(p−1). That is not a multiple of 9, a contradiction.
    Without the factor 9, an inner pair with D = 0 and one admissible V would collide for every spacer length q.
    The core count h = 9 rules out these infinite families.
  - For q ≥ 1 the right side does not depend on q. So, for a fixed inner pair, each spacer length that occurs makes
    2^(q+1)D + 3 a divisor of one of two fixed numbers.
  - Data: q = 0, 1, 2, 3 for 9,934, 4,754, 2,240 and 1,163 pairs, halving like j.

**Why the proof stops where THM-5141 does.**
- Collisions whose highest differing bit lies below about n/2 are few. There are at most about 2^n pairs of low
  parts there, and each allows at most two tops V (wave 21).
- The collision equation forces balance. When the windows have no runs longer than R, V/2^j lies within a factor
  2^(R + O(1)) of 102, as THM-5141 (C) balances its prefix against its suffix. So pairs that share their first n/4
  bits have windows of at most n/2 + R + O(1) bits. That is how THM-5141 gets 2^(3n/4) distinct values from one
  prefix family, and the same mechanism should give the same exponent for F (not written out here).
- Positive density needs the wide windows: collisions between different prefix families. There the window pairs
  outnumber the arguments, so almost every pair must be shown to have no lattice point in its box. In the data, all
  but a Poisson fraction have none.
- The missing input is equidistribution: the reversal-graph points (M, rev M) must meet these lattice lines at the
  random rate for all but a sparse set of outer data. That input has to use the particular slopes. The bit-reversal
  net puts 2^(k/2) points, the binary palindromes, on the slope-one diagonal, so no bound uniform over all lines can
  hold.

**Conjecture (stated here; no HYP number).** #{u < u′ ≤ U : F(u) = F(u′)} ≤ 0.0012·U for every U ≥ 2^18. With the
proved window this gives at least (1 − 0.0012)·X/917 − O(1) Hamiltonian-path counts up to X. That is positive lower
density, at least 2/917 of the odd numbers (0.87% when all four window values are counted, as observed). A bound
O(U^(5/4 − ε)) would already beat THM-5141.

**Pre-registered.**
- **P43 (registered in the session draft at 12:18 MST before the run; CONFIRMED).** At U = 2^28:
  - colliding pairs / U in [0.00105, 0.00120]: 0.001120;
  - between 140 and 260 values hit three times: 225;
  - no value hit four times: none.
- **P44 (registered at 12:20 MST before the program was written; CONFIRMED).** At U = 2^30:
  - colliding pairs / U in [0.00108, 0.00117]: 0.001104;
  - between 700 and 1,050 values hit three times: 910;
  - at most two values hit four times: none.
- **P45.** At U = 2^32 (one residue class modulo 9 at a time, split further by value range):
  - colliding pairs / U lies in [0.00107, 0.00113];
  - at least one value is hit four times (model: 1.9 expected);
  - no value is hit five times.

## Wave 23: thresholds, rank and reciprocity — a proof by contradiction, mapped

The owner asked for an attempt at a proof, possibly by contradiction, using analogies, in particular the Jacobian
conjecture: counterexamples exist from dimension three on, but none is known in the plane. Why do such thresholds
occur, and how do they relate to the symmetry structure of numbers? **Chan–Pak OPEN.** The result of the wave is that
Chan–Pak sits on the flexible side of every threshold that can be measured:
- its complete families have rank-two weights (2 and 3) or a transfer group that is Zariski-dense in GL_3;
- no reciprocity invariant appears in any test;
- square classes are automatically products.

A proof by contradiction therefore comes down to excluding sporadic exceptions above 2^32, and that needs one missing
theorem: a pointwise local–global principle for a thin orbit of higher rank.

**Prior art first.** The following are already in the repository:
- the atlas of dimension thresholds (HYP-9343, THM-4727);
- the thin-orbit analogy and the Apollonian counter-indication (HYP-9661 (D), the S82 note);
- the Heegner and roughness profile of F1 (HYP-9657);
- from session 142 ([minimal gap](../session142_tournament_odd_spectrum_20261102/minimal_gap/PROOF.md),
  [Jacobian bridge](../session142_tournament_odd_spectrum_20261102/jacobian_bridge/PROOFS.md)): the atoms of a
  minimal counterexample, the asymmetry of prime-count tournaments, and the Jacobian bridge;
- the tribonacci polynomial of the corridors (wave 16);
- the frame identity H = Σ 2^α 3^β (waves 6, 8);
- the density 1/2 of the one-hub value set (HYP-9469 (3)).

Sections 2–5 are new.

### 1. Why thresholds happen (ANALOGY; every example is a cited theorem)

Three mechanisms recur.
- **Room.** Below the threshold, the relevant boundary or symmetry is one-dimensional: a tree, an abelian group, a
  single invariant. A counterexample needs a second independent direction to hide in.
  - *Jacobian.* Every known Keller collision is driven by a C*-orbit over a two-dimensional quotient. In the plane,
    every C*-quotient is a line, where integrality kills the twisted equation (THM-4727, THM-1345).
  - *Automorphisms.* Aut(A²) acts on a tree (Jung–van der Kulk). Tame(A³) acts on a 2-complex, and Nagata's map is
    wild (Shestakov–Umirbaev).
  - *Banach–Tarski.* The isometries of the plane form an amenable group; SO(3) contains free subgroups.
  - *Poincaré–Bendixson.* Orbits of a planar flow separate the plane.
  - *Dehn.* Area is the only scissors invariant in the plane. In space a second one appears: length ⊗ angle in
    R ⊗ R/πQ.
- **Rank (the reverse direction).** When the symmetry group itself gains a second independent direction, objects
  become rigid.
  - *Furstenberg (1967).* ×2 alone has many closed invariant subsets of the circle, for example the binary
    expansions avoiding 11. Under ×2 and ×3 together, the only ones are finite sets and the whole circle. The
    integer form is Cobham (1969): a set recognizable in base 2 and in base 3 is ultimately periodic.
  - *Moduli.* Hyperbolic surfaces have moduli; hyperbolic 3-manifolds do not (Mostow).
  - *Conformal maps.* The Riemann mapping theorem holds in the plane; in space only Möbius maps are conformal
    (Liouville).
  - *Congruence subgroups.* SL_2(Z) has non-congruence subgroups; SL_n(Z) for n ≥ 3 has none
    (Bass–Milnor–Serre).
- **Reciprocity.** Where a local–global principle for representing integers fails along an infinite family that no
  congruence detects, the cause in every case below is a quadratic-character invariant: every prime separately
  allows the value, but the product formula forbids it.
  - *Ternary quadratic forms.* Infinitely many exceptions can occur, inside finitely many square classes t·m²
    (spinor exceptions; Duke–Schulze-Pillot). Colliot-Thélène and Xu identify them as a Brauer–Manin obstruction
    for integral points. With five or more variables, every large locally represented integer is represented
    (Tartakowsky).
  - *Apollonian packings.* Squares are missing (Haag–Kertzer–Rickards–Stange).
  - *Thin semigroups of SL_2(Z).* Squares are missing (Rickards–Stange).
  - The common feature: all three settings have real rank one (SO(3), SO(3,1), SL_2). There the congruence subgroup
    property fails, and the spinor or Kubota symbol is not a congruence character.
- **Sporadic exceptions.** These occur where there is not yet room:
  - the critical numbers of the 15 and 290 theorems;
  - Ramanujan's form x² + y² + 10z², which misses 3, 7, 21, 31, … (the list is complete under GRH,
    Ono–Soundararajan);
  - exceptional isomorphisms in small rank.

**Hypothesis (where thresholds come from).** The threshold is the dimension at which the governing symmetry stops
being commutative or one-dimensional. Such a symmetry carries one extra invariant: an invariant mean, a length function
on a tree, a winding number, or a character such as the spinor norm or the Kubota symbol. On one side that invariant
protects: no Keller fold, no paradoxical decomposition, no planar chaos. On the other side it obstructs: spinor
classes, missing squares, half-dense numerations. The invariant is usually a symbol from number theory:
- Frobenius in characteristic p, where the Jacobian conjecture already fails in dimension one (x + x^p has derivative
  1 but is not injective);
- whether angles are rational multiples of π (Dehn);
- quadratic reciprocity and the Hilbert product formula (spinor norm, Kubota symbol);
- the multiplicative independence of 2 and 3 (Furstenberg, Cobham, Baker).

**Sub-hypothesis (rank-one reciprocity).** Infinite reciprocity families occur only when the governing group has real
rank one. In higher rank, the congruence subgroup property makes every such symbol a congruence character, so only
congruence and sporadic exceptions remain.

### 2. Chan–Pak on the room and rank side: the rank ladder (PROVED algebra; KNOWN coverage)

For each family computed by a finite linear transfer (one matrix per letter), take:
- the associative algebra spanned by the letter words;
- its semisimple quotient, via the rank of the trace form tr(ab), computed exactly modulo two 61-bit primes;
- the letters' eigenvalues.

The scripts are [`scripts/rank_ladder.py`](scripts/rank_ladder.py) and
[`scripts/eigen_groups.py`](scripts/eigen_groups.py). Each transfer is checked against brute force first, and the
two-hub port of [`scripts/twohub.c`](scripts/twohub.c) against a subset DP on 300 random words.

| family | transfer | algebra | semisimple part | nonzero eigenvalues of letters | coverage (KNOWN) |
|---|---|---|---|---|---|
| one apex (O, Z) | 4-dim (H,E,L,h) | dim 8 | commutative, dim 3 | 1, 2 | 0.498 of the odd numbers (wave 21) |
| one hub | 6 states | dim 8 | commutative, dim 3 | 1, 2 | 1/2, numerically (HYP-9469 (3)) |
| two hubs | 36 states | dim 48 | commutative, dim 4 | 1, 2, 3 | every odd N ≤ 2.26·10^8 except 7, 21 (wave 17) |
| root changes (O, Z, J) | 4-dim | all of M_4(Q) | simple | O, Z: 1, 2; J: roots of x³ − x² − x − 1 | every odd N in (2647, 10^9] (HYP-9494 (D)) |

The two complete families are complete for different reasons:
- The two-hub transfer is triangularizable. Its counts are sums of products of the weights 1, 2, 3, which is the frame
  identity H = Σ 2^α 3^β of waves 6 and 8.
- The root-change transfer is not triangularizable at all.

The two half-dense families have weights generating only ⟨2⟩.

**P46 (registered before the computation): PARTLY REFUTED.**
- (a) One apex and one hub are triangularizable: confirmed.
- (b) Root changes and two hubs are not triangularizable: confirmed for root changes, **refuted for two hubs**.

So solvability does not separate complete families from incomplete ones; the rank of the weight group does. Rank one
(⟨2⟩) leaves room for a proper invariant set, like a ×2-invariant Cantor set. Rank two (⟨2, 3⟩, as in Furstenberg's
×2×3) leaves none, and neither does a non-commutative transfer (GL_3). The growth exponents agree:
- rank-one families have δ = 1 (2^m words, values up to about 2^m);
- two hubs have δ = log 4 / log 3 = 1.26.

**Theorem (PROVED; new).** Let Γ be the group generated by the invertible root-change blocks J, JOZ, JOOZ, JOZZ,
JOOOZ, JOZOZ, JOZZZ acting on (H, E, λ) (determinants 1, 4, 8, 8, 16, 16, 16;
[`scripts/blocks11.py`](scripts/blocks11.py)). The Zariski closure of Γ is GL_3.

*Proof.*
1. Let G be the Zariski closure; it is defined over Q.
2. The characteristic polynomial of J is x³ − x² − x − 1, the tribonacci polynomial that wave 16 met in the
   corridors. It is irreducible with discriminant −44, so its Galois group is S_3. J is regular semisimple, with
   eigenvalues τ_1 ≈ 1.839 and a complex pair of modulus τ_1^(−1/2).
3. The only multiplicative relation among the eigenvalues is τ_1τ_2τ_3 = 1. Suppose a relation had unequal exponents.
   Dividing it by its image under a transposition in the Galois group gives (τ_i/τ_j)^k = 1 for some i ≠ j and k ≠ 0.
   Applying a further Galois element turns this into a ratio of modulus different from 1, which is impossible.
4. So the Zariski closure of ⟨J⟩ is the norm-one torus T_1, a maximal torus of SL_3, and G° ⊇ T_1.
5. The Lie algebra of G° is its Cartan part plus the root spaces of T_1 for some set Φ′ of roots. Φ′ is stable under
   the Galois group, because G is defined over Q and the Galois group permutes the eigenlines of J. S_3 acts
   transitively on the six roots, so Φ′ is either empty or all of them.
6. If Φ′ is all of them, then G° ⊇ SL_3, and det(JOZ) = 4 makes G = GL_3.
7. If Φ′ is empty, G° is a torus between T_1 and the centralizer T of J. Then G normalizes T. Since the Weyl group has
   order 6, g^6 ∈ T for every g ∈ G, so g^6 commutes with J. But B^6 J ≠ J B^6 for each of the six other blocks B
   (exact check in [`scripts/rank_ladder.py`](scripts/rank_ladder.py)). ∎

**FINITE-EXACT: strong approximation, visible.** The seven blocks generate all of GL_3(F_3) (11,232 elements) and
all of GL_3(F_5) (1,488,000 elements) ([`scripts/modp_image.py`](scripts/modp_image.py)). The commutator subgroup
lies in SL_3(Z[1/2]) and is Zariski-dense in SL_3. By strong approximation (Matthews–Vaserstein–Weisfeiler, CITED),
it surjects onto SL_3(F_p) for all but finitely many odd primes p.

So the root-change orbit is a thin orbit of higher rank (real rank two), not of rank one like the Zaremba and
Apollonian orbits. By the sub-hypothesis of section 1 it can carry congruence and sporadic exceptions, but no
reciprocity families.

### 3. Chan–Pak on the reciprocity side: three tests

**P47 (registered before the run): CONFIRMED.** The root-change orbit has no Jacobi-symbol invariant.
- Data: the 262,240 states reached by O/Z/J words of length at most 14, and every pair of primitive linear forms
  a·v, b·v with coefficients in {−1, 0, 1}.
- There are 827 pairs with at least 2,000 samples once the degenerate values |a·v| ≤ 2 and b·v = 1 are excluded. For
  each of them, (a·v / b·v) = +1 with frequency between 0.459 and 0.591.
- Including those values, 25 of 830 pairs reach 0.60–0.705, all because of them
  ([`scripts/jacobi_inv.py`](scripts/jacobi_inv.py)).
- On the triangular one-apex orbit, the only exact "invariants" are constant coordinates (h = 1, E ∈ {0, 1}), which
  are congruences.
- The one-apex pair (g(h), 2h + 1) carries none: (g(h)/(2h+1)) = +1 for 50.08% of h < 2^20, and within 0.7
  percentage points of 1/2 in every class of h mod 8.

**P48 (registered): CONFIRMED for squares, with one class flagged.** This is an Apollonian-type test on F1, the
one-arc joins of two one-apex tournaments, in the coordinates of HYP-9657. The census of r(N) runs to 2^27
([`scripts/f1sq.c`](scripts/f1sq.c)).
- Odd squares m² in [2^22, 2^27): 22 misses among 4,769.
  - That is 0.86 times the rate of controls matched modulo 8·3·5·7·11·13.
  - It is 0.98 times the expectation of an empirical model whose strata are N mod 105, N mod 16, and the numbers of
    primes 11 ≤ p ≤ 97 dividing N and dividing N + 1 ([`scripts/f1strata105.c`](scripts/f1strata105.c)).
- Square classes t·m² for t = 5, 11, 13, 15, 19, 21, 35: within noise.
- But 3m² and 7m² (with t ∤ m) had 10 and 17 misses against 3.8 and 5.6 expected (z = +3.1 and +4.8).

**P49 (registered after that run, before any run on the disjoint window): REFUTED.** On [2^27, 2^29):
- 3m²: 2 misses against 2.1 expected;
- 7m²: 4 against 3.4;
- squares: 12 against 14.5.

So the in-sample excess was noise or a small-N effect. No square-class obstruction appears in F1, although F1 is the
critical rank-one family. A tempting mechanism existed: 7m² with 7 ∤ m is never of the form x² + 2^s y², which is
exactly the form of F1's self-joins g(h)² + 2^s(2h + 1)². But self-joins number about √X among the X log X
representations, too few to matter.

**F1's weakness on primes is roughness, not primality (FINITE-EXACT, [`scripts/f1prime.c`](scripts/f1prime.c)).**
- Take the odd N in [2^22, 2^27) with no prime factor up to 97. Primes and composites among them miss at the same rate:
  matched by N mod 16 and by the number of primes up to 97 dividing N + 1, observed/expected is 1.03.
- The miss rate falls with that number. For 0, 1, 2, 3 and at least 4 such primes it is 1.27%, 0.70%, 0.42%, 0.28%
  and 0.17%, so each small prime dividing N + 1 roughly halves it.
- This makes the Rabinowitsch stratum of HYP-9657 (B) quantitative.
- A model of independent class factors over the primes up to 97 underpredicts the misses of primes by a factor 1.83
  ([`scripts/f1local.c`](scripts/f1local.c)), so the local factors interact strongly.

For the whole spectrum, the Apollonian mechanism is unavailable for a structural reason (PROVED, elementary). H is
multiplicative over strong components, so t·m² is a count whenever t and m are. Every odd m < 2^32 is a count or one
of 7, 21, and 49 and 441 = 9·49 are counts. So every odd square below 2^64 is a Hamiltonian-path count, and a
square-class obstruction would need the roots m to be exceptions themselves.

### 4. The proof by contradiction, mapped

Suppose the exceptional set E is larger than {7, 21}, and let N* be its least element above 21. KNOWN: N* > 2^32,
and N* is a prime or 7p (THM-4094; minimal gap). Known mechanisms for exceptions come in three types.
1. **Congruence**: a residue class where the local density vanishes. No such class exists. F1 alone has positive
   local density in every class modulo each odd prime and each squarefree odd modulus, and in every odd class modulo
   every power of 2 (THM-5128).
2. **Reciprocity**: a square class, or a set of primes cut out by a splitting condition, through a quadratic-character
   invariant. Square classes are products (above). The tests of section 3 find no invariant. By the sub-hypothesis of
   section 1, the higher-rank transfer of section 2 cannot carry one.
3. **Sporadic**: values the families miss by accident. This is the type of 7 and 21.
   - Busch's floor (f(3k+1) = 5^k, f(3k+2) = 9·5^(k−1)) bounds the order of the strong core: at most 4 for 7
     (f(5) = 9 > 7) and at most 6 for 21 (f(7) = 25 > 21).
   - The strong spectra of orders 3 to 6 are {3}, {5}, {9, 11, 13, 15} and {15, 17, 19, 23, 25, 27, 29, 31, 33, 37, 41, 43, 45}, which skips 21
     ([`scripts/strong_small.py`](scripts/strong_small.py)).
   - Beyond 2^32 the least possible order is at least 42 (busch §5). There, about 2^(n²/2)/n! strong tournaments
     face about 0.36·b(n) targets. This is a heuristic, not a proof.

So a counterexample would have to be sporadic: a value missed by every family at a size where the families average
astronomically many representations. Turning that into a contradiction needs exactly one theorem:

> **The missing theorem.** Every sufficiently large odd N with positive local density is a value of the root-change
> orbit, or of some family that is supercritical and Zariski-dense in a group of real rank at least two.

The known theory gives at best density one, and only in rank-one settings (Bourgain–Kontorovich). The Apollonian and
Rickards–Stange counterexamples show that in rank one, density one cannot be upgraded. In higher rank there is
neither a counterexample nor a theorem. That is the precise gap.

**Jacobian parallel (ANALOGY).** In both problems the hard instances are exactly the symmetry-free ones, while
symmetric constructions (torus quotients; products and padding) handle everything else.
- A planar Keller counterexample would have to be torus-rigid: it admits no C*-equivariance (THM-1345,
  THM-4727 (D)).
- A tournament whose count is a prime p > 3 must be asymmetric (minimal gap, KNOWN).

### 5. Numerology of 7 and 21, checked

- 7 = 1 + 2 + 4 and 21 = 1 + 4 + 16 are 1 + q + q² for q = 2, 4: the point counts of the projective planes over F_2
  and F_4.
- 7 = |Z/7| and 21 = |Aut(QR_7)| = |Z/7 ⋊ Gal(F_8/F_2)|, the normalizer of a Singer cycle in GL_3(F_2).
- 7 is a Heegner number, and Ramanujan's form x² + y² + 10z² misses both 7 and 21.

Only the first reading has a mechanism, and it is the KNOWN reading of busch §2:
- 1 + q + q² = (1 + q) + q² is a Fermat cut only through a power-of-two port at count 1 + q.
- The two smallest Fermat blocks have none: Λ(3) = {3·2^t} and Λ(5) ∩ [5, 16] = {5, 6, 7, 9, 10, 12, 14}.
- The next member, 73 = 9 + 64, is attained.

The automorphism and Heegner coincidences constrain no count: a realizer of 7 would have to be asymmetric, and one of
21 would have an automorphism group of order dividing 21. Verdict: numerology.

### Pre-registered

- **P50.** For k hubs over a transitive spine, the transfer is triangularizable with letter eigenvalues
  {0, 1, …, k + 1}. So the weight group has rank π(k + 1): rank 2 for k = 2, 3 and rank 3 for k = 4, 5.

## Wave 24: the invertible root-change orbit — reset-free, Haar-distributed, and complete beyond 6223

Owner prompt (2026-11-14): attempt a density-one theorem for the root-change orbit using its GL₃ structure, and
keep looking for creative routes to Chan–Pak. **Chan–Pak OPEN. Density one OPEN; its local half is now PROVED, and
what remains is the archimedean / minor-arc half of an L² local limit law.**

The short answer of this wave:
- **The object.** The right object is the *invertible* sub-language: the one-apex source followed by invertible
  blocks only. It has no singular block, so no reset, and it attains every odd N > 6223 through 2^26.
- **The local half (proved).** Its local behaviour is settled: the blocks generate {g : det g ∈ ⟨2⟩} modulo every
  odd prime ([THM-5162](../../../01-canon/theorems/THM-5162-strong-approximation-for-the-invertible-root-change-blocks-heisenberg-unipotents-and-the-haar-local-law.md)),
  and the counts follow the Haar law exactly. The census matches this law to three decimals.
- **Primes.** Primes are not disadvantaged: their singular series is (7/8)ζ(3) > 1.
- **What remains.** Density one, and with it a proof by contradiction of Chan–Pak, now needs only the archimedean
  half of an L² local limit law. The data give that half with V(X) ≈ 8/mean(X) → 0.

**Prior art first (KNOWN, credited).**
- THM-5120 (Session 246) proves that rank-one resets force a factorization and that prime-count words avoid them.
- The local densities of the one- and two-root-change families come from p-divisible ("confined") states:
  THM-5128, THM-5134, THM-4916, HYP-9506, HYP-9660. Their uniform-gap conjecture HYP-9660 (B) concerns the rank-one
  digit operators.
- Wave 20 found that invertible blocks alone miss 34% of the odd numbers below 2^21. That census capped the source
  at four letters (`runterm.c`), and that cap turns out to be the whole effect.
- Wave 23 proved Zariski density in GL₃ and computed the images mod 3 and mod 5.

The shape formula, the Heisenberg unipotents, the surjection at every odd prime, the Haar law and the
free-source census are new here. No strong-approximation statement was in canon before THM-5162.

### 1. The invertible language attains every odd N > 6223 (FINITE-EXACT)

Letters act on (H, E, L, h) as usual, starting from (1, 1, 1, 1). Once the first J has occurred, two rules apply:
- no Z directly after a J;
- no J directly after an O.

So every block J w has w empty or w = O v Z, and these are exactly the invertible blocks (§2). The census
[`scripts/invorb.c`](scripts/invorb.c) counts the words ending in J or Z for each count.

| language (mode) | misses below 2^20 | largest | misses checked to | per-window misses from |
|---|---:|---:|---|---|
| one-apex source + invertible blocks, any length (1) | 43 | 6223 | **2^26: still 43** | 0 from [2^13, 2^14) on |
| no source: orbit of (1,1,1,1) (2) | 288 | 179,367 | **2^26: still 288** | 0 from [2^18, 2^19) on |
| mode 1, blocks with \|w\| ≤ 6 | 43 | 6223 | 2^20 | — |
| mode 1, \|w\| ≤ 5 (16 blocks) | 48 | 12,811 | **2^26: still 48** (P62) | — |
| mode 1, \|w\| ≤ 4 (8 blocks) | 71 | 968,947 | 2^20 | — |
| mode 1, \|w\| ≤ 3 | 721 | 1,048,507 | 2^20 | — |
| free source, J and J O^k Z only (one Z per block) | 350 | 1,041,965 | 2^26: 3823, the largest 67,076,437 (P63) | 0.01% of [2^25, 2^26) |
| free source, J and blocks with exactly two Z's | 355 | 1,028,217 | 2^20 | 0.03% of [2^19, 2^20) |
| free source, J and the 2^(n−2) blocks of one length n = 4, 5, 6, 7, 8 | 2324, 967, 494, 310, 257 | above 10^6 | 2^20 | positive in every window |

The 43 misses of mode 1 are 7, 13, 21, 41, 43, 55, 67, 77, 83, 163, 173, 231, 233, 247, 251, 277, 281, 299, 313,
335, 421, 449, 503, 533, 551, 569, 583, 637, 667, 821, 1067, 1339, 1387, 1661, 1811, 1861, 1877, 1957, 2647, 3523,
3551, 5231 and 6223. The full grammar's 29 are among them.

Multiplicities:

| window | mode 1 mean | primes | minimum (at) | mode 2 mean | primes | minimum (at) |
|---|---:|---:|---|---:|---:|---|
| [2^17, 2^18) | 116.5 | 107.7 | 18 (135,049) | 59.2 | 62.6 | 0 (169,935) |
| [2^20, 2^21) | 207.8 | 191.6 | 29 (1,142,249) | 105.1 | 110.8 | 6 (1,142,249) |
| [2^23, 2^24) | 369.1 | 339.8 | 99 (9,193,369) | 186.0 | 195.8 | 30 (8,560,341) |
| [2^25, 2^26) | 540.9 | 498.1 | 196 (33,727,757) | 272.1 | 286.4 | 68 (35,979,543) |

Fitted growth exponents of the mean multiplicity over [2^16, 2^24):

| | all N | primes |
|---|---:|---:|
| mode 1 | 0.2777 | 0.2771 |
| mode 2 | 0.2753 | 0.2748 |

The 11-block orbit of wave 21, which has resets, gives 0.38 for all N and 0.25 for primes. So the language is
supercritical with critical exponent δ ≈ 1.28 (#words with count ≤ X ≈ X^δ), and primes keep the full exponent.

**Finitely many generators, no source: a critical-exponent ladder.** Take the orbit of the single vector
(1, 1, 1, 1) under the blocks with |w| ≤ c (mode 2 with a cap,
[`scripts/m2cap4_24.out`](scripts/m2cap4_24.out) and its neighbours). This is exactly the Bourgain–Kontorovich
setting: a finitely generated semigroup, Zariski dense in GL₃, acting on one vector. Share of odd numbers missed
per window, with the mean multiplicity in brackets:

| cap | generators | [2^18, 2^19) | [2^20, 2^21) | [2^23, 2^24) | regime |
|---|---:|---|---|---|---|
| \|w\| ≤ 4 | 8 | 48.1% (3.57) | 49.4% (3.40) | 49.5% (3.30) | critical, about half missed |
| \|w\| ≤ 5 | 16 | 15.8% (9.56) | 12.8% (10.64) | 9.0% (12.23) | slowly supercritical, multiplicity ∝ X^0.07 |
| \|w\| ≤ 6 | 32 | 3.6% (17.71) | 1.8% (21.45) | 0.51% (28.04) | supercritical, multiplicity ∝ X^0.13 |
| unbounded | ∞ | 0 | 0 | 0 | multiplicity ∝ X^0.275 |

It mirrors Zaremba's ladder in the bound A of the partial quotients. With the free one-apex source, which is a
critical family on its own, already |w| ≤ 5 is complete beyond 12,811 (P62).

Mixing scales matters. A single block length leaves hundreds of misses up to the top of the range
([`scripts/invlen.c`](scripts/invlen.c)), while lengths up to 5 together already look cofinite.

**Correction to wave 20.** "The singular blocks are load-bearing" holds only with a capped source. With a free
one-apex source, the invertible blocks alone are complete beyond 6223, and blocks with |w| ≤ 5 are complete beyond
12,811 (to 2^20).

### 2. The shape of a block, and its Heisenberg skeleton (PROVED, THM-5162 (A), (B))

- **Shape.** Let w = O v Z, n = |w|, and let S be its set of Z-positions. Put P = Σ_{s∈S} 2^(s−1),
  R = Σ_{s∈S} 2^(n−s) (the n-bit reversal of P), and let Q be the pair energy of S. Then J w acts on (H, E, λ) by

      M_w = [[b − Q, −1, b], [1, 0, 0], [2^n − R, 0, 2^n]],   b = 1 + P,   det = 2^n.

  As v varies, b runs over the odd numbers of (2^(n−1), 2^n), once each.
- **Invertibility.** J w is invertible iff w is empty or w = O v Z. A final O kills the E-row, and J Z ignores E.
  Checked exactly for all 2047 blocks with |w| ≤ 12 ([`scripts/shape_lemma.py`](scripts/shape_lemma.py)).
- **Reading.** The last block reads N = b(H + λ) − E − QH, and E′ = H. Counts move through a shift register
  (H, E) = (current, previous), as in a negative continued fraction q_(k+1) = a q_k − q_(k−1), with λ as a carry.
- **Unipotents.** For blocks A, B of equal length, A⁻¹B = [[1, 0, 0], [α, 1, β], [γ, 0, 1]] with β = b_A − b_B and
  γ = (c_B − c_A)/2^n. In the basis (e_H, e_λ, e_E) this is lower unitriangular, so every such quotient lies in one
  Heisenberg group U, the unipotent radical of the flag ⟨e_E⟩ ⊂ ⟨e_λ, e_E⟩
  ([`scripts/unipotent.py`](scripts/unipotent.py), all 155 pairs with |w| ≤ 6).
- **Hammersley digits.** Relative to the anchor J O^(n−1) Z, the abelianized digits (γ, β) of the 2^(n−2) blocks of
  length n are, up to scaling, exactly the Hammersley set {(k, rev_(n−2)(k))}, a (0, n−2, 2)-net in base 2. So the
  invertible orbit is a walk whose step at scale 2^n is an anchor times a perfectly equidistributed Heisenberg digit.
  This is THM-5141's reversal graph, now seen as the digit set of a group walk.

### 3. Strong approximation at every odd prime (PROVED, THM-5162 (C))

The invertible blocks generate G_p = {g ∈ GL₃(F_p) : det g ∈ ⟨2⟩} for every odd prime p, and the analogous group
modulo every odd squarefree q. Proof in four lines:
1. u₁ = JOOZ⁻¹·JOZZ and u₂ = JOOOZ⁻¹·JOZZZ have abelianized coordinates (−1/4, −2) and (−3/8, −6), with determinant
   3/4. For p ≥ 5 they therefore generate the Heisenberg group U(F_p) (Burnside basis theorem); their commutator is
   the transvection I − (3/4)E_(E,H).
2. J_(H,E) = −1 and (J⁻¹)_(H,E) = 1, so J moves the flag of U to an opposite flag in every characteristic. U and an
   opposite unipotent radical generate SL₃(F_p).
3. The determinants 4 and 8 generate ⟨2⟩.
4. p = 3 is wave 23's exhaustive search. Squarefree q follows by Goursat.

The Heisenberg route genuinely fails at p = 3: over all 155 same-length pairs with |w| ≤ 6, the odd part of the gcd
of the 2 × 2 determinants is exactly 3. P52 is thus PROVED, beyond the exact search at p = 7 that it named.

**Every odd modulus (PROVED, THM-5162 (C′)).** The result lifts to all prime powers. Write u₁ = I + N. For odd p,
u₁^p = I + pN′ lies in the first congruence subgroup K₁, with N′ ≡ N + ((p−1)/2)N² ≠ 0 mod p. The adjoint
submodule of sl₃(F_p) that N′ generates is everything:
- for p ≥ 5 because the adjoint module is irreducible;
- at p = 3, and as a check at 5, 7, 11 and 13, by direct computation ([`scripts/padic_lift.py`](scripts/padic_lift.py)).

K₁ is uniform with Frattini subgroup K₂ (p odd), so the closure of Γ ∩ SL₃ is all of SL₃(Z_p). Goursat across
primes then gives Γ mod q = {g ∈ GL₃(Z/q) : det g ∈ ⟨2⟩} for every odd q. This is strong approximation with
exceptional set exactly {2}. The resulting singular series is the classical one for the first coordinate of a
primitive vector of Z³:

    S(N) = (7/8)ζ(3) ∏_{p | N, p odd} (1 − p^(−2)),

which depends on which primes divide N but not on their exponents.

**Corollary (Haar local law, PROVED for random words).**
- *From a primitive state:* from any state ≢ 0 mod p, long random invertible words give counts distributed mod p
  as the first coordinate of a uniform nonzero vector of F_p³. So ν_p(0) = (p² − 1)/(p³ − 1), slightly below 1/p.
- *With a source:* with a one-apex source, the only change is the source's p-divisible share π_p: the law is
  π_p δ₀ + (1 − π_p) ν_p.
- *Primes:* without p-divisible states, the prime singular series is ∏_{p odd} p³/(p³ − 1) = (7/8)ζ(3) = 1.0518.

### 4. The local law in the census (pre-registered P53, P56, P58; all CONFIRMED)

Ratio of the mean multiplicity over odd N prime to p to the overall mean, on [2^23, 2^24) (mode 2, no source):

| p | 3 | 5 | 7 | 11 | 13 | 17 |
|---|---:|---:|---:|---:|---:|---:|
| measured | 1.040 | 1.008 | 1.003 | 1.001 | 1.000 | 1.000 |
| Haar p³/(p³−1) | 1.0385 | 1.0081 | 1.0029 | 1.0008 | 1.0005 | 1.0002 |

In mode 2 the prime mean is 1.0526 times the overall mean, against (7/8)ζ(3) = 1.0518.

Prime powers (P65, registered after the lift was proved and before any class mean was computed): in mode 2 on
[2^23, 2^24),
- the mean over N ≡ 0 mod 9 is 0.998 times the mean over N ≡ 3, 6 mod 9;
- the mean over N ≡ 0 mod 25 is 1.000 times the mean over the other multiples of 5;
- the ratios (p | N)/(p ∤ N) are 0.8845, 0.9603, 0.9791 at p = 3, 5, 7, against 1 − p⁻² = 0.8889, 0.9600, 0.9796. With a free source (mode 1,
[2^21, 2^22), [`scripts/srcdiv.c`](scripts/srcdiv.c)), π_p is measured directly. The predicted ratio
(1 − π_p)p³/(p³ − 1) and the measured ratio are:

| p | π_p | predicted | measured (words) | census (values) |
|---|---:|---:|---:|---:|
| 3 | 0.0805 | 0.9548 | 0.9562 | 0.957 |
| 5 | 0.0308 | 0.9770 | 0.9771 | 0.977 |
| 7 | 0.0039 | 0.9990 | 0.9992 | 0.999 |
| 11 | 0.0050 | 0.9957 | 0.9958 | 0.996 |
| 13 | 0.0012 | 0.9993 | 0.9994 | 0.999 |

The relatively large π_11 and the visible effect at 17 match THM-5128 (F): 3, 5, 11 and 17 are the primes with a
balanced divisible one-apex state at small h. Compare the 11-block orbit with resets: 0.405, 0.672, 0.930, 0.958,
0.973 at p = 3, …, 13, and 0.164 for primes (wave 21).

So the higher-rank blocks wash out every local bias except the one carried in by the rank-one source. This is the
GL₃ structure at work, and it bypasses HYP-9660 (B) (the uniform gap of the rank-one digit operators) for every
primitive source.

### 5. Density one: what is proved, and what remains

*Reduction (elementary).* Let model(N) be the window mean times the class factors. Put
V(X) = mean over odd N in [X, 2X) of (r(N)/model(N) − 1)². By Chebyshev, the share of odd N in [X, 2X) with
r(N) = 0 is at most V(X). So density one follows from V(X) → 0, an L² local limit law.

*Data (P55, CONFIRMED; [`scripts/l2test.py`](scripts/l2test.py); model mod 9, 5, 7):*

| window | 2^15 | 2^17 | 2^19 | 2^21 | 2^23 | 2^25 |
|---|---:|---:|---:|---:|---:|---:|
| mode 1: V | 0.0985 | 0.0662 | 0.0459 | 0.0323 | 0.0228 | 0.0166 |
| mode 1: V·mean | 7.77 | 7.71 | 7.86 | 8.12 | 8.43 | 8.96 |
| mode 1: share with r < model/2 | 0.043 | 0.016 | 0.0052 | 0.0010 | 0.00011 | 0.000006 |
| mode 2: V | 0.171 | 0.116 | 0.0805 | 0.0560 | 0.0390 | — |
| mode 2: V·mean | 6.89 | 6.90 | 6.99 | 7.13 | 7.26 | — |

The multiplicities behave like overdispersed Poisson counts around the Haar model, with variance about 8 times the
mean. So V ≈ 8/mean ∝ X^(−0.28). The slow drift of V·mean (+0.15 per doubling in mode 1) is what a residual
structured variance of about 0.003 would give (unmodelled local factors, such as powers of 3, and the smooth
archimedean slope within a window).

*Where a proof stands, in the Bourgain–Kontorovich template.*
- **Local / major arcs.** PROVED at each fixed squarefree q (THM-5162 with the Markov chain on G_q).
  - *Uniformity in q:* would follow from expansion. Salehi Golsefidy–Varjú (2012) applies to the commutator
    subgroup, whose Zariski closure SL₃ is perfect; the transfer to the semigroup walk is the usual
    μ̌ * μ argument. CITED, not written out.
  - *Determinant part:* it mixes slowly (a walk on ⟨2⟩ mod q), but it does not bias the count, because G_q acts on
    nonzero vectors through SL₃ orbits.
  - *Numerically the gap is uniform (P64, CONFIRMED;
    [`scripts/gapwalk.py`](scripts/gapwalk.py)).* The walk applies one of the 16 invertible blocks with |w| ≤ 5,
    uniformly at random, to a nonzero vector of F_p³. Its second-largest eigenvalue modulus is:

    | p | 3 | 5 | 7 | 11 | 13 | 17 | 19 | 23 | 29 | 31 | 37 | 41 | 43 |
    |---|---|---|---|---|---|---|---|---|---|---|---|---|---|
    | λ₂ | 0.565 | 0.409 | 0.437 | 0.327 | 0.335 | 0.331 | 0.317 | 0.296 | 0.288 | 0.289 | 0.288 | 0.285 | 0.282 |

    It falls toward about 0.28 as p grows. A block of this walk averages 5 letters, so per letter the rate is about
    0.28^(1/5) ≈ 0.78 at large p and at most 0.565^(1/5) ≈ 0.89 at p = 3. The rank-one digit operators of
    HYP-9660 (B) run at 0.932 per digit, attained at modulus 7. More important than the size is the trend: λ₂ does
    not creep toward 1 as p grows. This is the numerical face of expansion: the higher-rank walk supplies the
    uniformity in the modulus that the rank-one theory lacks.
- **Archimedean.** Not done. It needs a renewal theorem for the norm cocycle of an infinitely generated semigroup,
  with congruence twists.
- **Minor arcs.** OPEN. The natural bilinear form is the one above: N = b(H + λ) − E − QH, a reader digit b
  (Hammersley-distributed) against a source state. The B–K minor-arc estimates lose too much at δ ≈ 1.28. Their
  Zaremba method needs δ near 2 in these units; Huang's improvement needs about 5/3.
- **The closest proved case is the Apollonian packing.** Its critical exponent is 1.3056 against our 1.28, and
  density one was proved there by Bourgain–Kontorovich (2014). That proof uses arithmetic subgroups inside the
  symmetry group: the stabilizer of a circle, whose orbits are binary quadratic forms. Here the analogue is the
  Heisenberg group U of §2. Its orbits are 2-parameter polynomial families, but U lives in the group and not in the
  semigroup: B = A·u means that replacing a block by another of the same length applies u, while u itself needs A⁻¹.
  Whether some semigroup family realises a full U-orbit is the concrete open question this wave leaves.

### 6. Chan–Pak by contradiction, re-mapped on the invertible orbit

Suppose N > 6223 is the least odd number that is not a count. By THM-4094 it is a prime or 7p. In the invertible
orbit:
- no congruence class is obstructed (THM-5162);
- no reset or confinement can force a factor (THM-5120, THM-5162 (E));
- primes even get a favourable singular series ((7/8)ζ(3) times the source factor).

So N would be a purely archimedean, sporadic hole in a supercritical, Haar-distributed orbit whose multiplicity
there is about 0.92·mean(N). In the overdispersed-Poisson model of §5, a hole at multiplicity m has probability
about 8^(−m/7). That is e^(−110) at m = 369, which is the scale of 2^23, so the expected number of holes above 2^13
is far below 1. This is the heuristic, not a proof. The single missing theorem of wave 23 now takes its cleanest
form:

> **Missing theorem (invertible form).** For the invertible root-change semigroup, V(X) → 0. That is, the counts of
> a Zariski-dense, reset-free, supercritical semigroup orbit in GL₃(Z[1/2]) obey an L² local limit law with the
> Haar singular series.

That gives density one. Pointwise, with an explicit rate for the minimum multiplicity, it would give HYP-9682 and
hence Chan–Pak.

### Verdicts and new predictions

[Pre-registration file](scripts/wave24_prereg.md) (2026-11-14T23:10:34Z; addendum before any mode-2 window
statistics were printed).
- **P51 CONFIRMED.** Mode 1 per-window misses fall to 0 from [2^13, 2^14) on; [2^21, 2^22) misses none.
- **P52 PROVED** (THM-5162 (C)), for every odd prime, not only the registered p = 7.
- **P53 CONFIRMED.** The ratio at p = 3 is 0.957 (registered interval [0.80, 1.00]). At 5, 7, 11, 13 the ratios are
  0.977, 0.999, 0.996, 0.999 (registered [0.93, 1.02]). The prime mean is 0.92 of the mean (registered ≥ 0.5).
- **P54 CONFIRMED.** The fitted exponents are 0.2777 for all N and 0.2771 for primes (registered: a difference of at
  most 0.03).
- **P55 CONFIRMED.** V falls from 0.0662 to 0.0323 between 2^17 and 2^21.
- **P56 CONFIRMED.** On [2^21, 2^22) the ratios are 1.042, 1.008, 1.003, 1.001, 1.000, against 1.0385, 1.0081,
  1.0029 (each within 0.01) and 1 (within 0.005).
- **P57 CONFIRMED.** Mode 2 misses nothing in [2^21, 2^24), and nothing in (179,367, 2^26] (minimum multiplicity
  68 on [2^25, 2^26), prime mean 1.0525 times the mean).
- **P58 CONFIRMED.** At p = 3 the ratio lies in [0.955, 0.960] on every window from 2^16 to 2^23. The law
  (1 − π_p)p³/(p³ − 1) holds within 0.0015 at the five primes measured.
- **P65 CONFIRMED** (fourth addendum). The ratios are 0.998 (mod 9; registered [0.99, 1.01]) and 1.000 (mod 25;
  registered [0.98, 1.02]). The ratio (3 | N)/(3 ∤ N) is 0.8845, within 0.0044 of 8/9.
- **P64 CONFIRMED** (second addendum, registered before the program was written). λ₂ ≤ 0.565 for every odd prime
  p ≤ 43; the bound registered was 0.60.
- **P62 CONFIRMED** (registered in the third addendum before the run). Blocks with |w| ≤ 5 and a free source miss
  exactly 48 odd numbers below 2^26, the largest 12,811. So 16 invertible matrices look cofinite: compare the 63
  matrices, singular ones included, of HYP-9677.
- **P63 CONFIRMED** (same addendum). The pure-block family (J and J O^k Z) misses 3823 odd numbers below 2^26, the
  largest 67,076,437, and every window from 2^12 to 2^25 holds more than 10 misses. One block per scale is not
  enough; the Hammersley digits are.

Open predictions (third addendum of the pre-registration file, registered before these runs):
- **P59.** Mode 1 misses no odd N in (6223, 2^30]. So far it holds to 2^26.
- **P60.** Mode 2 misses no odd N in (179,367, 2^28].
- **P61.** In mode 1, V(X)·mean(X) stays in [7, 11] on every window up to [2^27, 2^28). So far it holds to
  [2^25, 2^26), with 8.96. The run to 2^28 was stopped on an overloaded machine (load 12).

## Wave 25: hosts and cores — no Heisenberg orbit, the cyclic join in Z[sqrt(-2)], and band fullness

Owner prompt (2026-11-14):
- look for a full Heisenberg orbit inside the semigroup;
- choose a host from the absolute target to close Chan–Pak, thinking of cores with different ratios and of more
  creative constructions;
- once the target and the first core are fixed, the second core's two counts lie in a bounded arithmetic
  progression, and what remains is to show that this progression contains an actual tournament profile.

**Chan–Pak OPEN.** The short answer of this wave:
1. **No Heisenberg orbit (PROVED, negative).** At a fixed determinant the semigroup contains only J-insertions,
   which grow at least threefold per three letters. So no unipotent (polynomial) family fits inside it
   ([THM-5172](../../../01-canon/theorems/THM-5172-hosts-and-cores-no-heisenberg-orbit-profile-congruences-and-the-cyclic-join-as-a-z-sqrt-minus-two-inner-product.md) (A)).
   Polynomial freedom has to come from the cores' profiles instead.
2. **The host in Z[√−2].** The host is the other session's cyclic join
   ([chanpak_absolute_reflection_20261114](../chanpak_absolute_reflection_20261114/README.md)). It is an inner
   product in Z[√−2]: 3N = Re(α_A ᾱ_B) with α = (f − 3h) + f√−2. The first core's ratio h/f fixes the angle of α_A.
   The second core's bounded progression is the lattice line perpendicular to α_A.
3. **Profile congruences (PROVED).** f ≡ (−1)^n h (mod 3), f ≡ 0 (mod 4), and f ≡ 0 (mod 8) at even order. So a
   prime target needs a core of even order.
4. **The remaining task is band fullness.** The progression meets an actual profile as soon as the (h, f) cloud of
   the second core's order fills its ridge band. Occupancy jumps from 0.30 (order 9) to 0.71 (order 10); a
   calibrated model predicts 0.89 at order 11 and complete fullness at order 12.
5. **Coverage.** With all cores of order ≤ 10, every odd N in (36,833, 2,398,737] is a cyclic-join count.
6. **Reduction (argument).** Band fullness at every order ≥ 12, together with the census below 2^32, would close
   Chan–Pak ([HYP-9692](../../hypotheses/HYP-9692-band-fullness-two-core-hosts-close-chan-pak-if-every-order-from-twelve-fills-its-ridge-band.md)).

### 1. The Heisenberg question answered (PROVED, THM-5172 (A))

On every reachable state, 0 ≤ E ≤ H, 1 ≤ h ≤ H and L ≥ max(h, E). These are preserved by each letter (induction),
and checked on all 7,478,943,733 reachable states with H ≤ 2^18 ([`scripts/invariants.c`](scripts/invariants.c)).
Hence:
- H never decreases;
- H(J³s) = 3H + 3L − E ≥ 3H;
- a fixed block skeleton yields O((log X)^(r+1)) counts up to X. The skeleton J^a JOZ J^b gives 45, 211, 528,
  946, 1486 words with count ≤ 10^3, 10^6, 10^9, 10^12, 10^15 ([`scripts/polylog_skeleton.py`](scripts/polylog_skeleton.py)).

So a family P·u_t·Q with u_t unipotent cannot lie in the semigroup along any infinite arithmetic progression of t.
Its counts would grow polynomially, while the fixed-determinant class grows polylogarithmically. Hence there is no
full Heisenberg orbit, nor any orbit of a unipotent group: the arithmetic-stabilizer mechanism that
Bourgain–Kontorovich used for Apollonian packings has no counterpart here. The Heisenberg group of THM-5162 lives
only in the group generated by the blocks; inside the semigroup, its digits are the finite Hammersley nets.

### 2. The host and its Z[√−2] geometry (KNOWN host, NEW reading; THM-5172 (B), (C))

The cyclic join T = A → B → c → A has the count N = (f_A − h_A)(f_B − h_B) + 2h_A h_B. Here h = P(1) = H and
f = P(2) = F₂, and the identity is that note's; brute force confirms it on 300 random joins
([`scripts/join_check.py`](scripts/join_check.py)). With u = f − 3h, its reflection identity becomes

    3N = u_A u_B + 2 f_A f_B = Re(α_A ᾱ_B),   α = u + f√−2 ∈ Z[√−2],

which held on 100,000 random census pairs.
- **The first core sets the line.** Fixing the target N and the first core A puts the second core on the lattice
  line Re(α_A β̄) = 3N, that is a f_B − c h_B = N with a = f_A − h_A and c = f_A − 3h_A. It is bounded by h_B ≥ 1
  and f_B ≥ (n_B + 1)h_B.
- **Ratios set the slope.** The ratio r = h_A/f_A sets the slope c/a = (1 − 3r)/(1 − r), which lies between 1/3
  (r = 1/4, the 3-cycle) and 1 (transitive cores). The scale a decides which order's cloud the line crosses.
- **Congruences (PROVED).** P(−q) = (−1)^n P(q), with even odd-cycle weights, gives f ≡ (−1)^n h (mod 3),
  f ≡ 0 (mod 4) and f ≡ 0 (mod 8) at even order. Every profile of every order ≤ 10 obeys them (1,602,452 order–profile pairs); that is why the
  profile set has no step (2, 0) or (2, 12). Consequences:
  - N ≡ 2h_A h_B (mod 3), or N ≡ 0 when both cores have odd order;
  - 3N ≡ u_A u_B (mod 32);
  - a core joined to itself gives the norm 3N = u² + 2f², hence N ≡ 3 (mod 8): 3, 11, 99, 779.

### 3. The profile census and the ridge band (FINITE-EXACT)

[`scripts/prof.c`](scripts/prof.c) and [`scripts/profhf.c`](scripts/profhf.c) compute (h, f), and F₃ through
order 9, for every class that `gentourng` lists. Distinct profiles by order:

| order | 3 | 4 | 5 | 6 | 7 | 8 | 9 | 10 |
|---|---:|---:|---:|---:|---:|---:|---:|---:|
| distinct (h, f) | 2 | 3 | 8 | 24 | 156 | 1,727 | 50,173 | 1,550,357 |

The cloud of order n is a thin band around a ridge f ≈ α_n + β_n h ([`scripts/band.py`](scripts/band.py)):

| order | ridge | residual sd | sd(f given h) / sd(f) | occupancy B_n(0.25) | B_n(0.5) | B_n(1) |
|---|---|---:|---:|---:|---:|---:|
| 8 | f = 1351 + 11.70 h | 185 | 0.102 | 0.300 | 0.337 | 0.312 |
| 9 | f = 6115 + 13.68 h | 766 | 0.085 | 0.297 | 0.302 | 0.272 |
| 10 | f = 30574 + 15.69 h | 3545 | 0.072 | 0.751 | 0.714 | 0.664 |

Occupancy is the share of admissible points of the band that are actual profiles. Profiles are much rarer than
classes: converse classes share P_T, among other coincidences. A Poisson model with that collision factor (4.60,
fitted on the exact order-10 census) reproduces order 10 from 100,000 random tournaments (0.713). It then predicts
0.892 at order 11 and 1.0000 at order 12, where there are 112 classes per admissible band point
([`scripts/fullness_model.py`](scripts/fullness_model.py), [`scripts/sample.c`](scripts/sample.c)). By counting,
band fullness arrives at order 12.

### 4. What the two-core host covers (FINITE-EXACT; [`scripts/cover.py`](scripts/cover.py))

All ordered pairs of census profiles are used, with padding by transitive blocks, (h, f) → (h, 2^k f).

| cores | missed odd N ≤ X | largest complete stretch | primes |
|---|---|---|---|
| order ≤ 9 | 10,330 below 2·10^5 | — (misses throughout) | min 0 representations on [10^5, 2·10^5], mean 3.4 |
| order ≤ 10 | 3,952 below 2·10^6, all ≤ 36,833 | (36,833, 2,398,737] | min 24, mean 48 on [10^5, 2·10^5] |

- **The bottleneck at order ≤ 9** was the residue class. N ≡ 0 (mod 3) was never missed from 2^13 on, because two
  odd-order cores land there. The other classes need an even-order core, and order 8 is poor.
- **Order 10 fixes it.** Even-order cores with a 71%-full band appear, and the misses in [2^16, 2·10^5] fall from
  5,451 to 0.
- **The reach matches the window theory.** Beyond 2,398,737 misses return (10,491 in [2.4·10^6, 10^7], all with
  N ≢ 0 mod 3), exactly where order 10's largest usable reader scale runs out:

  | order | a*(n) = s_n/(12β_n) | window covered |
  |---|---:|---|
  | 10 | 18.8 | [8.9·10^4, 3.6·10^6] |
  | 11 | 85.2 | [5.4·10^5, 9.8·10^7] |
  | 12 | 428.8 | [3.6·10^6, 3.2·10^9] |

  These come from [`scripts/windows.py`](scripts/windows.py).

**The absolute target chooses its host.** For the prime 100,003, the reader TT₂ ((h, f) = (1, 4)) fixes the line
3f_B − h_B = 100,003. The census supplies the order-9 profile (1997, 34000) on it. The resulting 12-vertex
tournament TT₂ → B → c → TT₂ has exactly 100,003 Hamiltonian paths by direct count
([`scripts/cert_100003.py`](scripts/cert_100003.py)). Its outgoing masks are [2046, 2044, 3192, 2800, 4064, 3520,
2944, 3844, 2572, 3108, 2376, 3].

### 5. The reduction: band fullness closes Chan–Pak (argument; HYP-9692)

1. **A guaranteed point.** A reader of scale a crosses the ridge band of order n in an h-interval of length at
   least 2s_n/β_n. Admissible points on its line recur at most every 24a in h. So if B_n(1) is full and
   a ≤ a*(n), the progression contains an actual profile. This is the owner's remaining task, solved by fullness.
2. **Readers.** Readers TT_k step by a factor of about 2. A window spans U_n/L_n ≈ 2.1. So readers, plus a few
   small cores for margin, cover [L_n, a*(n)U_n].
3. **Overlap.** The windows of consecutive orders overlap with growing margin: a*(n) grows about fivefold per
   order, L_n about sixfold. Each N lies in two consecutive windows, one of even order, which takes care of the
   mod-3 condition.
4. **Below.** The census covers every odd N < 2^32, and order 12's window starts at 3.6·10^6.

Hence band fullness at every order ≥ 12, together with the persistence of the measured band geometry, gives
Chan–Pak.

The second count is exactly what makes this work. A one-dimensional local limit law for H alone would leave gaps:
the bulks of consecutive orders are a factor of about (n + 1)/2 apart, while each is only about twice as wide. The
reader's scale a slides the second core's line across orders.

What a proof still needs is a two-dimensional local limit theorem for (H, F₂) of random tournaments: at scale 1 on
the admissible sublattice, in the ridge band, from order 12 on. That is a statement about random tournaments. The
data say it holds already at order 12 (by counting) and is 71% true at order 10.

### Verdicts and new predictions

[Pre-registration file](scripts/wave25_prereg.md) (2026-11-15T01:02:30Z, with addenda before each run).
- **P66 REFUTED.** Cores of order ≤ 9 miss 9,981 odd N in (1000, 2·10^5] (registered: none).
- **P67 REFUTED.** The prime 100,291 has no representation with cores of order ≤ 9.
- **P68 CONFIRMED.** Band(0.5, 0.5) occupancy at order 10 is 0.714 (registered [0.40, 0.75]). No profile violates
  the congruences.
- **P69 REFUTED.** Bulk-only second cores of order ≤ 9 miss 36,774 odd N above 1000.
- **P70 CONFIRMED.** With order-10 cores the misses in [2^16, 2·10^5] fall from 5,451 to 0 (registered: at least
  tenfold).
- **P71 CONFIRMED.** Below 2^16 the misses change only by a factor of 1.23.
- **P72 CONFIRMED.** With cores of order ≤ 10 there is no miss in (36,833, 10^6], and none in [10^6, 2·10^6]
  either (registered: fewer than 1%).
- **P73 HALF CONFIRMED.** The model gives 1.0000 at order 12 (registered ≥ 0.999) but 0.892 at order 11
  (registered ≥ 0.95). The registered class count of order 12 was misquoted and corrected before the fit
  (A000568(12) = 154,108,311,168).

Open predictions:
- **P74.** At order 11, the band(0.5, 0.5) occupancy lies in [0.80, 0.95]. Testable by sampling classes, since the
  exact census of 903,753,248 classes is out of reach here.
- **P75.** Adding the order-11 cores pushes the complete stretch of cyclic-join counts beyond 5·10^7.
- **P76.** Every prime in (36,833, 2,398,737] has a representation whose first core has order ≤ 4.

## Wave 26: fans, slides and reservoirs — unit extensions never branch twice, a mod-3 law for passive channels, and marked two-block completeness

Owner prompt (2026-11-15), in short:
- caterpillars collapse into products, so the useful lead is fan invariance itself;
- ask whether two armed fan vertices, one at +1 and one at −1, stay stable under attachment, since they would
  generate the full binary tree, and search for such pairs;
- use the growing window (the 13/29 sliding cores of
  [chanpak_window_descent_20261115](../chanpak_window_descent_20261115/README.md), half-steps 145 and 39) as a
  repair, with a different left reservoir whose ordinary count scales the repair radius and whose marked endpoints
  set the residue class and the interval;
- find other constructions that slide along the number line without changing what the surrounding tournament
  reads, using number theory, primality and invariants.

Prior art read first:
- [chanpak_leaf_fans_20261115](../chanpak_leaf_fans_20261115/README.md): fan law, cross-port law, movable-port
  census and parallel fans, all in the leaf model where a newborn beats exactly one old vertex;
- the window-descent notes: matched slides, the CRT amplifier, the exposure theorem and primitive lanes;
- [CUT_STATE](../chanpak_selective_target_20261114/CUT_STATE.md): the exact transfer used throughout.

**Chan–Pak OPEN.** The short answer of this wave:
1. **No stable plus/minus fan pair, even in the widest one-vertex model (FINITE-EXACT, NEGATIVE).** Allow
   arbitrary out-sets, re-chosen at every node. Still no tournament of order ≤ 8 has both an armed +1 child and an
   armed −1 child. The children freeze their sign
   ([THM-5182](../../../01-canon/theorems/THM-5182-passive-channels-of-a-marked-one-apex-spine-mod-three-law-seam-lemma-and-the-fate-of-plus-minus-fans.md) (D),
   [HYP-9712](../../hypotheses/HYP-9712-unit-extensions-never-branch-twice-no-tournament-has-armed-plus-and-minus-children.md)).
   Two-vertex base-4 steps die just as early.
2. **Which channels can slide (PROVED).** In a marked spine a block's neighbour reads only (a, q, s, e).
   - The n-channel slides at flat weights (the CRT amplifier).
   - The t-channel slides at positional weights.
   - The count is an explicit pair expansion (THM-5182 (A)).
3. **A mod-3 law (PROVED).** T_v ≡ (−1)^(m+1) S_v (mod 3) for every vertex of every tournament. Hence the passive
   t of any mark set is fixed mod 3 by (a, s) (THM-5182 (B)).
4. **Slides are rare, and a reservoir cannot place them (census + seam lemma).**
   - Only cores of count 13 (order 5), 29 and 115 (order 7) have nonzero matched slides; there is no unit pump.
   - A slide that matches the right-facing data moves t geometrically.
   - The seam term s_L t_B ≥ 2^ℓ g a_B forces any reservoir to carry an exponentially wide palette.
5. **The reservoir works as an interface, not as a window (FINITE-EXACT; conjecture HYP-9702).** Marked reservoir
   L ⇒ marked repair B, one apex, no slides. With blocks of order ≤ 7 this attains every odd N in
   [23, 1,542,343); with orders (8,7), every odd N in [23, 8,667,349)
   ([HYP-9702](../../hypotheses/HYP-9702-marked-two-block-completeness-every-odd-count-except-7-and-21-is-a-reservoir-repair-pairing.md)).
   Fixed-interface palettes are trivial (mean fiber 1.05). The coverage comes from the variety of interfaces.
6. **Primality is invisible; only 3 is special.** Primes have exactly the multiplicity of their residue class.
   N ≡ 0 (mod 3) is favoured by 26%. Moduli 5 and 7 differ by under 2%.

### 1. Stable plus/minus fans with arbitrary out-sets (THM-5182 (D), HYP-9712)

Adding a vertex x that beats exactly X gives, by the odd-cycle expansion,

    H(T + x) = H + 2W(X),   W(X) = Σ_{u∈X, v∉X} Π(u,v),

where Π(u,v) sums H(T − P) over the even-order paths u → v. [`scripts/armed.c`](scripts/armed.c) builds Π by
subset DP and enumerates all 2^n cuts by Gray code. The formula is checked against direct DP on 1,000 random
(T, X).
- A unit extension is a cut of value (H ± 1)/2.
- T is armed when it has unit cuts of both values.
- The leaf model is |X| = 1. Here about half the unit out-sets are leaves or reverse leaves (750 of 1,391 +1 units
  at order 8), and the rest have intermediate sizes ([`scripts/xsize.c`](scripts/xsize.c)).

| order | eligible classes (H ≠ 3, 11) | armed | fraction |
|---:|---:|---:|---:|
| 5 | 8 | 5 | 0.625 |
| 6 | 50 | 18 | 0.360 |
| 7 | 448 | 50 | 0.112 |

H = 13 is never armed at orders 5–7, and neither are 17, 25, 27, 29 or 31 at order 7.

**Adaptive trees** ([`scripts/tree2.c`](scripts/tree2.c)). The out-set is re-chosen freely at every node, and a
branch whose target is 7 or 21 is not required. All classes of orders 3–8 were used as seeds (6,879 at order 8).
- Every seed with H ∉ {3, 5, 11} has depth ≤ 1.
- Depth 2 occurs only at H ∈ {5, 11}, where one branch is excused by 21.
- Depth 3 occurs only at H = 3: 3 → 5 → {9, 11} → {17, 19, 23}.
- Equivalently, no tournament of order ≤ 8 has an armed +1 child and an armed −1 child. Order 9 confirms it
  ([`scripts/childpair.c`](scripts/childpair.c)): 308 of the 191,536 classes are armed, 10 have an armed +1 child,
  27 an armed −1 child, and none has both. Arming collapses with the order: 11%, 1.9% and 0.16% of all classes at
  orders 7, 8 and 9.

[`scripts/childpair.c`](scripts/childpair.c) shows the mechanism, **sign freezing**. Typical order-6 seeds:

| seed H | +1 children | −1 children |
|---:|---|---|
| 19 | 39, only +1 extensions | 37, only −1 extensions |
| 9 | 19, armed | 17, never armed |
| 33 | 67, no unit extensions | 65, armed |
| 65 (order 7) | 131, one sign only | 129, armed |

The owner's binary tree fails at its second level. Branching survives on one side only.

**Two vertices per step** ([`scripts/base4.c`](scripts/base4.c)). Steps H → 4H + ε with ε ∈ {±1, ±3}, from all
seeds of orders 3–6: depth ≤ 1, except H = 3, which reaches depth 2 (every odd number from 33 to 63).

This extends the leaf-model obstructions of chanpak_leaf_fans_20261115 (the cross-port law and the repeat-stable
classification) to arbitrary out-sets and adaptive choices. Digit trees of exact affine moves need multi-vertex
steps whose intermediate defects are not units.

### 2. Passive channels: what a block can change without changing what its neighbour reads (THM-5182 (A))

The exact cut transfer is (h, N, S, d) → (ah, aN + (n−a)h + tS + fd, qS + hs, he). Unrolled, the apex count of a
marked spine B_1 ⇒ ... ⇒ B_r is

    N = A + Σ_i (n_i − a_i)A/a_i + Σ_{i<j} s_i t_j A_{<i} Q_{ij} A_{>j} + Σ_i e_i f_{i+1} A_{<i} A_{>i+1},

with Q_{ij} = Π_{i<k<j} q_k. Each term is a class of odd apex cycles.

- *Right-facing data* (a, q, s, e) compose as (a_Y a_Z, q_Y q_Z, q_Z s_Y + a_Y s_Z, a_Y e_Z).
- *Left-facing data* compose by the mirror rule (a_Y a_Z, q_Y q_Z, q_Y t_Z + a_Z t_Y, a_Z f_Y). This follows from
  converse duality and is checked on 400 random spines ([`scripts/compcheck.c`](scripts/compcheck.c)).

A neighbour reads only the right-facing data, so any block family with equal (a, q, s, e) is interchangeable.
Switching block i changes N by (Δn·h + Δt·S + Δf·d)·A_{>i}:

| channel | weight | what it gives |
|---|---|---|
| n | A/a_i, flat | the CRT amplifier: 13·29 makes 5 → 145 and 3 → 39 |
| t | S_{<i}A_{>i}, grows ≈ q/a per block | positional digits |
| f | ≈ h, flat | a small correction |

### 3. Slides: census, geometric t-slides and the seam lemma

**Census** ([`scripts/slides.c`](scripts/slides.c), [`scripts/slides2.c`](scripts/slides2.c)). Search every
strong marked core for transitive tiles z, w with Cz ≡ wC in (a, q, s, t), using the guards of
chanpak_window_descent_20261115:
- **Orders 3–4:** nothing.
- **Order 5, width ≤ 14:** only the count-13 core, at width 12, with Δ ∈ {±10, ±20, ±30}. All 48 solutions are
  confirmed by direct DP on the 20-vertex spines.
- **Order 6, width ≤ 30:** nothing (7.7 min exhaustive).
- **Order 7, width ≤ 28:** 16 slides:
  - two count-29 classes at width 28, with Δ ∈ {±6, ±32}; the Δ = 6 one agrees with the literal replay of the
    window-descent notes;
  - a count-115 core at width 24, with Δ = ±46. Its half-step 23 divides 115, so it alone meets only 5 residues.

So there is no unit pump (|Δ| = 2), and matched slides are rare objects.

**t-slides are geometric.** If Cz ≡ wC only in the right-facing data, the mirror rule gives
t(w^j C z^(k−j)) − t(w^(j+1) C z^(k−j−1)) = 2^(dj)(t(Cz) − t(wC)). A sliding core therefore produces a geometric
t-palette, never an arithmetic one.

**Seam lemma.** S_v ≥ H and T_v ≥ H: cut each Hamiltonian path at v. Put a reservoir L before a bank B whose last
unmarked vertex follows a transitive run of length ℓ. Then the seam term s_L t_B is at least 2^ℓ g a_B, while the
bank's window is O(ℓ). For fixed s_L, the attainable counts a_B n_L + f_B d_L + g(n_B − a_B) must bridge a gap of
t_B ≥ 2^ℓ a_B. So the reservoir needs an n-palette at least about 2^(ℓ−1) wide at fixed (g, s_L, d_L).

This is the precise reason the growing window cannot be placed by a reservoir of bounded size. It matches the
density bound of GROWING_WINDOWS. A reservoir works only if it is itself exponentially rich.

### 4. The reservoir as an interface: marked two-block completeness (THM-5182 (E), HYP-9702)

Drop the slides and let both blocks vary over all marked tournaments. For L ⇒ B with one apex,

    N = a_B(n_L − g) + g n_B + s_L t_B + d_L f_B,

a bilinear pairing of the reservoir's right-profile (g, n_L, s_L, d_L) with the repair's left-profile
(a_B, n_B, t_B, f_B). Wave 25's cyclic join is the case L all marked, B all unmarked.

[`scripts/palette.c`](scripts/palette.c) computes all profiles. It is checked against direct DP on 300 random
spines. Then [`scripts/twoblock.c`](scripts/twoblock.c) takes all pairs:

| block orders | profiles | pairs evaluated | every odd N in |
|---|---:|---:|---|
| (≤ 6, ≤ 6) | 2,831 | 8.0·10^6 | [23, 67,967) |
| (≤ 7, ≤ 7) | 52,313 | 2.7·10^9 | [23, 1,542,343) |
| (≤ 8, ≤ 7) | 1,693,156 × 52,313 | 8.9·10^10 | [23, 8,667,349) |

So every odd N < 8,667,349 except 7 and 21 is the count of a tournament of order ≤ 16 that has a vertex whose
removal leaves a two-block ordinal sum. 160 sampled witnesses, including the primes 999,983, 1,000,003 and
1,500,007, were rebuilt and re-counted by direct DP ([`scripts/twowit.c`](scripts/twowit.c),
[`scripts/verify_two.py`](scripts/verify_two.py)). Compare wave 25's unmarked cyclic join: (36,833, 2,398,737] with
cores of order ≤ 10.

**The repair radius is about 1.** A last block's palette is its fiber over (a, t, f). These fibers have mean size
1.05–1.09 at orders 5–8, and their longest run of consecutive odd apex counts is 4. So the owner's picture works,
but in a different form:
- the reservoir's count g does scale the repair, and its marked endpoints do set the residue class and the
  interval;
- the repair is the whole 4D left-profile cloud, not a window inside one interface.

Proving HYP-9702 is again a fullness problem, as HYP-9692 was. But the cloud is far denser. The marks enter (t, f)
linearly and n through a quadratic cut function, which points to Fourier analysis over the mark cube.

### 5. Number theory: the mod-3 law and the multiplicity pattern (THM-5182 (B))

**Mod-3 law.** For every tournament of order m and every vertex v, T_v ≡ (−1)^(m+1) S_v (mod 3).

*Proof.* In the squarefree algebra, H(X)² ≡ H(−X) (mod 3) coefficientwise. Indeed
[x^U]H(X)² = P_{T[U]}(2) ≡ P_{T[U]}(−1) = (−1)^|U| H(U), by the odd-cycle form (THM-5172 (B)) applied to every
induced subtournament. Combine this with the path identity H(X)P_v(−X) = −E_v(X) of the leaf-fan notes. ∎

Summing over the unmarked vertices gives t ≡ (−1)^m s + [m odd]·a (mod 3). This matches all 54 residue classes
in the census of orders 3–8, and explains why the t-values of every fiber differ by multiples of 3. The apex count
n itself is not fixed mod 3, 4 or 8 by the exterior.

**Multiplicities** ([`scripts/twomult.c`](scripts/twomult.c)). Odd N in [10^6, 1.4·10^6], family (7,7):

| class | mean representations |
|---|---:|
| N ≡ 0 (mod 3) | 215.9 |
| N ≡ ±1 (mod 3) | 171.0 |
| primes | 171.1 |
| composites | 188.5 |
| N mod 5 or mod 7 | 185.2–188.5 |

The minimum is 16. Primes are not harder than their residue class. The only local bias is the factor about 1.26
at 0 mod 3, the same prime that is special in THM-5162 (D) and THM-5172 (B).

### Reproduction

Profiles are regenerated with `gentourng m | ./palette run m > pal_m.txt` for m = 1..8. Order 8 has 1,761,280
lines and is not committed. [`scripts/twoblock.c`](scripts/twoblock.c) reads them from the working directory.
The unit-extension programs include [`scripts/armed_lib.h`](scripts/armed_lib.h). Outputs:
- `unit_census_order*.out`, `unit_trees_order*.out` (the order-8 trees as a summary) and
  `unit_children_order{8,9}.out`;
- `base4_trees_order*.out`;
- `slides_order5_w14.out`, `slides_order7_w28.out` and `slides_orders3_4_6.out`;
- `twoblock_{6x6,7x7,8x7}.out`, `twoblock_witnesses_7x7.out` and `twomult_7x7_summary.out`.

### Verdicts and new predictions

[Pre-registration file](scripts/wave26_prereg.md) (registered 2026-11-14 20:51 MST, addendum 21:25 MST, before
the runs they govern; the exploratory runs are listed there).
- **P77 REFUTED.** 36% of eligible order-6 classes are armed (registered ≥ 50%).
- **P78 REFUTED.** No seed of order ≤ 6 with H ∉ {1, 3, 5, 11} has adaptive depth ≥ 2 (registered: depth 5).
- **P79 CONFIRMED.** No depth 9, indeed no depth 2 beyond the 7/21-excused seeds.
- **P80 CONFIRMED.** Fixed port rules are a special case of the adaptive search, so depth ≤ 3.
- **P81 REFUTED.** (8,7) covers [23, 8,667,349) (registered ≥ 3·10^7).
- **P82 CONFIRMED.** No unit pump at order 7, width ≤ 28.
- **P86 CONFIRMED** (registered in HYP-9712 before its run). No order-9 tournament has armed children of both
  signs.
- **P83 VOID.** Its first half is automatic. The second half (no order-6 slide at any width) is confirmed only to
  width 30.

Open predictions:
- **P84.** Block orders (8, 8) cover every odd N in [23, X) with X ≥ 5·10^7 (HYP-9702). This needs about 3·10^12
  profile pairs, several hours here, and was not run.
- **P85.** The mod-3 multiplicity ratio tends to a limit in [1.2, 1.35]; other moduli are asymptotically unbiased.
- **P87.** The fractions of armed tournaments whose +1, respectively −1, children are armed both tend to 0, and
  their sum stays at most 1.
- **P88.** No marked core of order 8 has a matched slide with |Δ| = 2 at tile width ≤ 24.

## Wave 27: carriers, crossings and crowding — a finite list of digit-quadratic forms covers 333 to 10^8

Owner prompt (2026-11-15), in short:
- attack a fullness proof for the marked profiles, and try to prove sign freezing;
- work toward universal placement, raising the level of abstraction of the constructions until they become
  infinite;
- when something fails, ask why and what larger change gets around it, and pursue that by experiment, many times;
- find what the constructions have in common, an elegant form that also explains why 7 and 21 are impossible;
- place the other problems of the Hamiltonian-spectrum sessions around ours;
- take as context the other session's two-tape steering (two-sided O(N^(3/8)) brackets, O(N^(1/4)) on dyadic
  bands by fixing the carry), and think of the quasi-Riemann hypothesis.

Prior art read first:
- [chanpak_absolute_steering_20261115](../chanpak_absolute_steering_20261115/README.md): the two-tape carrier and
  its brackets;
- [chanpak_target_carry_20261115](../chanpak_target_carry_20261115/README.md): carry selectors and quartic shells;
- [chanpak_reservoir_slide_20261115](../chanpak_reservoir_slide_20261115/README.md): the count-45 unit slide,
  raw-unit components and slide entropy;
- [chanpak_viability_shape_20261115](../chanpak_viability_shape_20261115/LOCAL_DENSITY.md): the first moment of
  the fixed carrier;
- my waves 6–10 here (the two-hub frame identity H = Σ 2^α 3^β).

**Chan–Pak OPEN.** The short answer of this wave:
1. **Every one-apex carrier has an exact two-variable form (PROVED).** Guard, X-tape, block C, Y-tape, guard and
   one apex give a strong tournament with exactly

       F_C(X, Y) = qXY + tX + sY + a(p(X) + p(Y)) + n + f(X mod 2) + e(Y mod 2)

   Hamiltonian paths, for all X, Y ≥ 1. Here (a, q, s, e, t, f, n) are C's transfer data and p is the binary
   inversion surcharge (p(2u + 1) = p(u), p(2u) = p(u) + u + (u mod 2))
   ([THM-5192](../../../01-canon/theorems/THM-5192-carrier-forms-the-insertion-identity-monotone-unit-cuts-and-odd-cycle-crowding-at-7-and-21.md) (A)).
   This generalizes the other session's single count-45 carrier to every block.
2. **A finite list of forms covers everything tested (FINITE-EXACT; HYP-9722).** The 168 distinct forms of the
   strong marked cores of orders 3–5 attain every odd N in **[333, 10^9]**, and so does a 48-form basis chosen
   on [333, 10^7]. The order-6 forms alone cover [615, 10^8]. With the census below 333, the
   [Finite Carrier Conjecture](../../hypotheses/HYP-9722-finite-carrier-conjecture-every-odd-n-from-333-is-a-value-of-one-of-168-digit-quadratic-carrier-forms.md)
   implies Chan–Pak.
3. **Sign freezing has a mechanism (PROVED identity, FINITE-EXACT lemma).** Adding a vertex that beats X changes
   the count by H(T+x) − 2H(T) = E + R − Z. Through order 8 every cut with a double X→X^c crossing has defect ≥ 3.
   So unit cuts are monotone, and unit children have at most two types per sign
   ([HYP-9732](../../hypotheses/HYP-9732-crossing-lemma-a-double-crossing-cut-has-defect-at-least-three-so-unit-cuts-are-monotone.md)).
4. **Why 7 and 21 (FINITE-EXACT with classical inputs).** Multiplicativity, Busch's bound and odd-cycle crowding
   (Moon) leave exactly these two values unrealizable at the only orders allowed to them.
5. **Quasi-Riemann.** Each form is a binary problem: a hyperbola perturbed by a 2-regular digital term. Power-saving
   Fourier bounds would give density one form by form; every N needs the finite union plus a lower-tail bound on
   the clumping of representations.

### 1. Sign freezing: why unit moves cannot branch twice (THM-5192 (B), (C); HYP-9732)

**Insertion identity (PROVED).** Let x beat exactly X. A Hamiltonian path P of T has exactly 1 + #XX^c(P)
insertion slots: front, end, and every X^c→X junction. For binary words, [first ∈ X] + [last ∉ X] + #X^cX =
1 + #XX^c. All other Hamiltonian paths of T + x place x between u ∉ X and w ∈ X with w → u; call their number R.
Therefore

    δ(X) = H(T + x) − 2H(T) = E + R − Z,

where Z counts the paths of shape [X^c][X] and E the excess crossings.

**Crossing minimum (FINITE-EXACT, [`scripts/zer_all.c`](scripts/zer_all.c)).** Over every cut of every tournament
of orders 4–8, a double crossing (E ≥ 1) forces δ ≥ 3. At order 8 that is 1,620,135 of 1,747,520 cuts, and the
minimum 3 is attained. On 1,000 random order-9 classes the minimum is 15. So every cut with δ ≤ 1 is monotone:
each Hamiltonian path meets X in at most two blocks, the second at the end. On all 2,783 unit cuts of order 8,
E = 0 and R = Z ± 1 ([`scripts/zer.c`](scripts/zer.c)).

**Rigidity of unit children** ([`scripts/childtypes.c`](scripts/childtypes.c), [`scripts/sigpat.c`](scripts/sigpat.c)):
- at most two child types per sign through order 9 (307 armed order-9 classes: 281 have one +type and 26 two;
  289 have one −type and 18 two);
- all same-sign children of a seed share one continuation signature;
- the pair (armed, armed) never occurs.

**Armed does not imply strong** ([`scripts/strongarm.c`](scripts/strongarm.c)). 81 of the 133 armed order-8
classes are reducible. So the leaf-model raw-unit law of the other session does not survive arbitrary cuts, and
the rigidity has to come from monotonicity, not from strong components.

**Why it fails, and what larger change works.** A unit move needs a monotone cut whose reverse-junction covers
balance the [X^c][X] paths to within one. Such cuts are nearly unique. After the first move the newborn x takes
part in every later cut's [X^c][X] count, and the balance survives on one side only.
- *One vertex per step:* dead (wave 26).
- *Two vertices per step:* dead.
- *Affine trees of any small step:* blocked by this rigidity.

The larger change is to stop asking the count to move by exact affine steps, and instead let the core vary over a
finite list while two free tapes supply the scale. That is §2.

### 2. Carrier forms and the Finite Carrier Conjecture (THM-5192 (A), (D); HYP-9722)

**The form (PROVED).** The construction 𝒯_C(X, Y) is:
- a marked guard g₁;
- the X-tape: the binary digits of X below its leading 1, most significant first, mark = digit;
- the block C;
- the Y-tape: vertex j unmarked iff bit j−1 of Y is 1;
- an unmarked guard g₂;
- one apex.

Then H(𝒯_C(X, Y)) = F_C(X, Y). The proof has three steps:
1. *Prefix:* after g₁ and the X-tape the cut state is (1, p(X) + 1, X, X mod 2).
2. *Core:* C maps it to (a, a·p(X) + n + tX + f(X mod 2), qX + s, e).
3. *Suffix:* the Y-tape and g₂ send (h, N, S, d) to N + SY + h·p(Y) + d(Y mod 2).

Steps 1 and 3 are inductions on the bits, using p's two recursions.

Checks:
- the closed form against the transfer on 1,058,400 (core, X, Y) triples ([`scripts/closedform.c`](scripts/closedform.c));
- 48 random carrier tournaments of up to 18 vertices, re-counted by direct DP ([`scripts/verify_union.c`](scripts/verify_union.c)).

For the 3-cycle with one marked vertex the form is

    12XY + 6X + 3Y + 3(p(X) + p(Y)) + 5 + 2(X mod 2) + (Y mod 2).

**Coverage (FINITE-EXACT, [`scripts/unioncar.c`](scripts/unioncar.c)).**

| core set | distinct forms | Σ 2/q | outcome |
|---|---:|---:|---|
| order 3 | 4 | 0.667 | about 11% of odd N missed (to 2·10^7) |
| orders 3–4 | 20 | 1.667 | 10,853 misses to 10^8, density about N^(−0.58), unbounded |
| **orders 3–5** | **168** | **4.780** | **every odd N in [333, 10^9]** (the 32 misses above 21 are all ≤ 331) |
| order 6 alone | 1,950 | 12.597 | every odd N in [615, 10^8] |
| greedy basis (48 forms) | 48 | — | every odd N in [333, 10^9] |

The greedy basis was chosen on [333, 10^7] only ([`scripts/greedy_cores.c`](scripts/greedy_cores.c)):
- the 3-cycle form covers 84% alone;
- ten count-5 forms bring it to 99.9%;
- 36 order-5 forms (counts 9, 11, 13, 15) take sporadic leftovers.

**Multiplicity.** The representation count averages 5.03, 4.43, 4.48, 4.70 and 4.89 times ln N at 10^4 … 10^8.
That is the first-moment value Σ_C (2/q_C)·ln N = 4.78·ln N, since each form has about (N/q)·ln N pairs with
qXY ≤ N. The window minimum rises 4, 4, 8, 14, 18.

**Local structure.** Only 3 matters: near 10^7 the classes N ≡ 0, 1, 2 (mod 3) average 87.9, 63.0 and 76.2
representations. Primes average 69.1, the mean of their two classes; there is no primality effect.

**Failure, why, change (the threshold).** The orders-3–4 forms have Σ 2/q = 1.667, and every class mod 15 has mean
1.3–2.4·ln N. A Poisson/Borel–Cantelli model would then give finitely many misses. Yet zeros fall only like
N^(−0.58): 3.1·10^(−2), 8.0·10^(−3), 2.1·10^(−3), 5.7·10^(−4) and 1.5·10^(−4) per odd number in the decades
from 10^3 to 10^8. The misses cluster mod 3, mod 5 and mod 9. So the representations clump, and the threshold
is not the naive one.

Adding the order-5 forms brings counts 9, 11, 13, 15 and moduli q = 84–108 with new prime factors. That breaks the
clumping, and the misses stop at 331.

**Why this is the elegant form.** Every construction of the campaign is a word in the cut-transfer semigroup. One
apex over an ordinal spine gives the pair expansion of THM-5182. With two tapes around one block it collapses to
F_C: a hyperbola (the seam term s_X t_Y q_C), a 2-regular digital surcharge (the in-tape pairs) and two parity
bits (the signed channels at the guards). In this form:
- the slides of the window-descent notes, the unit pump and the reservoirs are special moves on (C, X, Y);
- the O(N^(3/8)) and O(N^(1/4)) brackets are tangent moves along the hyperbola;
- the reservoir's "ordinary count scales the radius" is the role of a and q;
- "fix the carry" is a statement about p.

### 3. Why 7 and 21 (THM-5192 (E))

1. **Multiplicativity.** H is multiplicative over strong components. 7 is prime and 21 = 3·7, so both need strong
   witnesses.
2. **Busch.** A strong n-tournament has H ≥ b(n) = 3, 5, 9, 15, 25, …. So a witness has order ≤ 4 for 7 and
   ≤ 6 for 21.
3. **Odd-cycle crowding.** H = 1 + 2c + 4d + … with c odd cycles and d disjoint pairs. Strong tournaments are
   pancyclic with at least n − 2 triangles (Moon).
   - *Census of orders 4–6* ([`scripts/cycles.py`](scripts/cycles.py)). Order 4 has c ≤ 2. Strong order 5 already
     has c ≥ 4 (at least 3 triangles and a 5-cycle). So c = 3, that is H = 7, never occurs.
   - *Order 6:* the strong classes have c + 2d ∈ {7, 8, 9, 11, 12, 13, 14, 15, 16, 18, 20, 21, 22}. That skips
     exactly 10, which is H = 21. Its neighbours H = 19 and H = 23 occur.

So 7 and 21 are the two values that crowding skips below the point where the Busch staircase starts producing
long runs. Above 333 the carrier forms take over. In this light the Fermat remark of the Busch note (7 = Φ₃(2)
divides no 2^j + 1) explains why Fermat-block constructions avoid 7. Crowding explains why nothing else reaches it.

### 4. The neighbouring problems, placed

| problem | where it sits | status here |
|---|---|---|
| Rédei, Grinberg–Stanley odd-cycle expansion | the positivity behind every transfer; H = 1 + 2c + 4d + … | used throughout |
| Moon pancyclicity, Busch's bound, Fermat primes | the small-order regime; why 7 and 21 fail | §3 |
| divisor problems, the hyperbola method | the multiplicative skeleton qXY + tX + sY of each form | §2 |
| Mauduit–Rivat digital functions (carry and Fourier properties) | the surcharge a·p(X); "fix the carry" | proposed tool |
| circle method, quasi-Riemann power savings | density one per form; every N needs more | HYP-9722 |
| Borel–Cantelli, Erdős–Rényi coverage | the log-density Σ 2/q and its clumping correction | §2 |
| Zaremba, Bourgain–Kontorovich, Apollonian orbits | thin-orbit local-global; density one versus every N | waves 24–25 |
| strong approximation (Heisenberg unipotents) | the local half for the invertible orbit | THM-5162 |
| Frobenius conductors, CRT | coprime slide steps 145/39 and the amplifier | window-descent notes |
| Erdős 3-smooth antichain sums | two-hub frames 2^α 3^β | waves 6–10 |
| Collatz functional graph and carries | binary addresses, carry decoders with terminating rank | target-carry notes |
| β-expansions, redundant digits | fat families (4 marks per vertex against growth 3) | waves 26–27 |
| Hart–Iosevich dot-product sets | the bilinear profile pairing of HYP-9702 | wave 26 |
| Rule 30, XOR gates | guarded bit forcing in carry selectors | other session |
| Jacobian conjecture / planar non-properness | whether locally compatible branches escape absolute bounds | other session |

The common thread is that a count is a bilinear reading of interface data, plus a carry-type digital correction.
Placement problems are the hyperbola part. Repair problems are the digital part. The impossible values live where
neither part has room yet.

### Reproduction

[`scripts/carrier_forms_orders3_5.txt`](scripts/carrier_forms_orders3_5.txt) lists the 168 forms as tuples
(a, q, s, e, t, f, n), with the 48 basis forms tagged. [`scripts/check_carrier_forms.py`](scripts/check_carrier_forms.py)
checks any range from that list and p alone, independently of the tournament code; it finds no miss in [333, 10^6]
for both lists.

The C enumerator [`scripts/unioncar.c`](scripts/unioncar.c) produced the 10^9 run, which is memory-bound and took
about 18 minutes. Its receipts are `unioncar_orders3_5_1e9.out`, `greedy_basis_1e9.out` and
`unioncar_order6_1e8.out`.

### Verdicts and new predictions

[Pre-registration file](scripts/wave27_prereg.md) (registered 2026-11-14 22:18 MST, with the exploratory runs listed).
- **P89 CONFIRMED.** Every odd N in [333, 10^9] is a value of the 168 forms (4.9·10^10 configurations).
- **P90 CONFIRMED.** On 1,000 random order-9 classes every double-crossing cut has δ ≥ 3 (minimum 15).
- **P91 CONFIRMED.** At most two unit-child types per sign: on the sample, and on all 307 armed order-9 classes.
- **P92 CONFIRMED.** The order-6 forms alone cover [615, 10^8] (registered X0 < 10^4).
- **P93 CONFIRMED.** Mean multiplicity near 10^8 is 4.89·ln N (registered [4.4, 5.2]).
- **P94 REFUTED.** The window minimum near 10^8 is 18 (registered ≥ 20).
- **P96 CONFIRMED** (registered in HYP-9722 before its run). The 48-form basis covers [333, 10^8] with no miss,
  and in fact [333, 10^9].

Open predictions:
- **P95.** The 168 forms cover [333, 10^10].
- **P97.** The minimum multiplicity on [10^9, 10^9 + 2·10^5] is at least 20.
- **P98.** The crossing lemma holds for every order-9 class.
- **P99.** Some list of at most 20 forms covers every odd N in [10^4, 10^9].

## Wave 28: bottlenecks — floor powers have rank one, so near-floor primes pay 3.5% per vertex

Session mac-mini-2026-11-15-bottlenecks. Owner prompt: "Attack the near-floor conjecture next through ribbon transfer
matrices and attack density one for the carrier forms, or try to prove the crossing lemma. Work to close chan-pak."

Every run was exploratory; none was pre-registered. New predictions P100 and P101 are filed in HYP-9742.

**Naming.** Here "near-floor conjecture" means the 7-21 repository's coverage statement:
- NF2(c): every odd N ≥ 25 is the count of a strong n-tournament with H < c·f(n), where f is the Moon–Busch floor.

The lower bound H ≥ f(n) + 3(hc − 1) is a different statement. It is proved in
[chanpak_near_floor_crossing_20261115](../chanpak_near_floor_crossing_20261115/README.md).

### §1. The ribbon mechanism (THM-5202)

- **Coalescence inequality (PROVED).** Glue A and B at a vertex v with A − v ⇒ B − v. Every Hamiltonian path
  decomposes uniquely as α β v γ ζ. Hence H = H(A)H(B) + Extra, with Extra ≥ 0 counted by two explicit bad junctions.
- **Floor powers have rank one (PROVED, exact).** In the minimal ribbon representations of THM-4880:
  - the span-2 unit 1.2 has characteristic polynomial x^13 (x − 3), and its fourth power is c·r^T;
  - the span-3 5-units 0.1.3 and 0.3.1 have rank-one cubes.

  Hence H(p u s) = L(p)·R(s) with integer L and R whenever a floor run u sits in the word. This was checked on 300
  of 300 random words, e.g. 2,130,237 = 753 × 2,829.
- **Prime avoidance.** A prime count must break every floor run in its bulk. Bulk runs also force factors 3^j or 5^j.
  This is the exact reason behind the linear growth of 3- and 5-adic valuations that the 7-21 session measured
  (notes/06 E3).
- **The premium (FINITE-EXACT on periodic words).** Every bulk-strong periodic block cheaper than 1.7700 per vertex is
  a bottleneck. The checks cover span 2 with period ≤ 10 (591 classes), span 3 with period ≤ 6 (408), and 1,926
  floor-plus-spoiler words. The cheapest bottleneck-free block, 0.2.1.2, has characteristic polynomial
  x^12 (x² − 9x − 8). So

      ρ* = ((9 + √113)/2)^(1/4) = 1.7700005 = 1.035104 · 5^(1/3).

### §2. What it does to NF2 (HYP-9742)

- **Census.** The span-3 DFS to order 14 reproduces the 7-21 window misses exactly: 637; then 1067, 1109, 1339. It adds
  seven misses at order 14: 1879, 1903, 1943, 1949, 1993, 2051, 2089. Every miss is prime or has a prime factor ≥ 13.
- **Smallest prime ratio.** The smallest prime count over f(k) rises 1.336 → 1.713 from k = 10 to k = 18. That is
  3.2% per vertex, against the predicted 3.51%.
- **Conjecture.** NF2(c) fails for every fixed c at all large primes. Extrapolation puts the failure near 1.5·10^8 for
  c = 3 and near 1.5·10^10 for c = 4.
- **Corrected form.** Primes live near the bottleneck-free floor: H < C·ρ*^n.
- **What survives.** Composites near the floor are products, and Chan–Pak is untouched. Primes are realized above the
  floor by additive families: recursive trees and carrier forms.

### §3. Density one for carrier forms (HYP-9722 update)

- **Reduction (PROVED).** Chebyshev's inequality reduces density one to a collision count.
- **Data.** The surcharge makes the counts Poisson: variance/mean is 1.03 at 10^5, 10^6 and 10^7 for the 3-cycle form.
  The same form with a = 0 (the pure hyperbola) has variance/mean 1.67 → 2.76, divisor clumping, and its misses
  barely fall.
- **Concurrent work, credited.** The 7-21 session (notes/06 E5) found the same Poisson law, and
  chanpak_near_floor_crossing_20261115/DENSITY.md computed the exact collision constant 131/121.

### §4. Crossing lemma (HYP-9732 update)

**REFUTED** by another session: an order-11 cut with H = 815, δ = −1 and E = 156; a two-shelf family makes δ
arbitrarily negative.

This session adds:
- a second form of the identity, δ(X) = Π(X) − s_{X^c} − e_X;
- the small-order picture: through order 8, δ = 3 occurs only at H = 5 (the alternating strong-4 core), and the
  minimum δ rises with H there. That rise does not persist.

### Reproduction

All scripts are in `scripts/` and share the `w28_` prefix. `w28_receipts.txt` has the console outputs.

```bash
cd 05-knowledge/results/chanpak_waves_20261113/scripts && python3 w28_bottleneck_factorization_check.py
```

```bash
cd 05-knowledge/results/chanpak_waves_20261113/scripts && python3 w28_bottleneck_free_charpoly.py
```

```bash
cd 05-knowledge/results/chanpak_waves_20261113/scripts && python3 w28_bottleneck_free_premium.py 3 6
```

```bash
cd 05-knowledge/results/chanpak_waves_20261113/scripts && gcc -O3 -o nfdfs w28_nearfloor_ribbon_dfs.c && R=../../../../04-computation/experiments/ribbon_spectrum_20261104 && ./nfdfs $R/span2_matrices.txt $R/span2_condensation.txt 2 18 4 3 nf_s2_18
```

### Verdicts

| claim | status |
|---|---|
| Coalescence inequality | PROVED |
| Rank-one floor powers; H = L·R through floor runs | PROVED (exact linear algebra) |
| Bottleneck-free premium ρ* = ((9 + √113)/2)^(1/4) | FINITE-EXACT on periodic words |
| NF2(c) fails at large primes | CONJECTURE (HYP-9742) |
| Density one ⇐ Poisson collision count | PROVED reduction; Poisson behaviour NUMERICAL |
| Crossing lemma | REFUTED (another session); small-order structure recorded |

## Wave 29: the growth bound — floor units compose by exact Hankel laws, and bottleneck-free words pay 1.7699 per vertex

Session mac-mini-2026-11-15-growthbound. Owner prompt: "Work toward proving chan-pak by finding more constructions and
analogies that apply to all sizes, come to understand the recursive patterns of tournament symmetry composition. try
many possible approaches to the proof. prove the bottleneck-free growth bound and attack the carrier collision count."

Exploratory throughout. Receipts are in `scripts/w29_receipts.txt`; prediction P102 is filed in HYP-9752.

### §1. The composition law of the floor units (THM-5212 (A))

The span-2 units A = 0.2 and B = 1.2 generate a language with an exact 3-dimensional representation in the Hankel
coordinates z = (H(W), H(WA), H(WB)):

    H(WAA) = 3 H(WA),                       H(WAB) = 4 H(WA) − 3 H(W),
    H(WBA) = (7 H(WB) + 4 H(WA) − 4 H(W))/3,  H(WBB) = (5/3)(H(WA) + H(WB) − H(W)).

- Compatible repetitions AA and BB have rank one, so they are exact products.
- The incompatible junctions AB and BA have characteristic polynomial x(x² − 9x − 8).
- So the bottleneck-free words are exactly the alternating ones, with growth ρ* = ((9 + √113)/2)^(1/4).

This is the recursive pattern behind THM-5202: compatible composition multiplies, and incompatible composition is
affine. The affine room is all a prime count can use.

### §2. The growth bound, proved in the span-3 five-unit language (THM-5212 (B))

The units 0.2, 1.2, 0.1.3, 0.3.1 and 1.2.3 generate a language with an exact 4-dimensional representation.
- Bottleneck-free words are exactly the paths of a 9-state automaton on image planes.
- Simple cycles have minimum growth ρ* at (0.2, 1.2).
- An exact rational certificate gives H ≥ C·1.7699^n for every bottleneck-free word: invariant sectors (with snapped
  exact boundary rays) and edge-expanding functionals.
- So inside this language the premium over the floor 5^(1/3) is a theorem.

Open: the other eight floor units, non-unit letters, and non-ribbon tournaments.

### §3. The carrier collision count (HYP-9752)

- **Rows are counts.** Setting Y = 1 in a carrier gives a row G(U) = cU + k p(U) + φ[U odd] of tournament counts.
- **Top-split reduction (proved).** A collision U ≠ U' splits at its highest differing bit into prefix w and suffixes
  b, b'. Then k d w + k (w mod 2) = c(b − b') + k(p(b) − p(b')) + φ(b₀ − b'₀), with d = rev b − rev b' odd. So the
  suffix pair determines the prefix.
- **Measured energy.** Ordered collisions per 2^L converge (exact to L = 24):
  - 0.2825 for 18U + 3p + 2[odd];
  - 0.0053 for 189U + 11p;
  - 1.267 for 2U + p.
- **The lemma.** The collisions with split level l and prefix w number at most C·2^l/(1 + w); the measured C is 14 for
  the 3-cycle row. This lemma gives linear energy, and Cauchy–Schwarz then gives positive lower density for the whole
  Hamiltonian-path spectrum.
- **Concurrent work, credited.** The 7-21 session measured one-hub energy ≈ 10.08·2^m and noted the same density
  consequence. chanpak_dyadic_collision_ports_20261115 proved the split identity and short-displacement bounds.

### §4. Constructions for all sizes

- **Concurrent.** chanpak_symmetry_interfaces_20261115 studies which data substitution interfaces pass and how
  symmetry composes recursively.
- **This wave.** It adds the floor-unit composition algebra of §1 and §2: exact Hankel laws whose compatible products
  are rank one and whose incompatible products carry a second eigenvalue.
- **Reading.** Constructions that work at every size must live on the incompatible side: carriers, recursive trees, and
  alternating unit chains.

### Reproduction

```bash
cd 05-knowledge/results/chanpak_waves_20261113/scripts && python3 w29_unit_language_rep.py 2 0.2,1.2
```

```bash
cd 05-knowledge/results/chanpak_waves_20261113/scripts && python3 w29_exact_growth_certificate.py w29_unitrep_w3_5units.pkl
```

```bash
cd 05-knowledge/results/chanpak_waves_20261113/scripts && python3 w29_row_collisions.py 22
```

### Verdicts

| claim | status |
|---|---|
| Floor-unit Hankel laws; alternating words are the only bottleneck-free span-2 unit words | PROVED |
| Five-unit span-3 growth bound H ≥ C·1.7699^n | PROVED (exact rational certificate) |
| Full span-3 and general growth bound | OPEN (HYP-9742) |
| Top-split reduction of row collisions | PROVED |
| Row collision lemma, and with it positive density of the spectrum | CONJECTURE (HYP-9752); implication PROVED |

## Wave 30: thirteen floor units, the twisted collision, and the local-limit program

Session mac-mini-2026-11-15-units13. Owner prompt: "Work toward a local limit theorem for chan-pak and try to prove
the row collision lemma, or extend the growth bound to all span-3 letters."

Exploratory throughout. Receipts are in `scripts/w30_receipts.txt`.

### §1. All 13 span-3 floor units (THM-5212 (D))

- **Exact representation (PROVED).** Hankel matrices of exact integer counts H(prefix·suffix) give an exact
  8-dimensional representation of the full 13-unit floor language. It reproduces 300 of 300 random words, and the
  same method recovers the five-unit rank 4.
- **Spectra.** Every unit has characteristic polynomial x^7(x − λ) with λ ∈ {3, 5}. Of the 169 pairs, 23 are rank one,
  73 compatible and 96 incompatible.
- **What failed.** Rank-3 images persist, so the wave-29 plane automaton grows past 3,000 states. The exact
  certificate does not extend this way.
- **Why, and what to change.** The five-unit language closed up because every recurrent product had rank exactly 2.
  With 13 units, rank-3 products persist and their planes rotate continuously. A proof needs a certificate that does
  not enumerate image planes: a cone family on the 8-dimensional reachable cone, together with a separate rank-two
  invariant.
- **Evidence (NUMERICAL).** Over all periodic words of at most five units, the cheapest bottleneck-free growth is still
  exactly ρ* = 1.7700005. It is attained by 0.2·1.2 and by the 12-letter cycle 0.1.3·0.1.123·0.1.23·0.3.1, with
  eigenvalue 9.815^3. The full-span-3 bound stays HYP-9742.

### §2. The row collision lemma as a twisted collision (HYP-9752 update)

Peeling the top suffix bit turns a collision at split level l and prefix w into

    A_w(x) − B_w(y) = (bounded),
    A_w(x) = cx + κp(x) − κ(2w+1) rev(x),   B_w(y) = cy + κp(y) − 2κw rev(y),

for (l − 1)-bit x and y. The lemma N(l, w) ≤ C 2^l/(1 + w) is therefore a statement about how often two twisted
digital value sets meet: about 2^l/w times, as for random sets of their size spread over a range of length κw·2^l.
The trivial window bound is weaker by a factor 2^l. A proof needs equidistribution of c·rev(ρ) + κp(rev(ρ)) under the
twist ρ ↦ −κwρ. **Not proved.**

### §3. The local-limit program

**Target.** A family 𝔉 of tournaments with representation counts r_𝔉(N) = #{T ∈ 𝔉 : H(T) = N} and an explicit
main term μ(N)κ(N) (a smooth density times local factors) such that

    r_𝔉(N) = μ(N)κ(N)(1 + o(1))   uniformly in the bulk.

A pointwise version with min_N μκ → ∞ gives every large N. What the waves establish:

| level | statement | status |
|---|---|---|
| local factors | κ(N) for carrier forms, explicit and positive | proved elsewhere (dyadic mixing, voltage profile) |
| L¹ (positive density) | E = O(Σ μ) for one carrier row | ⟸ HYP-9752 (Cauchy–Schwarz, proved) |
| L² (density one) | Σ (r − μκ)² = O(Σ μκ) for one form | ⟸ Poisson collision count (Chebyshev, proved); numerical variance/mean 1.03 |
| pointwise | min r ≥ 1 | needs μ → ∞ faster than log N |

**Where a local limit theorem cannot hold.** Near the Moon–Busch floor it fails, because counts factor through floor
runs (THM-5202) and are p-adically biased. The bottleneck-free growth (THM-5212) measures how far above the floor the
additive regime starts.

**Where it should hold.** In families with the "incompatible" composition of §1 everywhere and multiplicity growing
like a power of N: the union of carrier forms over growing cores, recursive trees, and many-hub words. The carrier
family is binary (μ ~ log N), so for one form only the Poisson L² version is available.

### Verdicts

| claim | status |
|---|---|
| Exact 8-dimensional representation of the 13 floor units | PROVED |
| Growth bound for 13 units | CONJECTURE; numerically ρ* (periods ≤ 5 units); certificate does not extend |
| Row collision lemma | OPEN; twisted reformulation |
| Local limit theorem | program stated; L¹ and L² reductions proved; pointwise needs super-logarithmic multiplicity |

## Wave 31: hosts for a prescribed count — no linear tape, and long prefixes never collide

Session mac-mini-2026-11-15-hosts. Owner prompt: "Work to choose a host that hits the prescribed absolute count. build
a cone certificate for the 13 units and attack the twisted-collision lemma. Aim for generating creative new methods to
prove chan-pak, not just extend the range for which we have verified it."

Exploratory throughout. Receipts are in `scripts/w31_receipts.txt`.

### §1. A host for a prescribed count: the linear-tape idea, tested (HYP-9762)

**The method proposed.** The only thing between a carrier form and an exact arithmetic progression is the quadratic
surcharge a·p(X) of its tape. Suppose some pair of digit gadgets made the host count exactly affine,
H = D·X + E. Then choosing the host for a prescribed N would be division with remainder, and finitely many
progressions covering the odd residues would prove Chan–Pak.

**What was tested.** A direct search over one-apex hosts: gadgets of up to 3 vertices, all marks, a 3-cycle core.
- 442,598 hosts with no affine one.
- 563,500 hosts with no affine-plus-digit-sum one.

**Why it fails.**
- By the cut formula H(T + x) = H(T) + 2W(X), the tape dependence is a cut function of the nonnegative path-weight
  matrix Π.
- Its quadratic part, −Σ_{X×X} Π, cannot cancel by positivity.
- Intermediate gadgets multiply long-range weights, so every pair of digits interacts.

**What a circumvention would need.** A counter v ↦ 2v + b, i.e. a unipotent element in some port algebra. Floor units
have single nonzero eigenvalues, and span-2 blocks have eigenvalue 1 at most simple. But the Hankel laws do subtract
(H(WAB) = 4H(WA) − 3H(W)). So the place to look is multi-port algebras with negative coefficients, not single tapes.

### §2. The twisted collision: long prefixes never collide (HYP-9752 update, PROVED)

- **Statement.** For prefixes w ≥ ((2c + 3κ)/(2κ))·2^l + O(1), the reversed difference is forced to vanish, and the
  remaining equation has no solution. So N(l, w) = 0 there.
- **Meaning.** Every collision of a carrier row is balanced: the two inputs share at most about half their leading
  bits, and the row is injective on blocks of length about √U.
- **Check.** No collision has w ≥ 8·2^l at L = 20.
- **What is left.** The lemma is now reduced to the balanced range w < W₀(l). There the bound N(l, w) ≤ C 2^l/(1 + w)
  is still open.

### §3. The 13-unit cone certificate: not built

- **Why the wave-29 certificate does not extend.** It tracked image planes, and with 13 units these rotate
  continuously.
- **The planned replacement, not carried out.** A syllable certificate: cut words into k-unit syllables that are not
  rank one, and find a piecewise-linear antinorm on the 8-dimensional reachable cone expanded by every such syllable.
- **The difficulty.** Compatible runs that are not yet rank one have floor-like growth for up to 3 units. So
  syllables must be long enough to exclude them, which makes the LP large.

### Verdicts

| claim | status |
|---|---|
| No affine or affine-plus-digit-sum one-apex tape (gadgets ≤ 3 vertices) | FINITE-EXACT; general form CONJECTURE (HYP-9762) |
| Long-prefix collisions vanish (N(l, w) = 0 for w ≥ W₀(l)) | PROVED |
| Row collision lemma in the balanced range | OPEN |
| 13-unit growth certificate | not built; plan stated |

## Wave 32: the decimal triple, sparse linear tapes, and no counter in span 3

Session mac-mini-2026-11-15-counters. Owner prompt: "Merge in prior work on the triple {1, 1001, 1000001} and how they
relate to 1000000001 and search multi-port algebras for a unipotent counter and build the 13-unit syllable
certificate, think prior work on 62 and arbitrary span-3 words."

Exploratory throughout. Receipts are in `scripts/w32_receipts.txt`.

### §1. The triple and 1000000001

**Decimal reading.**
- *Prior work.* Session 82 ([DECIMAL_TRIPLE](../session82_decimal_aliquot_square_sum_20261028/DECIMAL_TRIPLE.md))
  characterized the triple A = t² + 1, B = t + 1, C = t(t + 1) by A + B = C + 2 and B + C = B², at t = 1000.
- *Continuation.* The requested list continues the same pattern, with t = 10³:

      1, 1001 = t + 1, 1000001 = t² + 1, 1000000001 = t³ + 1 = (t + 1)(t² − t + 1) = 1001 · 999001.

- *For Chan–Pak.* Counts are multiplicative over ordinal sums, so 10^9 + 1 is realized by stacking witnesses for
  1001 and 999001.
- *The factor 7.* 1001 = 7·11·13 has a factor of 7, so it must use the seven lane: 1001 = 11 · 91, with 91 = 7·13
  realized directly.
- *A direct host.* The concurrent placement packet
  ([chanpak_cone_collision_placement_20261115](../chanpak_cone_collision_placement_20261115/README.md)) builds a
  31-vertex host. It is the 3-cycle carrier at addresses (804, 103594), and F(804, 103594) = 1,000,000,001
  (rechecked here).

**Binary reading.** Read as tape addresses, 1, 1001, 1000001, 1000000001 are 2^{3k} + 1. They are exactly where the
surcharge is linear: p(2^k + 1) = 2^(k−1) = (X − 1)/2. More generally, the padding law
p(2^k u + 1) = p(u) + (2^(k−1) − 1)u + (u mod 2) (k ≥ 2) confines the surcharge to the short payload u. So the linear
tapes that wave 31 ruled out exist only on geometrically sparse address families.

### §2. Multi-port counters in the full span-3 representation

In the 62-dimensional representation, checked exactly modulo primes:
- of the 4,680 blocks of length ≤ 4, 233 have eigenvalue 1, and none has a unipotent (Jordan) part;
- 141 of them are strong blocks, where eigenvalue 1 is a semisimple constant offset;
- only 8 fixed vectors are shared by distinct strong blocks, and every such pair is an H-twin with zero digit weight.

So no affine counter exists among span-3 blocks of length ≤ 4 (HYP-9762 update).

### §3. The 13-unit certificate and arbitrary span-3 words

- **Already done.** The 13-unit growth certificate was built by the concurrent packet: cones indexed by three-unit
  suffixes, with 12 promoted triples, give growth ≥ 1.7699 across all thirteen units. It is credited here and not
  rebuilt.
- **Still open.** Arbitrary span-3 letters (the full 62-dimensional language). Periodic blocks up to period 6 are
  covered (THM-5202 (D)). A suffix-cone certificate there would index cones by letter suffixes in 62 dimensions.

### Verdicts

| claim | status |
|---|---|
| Triple ↔ 10^9 + 1 relations; host F(804, 103594) = 10^9 + 1 | PROVED arithmetic; host credited to the concurrent packet |
| Padding law; linear surcharge exactly on 2^k + 1 addresses | PROVED |
| No unipotent block, and no affine counter, among span-3 blocks of length ≤ 4 | FINITE-EXACT |
| 13-unit certificate | done elsewhere (credited) |
| Arbitrary span-3 growth bound | OPEN |

## Wave 33: a survey of tangential prior work, and what it points to

Session mac-mini-2026-11-15-survey. Owner prompt: "spend this session surveying past work that may relate
tangentially to chan-pak for inspiration then continuing creatively toward its proof."

Three read-only surveys were run, each over 01-canon, 05-knowledge and 00-navigation:
- proved "every value" analogues;
- averaging and local-limit proofs on any problem;
- exact composition laws and obstructions.

Their findings, condensed, follow with paths, then one new test and the program they suggest.

### §1. Proved "every value" analogues and why they do not transfer

- **The solved analogues.**
  - BCS threshold graphs: every odd k is a dominating-set count. The proof is binary numeration with the cone move
    d(K1 ∨ G) = 2^|G| + d(G), which adds a power of two whatever G is.
  - Even graphs (THM-5040, THM-5055) and odd graphs (THM-5056, THM-5042): every odd value except {3, 7, 13}, resp.
    {1, 5, 37}. Three ingredients: padded cones write a leading "10"; offset moves +3·2^(n−1) (order +2) and
    +7·2^(n−1) (order +3) carry an offset across all larger powers of two; and a seed lemma, proved by strong
    induction on the offset.
  - Spanning trees (THM-4986): τ(C_N) = N.
- **The tournament analogues that exist.**
  - Rays c₀ + M·2^t with power-of-two response M (THM-5126, THM-5130, THM-5145) cover every odd N within about
    [−791, 647] of a power of two.
  - Endpoint counts are universal (THM-5001, THM-5113), because they are linear in the in-set.
- **What is missing.**
  - The tournament "cone" (near-join h + 2^k·B with B ≥ h) only reaches (2^K, 3·2^(K−1)] and does not chain
    (THM-5109).
  - So the precise missing pieces are a move writing a leading "11" and an inductive seed lemma. The type conjecture
    HYP-9675 is the closest statement.

### §2. Exact laws never used as generators (from the composition survey)

1. **Hyperbolic two-port kernels** (THM-5133): N₂₂H = LL′ − |det G|·hh′, an indefinite form. It lets a count be a
   small difference of two large, correlated core counts. The cyclic join used so far is definite, which is why that
   route needed band fullness (HYP-9692).
2. **Prime quotients of order ≥ 5 with three or more free modules** (THM-3121).
3. **Recursive binary cyclic joins**, with both inputs of C3[X, B, TT₁] built recursively. Analysed only mod 15.
4. **Iterated vertex substitution** with the negative k(E − O) term (symmetry interfaces, PORTS2).
5. **Exact switchboards.** The one-visit swaps TT₃ ↔ C₃ change H by exactly 2× the slot weight but were never
   combined across slots. THM-5124 forbids only independent parallel switches.
6. **The coalescence Extra term** (THM-5202 (A)), never parametrised.

### §3. Averaging routes and the one that could be pointwise (from the local-limit survey)

- **What the levels give.** L¹ energy gives positive density (HYP-9752), and L² gives density one (Chebyshev).
  Pointwise needs the minor-arc integral below one main-term coefficient. Sup bounds per member cannot do this
  (session205 and session193 cautionary notes).
- **The tilted all-tournament ensemble** has a proved reduction from a tilted local limit law to every large odd h
  (THM-5127, Proposition 5.1). The missing piece is an integral bound uniform in the tilt.
- **Twisted transfer operators for ribbon words** (HYP-9474). The multiplicity is about X^1.6 > X, so a sup-norm
  minor-arc bound would already be pointwise. The needed per-letter contraction is < e^(−λ) ≈ 0.59 away from the
  arcs a/(2·3^i·5^j); measured values are 0.7–0.8.

### §4. A test this wave: collective sign freedom fails

**Idea.** A single tournament cannot keep both unit moves H ↦ 2H ± 1 (sign freezing). But many recursive-tree
tournaments share each count. If for every large M some tree with count M admitted a +1 leaf and some admitted a −1
leaf, binary descent would give every odd N.

**Test** (`scripts/w33_recursive_tree_unit_offsets.c`, all trees):

| order | counts M | both ±1 | only +1 | only −1 | neither |
|---:|---:|---:|---:|---:|---:|
| 6 | 16 | 6 | 4 | 4 | 2 |
| 7 | 55 | 9 | 9 | 11 | 26 |
| 8 | 192 | 13 | 30 | 27 | 122 |

Unit offsets are rare even collectively at a fixed order. Binary descent inside recursive trees is dead.

### §5. The creative link: the major arcs are the floor bottlenecks

The twisted-operator route fails exactly at the arcs a/(2·3^i·5^j). By THM-5202 those denominators come from floor
runs: every bulk floor run forces an exact factor 3^j or 5^j and a factorization H = L·R. This gives a new program:

- **(P1) Bottleneck-free twisted operators.**
  - *Restriction.* Restrict the ribbon family to bottleneck-free words (THM-5212, and the 13-unit suffix cones of
    chanpak_cone_collision_placement_20261115). This removes the source of the 3^i·5^j arcs.
  - *Trade-off.* The minimum growth rises to ρ* = 1.7700, so the multiplicity exponent drops, and the contraction
    target changes accordingly.
  - *Test.* Measure the twisted contraction on the bottleneck-free automaton at θ = 1/3, 1/5, 1/8. If it beats the new
    target uniformly, a pointwise theorem follows.
- **(P2) A leading-"11" move from a hyperbolic kernel.** In N₂₂H = LL′ − |det G|hh′, fix one core. The count becomes
  a linear form αL − βh in the other core's ports, with a negative coefficient. Values in the "11" band then come
  from cores whose port ratio L/h lies near β/α. Test: a census of (L, h) for cores up to order 9, and the band of
  αL − βh values each fixed core reaches.
- **(P3) A seed lemma by strong induction on offsets**, modelled on THM-5055. Tournament rays already carry offsets
  across powers of two (THM-5126). What is needed is an operation turning a realized count c into an offset with
  power-of-two response. P2 or the switchboard (§2.5) are the candidates.

### Verdicts

| claim | status |
|---|---|
| Survey map (§1–§3) | collected with paths; claims are the cited files' own |
| Collective sign freedom (both unit leaves for every count) | REFUTED at orders 6–8 (FINITE-EXACT) |
| Major arcs at 3^i·5^j = floor bottlenecks | ARGUMENT (THM-5202 mechanism); P1 untested |
| P1–P3 | proposed program |

## Wave 34: norm forms for prime targets, and what bottleneck-free twisting does

Session mac-mini-2026-11-15-norms. Owner prompt: "work both the bottleneck-free twisting and the hyperbolic kernel
move, think about turning some prime targets into quadratic norm equations, and the number theoretic implications of
tournaments."

Exploratory throughout. Scripts are `scripts/w34_*`.

### §1. Cyclic joins are Z[√−2] arithmetic (THM-5222, PROVED)

For A ⇒ B ⇒ c ⇒ A,

    H = (h_A + P_A)(h_B + P_B) + 2 h_A h_B = Re(β_A·conj(β_B)),   β = (h + P₂) + h√−2,

where P₂ counts ordered covers by two disjoint nonempty paths. The formula matches on 200 of 200 random core pairs.

**Diagonal joins.** These are norms x² + 2y², always ≡ 3 (mod 8).
- A prime target p ≡ 3 (mod 8) becomes the two-port equation H(A) = y, P₂(A) = x − y for the unique representation
  p = x² + 2y².
- Primes ≡ 1, 5, 7 (mod 8) are unreachable by diagonal joins.
- Inert primes (≡ 5, 7 mod 8) divide diagonal counts to even powers.

**General joins.** With B fixed, a prime target is a point of the pair spectrum {(H(A), P₂(A))} on a line. So
Chan–Pak for primes reduces to a two-dimensional covering statement for the pair spectrum.

### §2. The hyperbolic kernel: an indefinite norm with a Pell supply (ARGUMENT)

THM-5133's first hyperbolic kernel at order 5 reads H = 3Π₂(A)Π₂(B) − H(A)H(B). Take B a copy of A. Then

    H = 3Π₂² − H(A)² = −Norm(H(A) + Π₂√3),

a norm from the real quadratic order Z[√3].
- *Which primes.* Its fundamental unit 2 + √3 has norm +1, and the class number is 1. So x² − 3y² = −p is solvable
  exactly for primes p ≡ 11 (mod 12).
- *Infinitely many representations.* Each solution generates an infinite orbit under multiplication by (2 + √3)^k.
- *What that buys.* Compared with the definite join of §1, a prime target ≡ 11 (mod 12) has infinitely many admissible
  port pairs (H(A), Π₂(A)), one at each geometric scale.
- *What is needed.* Some tournament A must land exactly on one of these pairs: a two-dimensional local statement along
  a Pell orbit instead of at one point.

The exact definition of Π₂ and of the kernel are THM-5133's. The diagonal substitution has to be rechecked against
that kernel before this becomes a theorem.

### §3. Bottleneck-free twisting: what causes the 3-adic arc (NUMERICAL)

Random span-2 words of length n, 40,000 samples per length; P(3 | H) (uniform 1/3):

| class | n = 16 | n = 28 |
|---|---:|---:|
| all words | 0.606 | 0.769 |
| strong words | 0.463 | 0.522 |
| strong, with an exact rank-one factor (≤ 10 letters) | 0.702 | 0.711 |
| strong, bottleneck-free (≤ 10 letters) | 0.439 | 0.470 |
| strong, with a factor of rank ≤ 1 mod 3 | 0.657 | 0.619 |
| strong, no rank-≤1-mod-3 factor | 0.418 | 0.452 |

P(5 | H) behaves the same way: 0.422 / 0.574 for all words, 0.231 / 0.281 for bottleneck-free strong words.

**Reading.** The 3- and 5-adic major arcs have three layers:
- non-strongness (products of small components);
- exact bottlenecks (THM-5202);
- *p-adic bottlenecks*: factors whose transfer collapses to rank ≤ 1 modulo p.

Removing all three still leaves a bias that grows with n. That remainder is a genuine local factor, and every class
keeps positive mass. So it belongs to the main term of a local limit theorem with explicit local densities, as in
THM-5162's singular series, not to an obstruction. Restricting to bottleneck-free words does not remove the arc; the
open problem is still the minor-arc contraction.

### §4. Number-theoretic reading

Tournaments carry quadratic arithmetic: cyclic joins pair port vectors as elements of Z[√−2] (definite; finitely many
representations). The hyperbolic kernel pairs them in Z[√3] (indefinite; Pell orbits). Prime targets become norm
equations whose solvability is governed by splitting:
- p ≡ 3 (mod 8) for the diagonal join;
- p ≡ 11 (mod 12) for the diagonal hyperbolic kernel.

Different kernels (different discriminants) cover different classes of primes. A natural program is a finite set of
kernels whose norm forms jointly represent every prime class, each with a two-port targeting problem. The indefinite
ones are the promising ones, because of their Pell supply.

### Verdicts

| claim | status |
|---|---|
| Cyclic join = Re(β_A·conj(β_B)) in Z[√−2]; diagonal joins ≡ 3 mod 8 and only reach primes ≡ 3 mod 8 | PROVED |
| Diagonal hyperbolic kernel = −Norm from Z[√3], with Pell orbits of targets for p ≡ 11 mod 12 | the norm identity holds; the Pell supply is REFUTED in wave 35 (positivity confinement) |
| 3-adic arc layers: non-strong, exact bottleneck, mod-p bottleneck, residual local factor | NUMERICAL |
| Bottleneck-free twisting removes the 3-adic arc | REFUTED as stated (residual bias persists) |

## Wave 35: positivity confinement — no kernel gives cancellation, so prime targets have finitely many port pairs

Session mac-mini-2026-11-15-positivity. Owner prompt: "spend another similar long session toward proving chan-pak."

### §1. The correction (PROVED)

THM-5133's visit-matrix form is

    H(K[A, B]) = Σ N_ij·P_i(A)·P_j(B),

where N_ij counts walks of the kernel visiting the ports i and j times, and P_i counts covers by i nonempty paths.
Every coefficient and every variable is nonnegative.
- *Where the sign comes from.* The elliptic / parabolic / hyperbolic trichotomy refers to the same form written in
  (H, Π₂), with Π₂ = P₂ + 2H. The shift makes det G negative; it cannot create cancellation.
- *The realizable cone.* Π₂ ≥ 2H always, because (path, ∅) and (∅, path) are among the covers. Checked: the minimum of
  Π₂/H over all 7,412 classes of order ≤ 8 is exactly 2. For a symmetric kernel the isotropic direction of the
  diagonal form has Π₂/H < 2, because N₁₁N₂₂ ≥ 0 gives |det G| < N₁₂². So every isotropic direction lies outside the
  realizable cone.
- *Consequence.* Each term of the sum is at most H, so a prime target p has only finitely many realizable port
  representations, all with ports ≤ p. Wave 34's "Pell supply" from the order-5 kernel 3Π₂Π₂′ − HH′ does not exist on
  tournaments: on Π₂ ≥ 2H that form is ≥ 11HH′.

This is the kernel-level version of HYP-9762: path counts never cancel.

### §2. What remains true: prime targets are linear equations on the pair spectrum

The pair (H, Π₂) is multiplicative under ordered joins (THM-5133). Every two-port kernel is a positive bilinear form
in the port vectors. With one core fixed, a prime target becomes a linear equation αH(A) + βP₂(A) = p on the pair
spectrum V = {(H(A), P₂(A))}. THM-5222's cyclic join is the case α = 3, β = 1 with a one-vertex B.

So for some fixed small kernel, Chan–Pak for primes follows from **pair-spectrum fullness**: V meets every line
αx + βy = p in its realizable cone.
- *How strong.* This is a two-dimensional covering statement. It needs more tournaments than the square of the value
  range (σ > 2 in the 7-21 notes' terms). Counting classes (~2^(n²/2)/n!) against range² (~(n!/2^n)²), that first
  happens near order 12.
- *Where it should hold.* Not near the floor, where counts factor (THM-5202). It should hold in the bulk of the
  all-tournament ensemble, which is exactly where the tilted local-limit reduction of THM-5127 lives.

### Verdicts

| claim | status |
|---|---|
| Two-port kernels are positive bilinear forms; isotropic directions lie outside Π₂ ≥ 2H; finitely many representations per prime | PROVED |
| Pell supply from hyperbolic kernels | REFUTED |
| Prime targets ⇐ pair-spectrum fullness for a fixed kernel | PROVED reduction; fullness OPEN (σ > 2 regime) |

## Wave 36: the pair spectrum is sparse in two dimensions but its lines are full

Session mac-mini-2026-11-15-pairs. Owner prompt: "test whether the pair set (H, P₂) covers every line at the largest
orders and other strategies to prove chan-pak."

Exact (H, P₂) for every tournament class of orders 8, 9 and 10 (`scripts/w36_pair_spectrum.c`; order 10, all
9,733,056 classes, in about 3 minutes).

### §1. Two-dimensional fullness: no

| order | classes | distinct pairs | distinct H | mean fill of the P₂-fibres of fixed H |
|---:|---:|---:|---:|---:|
| 8 | 6,880 | 1,727 | 320 | 0.038 |
| 9 | 191,536 | 50,173 | 1,520 | 0.038 |
| 10 | 9,733,056 | 1,550,357 | 7,663 | 0.032 |

The fill is the fraction of even P₂ values between a fibre's minimum and maximum that occur, over fibres with at least
20 slots. It does not grow with the order, so wave 35's pair-spectrum fullness is out of reach of computation. That
matches its σ > 2 threshold near order 12.

### §2. Line fullness: essentially yes

Fraction of the odd numbers in the central 10–90% window of each line α·H + β·P₂ that are attained:

| line | order 8 | order 9 | order 10 |
|---|---:|---:|---:|
| 3H + P₂ | 0.417 | 0.332 | 0.9993 |
| 5H + 3P₂ | 0.199 | 0.785 | 0.9992 |
| 7H + 2P₂ | 0.252 | 0.876 | 0.9998 |
| H + P₂ | 0.309 | 0.968 | 0.3333 |

**The one third is exact.** THM-4756 gives Π₂ ≡ (−1)ⁿH (mod 3), with Π₂ = P₂ + 2H.
- At even order P₂ ≡ −H, so H + P₂ ≡ 0 (mod 3).
- At odd order P₂ ≡ 0, so 3H + P₂ ≡ 0 (mod 3).

These are the only failures seen. On lines free of that congruence, the order-10 projection fills its central window
to within 0.1%.

### §3. What this changes

Prime targets need much less than wave 35 asked for. With one core fixed, a two-port kernel makes the count a line
α·H(A) + β·P₂(A). What is needed is **line fullness**: the projection of the all-tournament pair set onto the kernel's
line covers the odd numbers in the classes the mod-3 law allows.
- *Status.* Numerically it holds at order 10. It is a one-dimensional covering statement of the same kind as Chan–Pak
  itself, so it is not easier in principle.
- *Why it is a better target.* Its statistic is linear in two path-cover counts of a random tournament, and its only
  local obstruction is the explicit mod-3 law. That makes it the natural input for the tilted local-limit reduction
  (THM-5127), applied to αH + βP₂ instead of H.

### §4. Other strategies, ranked by what the last nine waves left standing

1. **A local limit theorem for linear statistics αH + βP₂ of random tournaments.** Its local factors are explicit (mod
   3 from THM-4756; the 2-adic laws). Line fullness at order 10 is the numerical evidence, and the THM-5127 reduction
   transfers.
2. **Seed lemma plus offset moves** (dominating-set analogy, wave 33). It needs a move writing a leading "11"; none is
   known.
3. **Collision energy for carrier rows** (HYP-9752): long prefixes are settled, the balanced range is open. Proving it
   would give positive density, not every N.

Ruled out by exact results: unit-offset descent (wave 26, wave 33), linear tapes and counters (HYP-9762), kernel
cancellation and Pell supply (wave 35), and near-floor coverage at large primes (HYP-9742).

### Verdicts

| claim | status |
|---|---|
| Pair spectrum not two-dimensionally full at orders 8–10 (fibre fill 3–4%) | FINITE-EXACT |
| Lines αH + βP₂ fill their central windows at order 10, except exact one-third lines | FINITE-EXACT |
| The one-third lines are explained by Π₂ ≡ (−1)ⁿH (mod 3) | PROVED (from THM-4756) |
| Prime targets ⇐ line fullness for one kernel | PROVED reduction; line fullness OPEN for all orders |

## Wave 37: fibre swaps — exact, H-preserving repairs in steps of 24

Session mac-mini-2026-11-15-repairs. Owner prompt: exact repairs that help place the baseline from the target, after
two independent apexes place a count within O(N^0.8963) of the target. The six-vertex intercept hosts' repair steps
are too coarse, and padding fails: one singleton forces a factor 3, and two leave holes.

### §1. The repair move

A reader outputs r·H(A) + M·Q(A) from a source A, with Q = P₂ + 2H (chanpak_pair_line_census_20261115). Replace A by
any A′ with the same order and **the same H**. The output changes by exactly M·(P₂(A′) − P₂(A)). Everything that reads
only H(A) is unchanged.

### §2. The step is exactly 24 at even order (PROVED from known laws; FINITE-EXACT at order 10)

At even order, Π₂ = Q satisfies:
- Π₂ ≡ H (mod 3), by THM-4756;
- Π₂ ≡ 0 (mod 8), by THM-5172 (B).

So Π₂ mod 24 is fixed by H mod 3, and

    P₂ ≡ 2H + 12 (mod 24).

The order-10 census confirms it on all 9,733,056 classes: P₂ mod 24 is 14, 18, 22, 2, 6, 10 for H ≡ 1, 3, 5, 7, 9, 11
(mod 12). The smallest gap inside every H-fibre is 24 (1,079,310 occurrences); every gap is a multiple of 24. At odd
order the same argument gives steps of 12.

This is the exact form of the rule that two modules with equal H give host counts differing by a multiple of 4: within
one order, the difference is a multiple of 24M.

### §3. Fibre swaps fill windows when combined (NUMERICAL, order 10)

Measured on the largest H-fibres, with H ≈ 5,000–6,000, about 400 P₂ values each, and a range of about 1,000 units of
24:

| sources swapped | fill of the central half of the repair range, in units of 24 |
|---|---:|
| one fibre | 0.39–0.46 |
| two independent fibres | 0.9960–1.0000 |
| three independent fibres | 1.0000 (all three triples tested) |

So three sources, each read once (additively, as the two-apex placement combines its parts), give an **exact repair
set**: every multiple of 24M across a window of about 3,000·24M around the baseline. No source's H changes.

### §4. How this would be used

1. **Place.** Put the baseline within the window of the target. At order-10 sources the window is about 7·10⁴·M. It
   grows with the source order, since fibres widen; that growth is unverified.
2. **Steer.** Fix the residue modulo 24M with the existing residue tools (dyadic laws, mod-3 law).
3. **Repair.** Swap each of three sources within its H-fibre to remove the remaining multiple of 24M exactly.

**What is not shown.**
- Fibre widths at the source orders needed for N^0.8963 windows (orders well above 10).
- That the physical construction reads three sources additively with fixed slopes M.

Neither is verified here.

### Verdicts

| claim | status |
|---|---|
| H-preserving source swap changes a reader output by M·ΔP₂ | PROVED |
| P₂ ≡ 2H + 12 (mod 24) at even order, so swap steps are multiples of 24 | PROVED from THM-4756 and THM-5172 (B); FINITE-EXACT at order 10 |
| Three independent fibre swaps give a gap-free repair set in units of 24 (order 10) | NUMERICAL |
| Use inside the two-apex placement | proposed; scaling and additivity unverified. **Wave 38: modules inside a two-apex source are read multilinearly, not additively, and a module swap is a coarse move** |

## Wave 38: fibre widths, two-apex coupling, and a synthesis with Collatz and the arithmetic of 2 and 3

Session mac-mini-2026-11-15-synthesis. The owner's prompt asked for four things:
- measure fibre widths at higher orders;
- test additivity in two-apex hosts;
- merge the earlier work on multiplication tables mod odd numbers (the back and forth between 3x+1 and 3x−1 edges);
- reach a grand synthesis with other open problems, including how the metagraph and the norm equations of ordered
  pairs of disjoint paths (x² + 2y², diagonal joins) fit with mod 8 in Collatz.

The synthesis is in [GRAND_SYNTHESIS.md](GRAND_SYNTHESIS.md). It is a companion to the concurrent
[mod-8 dictionary](../chanpak_mod8_synthesis_span3_20261115/SYNTHESIS.md) of klein-2026-11-15-synthesis, which
answers the same prompt and is credited throughout.

**Credit correction to wave 37.** The fibre step law (P₂ ≡ 2H + 12 mod 24 at even order, steps of 12 at odd order)
was already proved in [PAIR_LINE_STRUCTURE.md](../chanpak_pair_line_census_20261115/PAIR_LINE_STRUCTURE.md) §4 via
the Grinberg–Stanley odd-cycle expansion. That section also gives Q ≡ 4·HC (mod 8) at odd order and
Q ≡ 8·c_(n−1) (mod 16) at even order. Wave 37 re-derived the law from THM-4756 and THM-5172 (B) without citing it.

### §1. Fibre widths at higher orders

**Exact, orders 4–10** (`w38_fibre_widths.py` over the complete class lists):
- At odd order P₂ ≡ 6 (mod 12) for every class, while at even order P₂ ≡ 2H + 12 (mod 24).
- The mean of P₂/H over classes is 2(n − 1) to within 0.5.
- In the table, "central" means H in the middle half of the class distribution.

| order | step | central fibres with ≥ 2 values | median width / H | median steps | median fill | widest fibre (H, values, steps, width/H) |
|---:|---:|---:|---:|---:|---:|---|
| 7 | 12 | 27 | 0.92 | 5 | 0.33 | 123, 4, 17, 1.66 |
| 8 | 24 | 119 | 1.59 | 19 | 0.35 | 153, 10, 37, 5.80 |
| 9 | 12 | 492 | 2.34 | 256 | 0.20 | 1479, 55, 453, 3.68 |
| 10 | 24 | 2312 | 3.03 | 829 | 0.39 | 3297, 334, 1342, 9.77 |

**Sampled, orders 10–16** (`w38_sample_pairs.c`, `w38_sampled_fibres.py`): 10⁶, 10⁶, 5·10⁵, 2·10⁵, 10⁵, 5·10⁴ and
2·10⁴ uniform labelled tournaments, restricted to H within ±1% of the median H. The order-10 sample matches the exact
classes (sd 0.435 vs 0.432).

| order | median H | sd of P₂/H | q01–q99 | sampled range |
|---:|---:|---:|---:|---:|
| 10 | 7,077 | 0.435 | 2.08 | 3.89 |
| 11 | 38,855 | 0.379 | 1.88 | 3.93 |
| 12 | 232,801 | 0.333 | 1.56 | 3.10 |
| 13 | 1,510,685 | 0.303 | 1.42 | 2.26 |
| 14 | 10,551,169 | 0.282 | 1.27 | 2.07 |
| 15 | 79,239,983 | 0.249 | 1.19 | 1.63 |
| 16 | 630,698,211 | 0.229 | 1.03 | 1.49 |

**Reading (NUMERICAL, with a heuristic extrapolation).**
- *Relative widths shrink.* The relative spread of P₂/H shrinks like about n^(−1.4) between orders 10 and 16.
- *Absolute widths grow.* A fibre at H spans at least about 1.5·H in absolute terms at order 16, a lower bound from
  the sampled range of P₂/H; it is 3.0·H exactly at order 10. That is at least H/16 steps of 24, so absolute widths
  grow superexponentially along with H.
- *Fill should reach 1.* Classes per fibre grow faster still: from about 1.5 per step at order 10 to an estimated tens
  per step at order 12 (2^(n(n−1)/2)/n! classes over about H_max/2 fibres). So central fibres should become gap-free
  from about order 12. Fill above order 10 is not measured.

### §2. Additivity in two-apex hosts: exact, multilinear, coupled and coarse

Inflate two spine positions of a two-apex source (TWO_APEX_COARSE_PLACEMENT.md) by modules X₁, X₂. Every count is then
multilinear in the modules' port vectors (P₁, P₂, P₃, …), and fibre swaps (equal H) are additive iff the mixed
difference D = F(X₁′, X₂′) − F(X₁′, X₂) − F(X₁, X₂′) + F(X₁, X₂) vanishes, with F = H + Q
(`w38_two_apex_additivity.py`).
- **Never additive.** D ≠ 0 in 240 of 240 random sources (m = 6–8, order-5 and order-6 fibre pairs). The median
  |D|/min|ΔF| is 0.087, and 80% of cases lie in 0.071–0.113.
- **Exactly multilinear.** Six coefficient tensors c_kl were extracted with transitive probe modules, using
  P_k(TT_s) = k!·S(s, k). The predicted D equals the observed D in 6 of 6 cases. The double-split coefficients c₂₂
  are of the same size as c₁₁ (for example 19,544 against 6,481).
- **Coarse.** Let both modules range over whole order-7 fibres (12 and 15 classes). The 180 module pairs give only 25
  distinct counts, spread over 2·10⁹ to 6·10⁹ steps (`w38_two_apex_coupled_fill.py`). The mixed term is about 5% of
  the largest single swap. A module swap moves F by a constant fraction of F.

**Verdict.** The wave-37 plan "three sources read additively inside the two-apex placement" fails as stated.
- Module swaps in a two-apex host are exact but coupled.
- Their granularity is set by the host's splitting multiplier, which is a product of the rest's counts.
- Fine repair needs a **thin reader**: covers that split the repair block must complete essentially uniquely, so the
  multiplier is O(1), as in the triangle reader F = 3H + P₂ (M = 1). See GRAND_SYNTHESIS.md §6.

### §3. The synthesis: new statements

**[THM-5242](../../01-canon/theorems/THM-5242-the-triangle-law-cyclic-joins-of-three-blocks-are-a-symmetric-trilinear-jacobi-form-diagonal-joins-read-the-metagraph-closure-mod-16-and-the-collatz-group-is-the-2-adic-norm-group-of-q-sqrt-minus-two.md) (PROVED; finite checks):**
- **(A) Triangle law.** H(cyc(A, B, C)) = Σ admissible mult·P_a(A)P_b(B)P_c(C). The admissible triples are (k,k,k)
  with weight 3 and the rotations of (k+1,k,k) and (k+1,k+1,k) with weight 1. This is a symmetric trilinear form with
  tridiagonal Gram matrix J(v_C): J_aa = p_(a−1) + 3p_a + p_(a+1), J_(a,a+1) = p_a + p_(a+1). With C a point it is the
  ℤ[√−2] form. Checked on 400 of 400 random triples.
- **(B) The metagraph in the norm.** A diagonal join has N(β_A) ≡ 3H² + 8·HC(A) (mod 16) at odd order and
  3H² + 16·c_(n−1)(A) (mod 32) at even order. At odd order HC has the parity of the closure number of A's metagraph
  node. Checked on 300 of 300.
- **(C) Paley primes.** For p ≡ 3 (mod 4): the 3x±1 tables alternate mod p (THM-4792) ⟺ halving reverses P_p ⟺
  p = x² + 2y² ⟺ p ≡ 3 (mod 8). Checked on 332,397 primes below 10^7. By Hasse, alternating primes have density 17/24;
  the count below 10^7 gives 0.70817.
- **(D) 2-adics.** The closure of ⟨2, 3⟩ in ℚ₂^× is the norm group of ℚ₂(√−2), and negation is the non-norm class.
  Every Collatz-type cycle denominator 2^K − 3^k with K ≥ 3 is a non-norm.

**[HYP-9792](../../05-knowledge/hypotheses/HYP-9792-norms-are-thin-but-the-bilinear-family-is-locally-free-diagonal-join-primes-have-density-zero-while-two-apex-sources-equidistribute-modulo-every-prime-at-least-five.md) (CONJECTURE):**
- Diagonal-join primes have density zero. They are complete below 4,190,211: 26 of 73,940 primes ≡ 3 (mod 8). The
  angle obstruction is x/y ≥ 2n − 5.
- The two-apex family equidistributes modulo every ℓ ≥ 5: at order 12 the total-variation distance is 0.005–0.013.
- P103 and P104 are filed.

### Verdicts

| claim | status |
|---|---|
| P₂ ≡ 6 (mod 12) at odd order; widths and fills in §1 | FINITE-EXACT (orders 4–10) |
| relative spread of P₂/H shrinks like about n^(−1.4); absolute fibres widen with H; fill → 1 from about order 12 | NUMERICAL (orders 10–16) / heuristic |
| fibre swaps in two modules of a two-apex source are additive | REFUTED (240/240 nonzero mixed terms) |
| two-apex counts are exactly multilinear in module port vectors | PROVED (module decomposition); 6/6 tensor checks |
| THM-5242 (A)–(D) | PROVED |
| HYP-9792 (i), (ii) | CONJECTURE, with exact and finite evidence |
| the grand synthesis (GRAND_SYNTHESIS.md) | proved bridges, cited KNOWN results and labelled analogies; no open problem solved |

## Wave 39: the literal span-3 certificate (two obstructions) and the frame behind our tools

Session mac-mini-2026-11-15-frame. The owner's prompt had three parts:
- build the syllable-adaptive cone certificate for literal span 3;
- seek Chan–Pak by more creative routes than more branches;
- ask what the new mathematics represents as a deeper unifying theorem, given seven named connections.

The answer to the third part is [UNIFYING_FRAME.md](UNIFYING_FRAME.md). Concurrent and credited:
klein-2026-11-15-ocf attempted the same certificate
([LITERAL_CONES.md](../chanpak_apex_ocf_literal_cones_20261115/LITERAL_CONES.md)).

### §1. The certificate: built pipeline, two obstructions, no certificate yet

**Pipeline (all exact where it matters; `scripts/w39_*`).**
- The 119-state literal transfer. Its minimal representation after the first vertex has dimension 61.
- A deterministic suffix automaton whose memory dictionary is closed under deleting the last letter (floor-unit
  syllables and/or factors of cheap cycles).
- Bottleneck edges removed by an exact rank ≤ 1 test: all entries are nonnegative integers, so rank ≤ 1 means one
  primitive column direction.
- A sound cone LP on raw generators. An earlier draft deduplicated generators by their 61-dimensional image, which is
  unsound because raw generators need not lie in the reachable span. The bug was found before any result.
- An exact integer verifier that re-checks every rank decision, transition and inequality.

**Obstruction A: slow bottlenecks** (HYP-9802 (i), (iii)).
- A three-letter window leaves cycles such as (∅, {1,3})^∞ (growth 1.000), (0,1,7,4) (1.316) and (0,3,2) (1.442) uncut.
- Every cheap cycle of period ≤ 6 letters (3,343 with growth below 1.72) has a rank-one factor of ≤ 9 letters, so a
  2,527-string factor dictionary cuts them all.
- Periods 7–12 need longer memory. klein's census finds up to 12 code units. Syllable-adaptive memory alone cannot
  work: about 99% of the slow cheap patterns are not syllable concatenations, which was klein's correction to the plan.

**Obstruction B: the readout cap (new)** (HYP-9802 (ii)).
- At the letter level, "potential ≤ H on raw generators" is infeasible at every rate, because one letter can split a
  path fragment. The LP on the 6,262-state automaton is infeasible at ρ = 1.2–1.705.
- The fix is a bounded-growth lemma H(W ℓ) ≤ K·H(W). The observed maximum ratio is 39 at word lengths 6 and 7. The
  potential's Hankel coordinates are H(W w_i), so the lemma gives |Φ| ≤ C·H without a cap.
- The growth-only LP did not terminate in 25 minutes at ρ = 1.2 (inconclusive).

**Status.** The literal certificate is not built. The scheme that should work combines klein's rank potential (bounded
loss per rank drop) with the bounded-growth lemma.

### §2. The frame: Chan–Pak as a local–global principle for a positive thin semigroup

Every count we build is a matrix coefficient h·A_W·v₀ of a semigroup of nonnegative integer matrices. Chan–Pak is the
local–global principle for its value set, with exceptions {7, 21}. That is Zaremba's shape, where Bourgain–Kontorovich
proved density one.

Our results are the components:
- local laws at 2 and 3;
- the bottleneck torus with multipliers 3 and 5;
- the growth gap 1.7699 versus 1.7100;
- Hankel and cyclic-join bilinearity;
- equidistribution;
- placement;
- exact repair.

The seven connections are the classical models:
- **Hasse/Artin:** torus images mod p, and the major arcs at 3^i 5^j;
- **Furstenberg ×2×3:** non-lacunary ⟨3, 5⟩ plus bulk expansion;
- **Friedlander–Iwaniec:** Type II structure from concatenation;
- **Ramanujan's form:** a finite exceptional set, with major arcs playing the role of exceptional characters;
- **Tunnell/BSD:** main term versus cuspidal error in r_n(N);
- **class numbers:** the exact mass identity Σ_(H = N) 1/|Aut T| = labelled count / n!, which makes the HYP-9474 local
  limit law a Siegel–Weil formula;
- **Hadamard:** the spectrum's top and optimal flatness.

**Parameters measured for span 3** (`w39_frame_params.py`, `w39_torus_census.py`):

| parameter | value |
|---|---|
| joint spectral radius | ≥ 3.32260 (the letter {1,2,3}) |
| typical growth | 2.2449 per letter |
| counting exponent α | ≥ 1.732, above the 5/3 that Zaremba density one needs (Huang) |
| H mod ℓ over all words of length 9 | TV 0.020–0.034 for ℓ = 5–13, roughly halving per letter; no residue missed |
| eigenvalues of minimal rank-one words (length ≤ 6) | 19 distinct values, including the primes 11–43 (cuts), so the "torus" is ⟨3, 5⟩ only in its cheap part |

**New conjecture (UNIFYING_FRAME.md §1).** Density-one Chan–Pak for span-3 ribbons with a power saving: the number of
missed odd N ≤ X is O(X^(1−η)). Four of its five ingredients are listed with their status in §4 there.

### Verdicts

| claim | status |
|---|---|
| every cheap cycle of period ≤ 6 letters (growth < 1.72) has a rank-one factor of ≤ 9 letters | FINITE-EXACT (3,343 cycles) |
| a letter-level cone certificate with a readout cap is infeasible at every rate | PROVED in outline; FINITE-EXACT (LP) |
| bounded one-vertex growth H(W ℓ) ≤ 39·H(W) for span-3 ribbons | CONJECTURE (exact maximum 39 at lengths 6 and 7) |
| literal span-3 growth certificate | NOT built; the obstructions and the fix are specified (HYP-9802) |
| Chan–Pak = local–global principle for a positive thin semigroup (UNIFYING_FRAME.md) | frame: component statuses as labelled; the density-one form is a CONJECTURE |

## Next waves (queued predictions)

- **P13.** Classify tournaments with no two disjoint triangles; prove the Ω-complete spectrum has density < 1.
- **P14.** P1 at n = 14 with the full prime set: (5,9)/(6,8), (5,9)/(7,7), (6,8)/(7,7).
- **P15.** The exponent pattern of Gap_b. For b ≥ 60, every m < b with 2m + 1 in the 22-list is a gap. Predict which
  of them survive for m ≥ b, extending THM-5070's 30,000 bound to b = 9–12.
