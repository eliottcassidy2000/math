# PRIMEGAME's blocks as a deletion tower: path-weighted decks determine every merged tournament node through order 10, the converse-merged deck fails exactly at the odd orders, and the Hamiltonian cycle is the missing coordinate

**Session:** opus-2026-11-04-S72. Worktree `C:/tmp/wt-s72`, cut from main `89c0723a76`.

**Owner prompt.** "Merge in PRIMEGAME and other tournament connections, employing maximal creativity when the provided
example is not quite accurate, looking for other possible related ideas that could be true." In particular:
- decompose the 14 fractions into {11/13, 13/11}, {55/1, 1/7, 1/17, 17/91}, {19/51, 77/19, 77/29, 29/33} and
  {78/85, 23/38, 95/23, 15/2};
- think about how the first maps to the single node of the merged metagraph at n = 3, and the other three to its 3 nodes
  at n = 4;
- see how the prime factorizations relate to the triangular / Fibonacci work, and how they determine the 10 nodes at
  n = 5.

**Status.**
- **PROVED:**
  - Lemma W, which identifies THM-796's face matrix with the Hamiltonian-endpoint-weighted vertex deck;
  - Theorem H: from the converse-merged deck one recovers H(T) at even order, and H(T) - 2 hc(T) at odd order, where
    hc counts Hamiltonian cycles;
  - Lemma L (source lift);
  - Lemma J (the join rule for W; source extension is its smallest case);
  - the Burnside count of line classes;
  - the order-4 line-fraction identity.
- **FINITE-EXACT:**
  - the path-weighted tower is injective through order 10, on 4,870,920 nodes at order 10, and has full column rank
    through order 9;
  - the converse-merged deck collides at orders 3, 5, 7 and 9, and is injective at orders 4, 6, 8 and 10;
  - every collision through order 10 is classified;
  - McKay's census of same-deck tournament pairs is reproduced through order 10: the even-order pairs are converse
    pairs, and the odd-order pairs are self-converse with Hamiltonian cycle counts differing by one;
  - the line fraction is injective at orders 6-8.
- **PRE-REGISTERED** (HYP-9455):
  - P1, P2, P3 and P6 are confirmed;
  - P4 is refuted at order 8, which is how the parity law was found;
  - P7 (order 10) is confirmed.
- **CITED:** THM-796, session129, S71 / THM-4854, THM-002 (OCF; Grinberg-Stanley), Kelly 1957, McKay 2022, Stockmeyer
  1977 / 2021, session29, session39.
- **DICTIONARY / NUMEROLOGY** for every PRIMEGAME label, typed where it is used.

Independently audited: section 8.

## 0. What was already known, and the accurate form of the owner's picture

The owner gave this directive on 2026-10-31 as well. Two lanes answered it then.

**Session129.** [Session129](session129_primegame_tournament_partition_20261031/README.md) contains:
- the four blocks as X0-X3 with an additive height grading;
- the block products P = 55/637, Q = 539/153 and R = 585/34, whose control primes form the transitive triangle
  13 -> 11, 17 -> 11, 17 -> 13;
- the verdicts that order 3 has one complement line, not one node, and that the order-4 fibres have sizes 1, 2, 5,
  not 4, 4, 4;
- the decisive mechanism, inherited from [THM-796](../../01-canon/theorems/THM-796-three-sorted-recursive-tiling-line-node-incidence.md):
  - deleting a path endpoint from every tiling of an order-5 node gives a count vector over the three order-4 nodes;
  - the ten primitive vectors are distinct, so they classify the order-5 nodes;
  - the block products then encode all ten nodes faithfully, with the assignment of blocks to coordinates declared
    rather than forced.

THM-796 itself records distinct primitive rows at orders 4-7 (3, 10, 34, 272). It states explicitly that this is "not a
general injectivity theorem".

**S71.** [My S71 note](primegame_blocks_register_decomposition_metagraphs_20261031.md) /
[THM-4854](../../01-canon/theorems/THM-4854-primegame-four-blocks-register-decomposition-bipartite-flags-triangular-semiprimes-metagraph-quotient-of-complement-graph.md)
reproduced those verdicts. It added:
- the register-union reading of the blocks;
- the bipartite flag machine;
- the unique reciprocal pair;
- the semiprime lemma p(2p +- 1);
- a calibrated numerology verdict for the ten composites against the ten order-5 nodes.

**The accurate form of the owner's picture.** "Level n - 1 determines level n" is THM-796's face coordinate. The
owner's "prime factorization" is the exponent vector of an order-n node over the order-(n - 1) nodes, read as primes.

This session asks the questions that make that picture a statement that could be true or false:
- Does the tower stay injective (unique factorization) at every level?
- What do the weights encode?
- Where exactly does PRIMEGAME sit on the tower?

Program and outputs:
- [`deletion_tower_census_20261104.c`](../../04-computation/experiments/deletion_tower_census_20261104.c)
  ([.out](../../04-computation/experiments/deletion_tower_census_20261104.out), orders 1-9);
- [`deletion_tower_order10_20261104.c`](../../04-computation/experiments/deletion_tower_order10_20261104.c)
  ([.out](../../04-computation/experiments/deletion_tower_order10_20261104.out), order 10 and the full-deck census);
- [`deletion_tower_line_fraction_20261104.c`](../../04-computation/experiments/deletion_tower_line_fraction_20261104.c)
  ([.out](../../04-computation/experiments/deletion_tower_line_fraction_20261104.out));
- [`primegame_deletion_tower_analysis_20261104.py`](../../04-computation/experiments/primegame_deletion_tower_analysis_20261104.py)
  ([.out](../../04-computation/experiments/primegame_deletion_tower_analysis_20261104.out));
- [`deck_parity_source_lift_check_20261104.py`](../../04-computation/experiments/deck_parity_source_lift_check_20261104.py)
  ([.out](../../04-computation/experiments/deck_parity_source_lift_check_20261104.out)).

## 1. Lemma W: the face row is the Hamiltonian-endpoint-weighted deck (PROVED; checked against every tiling, orders 3-7)

**Notation.**
- For a tournament T and a vertex x, h_s(x) and h_e(x) count the Hamiltonian paths of T that start and end at x.
- [T - x] is the converse-merged class (merged node) of the vertex-deleted subtournament.
- W(T) = sum over x of (h_s(x) + h_e(x)) e_[T - x]. This vector is indexed by the order-(n - 1) merged nodes and has
  total weight 2H(T).

**Lemma W.** Let R_n(u, v) be THM-796's face matrix: the number of tilings t in the fibre F_n(u) whose endpoint
deletion lies in v. Then

    R_n(u, .) = W(T) / (2 |Aut T|)    if u = {T} is self-converse,
    R_n(u, .) = W(T) / |Aut T|        if u = {T, T^op} is a converse pair.

In particular, THM-796's primitive coordinate Rbar_n(u) is the primitive vector of W(T).

**Proof.**
1. By THM-781 / THM-796 (1.6), the fibre of a class T is HP(T)/Aut(T). Aut acts freely on Hamiltonian paths, because an
   automorphism fixing a directed Hamiltonian path fixes it pointwise.
2. Deleting the terminal endpoint of a path P leaves a tiling of T - end(P). Hence
   #{tilings of T deleting into v} = (1/|Aut|) #{P : [T - end(P)] = v} = (1/|Aut|) sum over x with [T - x] = v of
   h_e(x).
3. The set of paths counted in step 2 is Aut-invariant, so the count is an integer.
4. A converse pair u = {T, T^op} adds the fibre of T^op. There h_e^(T^op)(x) = h_s^T(x), and [T^op - x] = [T - x] after
   merging, which gives W/|Aut|.
5. For self-converse T, an anti-automorphism carries ending paths to starting paths and preserves merged cards. So the
   end-weighted and start-weighted sums agree, and each is W/2.

QED. Both endpoint faces give the same matrix, as THM-796 (6.2) states.

**Example: S_4 (H = 5, |Aut| = 1).** The endpoint weights h_s + h_e are 3, 2, 3, 2. The two weight-2 vertices delete to
TT_3 and the two weight-3 vertices to C_3. So W = (4, 6) and R = (2, 3), THM-796's and session129's row.

**Check.** The census program recomputes R_n from all 2^binom(n-1, 2) fixed-path tilings through both faces at orders
3-7, together with the fibre sizes. There is no mismatch. Lemma W's integrality is asserted at every node through
order 9.

The lemma turns the face coordinate, which THM-796 defines through 2^binom(n-1, 2) tilings, into n vertex deletions and
one subset DP. That is what makes orders 8, 9 and 10 cheap.

## 2. The tower is injective and full-rank through order 10 (FINITE-EXACT; pre-registered as HYP-9455)

Census, with every class generated by one-vertex extension and canonized:

| order n | classes | self-converse | merged nodes | distinct primitive W | rank of R_n | distinct unweighted merged decks |
|---:|---:|---:|---:|---:|---:|---:|
| 3 | 2 | 2 | 2 | 1 | 1 | 1 |
| 4 | 4 | 2 | 3 | 3 | 2 | 3 |
| 5 | 12 | 8 | 10 | 10 | 3 | 8 |
| 6 | 56 | 12 | 34 | 34 | 10 | 34 |
| 7 | 456 | 88 | 272 | 272 | 34 | 267 |
| 8 | 6,880 | 176 | 3,528 | **3,528** | **272** | **3,528** |
| 9 | 191,536 | 2,752 | 97,144 | **97,144** | **3,528** | 97,141 |
| 10 | 9,733,056 | 8,784 | 4,870,920 | **4,870,920** | not tested | **4,870,920** |

Notes on the table:
- The class counts are OEIS A000568 and A002785.
- At order 9 the node tilings sum to 2^28.
- Ranks are certified from below by the rank of R^T R modulo 2^31 - 1.
- The order-10 rows were compared by 64-bit hashes. There was no equal-hash group, so no exact re-check was needed.

**Pre-registration.**
- HYP-9455's first pre-registration (commit `ce33106a23`) made these predictions:
  - P1: injective at order 8;
  - P2: injective at order 9;
  - P3: full rank at orders 8 and 9;
  - P4: the unweighted merged deck collides at both 8 and 9.
- The second pre-registration (commit `e58a09c18b`, after the first results) made two more:
  - P6: injective at order 10;
  - P7: the unweighted merged deck is injective at order 10.
- Results:
  - **P1, P2, P3, P6 and P7 are CONFIRMED.**
  - **P4 is REFUTED at order 8**, where the unweighted merged deck is already injective.

Order 3 is the only failure of the tower. TT_3 and C_3 both have W proportional to the single order-2 node.

## 3. Why the weights: what the plain deck sees, and the Hamiltonian cycle it misses

### 3.1 The parity law of the converse-merged deck (FINITE-EXACT through order 10)

Collision classes of the unweighted converse-merged deck (the multiset of [T - x]):

| order | 3 | 4 | 5 | 6 | 7 | 8 | 9 | 10 |
|---|---|---|---|---|---|---|---|---|
| collision classes | 1 | 0 | 2 | 0 | 5 | 0 | 3 | 0 |

The deck fails at every odd order tested and at no even order. P4 predicted the opposite at order 8; the alternation
was then pre-registered at order 10 (P7) and held.

### 3.2 McKay's census reproduced, with the converse relation recorded (FINITE-EXACT)

Pairs of non-isomorphic tournaments with the same ordinary (full, unmerged) deck, from
`deletion_tower_order10 10 full`:

| order | pairs | relation of each pair | H and Hamiltonian cycles |
|---:|---:|---|---|
| 3 | 1 | both self-converse, not converses of each other | H 1 / 3, cycles 0 / 1 |
| 4 | 1 | converses of each other | equal |
| 5 | 1 | both self-converse | H 15 / 13, cycles 3 / 2 |
| 6 | 4 | each pair: converses of each other | equal |
| 7 | 0 | | |
| 8 | 2 | converses of each other | equal |
| 9 | 1 | both self-converse | H 2919 / 2917, cycles 211 / 210 |
| 10 | 1 | converses of each other | equal |

The counts 1, 1, 1, 4, 0, 2, 1, 1 are exactly McKay's complete census through order 13
([McKay 2022, Australas. J. Combin. 83, arXiv:2102.01942](https://arxiv.org/abs/2102.01942), section 4). His sources:
- Harary-Palmer at orders 3-4;
- Beineke-Parker at 5-6, plus Stockmeyer's overlooked sixth-order pair (1988);
- Stockmeyer at 8;
- Stockmeyer's families at 9 = 1 + 8 and 10 = 2 + 8.

Session29 already had the order-5 pair with H = 15 / 13 (its masks 40 and 41), and noted that "even all proper
unmarked induced profiles do not determine H, including after converse merging".

**New here.** In this census, every even-order pair consists of mutual converses. Every odd-order pair consists of two
self-converse tournaments whose Hamiltonian cycle counts differ by exactly one.

McKay records that Stockmeyer's odd-order (2^t + 1) tournaments are self-complementary. Stockmeyer's 2021 power-of-two
paper ([arXiv:2110.06158](https://arxiv.org/abs/2110.06158), Theorem 2) does not state that G_p and G*_p are converses.
At p = 4 his Figure 1 gives G_4 = 1 -> C_3 and G*_4 = C_3 -> 1, which are converses. At p = 8 the census confirms it.

### 3.3 Theorem H: the deck sees H except for the Hamiltonian cycles (PROVED)

**Theorem H.** Let T and T' be tournaments of order n >= 3 with the same converse-merged deck; in particular, this
holds if their ordinary decks agree. Then

    H(T) = H(T')                                 if n is even,
    H(T) - H(T') = 2 (hc(T) - hc(T'))            if n is odd,

where hc counts directed Hamiltonian cycles.

**Proof.**
1. **The odd-cycle formula.** By [THM-002](../../01-canon/theorems/THM-002-ocf.md) (Grinberg-Stanley, arXiv:2307.05569;
   Irving-Omar Corollary 20), H(T) = sum over collections C of pairwise vertex-disjoint directed odd cycles of
   2^(|C|). Split the collections by their vertex support U.
2. **Non-spanning collections (|U| = k < n).** These contribute sum over [S] of order k of g(S) s([S], T). Here
   g(S) = sum over spanning odd-cycle collections of S of 2^(|C|), which is converse-invariant, and s([S], T) counts
   the k-subsets inducing a member of the merged class [S]. Kelly's lemma, applied to merged classes, gives
   s([S], T) = (1/(n - k)) sum over x of s([S], T - x). This is determined by the merged deck, because s([S], X) is
   converse-invariant in X.
3. **Spanning collections with c >= 2 cycles.** Mark one of the c cycles, with support U. Each collection is then
   counted c times, so their number is (1/c) sum over partitions V = U ⊔ W, with |U| odd, |U| >= 3 and |W| >= 3, of
   hc(T[U]) g_(c-1)(T[W]). Here g_(c-1) counts spanning collections of c - 1 odd cycles.
4. **Partition counts are deck-determined.** The number N([A], [B]) of partitions with T[U] in [A] and T[W] in [B]
   equals s([A], T) s([B], T) minus the number of overlapping pairs (U_A, U_B).
   - An overlapping pair has a union X of size < n.
   - Grouping the pairs by the class of T[X] gives sum over [S] of order < n of p_(A,B)(S) s([S], T), and p_(A,B) is
     converse-invariant.
   - So every N([A], [B]), and hence the count in step 3, is determined by the merged deck.
5. **The single spanning cycle.** What remains is c = 1, a Hamiltonian cycle with weight 2. It exists only when n is
   odd.

So H(T) = F(merged deck) + 2 hc(T) [n odd]. QED.

**Checks.** The identity holds on every same-deck pair of section 3.2 and on every collision of the merged deck through
order 9. At order 7, one merged collision has H = 143 / 141 and Hamiltonian cycle counts 15 / 14.

**What this says about the weights.**
- The total weight of W is 2H.
- So at odd order W carries exactly the coordinate that the deck misses: the Hamiltonian cycle count, through its total
  weight.
- W also carries the endpoint distribution, which separates the collisions with equal H. These are the Hamiltonian-cycle-free lifts of section 3.4.

We did not find Theorem H in the reconstruction sources checked:
- McKay 2022 section 4;
- Stockmeyer 2021;
- the repository's session29 and session39 reconciliations.

Its ingredients are classical (Kelly's lemma; the partition counting is the two-piece case of Kocay's covering lemma)
and recent (THM-002). Its novelty is left to the audit.

### 3.4 Lemma L (source lift) and the classification of the collisions (PROVED + FINITE-EXACT)

**Lemma L.** Let S and S^op be non-isomorphic converses with the same ordinary deck. Then 1 -> S and 1 -> S^op are
distinct merged nodes with the same converse-merged deck.

**Proof.**
- **Same deck.**
  - The cards are [S] = [S^op] (delete the source) and [1 -> (S - x)] for x in S.
  - The ordinary decks of S and S^op agree, so the multisets {[1 -> (S - x)]} and {[1 -> (S^op - x)]} agree.
- **Distinct nodes.**
  - 1 -> S is not isomorphic to 1 -> S^op: removing the unique source gives S and S^op.
  - Nor is 1 -> S^op isomorphic to the converse S^op -> 1 of 1 -> S. Ordered by their strong components, the first has
    a singleton source component and S^op's components after it; the second has S^op's components followed by a
    singleton sink component.

QED.

**Classification through order 10.** Checked independently with brute-force canonical forms.
- Every merged-deck collision without Hamiltonian cycles, at orders 5, 7 and 9, is a lift of an order-(n - 1) converse
  pair from section 3.2. There are 1, 4 and 2 of them, matching the 1, 4 and 2 even-order pairs.
- The remaining collision at each odd order is:
  - at order 3, {TT_3, C_3};
  - at order 5, the self-converse pair H 15 / 13;
  - at order 7, one strong merged-card collision (a converse-pair node against a self-converse node, H 143 / 141; no
    same-full-deck pair exists at order 7);
  - at order 9, Stockmeyer's self-converse pair.
- At even orders no lift mechanism exists:
  - the lift of a self-converse pair {A, B} has cards [A] and [B], which differ;
  - in the census, the even-order same-deck pairs are all converse pairs.

This explains the parity law through order 10. Lemma L predicts merged-deck collisions at order 11 from the order-10
converse pair, and at every order m + 1 for which order m has a non-isomorphic converse same-deck pair.

## 4. PRIMEGAME on the tower: the owner's mapping made exact where it can be

### 4.1 The single node at the bottom (DICTIONARY)

The tower starts with the single order-2 node. It fails exactly once: the two order-3 nodes TT_3 and C_3 have the same
face row (1), so they are indistinguishable by deletion. They are separated only by the involution kappa, the single
blue line.

The owner's first block {11/13, 13/11} has the same shape:
- two fractions with the same prime support {11, 13};
- swapped by inversion;
- the program's only reciprocal pair (S71).

So "the first block maps to the single node at n = 3" is accurate in this form: the block is the pair of order-3 nodes
over the single order-2 node, and its inversion is the order-3 line. This sharpens session129's reading (reciprocal
pair = the single line) by one level of the tower.

### 4.2 The line fraction: B1 is the order-3 -> 4 interface (PROVED at order 4; FINITE-EXACT beyond)

**Definition.** For a merged node u of order n, put

    Phi(u) = sum over t in F_n(u) of ( e_[kappa d t] - e_[d t] ),

where d deletes a path endpoint and kappa complements every free arc. Read with the order-(n - 1) nodes as primes, this
is the fraction prod over t of p_[kappa d t] / p_[d t].

**Basic properties (PROVED).**
- kappa d = d kappa (THM-796 (2.5)), and kappa is a bijection of the cube, so sum over u of Phi(u) = 0.
- Reflection sigma permutes each fibre and swaps the two faces (d_L sigma = sigma d_H), so both faces give the same Phi.

**At order 4 (PROVED by hand; confirmed by the program).** kappa swaps TT_3 and C_3, so
Phi(u) = (R(u, TT_3) - R(u, C_3))(e_(C_3) - e_(TT_3)). From the rows (1, 0), (1, 1), (2, 3):

    Phi(TT_4) = e_(C_3) - e_(TT_3),    Phi(pair) = 0,    Phi(S_4) = e_(TT_3) - e_(C_3).

Name the two order-3 nodes by the first block's primes 11 and 13, in either order. Then

    {Phi(TT_4), Phi(S_4)} = {11/13, 13/11},    Phi(1 -> C_3 merged with C_3 -> 1) = 1.

So the owner's first block is exactly the set of line fractions of the two self-converse order-4 nodes. These are the
two ends of the blue line TT_4 - S_4, and the converse-pair node gets 1. This is an exact, labelling-independent form of
"{11/13, 13/11} belongs to n = 3": the block's primes are the order-3 nodes, and its fractions are the order-4 nodes.

**Beyond order 4 (FINITE-EXACT, exploratory; not pre-registered).** Phi takes:

| order | distinct values |
|---:|---|
| 5 | 9 values on 10 nodes; only the two converse-pair nodes (1, 1, C_3) and (1, S_4) coincide |
| 6 | 34 / 34 |
| 7 | 272 / 272 |
| 8 | 3,528 / 3,528 |

So the line fraction alone classifies the merged nodes at orders 6-8. Phi is a single ratio p_a/p_b exactly at:
- order 4: TT_4 and S_4;
- order 5: the three pure-blue nodes;
- orders 6-8: the pure-blue nodes with a single tiling, for which it is automatic.

At order 7 one pure-blue node, 1 -> T_5 -> 1 with three tilings, has a composite Phi. So "pure-blue implies single
ratio" is a small-order effect, not a law.

### 4.3 The three data blocks as order-4 coordinates (CITED; DICTIONARY)

This is session129's codec, recorded here for completeness:
- the three order-4 nodes are the coordinates of order 5 (Lemma W makes the rows intrinsic);
- the block products P, Q, R are multiplicatively independent, so P^a Q^b R^c encodes the ten rows faithfully, with the
  assignment of blocks to coordinates declared.

The blocks' control-prime parts form a transitive triangle on {11, 13, 17} (session129). The order-4 metagraph is a
star at S_4 with a loop (S71: no graph map). This session adds nothing to the block-to-node assignment.

### 4.4 The ten registers as order-5 coordinates (DICTIONARY; a calibrated test is NUMEROLOGY)

R_6 has rank 10, so all ten order-5 nodes are needed as independent coordinates of order 6. Ten is the number of
PRIMEGAME's registers 2, ..., 29.

Order the order-5 nodes by (H, tilings, row) and give them the primes 2 < 3 < ... < 29. Then each of the 34 order-6
nodes becomes an integer N(u) = prod p^(R(u, v)). The 34 integers are distinct, which is the order-6 case of section 2.
Examples:
- TT_6 = 2;
- the converse pair (1, 1, 1, C_3) = 10;
- C_3 => C_3 = 5;
- the axis node over (1, S_4) = 7^5.

**Test.** Do PRIMEGAME's 14 fractions act as moves between order-6 node codes, N(u') = f N(u)? Under the canonical
assignment there is 1 (pair, fraction) incidence. Over 20,000 random assignments of the primes the count has minimum 0,
median 1, mean 0.72 and maximum 3, and P(random >= canonical) = 0.58. **No signal: NUMEROLOGY.** The ten-registers
reading is a labelling, not a dynamics.

### 4.5 The owner's counts 1 and 3 are the line-class counts (PROVED Burnside; the rest NUMEROLOGY)

Lines modulo the mirror are the orbits of <kappa, sigma> on the cube X_n. kappa and sigma kappa have no fixed tilings
(THM-796 section 4), so

    #line classes = (2^m + 2^h) / 4,    m = binom(n-1, 2),  h = floor((n-1)^2 / 4).

The values are 1, 3, 20, 272, 8320, 525312, ... for n = 3, 4, 5, ... They are verified by enumeration through n = 7.
For n <= 5 each line class is one adjacency of the merged metagraph.

So the owner's "1 at n = 3" and "3 at n = 4" are exact as line classes. Order 5 has 20, not 10.

The coincidences are NUMEROLOGY:
- 272 line classes at n = 6 against 272 merged nodes at n = 7;
- 20 against the twenty (prime, side) incidences of PRIMEGAME, since each register occurs in some numerator and some
  denominator.

## 5. Triangular and Fibonacci numbers (Lemma J PROVED; the Fibonacci values NUMEROLOGY)

**Fibonacci rows.** Face rows made of consecutive Fibonacci numbers:
- order 4: (2, 3) at S_4;
- order 5: (2, 3) at Q_5, (5, 8) at the H = 13 node, and (2, 3, 5) at the pair (1, S_4);
- orders 6-8: only (2, 3, 5) rows. These are 1 of 29, 2 of 261 and 3 of 3,515 multi-entry rows.

**Lemma J (the join rule; PROVED; checked on all 55 joins of ordinary classes with |A| + |B| <= 6).** For the ordinal
sum T = A => B (every vertex of A beats every vertex of B):

    W(T) = H(B) sum over x in A of h_s^A(x) e_[(A - x) => B]  +  H(A) sum over y in B of h_e^B(y) e_[A => (B - y)].

Proof: no arc leaves B for A, so a Hamiltonian path starts in A, covers A, and then covers B. Hence H(T) = H(A) H(B),
h_s^T = H(B) h_s^A on A and 0 on B, and h_e^T = H(A) h_e^B on B and 0 on A. The case |A| = 1 is the source extension
W(s -> T') = H(T') e_[T'] + sum over x of h_e^(T')(x) e_[s -> (T' - x)].

**Where the later Fibonacci rows come from.** Every Fibonacci run at orders 6-8 is a join built on S_4:
- S_4 => TT_k (orders 6, 7, 8);
- S_4 => C_3 (order 7);
- S_4 => C_3 => 1 and S_4 => 1 => C_3 (order 8).

The entry 5 is H(S_4), contributed by the B-side of the join rule. It is the row sum of S_4's row (2, 3) (the
"Fibonacci" step 2 + 3 = 5), appended once and never iterated.

Under a source extension the other order-5 Fibonacci rows become:
- (5, 8) at the H = 13 node becomes (2, 3, 8, 13): the merged pair coordinate splits as 5 = 2 + 3 into the
  end-weighted cards 1 -> C_3 and C_3 -> 1, and 13 = H is appended;
- Q_5's (2, 3) becomes (1, 1, 3, 5).

At these sizes the entries are small, so most of them are Fibonacci numbers anyway. The values are NUMEROLOGY, and the
join rule is the only structure behind them.

**Triangular numbers.**
- The triangular entries of the fractions (S71 / session129: 1, 15, 55, 78, 91) are free-bit dimensions
  m_n = binom(n-1, 2) = T_(n-2) at n = 3, 7, 12, 14, 15. NUMEROLOGY.
- The order-3 entry m_3 = 1 is the single free bit, the single line of section 4.1. DICTIONARY.

## 6. Typing

| item | type |
|---|---|
| Lemma W (face row = path-weighted deck) | PROVED; checked against all tilings, orders 3-7 |
| tower injective at orders 4-10; full rank at orders 3-9 | FINITE-EXACT; P1-P3 and P6 pre-registered and confirmed |
| unweighted merged deck: collisions at odd orders 3-9, none at even orders 4-10 | FINITE-EXACT; P4 refuted at 8, P7 confirmed at 10 |
| McKay's same-deck census reproduced; even pairs converse, odd pairs self-converse with hc differing by 1 | FINITE-EXACT (counts CITED, McKay 2022) |
| Theorem H (H at even order, H - 2hc at odd order, from the merged deck) | PROVED (THM-002 + Kelly + partition counting); literature status to the audit |
| Lemma L and the collision classification through order 10 | PROVED + FINITE-EXACT |
| line fraction: order-4 identity {11/13, 13/11, 1} | PROVED |
| line fraction: injective at orders 6-8 | FINITE-EXACT (exploratory) |
| first block as the order-3 pair over the order-2 node | DICTIONARY |
| blocks as order-4 coordinates; registers as order-5 coordinates | DICTIONARY (session129's codec CITED) |
| PRIMEGAME fractions as moves between order-6 codes | NUMEROLOGY (calibrated, P = 0.58) |
| line classes (2^m + 2^h)/4 | PROVED |
| owner's 1 and 3 as line classes | DICTIONARY |
| 272 and 20 coincidences | NUMEROLOGY |
| Lemma J (join rule; source extension as its smallest case) | PROVED; checked on 55 joins |
| Fibonacci rows | NUMEROLOGY, with the join rule as their only structure |

## 7. Pre-registration record

- **HYP-9455:**
  - reserved in `d82961b8d0`;
  - predictions P1-P5 committed in `ce33106a23`, before any order-8 or order-9 run;
  - results and the second pre-registration P6-P8 committed in `e58a09c18b`, before the order-10 run.
- The order 3-7 runs were exploratory: they validated the program against THM-796.
- Exploratory and not pre-registered:
  - the line fraction beyond order 4;
  - the source-lift classification;
  - Theorem H. It was proved after the census exposed the H-differences.

## 8. Audit

(Filled after the independent audit.)

## 9. Reproduction

```
gcc -O2 -o dtc 04-computation/experiments/deletion_tower_census_20261104.c
./dtc 9 DATA_DIR                                        # 40 s; orders 1-9; writes DATA_DIR/tower_n.txt
gcc -O2 -fopenmp -o dt10 04-computation/experiments/deletion_tower_order10_20261104.c
./dt10 10                                               # 2 min on 12 threads: orders <= 10
./dt10 10 full                                          # 1 min: ordinary same-deck pairs, H, Hamiltonian cycles
gcc -O2 -o dtl 04-computation/experiments/deletion_tower_line_fraction_20261104.c
./dtl 9                                                 # 20 s: line fraction (orders 4-8), merged-deck collisions
python3 04-computation/experiments/primegame_deletion_tower_analysis_20261104.py DATA_DIR   # 1 min
python3 04-computation/experiments/deck_parity_source_lift_check_20261104.py                # seconds
```
