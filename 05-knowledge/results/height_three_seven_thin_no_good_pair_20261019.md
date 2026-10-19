# A rigid ten-point seven-thin height-three poset with no good pair, and strict Gold

Date: 2026-10-19. Status: explicit counterexample to two structural hypotheses, independently verified by exact integer counting and full linear-extension enumeration. It is not a counterexample to balance or Gold. No minimum-order or novelty claim is made.

## 1. Main result and the hypotheses it refutes

There is a ten-point poset of height three and maximum incomparability degree seven with no very good pair and, more strongly, no good pair in either the poset or its dual. It is rigid, has no exact half-probability comparison, and nevertheless has a fixed two-comparison weak Gold strategy with strictly positive margins on all four paths.

Thus both proposed statements are false:

- Every seven-thin poset has a good pair.
- Every seven-thin poset of height at most three has a very good pair.

The [height-two seven-thin theorem](seven_thin_height_two_very_good_gold_20261019.md) remains valid and its height restriction is substantive. The [bounded K4 deformation result](k4_bounded_seven_thin_deformations_20261019.md) also remains valid within its explicitly restricted family. The general seven-thin balance/Gold questions are not resolved by the present example; the example itself satisfies both.

We use Zaguia's probability-sensitive definition: in P or its dual, an ordered incomparable pair (a,b) is good when D(a)⊆D(b), U(b)\U(a) is a chain and p(a<b)≤1/2. A very good pair has equal downsets and both directed upset differences chains. [Primary definitions 1 and 5](https://arxiv.org/pdf/1610.00809v3). The [attributed good-pair Gold synthesis](good_pair_gold_synthesis_20261019.md) gives sufficient conditions for Gold, not necessary ones; this example illustrates that distinction.

## 2. The exact poset

Label the vertices 0,...,9. Its strict predecessor masks, with bit i representing vertex i, are

    (0,0,0,6,0,19,2,18,86,151).

Equivalently, the complete list of cover relations is

    1<3, 2<3;
    0<5, 1<5, 4<5;
    1<6;
    1<7, 4<7;
    2<8, 4<8, 6<8;
    0<9, 2<9, 7<9.

Take the transitive closure of these relations. The displayed masks already include that closure. The chains 1<6<8 and 1<7<9 show height at least three, and every chain has at most three vertices. A valid three-level rank assignment is

    (0,0,0,1,0,1,1,1,2,2).

The incomparability degrees, in vertex order, are

    (7,3,6,7,5,6,7,6,5,4).

The poset is therefore seven-thin and not six-thin. It has exactly 11,971 linear extensions.

## 3. Elementary certificate of no very good pair

The only distinct vertices with equal strict downsets are pairs among the four minima {0,1,2,4}. The only distinct vertices with equal strict upsets are pairs among the four maxima {3,5,8,9}. All other predecessor and successor sets are distinct within their respective lists.

For each pair of minima, the following indicated directed upset difference contains two incomparable vertices:

| Pair | Directed upset difference | Incomparable members |
|---|---|---|
| 0,1 | U(1)\U(0) | 3,6 |
| 0,2 | U(2)\U(0) | 3,8 |
| 0,4 | U(4)\U(0) | 7,8 |
| 1,2 | U(1)\U(2) | 5,6 |
| 1,4 | U(1)\U(4) | 3,6 |
| 2,4 | U(4)\U(2) | 5,7 |

For each pair of maxima, a directed downset difference similarly contains an incomparable pair, ruling out a very good pair in the dual:

| Pair | Directed downset difference | Incomparable members |
|---|---|---|
| 3,5 | D(5)\D(3) | 0,4 |
| 3,8 | D(8)\D(3) | 4,6 |
| 3,9 | D(9)\D(3) | 0,4 |
| 5,8 | D(8)\D(5) | 2,6 |
| 5,9 | D(9)\D(5) | 2,7 |
| 8,9 | D(9)\D(8) | 0,7 |

These twelve checks, directly visible from the masks, exclude every very good pair without using probability computation or a numerical solver.

## 4. Exact probability certificate of no good pair

For a placed ideal S, define F(S) as the number of extensions of the unplaced remainder. The recurrence

    F(S)=sum F(S∪{v}), over v not in S with D(v)⊆S,
    F(P)=1

uses exact integers. Forward ideal counts multiplied by the suffix counts F give each numerator N(a<b). Every ordered incomparable pair is tested in both P and its dual for the structural conditions D(a)⊆D(b) and chain U(b)\U(a). Exactly 48 orientations are structurally eligible. Every eligible numerator is at least 6,547, whereas E/2=11,971/2. Thus all eligible probabilities strictly exceed one half, and no good pair exists.

The minimum is attained by (a,b)=(2,0) in P, at 6,547/11,971. The full receipt lists all 48 eligible orientations and their exact rational probabilities; the proof does not infer absence of good pairs merely from absence of very good pairs. The verifier checks eligibility directly from masks in both orientations and tests the exact integer inequality 2N>E.

## 5. Explicit strict Gold and balance

First compare vertices 3 and 5, then compare vertices 2 and 4, regardless of the first answer. The four extension cells, ordered by answers (true,true), (true,false), (false,true), (false,false), are

    (3281,3281,1899,3510).

They sum to 11,971. The first comparison leaves 6,562 extensions on its true branch and 5,409 on its false branch. Hence the Gold deficits E−t1−t2 on the four paths are

    (2128,2128,4663,3052),

all strictly positive. Full enumeration of all linear extensions independently verifies these cell counts. Thus the example has a fixed non-adaptive two-comparison Gold strategy, even without using the permission to repeat a determined comparison.

The first comparison is balanced: p(3<5)=6,562/11,971 lies strictly between 1/3 and 2/3. No pair has probability exactly one half because the extension total is odd.

Nor is this an automorphism-based explanation. The pairs (|D(v)|,|U(v)|), in vertex order, are

    (0,2),(0,6),(0,3),(2,0),(0,4),
    (3,0),(1,1),(2,1),(4,0),(5,0).

They are all distinct, so every order automorphism fixes every vertex. The poset is rigid.

## 6. Additional twelve-point witness and discovery scope

The same bounded search also returned the twelve-point predecessor tuple

    (0,0,0,0,7,13,14,12,206,5,685,535).

It too has height three, maximum incomparability degree seven, no very good pair, and no good pair. Its degree sequence is

    (6,7,3,6,7,7,7,7,6,7,5,6),

and it has 61,180 extensions. Here vertices 0 and 3 form a half pair: each orientation has 30,590 extensions. Comparing this pair and repeating gives weak Gold. The independent exact verifier checks both examples, including this half strategy by full enumeration. The ten-point example is the main result because it is smaller and requires no half-pair symmetry explanation. No order-minimality claim follows.

For discovery only, an installed HiGHS mixed-integer solver was given naturally labeled relation Booleans, transitivity, degree≤7, integer ranks 0,1,2 increasing along relations, non-totality, and exact Boolean conditions excluding very good pairs in P and its dual. Each size had a 120-second limit. Feasible models were found at sizes 10 and 12 in approximately 5.2 and 25.6 seconds, respectively. No large census or UNSAT theorem was attempted.

The no-very-good condition says for each incomparable pair: either the downsets differ, or at least one directed upset difference contains an incomparable pair; the dual condition is imposed separately. Auxiliary binary witnesses encode these disjunctions. The solver's floating-point output is only a source of candidate masks. The mathematical conclusion depends solely on the elementary mask checks and exact counting described above, not on trusting its feasibility status.

## 7. Arbitrarily large rigid no-good-pair families with strict Gold

Let Q be the ten-point witness. For every k≥1, take the ordinal sum of k copies of Q, optionally inserting finite chains before, after or between the copies. Every element of an earlier block is below every element of a later block.

The result remains seven-thin: all new cross-block pairs are comparable, so each vertex in a Q block keeps its original incomparability degree, and chain vertices have degree zero. It has no good pair. Indeed, only within-block incomparable pairs need consideration; for two vertices in one Q block, earlier vertices are added to both downsets and later vertices to both upsets. Downset inclusions and directed upset differences are therefore unchanged. Every extension is the concatenation of independent block extensions, so each within-block comparison probability is also unchanged. Chain blocks contain no incomparable pair.

The ordinal sum remains rigid. An automorphism preserves downset sizes. Vertices in a block of size m preceded by t vertices have downset sizes in the interval [t,t+m−1], disjoint from those of all other blocks. Every block is therefore preserved setwise; each Q block is rigid and each chain block is fixed pointwise.

The extension total is E=11,971^k, which is odd, so there is no exact half-probability pair. Make the two comparisons from Section 5 inside any one Q block. Every cell count and Gold margin is multiplied by 11,971^(k−1). Thus all Gold margins remain strictly positive, and the minimum normalized margin is exactly 2,128/11,971. The family has arbitrarily many vertices and unbounded height (height 3k before optional chain insertions), with no good pair, rigidity, no half pair, and strict Gold simultaneously.

This corollary refutes no additional balance or Gold statement and claims no novelty. It shows that the structural obstruction is not confined to one bounded solver instance; only the original witness has height three.

## 8. Replay and attribution

Run `python -O ../experiments/verify_height_three_seven_thin_no_good_20261019.py` from the results folder. This standard-library verifier independently checks natural labeling, transitivity, height, actual incomparability degrees, every very-good-pair candidate, every good-pair candidate and probability, and both explicit Gold strategies. Its dynamic-programming extension totals agree with full extension enumeration. Normal and optimized receipts are supplied; explicit guards remain active under `-O`.

The optional discovery script requires the already installed NumPy/SciPy HiGHS solver; those packages are unnecessary for the proof replay. Separate discovery receipts record the bounded runs. Zaguia's definitions and the inherited Gold framework are explicitly attributed above. General relation encoding, ideal counting, and certificate verification are standard methods; no historical priority claim is made for the examples or their analysis.
