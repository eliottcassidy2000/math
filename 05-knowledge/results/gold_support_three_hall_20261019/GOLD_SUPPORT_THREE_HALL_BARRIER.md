# A structural support-three barrier, and support-four repair

19 October 2026. **PROVED structural lemmas; FINITE-EXACT application and repair.**

This note does not prove or refute seven-thin Gold, H7, or the existence of a support-three certificate after a different first comparison. It concerns the capacity of all exchanges moving at most three labels after one specified first comparison. No historical-novelty claim is made.

## 1. Main result

Use the established rigid, seven-thin, no-good-pair seed

    Q = (0,0,0,6,0,19,2,18,86,151).

A tuple records strict predecessor bitmasks; labels are 0 through 9. The closest inherited sources are packet136, `conditional_gold_exchange_and_minimal_obstruction_20261019.md`, and packet140, `gold_module_reduction_and_triple_exchange_20261019.md`.

Choose the first comparison 5<9. Its larger child is

    R = (0,0,0,6,0,19,2,18,86,183).

The exact totals are

    E = e(Q) = 11971,
    b = e(R) = 7741,
    3/5 < b/E < 2/3,
    required matching capacity = 2b−E = 3511.

Among every possible second comparison, only the unordered pair {0,2} can complete this first comparison to weak Gold. Its heavy cells are 3738 and 4003. Both heavy Gold margins are strictly positive, namely 492 and 227. The smaller first branch is automatic, as in packet136.

Nevertheless, **every matching or fractional matching built from every legal label permutation of support at most three has capacity at most 3187 < 3511** for the pair {0,2}. This is a short structural Hall bound, proved below from four small extension counts. A separate complete graph certificate sharpens the optimum to 3017. There are no isolated vertices on either side: the obstruction is genuinely collective, not just an extension with no legal local reversal.

Support four repairs the same first comparison. A supplied matching of 3738 legal support-at-most-four exchanges saturates the smaller second cell. It proves the heavy Gold margin at least

    E−2b+3738 = 227.

Thus this first comparison admits strict Gold, admits no support-three matching/fractional-matching certificate for any second pair, and admits a support-four certificate. This is not a failure of Gold and does not exclude a different first comparison with a support-three certificate. The inherited work already supplies successful alternative first comparisons on Q.

## 2. Exchange graphs and a Hall upper bound

Let R be any finite poset and let x,y be incomparable vertices. Put

    A = {L in L(R): x precedes y},
    B = {L in L(R): y precedes x}.

The complete support-k exchange graph joins L in A to M in B exactly when a permutation of the labels moving at most k labels sends L to M. Since labels occur once, this is equivalent to the two full words having Hamming distance at most k. Both words must remain legal extensions; no permutation is presumed to be an automorphism.

For any S contained in A, a matching has at most

    |A\S| + |N(S)|

edges: edges from S must use distinct vertices of N(S), while the other edges use distinct vertices of A\S. Exactly the same bound holds for fractional matchings whose incident weight at each extension is at most one. It suffices to bound N(S) from above; no search for an optimum matching is needed.

## 3. Two chains identify all possible escape cycles

Suppose four distinct vertices satisfy x<X and y<Y in R. No assumption that these four induce precisely 2+2 is needed. Define

    S = {L: x<y and X<Y in L},
    T = {L: y<x and Y<X in L}.

Here and below inequalities inside an extension describe its word order, not necessarily an additional relation of R.

**Lemma 1 (two-chain crossing classification).** If a legal permutation of support at most three sends a word in S to a word in T, then it is one of the following three-cycles for a vertex z outside {x,X,y,Y}:

    sigma_A(z) = (x z X),
    sigma_B(z) = (z y Y).

The first cycle maps x to z, z to X, and X to x. The second maps z to y, y to Y, and Y to z. Moreover z is incomparable in R with every one of x,X,y,Y.

More precisely, the order of the five displayed vertices must be respectively

    source: x<y<X<Y<z       target: z<y<x<Y<X,

or

    source: z<x<y<X<Y       target: y<x<Y<X<z.

Other labels may occur between the displayed labels and do not change the statement.

**Proof.** The four-label source word is either x X y Y or x y X Y; the four-label target word is either y Y x X or y x Y X. If the exchange moves no label outside these four, each allowed source/target pair changes all four occupied positions. Support at most three is impossible.

Otherwise, at most two of the four labels move. To reverse both comparisons x/y and X/Y, one moved label must belong to each pair. The two crossed choices {x,Y} and {y,X} cannot work. For {x,Y}, their new positions would both have to lie strictly between the fixed y and X, while their old positions lie outside that interval; one outside label supplies only one new position there. For {y,X}, the new positions must lie on opposite sides of the fixed x and Y, again requiring two new positions outside the old positions. These statements also rule out transpositions.

For the remaining choice {x,X}, fixed y and Y force the new X to occupy z's position after Y, and the new x to occupy X's old position between y and Y. This gives exactly the first displayed cycle and five-label orders. The choice {y,Y} is its dual and gives the second. In either case z crosses all four other labels, so legality before and after forces it to be incomparable with all four. ∎

Let Z be the common incomparability neighborhood of the four vertices. Adjoin the source or target comparison relations and close transitively:

    R_S = R + (x<y) + (X<Y),
    R_T = R + (y<x) + (Y<X).

For a label permutation sigma, let sigma R_S be its relabeling. An inconsistent relation union has zero extensions. Define

    d_A(z) = e(R_T union sigma_A(z) R_S),
    d_B(z) = e(R_T union sigma_B(z) R_S).

All unions include transitive closure. These counts are exactly the target extensions in T whose inverse exchange lies in S.

**Theorem 2 (support-three Hall cap).** Every integral or fractional matching in the complete support-three exchange graph satisfies

    nu <= D + sum over z in Z of (d_A(z)+d_B(z)),                 (1)

where

    D = #{L: x<y and Y<X} + #{L: y<x and X<Y}.

In words, D counts extensions on which the lower-endpoint and upper-endpoint comparisons disagree. The extra terms are the only possible simultaneous reversals, and require a common incomparable escape label.

**Proof.** Choose the Hall set S above. Its complement in A is the first discordant cell. Every neighbor of S either lies in the second discordant cell or in T. Lemma 1 places the latter neighbors in the union of the stated cycle domains. The union bound and the Hall upper bound prove (1). Overlap of cycle domains only makes this an upper bound, which is the needed direction. ∎

If Z is empty, the cap simplifies to nu<=D. The theorem is a structural obstruction and not a lower bound: abundant escape domains do not, by themselves, supply a matching.

## 4. Four escape counts prove the rigid-seed obstruction

In the child R of Section 1, use

    (x,X,y,Y) = (0,5,2,3).

The relations 0<5 and 2<3 hold. The only labels incomparable with all four are 6 and 7. Thus Lemma 1 eliminates every other candidate escape before any exchange-graph construction.

The exact ideal counts are

    |A| = 3738,
    |S| = e(R + 0<2 + 5<3) = 2616,
    #{L: 2<0 and 5<3} = 1844.

Hence D=(3738−2616)+1844=2966. The only four possible escape terms are:

- cycle (0 6 5): 46;
- cycle (6 2 3): 72;
- cycle (0 7 5): 55;
- cycle (7 2 3): 48.

Their sum is 221, giving

    nu <= 2966+221 = 3187 < 3511.

For direct replay, the four transitive predecessor tuples of the relation unions are, in that order:

    (86,0,66,87,0,95,2,18,127,255),
    (6,0,2,7,0,31,63,18,127,191),
    (150,0,146,151,0,159,2,18,214,191),
    (22,0,18,23,0,31,2,63,95,191).

These are tiny ten-element ideal recurrences. With a predecessor tuple p, the recurrence starts with f(0)=1 and adds a vertex v precisely when all bits of p[v] are already in the ideal. Summing f over the available additions gives e(p). No floating-point optimization is used.

The exact union of the four exceptional target sets has size 184, rather than 221; its sharper Hall cover has size 3150. Neither refinement is needed to prove failure. The complete graph optimum 3017 is stronger again. It is important not to add these overlapping cycle domains as a lower-bound matching certificate; only the upper-bound union estimate is used here.

## 5. Why no second-pair choice repairs support three for this first comparison

Gold requires the affirmative heavy count c to lie in

    [2b−E,E−b] = [3511,4230].

The following list gives all genuinely undetermined unordered pairs in R, with the lower numeric label asked first. Comparable pairs have count 0 or 7741 and cannot help.

    01:2416  02:3738  03:6777  04:3132  06:5516  07:6434  08:7320
    12:5231  14:4724  24:3300  25:6835  26:5728  27:6471
    34:525   35:3281  36:1946  37:2799  38:4926  39:5853
    46:6284  56:2268  57:3217  58:5628  67:4939  69:7239
    78:6054  89:5002

Only 02 lies in the required interval. Reversing a pair replaces c by b−c, which preserves membership in this symmetric interval. This finite-exact check exhausts all second comparisons for this one first comparison. It is not a census of posets or first comparisons.

## 6. Support four is sufficient on the same child

The supporting certificate files index A and B by lexicographically sorted full extension words. The small verifier reconstructs those words from R.

- `rigid_support3_certificate.json` contains 3017 disjoint exchange edges, a left cover of size 1587, and a right cover of size 1430. Every possible support-at-most-three edge meets this cover. Therefore the optimum is exactly 3017. The graph has 30897 edges and no isolated vertices.
- `rigid_support4_certificate.json` contains 3738 disjoint exchange edges. Each word pair has Hamming distance at most four and opposite 0/2 orientation. Every word of A occurs once. Thus no optimal-matching claim or exhaustive support-four edge search is needed to prove the repair. For comparison, the full graph has 157835 edges.

The structural cap in Sections 3–4 is the primary obstruction proof. The optional full matching data sharpen the numerical boundary and verify the support-four survivor. At the exact support-three optimum the shortfall is 3511−3017=494. The support-four matching exceeds the threshold by 227.

## 7. Reproduction and scope ledger

Run from this directory:

    python verify_gold_support_hall.py
    python -O verify_gold_support_hall.py

The verifier is self-contained and uses only the standard library. It checks:

1. Q and R, the two relevant degrees/counts, and absence of good pairs in Q and its dual;
2. every second comparison in R and uniqueness of {0,2};
3. the general five-label crossing classification on every five-label order respecting the two chains and every support-at-most-three permutation;
4. all Hall counts, the common incomparability set, all four escape tuples, and the cap;
5. every supplied matching edge and endpoint capacity;
6. the support-three cover against every legal transposition/three-cycle;
7. independently, equality of that generated support-three graph with the full-word Hamming-distance graph.

The resulting compact receipt records the consequence objects: failed support-three threshold, the structural and optimal upper bounds, and a passing support-four Gold margin. Reproduction is finite-exact evidence for the specified seed. Lemma 1 and Theorem 2 are general proofs.

Inheritance and next frontier: packet140's positional cap obstructed all literal pair-swapping involutions while support three repaired its chosen tree. Here the missing coordinate is a joint ordering of two chains, and its only support-three escape is a common incomparable vertex. Hall capacity can fail even when every extension has local neighbors. The strongest demonstrated survivor is a support-four repair of this specific strict-Gold first comparison. A universal theorem selecting a first comparison and a sufficiently large bounded-support matching for seven-thin posets remains open.
