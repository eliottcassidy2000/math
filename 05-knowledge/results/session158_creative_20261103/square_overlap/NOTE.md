# Square-sum intervals and two-edge interference

Status: **proved by path classification**; independently checked on finite examples by two enumeration procedures. The constructions were inspired by square-sum edges and by binary/ternary digit splitting. Novelty relative to the literature and this repository is not asserted. They do not establish a complete tournament path-count spectrum.

## Strong crossing-edge gadget

Let A,B,C be positive integers. Put a=0, c=A, b=A+B, d=A+B+C. Take the increasing tournament on the integer interval [a,d], except reverse the two edges b→a and d→c. With H the number of directed Hamilton paths,

    H(A,B,C) = 1 + 2^(A+B−1) + 2^(B+C−1)
                 − [C=1] 2^(A−1) − [A=1] 2^(C−1)
                 + 2^(A+C−1) 3^(B−1).                    (1)

Here [condition] is 1 when true, 0 otherwise. The tournament is strong: the directed cycles a→a+1→⋯→b→a and c→c+1→⋯→d→c intersect and cover all vertices. It has a feedback vertex set of size at most two, namely {b,d}.

**Proof.** A Hamilton path uses neither, exactly one, or both reversed edges. With neither it is increasing, contributing 1. If only b→a were reversed, the number of paths using it would be 2^(A+B−1): each vertex strictly between its endpoints goes independently before b or after a, in increasing order. Reversing the second edge removes precisely those such paths using the formerly forward edge c→d. This happens only if C=1, because every vertex between b and d otherwise must appear after a but cannot be positioned around the adjacent c→d. In that case all B−1 middle vertices precede b, while A−1 left vertices independently choose a side, giving 2^(A−1) removals. Interchanging the two reversed edges gives the second correction.

For paths using both reversed edges, their two blocks occur in either order. Once an order is fixed, partition the path into increasing zones before, between, and after those blocks. Every one of the A−1 left-outer and C−1 right-outer vertices has exactly two available zones, while every one of the B−1 overlap-interior vertices has three. Assignments uniquely determine the path, since each zone is increasing. The mixed contribution is thus 2·2^(A+C−2)·3^(B−1). This proves (1).

**13 witness.** For (A,B,C)=(1,2,1), on labels 1 through 5 reverse 4→1 and 5→2. Equation (1) yields H=1+4+4−1−1+6=13. Deleting any one vertex leaves a directed triangle: after deleting 1,2,3,4,5, respectively, use the triples (2,3,5), (1,3,4), (1,2,4), (2,3,5), (1,2,4), listed in cycle order. Thus this witness has feedback vertex number exactly 2. This complements the S157 observation that no strong one-apex tournament has 13 paths, but it does not assert that the 13 witness is new.

**Square-sum witness.** On labels 1,...,10 reverse exactly 8→1 and 10→6. Their endpoint sums are 9 and 16, the gaps are (5,2,2), and (1) yields H=265. This specifies only the two selected square-sum edges: all other edges, including other square-sum pairs, remain oriented forward.

**Exact square-sum feasibility.** With a positive starting label a, the two endpoint sums are 2a+A+B and 2a+2A+B+C. Thus both are squares exactly when there exist integers 0<s<r with

    r²−s²=A+C,  s²≡A+B (mod 2),  s²≥A+B+2,

and then a=(s²−A−B)/2. In particular A+C≡2 (mod 4) is impossible. The 13 gadget has A+C=2, so it cannot itself acquire two square-sum reversed edges by translation. Further, for fixed A,C the inequality s≤(A+C−1)/2 shows only finitely many B can qualify; this gives an exact obstruction to reading the following unrestricted recurrence as a simultaneous-square-sum recurrence.

**Recurrence.** For fixed A,C≥2 and varying B≥1, (1) simplifies to

    H_B = 1 + (2^(A−1)+2^(C−1)) 2^B + 2^(A+C−1) 3^(B−1),
    H_(B+3) = 6 H_(B+2) − 11 H_(B+1) + 6 H_B.

This recurrence is for the unrestricted crossing-edge gadget. It makes no claim that the two reversed edges can have square-sum endpoints simultaneously for every B.

## Product construction from isolated square-sum edges

For arbitrary positive k_1,...,k_m, there is a tournament on consecutive positive labels whose only backward edges join pairs with square sum and whose Hamilton path count is ∏_i (1+2^k_i). Put d_i=k_i+1 and choose s_i arbitrarily large of the same parity as d_i so that

    x_i=(s_i²−d_i)/2, y_i=(s_i²+d_i)/2

are positive integers and the intervals [x_i,y_i] are successively disjoint. Reverse y_i→x_i and orient all other edges forward. Each interval's induced tournament has exactly 1+2^(d_i−1) Hamilton paths: either the increasing path, or a path traversing the reversed edge, with each intermediate vertex independently on either side of it. The disjoint intervals and gap singletons occur in their given order in every Hamilton path, giving the product. This is only a multiplicative subfamily, not a coverage theorem.

## Finite checks

Run `python3 check.py` and `python3 -O check.py`. The endpoint-subset DP checks 108 gap triples A,B,C≤6 with at most 11 vertices. For at most 8 vertices, a separate permutation enumerator checks the four edge-use categories, not just their sum. The code also checks the 13 example, its five vertex-deletion triangles, and the 265 example. These finite checks support but do not replace the general path-classification proof.
