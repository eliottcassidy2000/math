# A finite boundary atlas at every bounded joint response

**PROVED all-size reduction and formulas; FINITE-EXACT complete tail lists.**
Let a->b be two in-neighbours of a root r. The residual tournament is
the transitive chain u1->...->um; r beats that chain. Write each chain
vertex as the letter `(A_i,B_i)=([ui->a],[ui->b])`. Assume first11 and
last not00. These are exactly the strong words in this chart with q=s.
Set U=T-r and C=W_ab(U)+W_ba(U)=F-H, using ordered distinct path bins.

## A1. Positive energy bounds the tail

An exact closed formula, independent of the inherited automaton, is

    C = A_m+B_m
        + sum_(i=1..m-1) [A_i(1-B_(i+1))+B_i(1-A_(i+1))]
        + sum_(j<i) 2^(i-j-1)[(1-B_j)A_i+(1-A_j)B_i].      (A1)

Proof: in a W_ab cover assign every chain vertex to the path before a
or the path after b. The all-before assignment contributes A_m. For a
mixed assignment let i be the last before-a vertex and j the first
after-b vertex. If i<j, necessarily j=i+1, with a unique assignment. If
j<i, all vertices between them can be assigned freely, giving2^(i-j-1).
The two endpoint conditions are A_i=1 and B_j=0. First11 excludes an
all-after assignment. Interchange a,b and add. In particular this sum
does not use the arc between a and b.

Remove the maximal initial11 prefix and let L be the remaining length.
For fixed C, either L=0 or

    L <= 2 floor(log2 C)+4.                              (A2)

Let ell be the first non11 position. If its zero coordinate couples
to a1 in the final non00 letter, A1 has a term at least2^(L-2). Otherwise
the endpoints are10 and01, or the reverse. Choose a middle position i.
For endpoints10,01, its bit A_i contributes either
2^(i-ell-1) or2^(m-i-1). At a middle position the smaller of these two
quantities is at least2^floor((L-3)/2). This yields A2; lengths
one and two already satisfy it. The reversed case is identical.

Thus there are at most256*C^4 possible tails: the first non11 letter
has three choices and each subsequent letter at most four. This gives
O(C^4) possible tails, not an unbounded experimental list.
Initial11 padding has zero additional response cost and is the only
unbounded part of the language. This is a precise finite-boundary theorem.

## A2. An effective capped positive automaton

Use the inherited state `(G1,D1,G2,D2,A_last,B_last)` and recurrence (8)
of [READY.md](../two_entry_funnels_20261103/READY.md). Cap each of the four
counts at C+1. Truncation commutes with every nonnegative addition and
multiplication in the transition, so it preserves the test for exact C.
There are at most4(C+2)^4 states. Backward reachability from accepting
states removes states with no legal terminal continuation. The only
productive cycle is the initial11 loop: any other one would violate A2.

This supplies a finite algorithm for every even C, with a polynomial
state bound. The script uses backward BFS and enumerates productive
tails. For each fixed C the rooted-class generating function in this
chart is

    sum_(accepted tails w) x^(4+|w|)/(1-x).               (A3)

The ordered ports and the unique chain recover the word, so these are
isomorphism classes of rooted tournaments, without multiplicity or
automorphism weights. Every such rooted tournament is rigid.

## A3. The complete neutral four-defect layer

For any strong rooted tournament with q=s and F-H=4, the cut theorem
gives d(d-1)<=4 and a nonzero joint defect, hence indegree d=2. Moreover
H(T-{r,a,b})<=2. Hamilton counts are odd; therefore this count is1 and
the residual tournament is transitive. Consequently the chart above
captures **every** neutral four-defect root, at every size.

Its sixteen tails and H=A*2^k+B, k>=1, are:

|Tail|A|B|Tail|A|B|
|---|---:|---:|---|---:|---:|
|00,01|9|-1|00,10|11|1|
|01,01|9|-3|01,11|15|-1|
|10,10|12|1|10,11|16|1|
|01,00,10|25|5|01,01,10|34|5|
|01,10,01|30|5|01,10,10|30|7|
|01,11,10|34|5|10,00,01|22|5|
|10,01,01|25|5|10,01,10|33|9|
|10,10,01|33|7|10,11,01|31|5|

Thus the complete class series is `(6x^6+10x^7)/(1-x)`. Distinct tails
can have identical H functions: `01,01,10` and `01,11,10` both have
H=34*2^k+5. Keeping H and the growth clock would merge different actual
rooted classes. The full word retains their incidence.

At C=2 this recovers the inherited five families. The same exact compiler
gives31,44,64,84 tails at C=6,8,10,12, respectively. These higher lists
classify the **transitive-residual two-in-port chart**, not all neutral
roots: at C6 a third in-port or a nontransitive residual can occur.

## A4. Uniform count and transport formulas

Let W consist of a,b and the tail, J={a,b}. Prefixing by11^k gives

    H = 2^k E_J,W(2)+H(W)-C,
    s=q=H(W),  F=H+C,
    E_r(1)=E_J,W(1),
    S_r(2)=2^k P_W(2)-S_J,W(2).                         (A4)

The ordered two-bin formulas follow by deciding which bin contains each
prefix vertex; in the joint C term all prefix vertices must be before
the terminal in-port. They also follow from the inherited polarized
path-bin calculus. The small W quantities, not fitted samples, give
every A,B in the table. The original neutral-root birth R prepends11,
preserves C and gives H_next=2H-B. Therefore this is a closed growth
library at every fixed C.

The golden observer `(H,(S(2)-E(1))/2,(S(2)+E(1))/2)` is affine in2^k.
Its mod7 phase has period dividing3, including all five Q4 phase images.
This projection omits the word and the full joint polynomial carrier;
it is sufficient only for the specified phase queries. Actual native
complements still require their literal path guard and source relation.

## A5. The next layer already reaches prime ready-image obstructions

O3's complete inverse excludes the attained prime values41,43,67,83,101
from both ready lanes. The chart nevertheless supplies explicit infinite
families, with the same verified finite-core formula:

|C|Tail|H(k), k>=1|Examples|
|---:|---|---|---|
|6|00,00,10|19*2^k+3|41,79,155|
|6|10,00,10|20*2^k+3|43,83,163|
|6|01,00,00,10|45*2^k+11|101,191,371|
|4|01,10,10|30*2^k+7|67,127,247|

No infinitude of primes in these exponential families is asserted.
Among positive neutral defects, the minimum at H35 and H67 is4; the
minimum at H41,H43,H83,H101 is6. To check the latter statement, C2 is
excluded by the ready inverse, and none of the sixteen C4 formulas has
those four values. At k>=4 all C4 values exceed101, so that check has
only three k values per tail. The displayed C6 words attain the bound.
The word "positive" matters: zero-defect guards are a different class.

Thus enlarging the interface is relevant to the prime lane itself,
not just to composite counts such as35. A reusable nonzero-offset
operation on these specific words would be a concrete next advance.
