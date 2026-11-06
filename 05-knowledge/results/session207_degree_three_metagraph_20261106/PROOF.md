# Exact degree-labelled incidence at every H=3 source

Status: PROVED, using the inherited odd-cycle collection formula. No novelty or global spectrum-completeness claim. All tournaments are nonempty and H counts directed Hamilton vertex words.

## 1. Complete source classification

A tournament has H=3 if and only if it is

    T_ab = Tr_a -> C3 -> Tr_b,  a,b>=0.

The inherited odd-cycle collection formula (THM-002-ocf.md) states H(T)=sum_C 2^{|C|}, over collections of vertex-disjoint directed odd cycles of length at least three, including the empty collection. Thus H=3 forces exactly one directed odd cycle, and in particular exactly one directed triangle. Conversely a unique directed triangle forces the asserted form: any vertex with mixed orientation toward that triangle would form another triangle with two of its vertices. So each exterior vertex either beats the whole triangle or loses to it. An arc from the latter set to the former would form another triangle with a triangle vertex. Each exterior set is itself triangle-free, hence transitive. The displayed ordinal join follows, and its only directed cycle is the middle triangle, giving H=3.

The parameters a,b uniquely determine the ordinary isomorphism class, and |Aut(T_ab)|=3. All three Hamilton paths are in one automorphism orbit. Consequently each source has ordinary tiling-fibre degree one, and a unique incident line in the original class multigraph. This is a classification of sources with ordinary count3; it is not a classification of all fibre-degree-one sources.

The source classification and fibre identity are inherited elementary mechanisms. The explicit opposite-endpoint formula below is the present refinement.

## 2. The exact neighbour

Put n=a+b+3 and fix the path 0->1->...->n-1 in T_ab. Its only backward arc is (a+2)->a. Let R_n preserve consecutive forward arcs and point every longer arc backward, the known tribonacci partner of the transitive tournament. The path-complement neighbour U_ab is precisely R_n with (a+2)->a reversed to a->(a+2).

Let g(a,b)=H(U_ab). Then the complete set of ordinary-count pairs incident to H=3 vertices at order n is

    {(3,g(a,n-3-a)): 0<=a<=n-3}.

This is exact support, not a multiplicity-free class classification of the target tournaments: different parameters can give isomorphic targets or equal target counts. There are n-2 ordinary source classes, each contributing one incident line. No two such sources are adjacent: U_ab is transitive if n=3, and is strong if n>=4, whereas T_ab is nonstrong for n>=4. For n>=4 the retained path closes using (n-1)->0, which is not the exceptional chord; at n=3 the complement is transitive.

## 3. Rational formula and recurrences

Define D(z)=1-z-z^2-z^3. Then

    G(x,y):=sum_(a,b>=0) g(a,b)x^a y^b = N(x,y)/(D(x)D(y)),

where

    N=1+4x+4y+3x^2+3y^2+4xy+2x^2y+2xy^2+x^2y^2.

Equivalently g obeys a tribonacci recurrence independently in each parameter:

    g(a,b)=g(a-1,b)+g(a-2,b)+g(a-3,b) for a>=3,
    g(a,b)=g(a,b-1)+g(a,b-2)+g(a,b-3) for b>=3,

with seed matrix, rows a=0,1,2 and columns b=0,1,2,

    1   5   9
    5  13  23
    9  23  41.

In particular g(a,b)=g(b,a). These are all-order recurrences obtained from cycle classification, not extrapolation from a finite table.

### Cycle classification proof

Write h_m=H(R_m), with h_0=1 for the empty tournament. Every directed cycle of R_m is exactly an interval i,i+1,...,j,i of length at least three: starting at its minimum vertex forces successive increments until the only legal jump to an already visited vertex closes the cycle. Odd-cycle collections are therefore packings of disjoint odd intervals of length at least three, each carrying weight2. Uncovered vertices have weight1. Consequently

    H(z):=sum_(m>=0) h_m z^m
         =1/(1-z-2z^3/(1-z^2))=(1-z^2)/D(z).

This recovers the inherited tribonacci partner formula (THM-4731-tiling-line-metagraph-self-lines-by-twisted-burnside-sigma-twins-depth-and-the-tribonacci-partner.md).

Reversing (a+2)->a removes only the interval triangle [a,a+2] from the directed-cycle inventory. Every new cycle must use a->a+2. Starting at its minimum i, it proceeds consecutively to a, skips a+1, then proceeds consecutively to some j>=a+2. It then either:

A. closes directly to i, producing [i,j] minus {a+1}, where i<=a and j-i>=3 (at least three cycle vertices); or
B. jumps to the skipped a+1 and then closes to i, producing the full interval [i,j], where i<=a-1 and j>=a+3.

These are the only possibilities: any other downward jump hits an already visited vertex other than the minimum; after visiting a+1 its only possible unvisited continuation is closure to the minimum. Type A is odd exactly when j-i is odd, and type B exactly when j-i is even. Existing interval cycles on the same support as a type B cycle remain distinct directed cycles and are counted separately.

Every changed cycle contains both a and a+2, so a disjoint collection contains at most one changed cycle. Other cycles in such a collection lie wholly before i or after j; the single omitted a+1 in type A cannot support any cycle. The remaining weights therefore factor as h_i h_(n-1-j). Set p=a-i, q=b-(n-1-j). Subtracting collections using the removed triangle and adding those using the new cycles gives

    g(a,b)=h_(a+b+3)-2h_a h_b
       +2 sum_(0<=p<=a,0<=q<=b,p+q odd) h_(a-p)h_(b-q)
       +2 sum_(1<=p<=a,1<=q<=b,p+q even) h_(a-p)h_(b-q).

The generating functions of the two index constraints are respectively

    (x+y)/((1-x^2)(1-y^2)),
    xy(1+xy)/((1-x^2)(1-y^2)).

For the baseline term set F(z)=sum_(m>=0)h_(m+3)z^m=(3+2z+z^2)/D(z). Its bivariate generating function is (xF(x)-yF(y))/(x-y), interpreted as a formal divided difference (so no difficulty at x=y). Therefore

    G=(xF(x)-yF(y))/(x-y)
       +2H(x)H(y)[-1+(x+y+xy+x^2y^2)/((1-x^2)(1-y^2))].

Simplification gives the displayed numerator N. Multiplying by D(x)D(y) proves the two recurrences and the seed matrix.

## 4. Actual fibre degrees, colours and an aggregate recurrence

Every U_ab is rigid. For n>=5, the scores of R_n are (1,1,2,...,n-2,n-2), and the flip adds1 at a and subtracts1 at a+2. If 1<=a<=n-4, the resulting equal-score triple {a,a+1,a+2} is transitive, and all other tied cells are the remaining endpoint pairs, also transitive. For a=0 or a=n-3, the modified cells and the untouched endpoint pair are likewise transitive pairs or singletons. At n3 the target is transitive; at n4 each target's two score cells are directed pairs. An automorphism preserves scores and hence each such transitive cell, whose automorphism group is trivial. Thus it fixes every vertex.

It follows that g(a,b) is also the target's ordinary tiling-fibre degree H/|Aut|. The source's ordinary fibre degree is1. Hence this gives an explicit degree-labelled row of the original line multigraph at every H=3 source, not merely a Hamilton-count relation.

Converse sends T_ab to T_ba. Thus in the merged graph:

- a=b gives one self-converse source with one blue line;
- a<b gives one non-self-converse source with two black lines, both ending at the same merged target.

These use the inherited reflection/twin rules from THM-643 and THM-4731. The merged target's fibre is g(a,b) or2g(a,b), according to its own converse status; no classification of that status is claimed.

The sum of target ordinary degrees over all n-2 ordinary source classes, counted once per source even when targets coincide, has generating function

    sum_(m>=0) [sum_(a=0)^m g(a,m-a)] z^m
      =G(z,z)=(1+8z+10z^2+4z^3+z^4)/(1-z-z^2-z^3)^2.

This is an exact size-indexed incidence statistic, with a linear recurrence of order at most six. It must not be interpreted as a sum over distinct target classes.

## 5. Exact count-quotient edge multiplicity and failure of scalar lumping

Let E_n(3,k) be the number of unordered edges in the labelled path-complement graph joining a tournament with H=3 to one with H=k. Then for n>=3,

    E_n(3,k)=n! * #{a in {0,...,n-3}: g(a,n-3-a)=k}.

There is no double counting: no edge joins two H=3 vertices. Each source class T_ab has n!/3 labelled realizations, and each has three distinct neighbours, all isomorphic to U_ab. Therefore it contributes exactly n! labelled edges, independently of possible target-class coincidences. In the original class line multigraph the analogous edge count is simply the cardinality on the right without n!. This gives an exact degree-labelled incidence row, not merely possible pairs.

The ordinary-count partition is already not equitable at order5. The sources T_0,2 and T_1,1 both have H=3, but all three labelled neighbours of the former have H=9, while all three neighbours of the latter have H=13. Consequently the walk that chooses a uniformly random Hamilton path and complements its other arcs is not strongly lumpable under H alone. Knowing order and count, or the fixed-order count spectrum, does not determine the next count distribution. This is a concrete obstruction at a fixed size, distinct from the previous across-size unbounded-neighbour theorem. It does not preclude a specially initialized weakly lumpable process or a richer quotient retaining additional coordinates.

## 6. Scope and provenance

The ordinary-count pair support above is complete for first endpoint H=3, at every order n>=3. It does not determine all pairs with arbitrary H, the entire metagraph, or whether every positive odd value other than7 and21 is attained. Single-arc-flip graphs and vertex-extension graphs remain different objects.

Dependencies inspected from cached origin/main on 2026-11-06: THM-002-ocf.md (odd-cycle formula), THM-4731-tiling-line-metagraph-self-lines-by-twisted-burnside-sigma-twins-depth-and-the-tribonacci-partner.md (R_n and fibre conventions), and session151_spectrum_metagraph_20261103/metagraph/STRUCTURE.md (line conventions and mirror pairing). No claim of a refreshed complete literature search is made.
