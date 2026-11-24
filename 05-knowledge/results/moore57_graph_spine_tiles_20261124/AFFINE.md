# Why the affine reversal gives girth five, and produces graph spines

In the compact rule use U_(i,x), W_(j,y). The reversed cross arrows encode the graph contacts y=x+ij in F5. The within-cohort graph edges have steps ±1 on U and ±2 on W.

- Two vertices in different U cohorts can have a common W contact only at the unique slope j solving x+ij=x'+i'j. The symmetric statement holds for different W cohorts. This rules out a cross four-cycle; it is the same uniqueness as in an affine incidence table.
- Two vertices in the same U cohort, or same W cohort, get their adjacency or unique midpoint inside that five-cycle. They cannot share a cross contact (cross maps are matchings).
- For a mixed pair define δ=y−x−ij. If δ=0 it is a contact, with no common neighbor. Otherwise exactly one of δ∈{±1}, δ∈{±2} holds. In the first case the unique common neighbor is U_(i,x+δ); in the second it is W_(j,x+ij). Thus **the two sides partition the nonzero residuals**, with no collision or gap.

These checks give A²+A=(d−1)I+J. The pentagonal cohorts supply girth exactly five. Taking one cohort per side instead gives the degree-three Petersen graph. Taking all five gives d=5+2=7.

## Why not reuse the same recipe with a bigger finite field?

For the analogous symmetric scaffold over F_q with a cohorts per side, require regular symmetric Cayley steps within each cohort. On every pair of mixed cohorts the two step sets must partition F_q excluding zero for exact unique routes. Regularity across sides forces each half to have size (q−1)/2. Inside an individual cohort, the exact Moore condition must already be met internally, because two vertices in that cohort cannot share a cross neighbor. Thus q=1+((q−1)/2)², so q=1 or 5. If there are at least two distinct cohorts per side, the unique intersection slopes of different cohorts sweep F_q as residues vary, so all slopes must be present: a=q. Hence this specific symmetric scaffold (unmodified Cayley cohorts and complete affine matchings) gives precisely the Petersen and Hoffman–Singleton nontrivial cases, not a degree-57 recipe by changing 5 to a larger prime.

The full Hoffman–Singleton restriction comes from a different use of girth five. A diameter-two d-regular graph with no triangles or squares has n=d²+1 and nonprincipal eigenvalues (-1±√(4d−3))/2 by the Moore identity. If the discriminant root is irrational, its two eigenvalues have equal multiplicities and the trace forces d=2. If the root is an odd integer s, the trace forces s to divide d(d−2). Since 16d(d−2)=(s²+3)(s²−5), s divides 15. Its positive possibilities 1,3,5,15 give the trivial d=1 and the degrees 3,7,57. The local affine argument explains the construction at 7; the integrality argument explains why 57 is on the list despite that construction not extending to it. See the [original Hoffman–Singleton paper](https://doi.org/10.1147/rd.45.0497).

## The same cross rule builds an explicit Hamilton cycle

Visit the whole cohorts in order U_0,W_0,U_1,W_1,...,U_4,W_4. Enter U_i at a_i, traverse it by +1, and enter W_i at b_i=a_i−1+i². Traverse W_i by +2. The next affine contact gives

    a_(i+1) = a_i−3−i  (mod 5),

including the transition from i=4 back to i=0. The sum around the indices is −25=0, so this is a spanning cycle. Starting at a_0=0 gives entries 0,2,3,3,2. Breaking any edge gives a graph spine. In Petersen one can likewise visit both five-cycles in turn and cross at their matching; the resulting path does not close.

More generally, choose two lists of cohort indices u_i,v_i, each a permutation of F5, and visit U_(u_i),W_(v_i) alternately. Choose signs on the steps ±1, ±2 inside every cohort and choose the first residue. Every such choice gives an oriented Hamilton path that traverses each cohort in one block. There are exactly 2·(5!)²·5·2^10 = **147,456,000** oriented paths in this *restricted family*, counting either half first.

Closure of one starting in U is a finite-field area test:

    Σ_i (u_i−u_(i+1)) v_i = Σ_i (ε_i+2η_i) mod 5,

with ε_i,η_i=±1 and u_5=u_0. The phase sum on the right has 204 sign assignments at 0 and 205 at each nonzero residue (small convolution, or multiply the two cosine factors at a nontrivial fifth root to get −1). The area counts are 6400,2000,2000,2000,2000. To see the zero count without searching all order pairs, normalize the first two U indices to 0,1; the six permutations of the other three indices show that the cyclic edge differences are either constant or of multiplicities 2,2,1. There are 20 constant cases out of 120 U orders (the arithmetic progressions), each allowing all 120 W orders. Otherwise subtract the singleton coefficient from all five coefficients, without changing the area, and rescale: the two repeated coefficients are opposite because the coefficients sum to zero. Thus the area-zero test has coefficients 1,1,−1,−1,0. If the value assigned to coefficient zero is z, the equal-sum pairs must be {z−1,z+1} and {z−2,z+2}; there are 8 orders after the choice of z. That is 40 zero-area W orders in each of the other 100 U orders. Hence 20·120+100·40=6400, and scaling W values makes the nonzero counts equal. Divide by ten recorded starts (five U-block entry points in each of two directions). This yields **1,472,800 undirected Hamilton cycles in the whole-cohort family**, not the total for the graph.

For the displayed all-positive spine, reflection interchanges U and W cohorts and reverses their slot numbers. Within-cohort edges are preserved. On each pair of opposite-type cohorts, compare the old and reflected cross contact as equations in the two physical slot numbers t,s. Their coefficient rows are (-1,2) and (-2,1) modulo 5, so their determinant is 3. Exactly one of the five contacts agrees on every one of 25 blocks. Hence reflection changes precisely 2·(125−25)=200 graph-chord positions. This is an exact black-defect law for this spine, stronger than merely proving it is not blue.

## The entire alternating-block family stays outside the blue class locus

An exact finite extension checks all the alternating whole-cohort spines, not all possible Hoffman–Singleton graph Hamilton paths. Use the graph relabelings over F5

    i'=a i+b,    j'=c j+d,    s=ac=±1,
    x'=s x−a d i+e,           y'=s y+b c j+b d+e.

Indeed y'−x'−i'j'=s(y−x−ij), and the intra-cohort steps are preserved. Start in U (reversal handles the other start). Normalize the first two U cohort indices to 0,1, the first W index to 0, and the second W index to either 1 or 2 by the sign freedom; normalize the entry residue to zero. Then just 72 cohort-order pairs and 1024 direction choices remain.

In all 73,728 representatives, the sparse tournament has no anti-isomorphism. The short [exact checker](hs-block-full.cpp) uses singleton score classes: an anti-isomorphism must send a unique score z to the unique score 49−z, yet some edge on forced vertices does not reverse. Thus there is no incident blue line even at a tournament reroot of one of these classes. The exact reflection Hamming defects within this family range from **124 to 204** (bounds attained). Neither claim exhausts all Hoffman–Singleton Hamilton spines.
