# Physical ternary fork with a common strong source

Let W(w) consist of the transitive spine indexed left-to-right, plus apex A
with A beating precisely the 1-bits. If w begins 1 and ends 0 this is strong:
every spine vertex is reached from A through the first spine vertex, and
every spine vertex reaches A through the last. The weighted-inversion clock
(S302) starts with H=1 and (X,b)=(0,0); for each next bit e add (1-e)(X+b)
to H and replace (X,b) by (2X+e,e). Here X is the integer represented by the
bits read and b the last bit. For a repeated run, 1^k leaves H unchanged and
maps X to (X+1)2^k-1; 0^k maps X to X2^k and adds X(2^k-1)+b for k>=1.

The **same** source S_k=W(1010^k10), k>=2, has H=15*2^k+1 and offers:

* Digit 0, H'=3H-2=45*2^k+1. Insert a 0 immediately before the second
  1 and a 1 immediately after the last old 1. This gives W(10010^k110),
  with the displayed H' by the clock.
* Digit 1, H'=3H. Add one vertex to complete the first two vertices of
  the middle zero-run into a cyclic triangle; give it exactly the pair's
  orientations toward all exterior vertices. Any excursion that leaves
  the module and returns must visit A, since without A the exterior is
  ordered with the module. So at most two runs through the module are
  possible. The ordered cover profiles of TT2 and C3 have first two
  terms (1,2) and (3,6), respectively, giving exact multiplication.
* Digit 2, H'=3H+2=45*2^k+5. Add one vertex z losing to spine indices
  0 and 2 only (zero-based), beating all others including A.
  In the quotient with the middle run contracted to a vertex v, indices
  for k=1 are 0,1,2,v=3,4,5,A=6,z=7. The only out-to-in edges for v are
  A→0,A→2, which share A. A three-run path would need two disjoint
  out-to-in crossings; thus the whole output has at most two runs through
  TT_k. Its count is affine in 2^k, since the two ordered cover counts
  of TT_k are 1 and 2^k-2. Direct subset DP at k=1,2 gives 95,185,
  forcing the claimed formula.

All extensions preserve all old vertices/arcs and stay strong (the apex-word
witnesses as above, and each added vertex has an in/out neighbor in the
strong core). Thus m=(H-1)/2 maps to 3m,3m+1,3m+2 on one common class, not
merely three different number-theoretic charts. The destinations have not
been returned to a common ready class. Consequently this is a ternary
**fork**, not a ternary **compiler**, and does not prove cofinite coverage.

Companion −2 one-birth chart: U_k=W(1001^k00), H=15*2^k+3, k>=1.
A newborn losing only to spine index 1 and apex A gives H'=45*2^k+7.
Contract the middle run; at k=1 its index is 3, and the only out-to-in
edges (4,6),(5,6) share A=6. As above, affine interpolation from
97,187 verifies the formula. The k=2 source is the S302 count-63 seed,
up to the stated marked spine, and its two exact exits are 187 and 189.
A second realization of H'=45*2^k+7 inserts 0 before the repeated 1-run
and 1 between the final two zeros, giving W(10001^k010).
