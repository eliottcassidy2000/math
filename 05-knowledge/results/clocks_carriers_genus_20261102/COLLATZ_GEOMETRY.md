# A rigorous boundary for the Collatz honeycomb picture

Status: PROVED, with exact native-clock replay and same-agent proof review.
No independent-agent or formal audit.
Use the standard map C(n)=n/2 for even n, C(n)=3n+1 for odd n, on
positive integers. Its directed functional graph has exactly one outgoing
edge per vertex and at most two incoming edges. The predecessors of m are
2m and, only if m=4 mod 6, the odd positive integer (m-1)/3.

## B1. A binary inverse grammar with bounded native time

For odd n>1 with 3 not dividing n, let K(n) consist of those k in 1..6 for
which y_k=(2^k n-1)/3 is a positive odd integer not divisible by 3.

**Lemma.** K(n) has exactly two elements, and each y_k reaches n in exactly
k+1 standard Collatz steps, with no intermediate odd vertex.

Proof: integrality requires even k when n=1 mod 3 and odd k when n=2 mod 3.
There are three such k in 1..6. Modulo9, their values 2^k n are precisely
1,4,7. The value1 gives a multiple of 3 after division, and the other two
give units. All y_k are odd. The identity 3y_k+1=2^k n proves the clock.

Start at 85, which reaches1 via256. Apply both inverse branches repeatedly.
At depth r there are exactly 2^r distinct odd vertices, all known to converge
to 1, and each is at undirected native graph distance at most7r from 85.
Distinctness follows from the unique accelerated odd successor: a repeated
label would give the same parent chain. No nontrivial cycle occurs in this
constructed basin, since every node reaches85 and then1;85 is not periodic.
The macro edges can share even halving tails. We claim a binary grammar on
odd vertices, NOT disjoint native paths or a topological binary subdivision.

## B2. No uniform native-edge dilation into a planar periodic lattice

Give the honeycomb its usual vertex graph. A coordinate model has vertices
(a,b,e) in Z^2 x{0,1}, with neighbours of (a,b,0) at
(a,b,1),(a-1,b,1),(a,b-1,1). A radius R ball has at most2(2R+1)^2
vertices. This simple upper bound is sufficient; no exact sphere count is used.

Suppose an injective labelling of the Collatz component containing1 into
this graph sends every native edge's endpoints to distance at most L.
The 2^r vertices from B1 must all lie within radius7Lr of the image of 85.
Consequently

    2^r <= 2(14Lr+1)^2                               (B1)

for every r, impossible for fixed L. More quantitatively, if L_r is the
largest lattice distance assigned to any native edge in these depth-r
certificates, then

    L_r >= [2^((r-1)/2)-1]/(14r).                     (B2)

The argument also excludes bounded-dilation injective placements in every
fixed polynomial-growth graph, including the triangular lattice of hexagon
centres. It uses only an explicit convergent subtree, not the Collatz conjecture.
Even the trivial1->4->2->1 cycle already precludes a literal edge embedding
in a bipartite honeycomb, but the growth obstruction survives contracting
that cycle and allowing uniformly bounded subdivisions.

## What can replace the picture

A functional graph is not one simple forward space-filling curve:5 and 32
both map to 16, so a forward path cannot visit both and remain injective.
A contour traversal of an ordered inverse tree can revisit vertices and
encode child choices. A dyadic interval address can likewise retain an
inverse branch word. These are useful alternative representations, with
branch order, revisits and time change explicitly stored. They do not by
themselves prove that every positive integer enters the tree rooted at 1.

Precise next problem: construct a source-labelled, variable-scale layout
of a certified inverse grammar with explicit distortion and decoder bounds.
Any route to universal coverage must separately prove membership for every
source. Drawing a connected or space-filling surrogate is not that proof.

## B3. A positive replacement: a dyadic, hyperbolic macro-layout

Number binary addresses of length r by 0<=j<2^r. Place them at

    z_(r,j)=((j+1/2)2^(-r), 2^(-r))

in the upper half-plane. Decode each address arithmetically using B1,
ordering its two exponents increasingly. This preserves the exact integer
source and stores the native clock k+1 on every macro edge.

Children have vertical scale1/2 and horizontal displacement plus/minus1/4
of the parent's height. Thus their Euclidean lengths shrink with depth.
Using the standard [upper-half-plane distance formula](https://www.maths.dur.ac.uk/users/anna.felikson/Geometry/Geometry24/notes.pdf), every parent-child distance satisfies

    cosh d=1+[(1/4)^2+(1/2)^2]/[2*(1/2)]=21/16.

So the certified binary family has an injective constant-macro-length
layout in the hyperbolic plane, while its Euclidean resolution changes
at every depth. This is a proved layout of the odd macro-grammar, not a
claim that native Collatz edges are disjoint or are edges of a particular
hexagonal tiling. The time labels cannot be dropped when returning to C.
Universal entry of arbitrary integers remains a separate OPEN obligation.
