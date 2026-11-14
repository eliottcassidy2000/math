# Two exact boundaries on growing the seven-vertex exterior

The [portable skew criterion](PAINTED_DOMINO.md) gives a larger
exterior C with rows (14,28,104,112,37,67,19), tester masks U=9,
V=68, and order-three automorphism p=(2,3,5,6,4,0,1) in vertex-image
notation. In particular p fixes vertex 4 and sends the first of the two
tester masks to the second. Its bracket in equation (1) is 2.
Here are two ways to keep the symmetry but lose the unit step.

## Replacing the fixed point is an all-order obstruction

Replace vertex 4 by any nonempty tournament Z. Its vertices are
outside both tester masks. The permutation p extends by the
identity on Z. The skew of the inflated exterior is exactly

    2H(Z)+17 C_2(Z)+12 C_3(Z).

A pair of nonempty covering paths of the inflated exterior
projects to two paths through the old symbols, where Z can be
visited repeatedly. It can use at most eight runs in total:
each of six other vertices can separate runs, and each of the
two paths can have an initial run. For j runs there are C_j(Z)
ways to fill the ordered nonempty slots, with each internal
path's directions built literally inside Z. Thus the bracket
is some linear functional in C_1 through C_8, independently of
size of Z. To compute its coefficients one may evaluate it for
Z=TT_n, n=1,...,8: the values C_j=j! S(n,j) form an invertible
lower-triangular evaluation matrix. Exact DP returns the
coefficient list (2,17,12,0,0,0,0,0) in exterior_fixed_point.py.
This is complete by the run bound, not a fit extrapolated past
an unknown degree. For every nonsingleton Z, its C_2 is at least
2H(Z), so the skew is at least 36; a unit step is impossible.
The issue is not loss of the mask automorphism.

## Adding one genuinely new symmetry orbit

Preserve the old seven labelled vertices, and adjoin a directed
three-cycle on 7,8,9 on which p extends as 7->8->9->7.
There are exactly 256 labelled equivariant orientations:
the arcs between the new cycle and each of the two old
3-orbits have three independent directions, the arcs to the
old fixed point have one, and the new internal cycle has two
orientations. Exactly 252 are strong (every cross assignment
except the two uniform choices for each internal orientation).
For every one of these 252 and every choice of tester mask U,
test V=p(U) by the two-cover endpoint matrix. None gives
a skew of +2 or -2. See equivariant_extension.cpp. This is a
finite exact boundary for this extension prescription, not an
obstruction to other cores, other symmetries, or more than one
added orbit. It shows why simply keeping the old rotation
cannot serve as the recursive selector.
