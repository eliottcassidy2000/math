# Known bare K4,4 point: independent addendum

The newly checked primary source establishes that bare K4,4 existence is known, not an unresolved inference to be drawn from dimension. The earlier geometry proof and distinction from the decorated magic-square surface remain unchanged. The initial audit's statement that those calculations themselves do not produce a bare point must not be read as a literature-status claim or applied to this addendum.

Ajai Choudhry, *Two sets of integers such that all elements of the sumset of the two sets are perfect squares*, arXiv:2508.07806v1, 11 August 2025, introduction: https://arxiv.org/html/2508.07806v1 . Retrieved directly on 25 October 2026. The introduction records A=(324,54756,119716,264196), B=(0,79200,227205,1258560), citing Guy's book, page 267. In the same paragraph it identifies a parametric construction with Jean Lagrange. Without inspecting Guy directly, the careful attribution for the numerical witness is 'recorded by Choudhry from Guy, alongside Jean Lagrange's parametric family.' No originality is claimed.

Taking L=A-1 and R=B+1 yields

L=(323,54755,119715,264195),
R=(1,79201,227206,1258561).

These eight integers are positive and pairwise distinct. Their 16 cross-sum roots are

18, 282, 477, 1122;
234, 366, 531, 1146;
346, 446, 589, 1174;
514, 586, 701, 1234.

Exact checks verify all 16 squares, all nine X44 equations, nonzero entries, and distinct squared boundary entries. This is an explicit point of U44(Q), proving the corresponding bare positive-integer K4,4 exists.

All 16 squared entries are distinct. The specified decorated locus has at most nine distinct squared entries because of its repeated linear forms. Therefore this point is not in that decoration, even after arbitrary row and column permutations, transposition, gauge shift, or nonzero common scaling. It gives no magic square of squares or prize claim.

The separate known_point.py uses explicit exception checks and passes identically in normal and python -O modes; normal and optimized JSON outputs were compared byte-for-byte. The original geometry/check certificate scripts use assert and therefore are not themselves guaranteed under -O unless their producer hardens them. The independent audit.py was also run in both modes, with identical outputs and all 72,097 checks active. These are separate guarantees, not a claim that an optimized-mode run validates the original assertion-based scripts.
