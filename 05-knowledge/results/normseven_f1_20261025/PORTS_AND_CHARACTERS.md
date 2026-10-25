# The norm-seven repair changes a cubic observation into a pair observation

2026-10-25. **PROVED**, explicit eight-state calculation. This integrates the
incoming `1f804be357` projective-eight and complete-boundary results. Their
metric and palette conclusions are inherited, not strengthened here.

The D-side change in K3 replaces the character masks (1,2,4) by (1,2,5).
At a fixed branch-free source choose a root in each coordinate as reference,
and write the other choices by signs x=(x1,x2,x3) in F2^3. Normalize the new
reference roots using the actual formula K4, including b5*b6. The induced
map on sign choices is exactly

    y=(x1,x2,x1+x3),   x=(y1,y2,y1+y3).                      (P1)

This preserves the full eight-state fibre and its affine sign action. But
the old cubic character becomes a new pair character:

    (-1)^(x1+x2+x3)=(-1)^(y2+y3).                           (P2)

For a precise hostile, let R_epsilon={x:x1+x2+x3=epsilon}. Both parity
relations have identical counting marginals on every proper subset of the
old three coordinates. Every old pair assignment has exactly one completion.
After P1 they are the relations y2+y3=epsilon, with y1 free. Their counting
marginals on the new pair (y2,y3) are now different: two completions on the
required parity and zero on the other parity.

Nothing has recovered a missing old measurement for free. A new-coordinate
pair query may read the old combined bit x1+x3 (a product of signs); it is
a different query family.
To retain the original ports, transport their observables as

    x1=y1, x2=y2, x3=y1+y3,

rather than declaring the three y coordinates to be the same ports. The
cube metric exposes the same distinction: the old edge displacement100
becomes101, a distance-two displacement. P1 is an affine-design symmetry
but not a cube-graph automorphism. Incoming work proves that the marked
three weight-one modes are actual divisor relations and unit-cube coordinates;
the seven unmarked affine directions do not recover them.

The alternative all-rank Laurent repair in section3 of `KUMMER_AND_F1.md`
keeps every individual binary character: Q=C+2D, w'_j=w_j b^D_j. It changes
numerical coordinates by known source units, rather than changing which
character is called a port. Both repairs recover the full arithmetic cover,
but only the latter retains the original character names without a linear
relabeling. Neither preserves Euclidean lengths automatically.

Source: two marked root bases and the parity boundary tables. Map: K4 and
its induced shear P1. Preserved: the entire joint fibre, character action,
and point reconstruction. Destroyed unless carried: the old proper-coordinate
query class and cube adjacency. Sidecar: transported port observables and
the original metric. Decisive tests: the parity marginals and edge100 ->101.

More generally GL3(F2) is transitive on nonzero characters (extend one
nonzero vector to a basis). Thus interaction order1/2/3 is not intrinsic to
an unmarked seven-character set. The decomposition3+3+1 becomes meaningful
only after choosing the three coordinate observables. No new six-colour
forcing context follows from a freely chosen arithmetic basis.
