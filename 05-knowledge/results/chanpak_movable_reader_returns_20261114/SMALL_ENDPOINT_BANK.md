# Complete small-endpoint bank for strong ready sources

Let T be a strong tournament with at least two vertices, rooted at v,
and suppose its unsigned rooted response equals its ordinary count, R_+=H.
Let E be the number of ordinary Hamiltonian paths ending at v.
These are the low-ratio sources for the movable-reader chart. The exact
supply of small E is more restricted than the inequality H≥3E suggests.

## Full classification for E=1,3,5

Use the [ready-source classification](https://github.com/eliottcassidy2000/math/blob/2d70eb9cd31437bbfa8cf418c6285f79fc8afaaa/05-knowledge/results/chanpak_11_13_26_ready_barrier_20261114/READY_CLASSIFICATION.md).
The root has exactly one outneighbor u. Deleting v leaves A→C, where A
is the strong component of u and C is nonempty. There cannot be an earlier
component, because the root only beats u and must reach all vertices.
The endpoint count E is H(A)H(C), and
\[
 H=E+R_+(A,u)F_2(C).
\]
Indeed, every Hamiltonian path of A→C ends in C, so it can be followed
by v; any full path ending at v arises uniquely that way.

The only strong components of counts 1,3,5 are, respectively, the
singleton, C3, and the unique strongly connected four-vertex class.
An elementary check avoids a large-order census: the odd-cycle packing
formula bounds the number of directed triangles by (H−1)/2. In a
nonsingleton strong tournament every vertex lies on a triangle (if
there were no arc outneighbor→inneighbor then the vertex cannot reach
its inneighbors). One triangle thus forces order three. For two
disjoint triangles, each outside vertex sees each triangle uniformly,
since a mixed view creates another triangle, so all cross arcs between
the two triangles have one direction, contradicting strength. For
two triangles meeting in exactly one vertex, each of the two external
vertices of the second sees the first uniformly; the second triangle
then copies at all three vertices of the first, a contradiction.
With at most two triangles the strong order is therefore at most four.
At order four uniqueness follows by orienting the remaining vertex
against a triangle: it must beat exactly one or two vertices; these
two cases are isomorphic (the four-vertex degree profile is 1,1,2,2,
and the unique edge within each same-degree pair fixes the labels).
For the four-vertex class, use outgoing masks (8,1,3,6); it has
F_2=32 and root-start response menu {5,6,7,9}, by direct subset sums.
The 3-cycle has F_2=12 and each root-start response is 3.

As E is prime for E=3,5, either A is singleton or H(C)=1, and the
nontrivial count in C, if any, comes from one component with optional
ordinal singleton padding. Thus the possible full path counts are exactly:

| E | possible H |
|---|---|
| 1 | 1 + 2^r, r≥1 |
| 3 | 3 + 3·2^r, r≥1 (includes the C3-in-C case) |
| 5 | 5 + D·2^r, D∈{5,6,7,9}, r≥1; or 5 + 32·2^r, r≥0 |

Every listed count is physical by taking the appropriate core as A
or C and any transitive padding indicated by the exponent.

## Intersections: a finite supply, proved at all orders

* E=1 and E=3 intersect only at H=9. In 2^n−3·2^r=2, necessarily
  r=1 (n≥2), hence n=3.
* E=3 and E=5 intersect only at H=15. In 3·2^n−D·2^r=2, necessarily
  r=1 because n≥2 by the H≥15 bound. Then D+1 is 3·2^(n−1);
  among 6,7,8,10, only 6 works. The separate 32·2^r branch has
  no solution by reduction mod 4.
* E=1 and E=5 intersect only at H=17,33. Since n≥4,
  2^n−D·2^r=4 imposes r+v_2(D)=2; then D/2^v_2(D)+1 must be a
  power of two. D=6 or 7 yields precisely these witnesses. The
  separate 32·2^r branch has no solution mod 8.

So the smallest distinct source endpoints E=1 and E=3 can repair
two counts with the same older reader at H=9, but changing the
strong source order cannot make the repair available for other H.
The next unit shift, 3→5, is likewise confined to H=15. The natural
way to escape this obstruction is to let the endpoint intercepts
grow with the baseline, rather than search for more small-endpoint
sources at large order.

## For a growing E, a unit shift forces a singleton terminal component

More generally, take any two of these strong rooted sources with the
same H and end counts E, E+2. Exactly one of the even quantities
H−E and H−(E+2) has 2-adic valuation 1. In that source C must consist
of exactly one vertex, and R_+(A,u) must be odd. Indeed F_2(singleton)=2,
but F_2(C) is divisible by 4 for every |C|≥2: the odd-cycle packing
terms after two-color evaluation have weight 2^(|C|−s+2k) when k
nontrivial odd cycles occupy s vertices. For k=0 the exponent is
|C|≥2; otherwise it is at least 2. Thus at least one side of **every**
unit endpoint repair at fixed H is a single-terminal construction with
H(A)=E on that side and H=H(A)+2R_+(A,u). The other side must have
4 dividing the corresponding product. This constraint does not extend
to non-ready sources where R_+ differs from H.
