# Derived corollary: homotopy type and integral homology

20 October 2026. **PROVED**, as a consequence of the independently audited
finite-field normalization theorem. This is a separate audit supplement;
it does not alter or enlarge the exact six-file producer approval.
The corollary is derived from that normalization and audited separately.
No novelty or tournament interpretation is claimed.

Let q=p^r be odd. Retain the audited notation:

```
V=(q^2-1)/2,
m=q/p,
c=q(q^2-1)/(p(p^2-1)),
g_p=1+(p^2-1)(p-6)/24.
```

The normalization of K_q is a disjoint union of c connected, closed,
orientable genus-g_p surfaces. Each original vertex identifies m points,
one from each of m distinct normalization components. Edges and faces
undergo no identifications beyond their endpoints.

## The incidence graph

Form a bipartite graph B whose nodes are the c normalization surfaces and
the V original vertices. Connect a surface to an original vertex whenever
that surface contains its unique corresponding normalized vertex. There
are exactly I=mV edges. The graph is connected: any path in K_q can be
followed through its surfaces and their intersection vertices, and K_q
is connected. Its cycle rank is therefore

```
b=I-(c+V)+1=(m-1)V-c+1.
```

In particular b is nonnegative. When q=p it is zero.

## Homotopy decomposition

Replace each original vertex identification by a star: introduce a central
point labelled by the original vertex, and join it by an interval to each
of its m normalized lifts on the disjoint surfaces. Call the resulting
finite CW complex Y. The closed stars for distinct original vertices are
pairwise disjoint contractible subcomplexes. Collapsing each star to one
point recovers exactly K_q.

Collapsing a contractible CW subcomplex is a homotopy equivalence, since
its inclusion is a cofibration. Applying this finitely many times gives
Y homotopy equivalent to K_q.

Now regard Y as a graph of spaces over B, with the surfaces and points
as vertex spaces and single intervals as edge cylinders. Choose one
basepoint in each surface. Each interval endpoint on a connected surface
can be moved to that basepoint along a path without changing the homotopy
type of the attachment. After all endpoints are so moved, contract a
spanning tree of B. The remaining I-(c+V)+1=b edges are circles, and the
surface vertex spaces are wedged at their basepoints. Consequently

```
K_q is homotopy equivalent to
    (wedge of c genus-g_p closed orientable surfaces)
    wedge (wedge of b circles).
```

This is a homotopy equivalence, not a homeomorphism. In particular it does
not remove the singular vertices from K_q or turn K_q into one surface.

## Integral homology and examples

The reduced homology of a finite wedge is the direct sum of that of its
summands. Thus

```
H_0(K_q; Z)=Z,
H_1(K_q; Z)=Z^(2*c*g_p+b),
H_2(K_q; Z)=Z^c,
H_i(K_q; Z)=0 for i>=3.
```

There is no integral homology torsion. These are uniform consequences of
the topological decomposition, not inferences from finite-field rank tests.
The fundamental group is correspondingly the free product of c genus-g_p
surface groups and a free group of rank b.

For F9, V=40, m=3, c=30, g_p=0, so b=51. Hence H_1=Z^51 and H_2=Z^30.
For F25, V=312, m=5, c=130, g_p=0, so b=1119. Hence H_1=Z^1119 and
H_2=Z^130. Their Euler characteristics are respectively -20 and -988,
agreeing with the original simplicial counts.

As an independent finite corroboration, `independent_surfaces.py` computes
the actual simplicial boundary ranks over F3 at q=3,5,7,11,9. It obtains
the expected Betti vectors, including (1,51,30) at q=9. The proof above,
not that modular computation, establishes the integral and torsion claims.
