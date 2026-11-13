# A field-with-one-element style jet that can see a wedge

The older [F1 geometry session](../session59_f1_geometry_and_information_20261025/SYNTHESIS.md)
uses projective phases of a cyclic unit group on a labelled support;
[the carousel work](../module_shattering_lines_carousel_finite_field_20261026.md)
separately retains Frobenius and orbit data. Neither constructs a
literal field with one element. Here is a precise point-count
specialization using the first of those mechanisms.

Take an unordered spanning cover by nonempty directed paths. Decorate
the vertices of each path by units of F_q, up to common multiplication
by a unit **within that path**. A path on r vertices has (q−1)^(r−1)
classes; pick one named vertex to normalize to 1. Distinct disjoint
supports identify the paths, so passing from ordered to unordered
nonempty slots contributes exactly a division by k! for k paths.
There is no Frobenius quotient implicit in this division.

Mark each undirected path edge e by a variable y_e. For t=q−1 the
point-count polynomial of these configurations is
```
J_T(t;y) = sum_C t^(n−|C|) product_(e used in C) y_e.
```
The sum is over unordered nonempty path covers. The exponent is
literally the number of used edges. At q=1 this is formal
specialization, not an actual finite-field point set.

- The constant term is 1 (all singleton paths).
- The linear term is t sum_e y_e, independent of every orientation.
- For two distinct edges, the coefficient of t² y_e y_f is always 1
  if the edges are disjoint. If they share a vertex, it is 1 exactly
  when they form a directed path through that vertex, otherwise 0.
  All unmarked vertices must be singleton paths.

Therefore the second-order **labelled** jet already reconstructs
the tournament up to converse, by [wedge integration](WEDGE_RECONSTRUCTION.md).
The whole zeroth/first-order jet is universally blind. If the edge
variables are forgotten, even the full series is merely the number of
unordered covers by each number of paths, read from the slot polynomial,
and is blind to reversing any module. The field-style limit itself does
not do the recovery: retaining two labelled contacts does.

This construction uses only the unit group; an abstract cyclic phase
group of order q−1 gives the same count. The analogy stops at the
formal q→1 specialization and at quotienting phase symmetries. The
q in the unit field and the z counting optional slots are distinct
parameters; the F1 join-family name used in Chan–Pak work is distinct
from both.
