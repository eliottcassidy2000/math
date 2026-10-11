# The four-point interface changes with the predicate and the scale

2026-10-11. **PROVED:** the elementary geometry, palette interfaces and
restricted six-colour certificate below. **FINITE-EXACT:** the declared
controls. **OPEN:** a universal plane six-colouring or seven-chromatic
unit-distance witness. These small gadgets do not decide that alternative.

## Inheritance and the three kinds of missing-pair bit

The user-specified four-arc tournament is recovered in
[the partial-four packet](partial_four_tournament_20261011c.md), from the
earlier [four-state section](tournament_recursive_four_state_20261004.md).
Its underlying graph is C4, with its two diagonals missing. Completing each
orientation gives a two-bit section of the four tournament types.

The closest plane mechanism is
[boundary palette compression](boundary_palette_compression_20261011.md).
The hostile is assuming an orientation choice also gives a geometric edge
or a colour equality. The corrected near miss is a basis of the original
graph used after nonlinear deletion of differently labelled edges; the
[K4 hostile](plane_odd_cycle_channels_20261011b.md) rules this out. The
least-used sidecar here is the actual Gram matrix and the selected palette.
The board is **oriented completion / unit geometry / equality relation /
triangle propagation / lattice marking / geometric scale**.

There are three different two-bit questions on the same four vertices:

| Predicate | Allowed states | What the bit means |
|---|---|---|
| Partial-tournament completion | All four | Direction of each missing arc |
| Unit rhombus geometry | 00,01,10 | Whether each diagonal is a unit edge |
| Proper colouring of C4 with k=2 | 00 | Whether each diagonal has unequal colours |
| Proper colouring of C4 with k=3 | 00,01,10 | Same colour-difference question |
| Proper colouring of C4 with k>=4 | All four | Same colour-difference question |

The equality of two three-state sets is not a map identifying their meanings.
In a unit rhombus with one unit diagonal, that diagonal must have different
colours. For three colours the opposite diagonal must then have equal colours.
The resulting compatibility rule is geometric plus combinatorial, rather
than an identification of the two bitstrings.

## Unit geometry couples the diagonals

Let distinct A,B,C,D have AB=BC=CD=DA=1. Both B and D are intersections of
the two unit circles about A and C. Since B and D are distinct, they are
the two intersections and are reflections in the line AC. Their midpoint
is also the midpoint of AC, so the four points form a nondegenerate rhombus.
Its perpendicular diagonals have lengths p,q satisfying

\[
                         p^2+q^2=4.                             \tag{1}
\]

Thus both diagonals cannot be unit. A square realizes neither; rhombi with
angles60 or120 degrees realize the two single-unit cases, with other
diagonal length sqrt3. The script verifies all three using exact algebraic
coordinates, rejecting coincident vertices. The orientation completion
cube therefore cannot be transported as an independent edge-presence cube
inside this unit framework.

## Palette three is rigid; palette six is not

In a proper colouring of C4, if both opposite pairs have distinct colours,
all four vertices have distinct colours. Hence at palette three at least
one opposite pair must agree. At palette two both must agree, while four
colours realize all states.

For a diamond K4 minus one edge, let its two nonadjacent apices be A,C and
the shared edge be BD. Its exact labelled chromatic count is

\[
                         k(k-1)(k-2)^2.                          \tag{2}
\]

Choose distinct colours on B,D, then independently choose either apex
from the remaining k-2 colours. At k=3 the apices are forced equal. At
k>=4 they can agree or differ. In particular the elementary three-colour
equality gadget is not a six-colour equality gadget. This is a precise
micro-to-macro boundary: each additional triangle has one available third
colour at palette three, but four at palette six.

Also, a C4 alone never rejects a candidate in the previous ternary-label
six-colour algorithm: deleting edges between different labels leaves a
subgraph of a bipartite graph. Every ternary assignment therefore lifts.
Its apparent three-state proper-colour rule is not a constraint on those
ternary labels, because adjacent labels are allowed to agree when their
binary bits differ.

## Triangle patches, exact propagation, and a positive six-colour family

Consider a finite unit-distance graph whose vertices are covered by a
nonempty collection of triangles, connected under sharing an entire edge.
Fix one root triangle. Its vertices lie in a unit triangular lattice L.
Across a known unit edge, either possible third equilateral vertex also
lies in L: rotation by60 degrees preserves L. Induction along the connected
triangle collection puts every vertex in L.

In coordinates a u+b v, where |u|=|v|=1 and u dot v=1/2, squared distance is
Q(a,b)=a^2+ab+b^2. Since Q(a,b)=(a-b)^2 modulo3, every unit edge changes
(a-b) modulo3. This gives a proper colouring. Once the three colours on
the root triangle are chosen, each shared-edge extension forces the new
vertex's colour. Consequently the patch has exactly six labelled proper
three-colourings, one orbit under S3. Any additional unit edges among its
vertices remain proper under the same lattice formula.

This uses only the norm-one triangular-lattice argument of
[THM-418, lattice chromatic dichotomy](../../01-canon/theorems/THM-418-lattice-chromatic-dichotomy-of-the-plane.md).
It does not inherit the broader lattice claims suggested by that file's
historical title or consequences. Triangles touching only at a vertex can
rotate independently, and are outside the propagation hypothesis.

**Marked two-carrier theorem.** Suppose a finite point set is contained in
the union of two isometric copies of the unit triangular lattice. Supply
origins o_j, bases u_j,v_j with the above Gram matrix, and for each point
one address (j,a,b) of exact integers satisfying

\[
                 P=o_j+a u_j+b v_j.
\]

Then its entire induced unit-distance graph has the explicit colouring

\[
                 c(P)=3j+(a-b\bmod3),\qquad j\in\{0,1\}.        \tag{3}
\]

Same-carrier edges are proper by the norm congruence. Different-carrier
edges use disjoint palettes. Points in both carriers receive just one
supplied address; no uniqueness of address is assumed. This is a directly
checkable sufficient family for six colours, independent of girth.

The implementation authenticates the Gram matrices, every source coordinate
identity and every actual unit edge, then exports the full colouring. An
exact two-carrier representation of the Moser spindle supplies a positive
control. Its three-colour obstruction is checked independently; (3) is a
six-colouring certificate, not an optimal colouring claim.

The carrier hypothesis is substantial. The three points
(0,0),(1/3,0),(0,1/3) are pairwise closer than one, so at most one belongs
to any single isometric unit triangular lattice. They cannot be covered
by two such carriers, although their unit graph is empty. Thus a universal
two-unit-carrier cover is false, and this sufficient certificate is far
from necessary. Allowing scaled triangular lattices removes this particular
spacing obstruction but requires the appropriate norm/3-adic decoder;
the present implementation certifies only the explicitly stated unit scale.

## Next geometric target and verification

The exact positive family suggests combining geometric carrier certificates
with the prior complete boundary relations. A fragment can use its known
lattice colouring while unstructured fragments retain their full relation.
To prove a universal bound one must still show that the interfaces admit
compatible six-colour choices; counts or orientation types do not supply
those choices. The close-point obstruction also links to the previous
[high-girth crowding theorem](plane_girth_crowding_20261011b.md): it prevents
using a fixed collection of unit-scale carriers as a universal shortcut.

```text
python -B -X utf8 04-computation/experiments/unit_rhombus_interfaces_20261011c.py
python -O -B -X utf8 04-computation/experiments/unit_rhombus_interfaces_20261011c.py
```

The finite universe includes all C4 label assignments for palettes2..7,
all81 ternary assignments, every three-colouring of five triangular strips
of orders3..7, exact rhombus/spindle coordinates,121 signed integer address
controls, and malformed geometric/address inputs. There are **243 exact
checks**, including complete enumeration within those scopes. Normal and
optimized outputs agree; [saved output](unit_rhombus_interfaces_20261011c.out).
