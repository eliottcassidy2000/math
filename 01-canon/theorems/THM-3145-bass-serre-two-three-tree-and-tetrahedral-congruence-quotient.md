---
id: THM-3145
title: "Bass-Serre two-three tree and tetrahedral congruence quotient"
status: "PROVED + VERIFIED-EXACT + INDEPENDENTLY HOSTILE-AUDITED"
source: root/codex-thm3088-push-2026-08-02
audit: >
  An independent hostile audit rebuilt the coset graph, free rank-three
  kernel, local C3 orders, three binary sectors, regular V4 action,
  opposite-edge quotient, and Gamma(3)-fibre gate hostile.  It caught and
  repaired the initial false claim that the subdivided K4 was not a partial
  cube, plus four-point-orbit, nonabelian-kernel, primitive-cycle, and rooted
  binary scope.  The repaired companion independently verifies the isometric
  Q4 embedding and four Theta classes.  Fresh normal and optimized runs
  LF-byte-match the stored twelve-line transcript and both declared hashes.
depends_on:
  - THM-3141-quartic-v4-modular-congruence-shadow-and-gamma3-sidecar-boundary
  - THM-2056-kelvin-polar-farey-defect-certificate
related:
  - THM-2950-three-conjugate-pair-v-four-torsor-and-quartic-resolvent-frame
  - THM-2951-fifth-compound-reconstruction-and-v-four-phase-scalarization-boundary
script: 04-computation/modular_bass_serre_two_three_tetrahedral_quotient_thm3145.py
output: 05-knowledge/results/modular_bass_serre_two_three_tetrahedral_quotient_thm3145.out
hash_basis: LF-normalized bytes
script_sha256: b792294db2bbdfac1f0fad27f00516a747d21f413bf42f79701a2a03843235b4
output_sha256: 11735f4e7237f7f9bf93ecb60ce24b0ffec13d5bacbbf04ca43ce62e6dee449a
---

# THM-3145 -- Bass--Serre two-three tree and tetrahedral congruence quotient

**PROVED + VERIFIED-EXACT + INDEPENDENTLY HOSTILE-AUDITED.**

The binary and ternary tower grammars have a literal common object.  It is
not a homogeneous tree and not the finite four-point torsor: it is the
Bass--Serre `(2,3)`-biregular tree of

```text
G = C2 * C3 = <s,r | s^2=r^3=1>.                         (1)
```

The rooted binary tree and the ternary triple are two different views of
this one tree.  THM-3141's four-point `V4` torsor is the four-vertex orbit
inside its level-three finite quotient after cycle-lift data have been
forgotten.

## 1. The common `(2,3)` tree

Let `H2=<s>` and `H3=<r>`.  Define a bipartite graph `T` by

```text
V(T)=G/H2 disjoint_union G/H3,        E(T)=G,             (2)
```

where the edge labelled `g` joins `gH2` to `gH3`.  Left multiplication by
`G` acts on this graph.  Every `H2` vertex has the two incident edges `g,gs`,
and every `H3` vertex has the three incident edges `g,gr,gr^2`.  Thus the
degrees are exactly `2` and `3`.

The graph is connected by reduced-word traversal.  A nonbacktracking closed
path would give a nonempty alternating reduced word in `C2*C3` equal to the
identity, contradicting free-product normal form.  Hence `T` is a tree.  Its
edge labels are the left regular `G`-set, so the action is faithful.

## 2. Why the same object looks binary and ternary

Suppress every degree-two vertex of `T`.  The result is a `3`-regular tree.
After choosing a root of the `H3` type, the root has three outgoing branches,
while every other vertex has one parent and exactly two forward children.
Thus three rooted binary forward sectors meet at the root; choosing one
oriented root edge selects one ordinary rooted binary tree.  This is the
binary fraction-tree grammar, not a canonical identification with any one
continued-fraction labeling.

Conversely, at `gH3` the incident labels `g,gr,gr^2` have a cyclic order,
well-defined up to cyclic rotation of the representative `g`.  The conjugate
stabilizer `gH3g^-1`, not one global copy of `H3`, cycles that local star.
This is the ternary triple grammar.  Changing `r` to `r^-1` reverses the
orientation, so it is a sidecar and not part of the unpointed abstract tree.

## 3. The mod-three quotient is the tetrahedron

Use THM-3141's epimorphism

```text
pi:G -> PSL2(F3)=A4,          s |-> S,   r |-> R,          (3)
```

and put `N=ker(pi)`, the projective principal congruence subgroup
`barGamma(3)`.  Because `S` and `R` retain exact orders `2` and `3`, `N`
meets every conjugate of `H2,H3` trivially.  Therefore `N\T` is the ordinary
bipartite quotient graph

```text
vertices: A4/<S> (6)  disjoint_union  A4/<R> (4),
edges:    A4 (12),                 degrees 2 and 3.        (4)
```

Suppressing its six degree-two vertices gives `K4`.  More precisely, every
`A4/<S>` coset joins a different unordered pair of the four `A4/<R>` cosets,
and the six resulting pairs are

```text
01,02,03,12,13,23.                                      (5)
```

Thus the four `H3`-type quotient vertices are the tetrahedron's vertices and
the six `H2`-type quotient vertices are its edges.  The normal `V4` in `A4`
acts regularly on the four vertices, recovering exactly the THM-2950/3141
torsor.  The full tetrahedral automorphism group is `S4`; quotienting its
regular `V4` gives the cubic-resolvent shadow `S4/V4=S3`.

The same quotient appears directly on the six edge slots.  The regular `V4`
has three orbits, each an opposite-edge pair,

```text
{01,23}, {02,13}, {03,12}.                               (6)
```

The residual `S3` permutes these three pairs.  Thus the familiar
`4 vertices -> 6 pairs -> 3 opposite-pair classes` ladder is the literal
tetrahedral model of the quartic / cubic-resolvent quotient.

Both the subdivided graph and `K4` have first Betti number

```text
12-(6+4)+1 = 6-4+1 = 3.                                 (7)
```

Accordingly `N` is isomorphic to the free rank-three fundamental group of
this covering graph, while its abelianization is the rank-three `H1` cycle
lattice.  The four tetrahedral triangles span that three-dimensional cycle
space.  One primitive six-edge cycle before suppression is represented by
`(SR)^3`, a first `Gamma(3)` cusp/cycle coordinate that THM-3141 shows is
invisible on the four points; it is not the whole nonabelian kernel.

## 4. Partial-cube and gate boundaries

Every tree is a partial cube: fix a root and map a vertex to the set of edge
cuts on its root path; symmetric difference has size equal to tree distance.
The finite quotient `N\T`, the once-subdivided `K4`, is also a partial cube.
Identify its four `H3` vertices with the singleton subsets of a four-set and
its six `H2` vertices with the two-subsets.  Incidence changes one coordinate,
and the graph distance equals Hamming distance.  Its four
Djokovic--Winkler `Theta` classes are the four coordinates, each containing
three edges.  Only *suppressing* the intermediate two-subsets creates the
triangles of `K4`; that nonbipartite suppressed graph is not a partial cube.
Thus the quotient folds infinitely many tree cuts to four coordinates but
does not destroy partial-cube structure until the six edge slots are forgotten.

The lost cycle is already load-bearing for THM-2056's sufficient gate.  In
the one-tail chart put

```text
D(a,b)=max(13|b|,|a-12b|),
d=(43,2),                  T^3 d=(49,2).                  (8)
```

The two primitive slopes are identical in `P1(F3)`, but

```text
91D(d)=2366 > 1853=||d||^2,
91D(T^3d)=2366 <= 2405=||T^3d||^2.                       (9)
```

Thus even the truth value of this fixed sufficient certificate is not a
function of the tetrahedral quotient point.  This says nothing about physical
unsafety of `d`; it proves that a `Gamma(3)` lift is indispensable before the
certificate can be transported.

## 5. Evidence and scope

The companion pins THM-3141 and uses finite permutation arithmetic only.  It
reconstructs `A4`, all `6+4` cosets and `12` incidence edges, the complete
tetrahedral edge list, the regular `V4` action and its three opposite-edge
orbits, the cyclic `C3` stabilizer, four triangles, both Betti-number invoices,
the isometric `Q4` embedding and four `Theta` classes, and the gate-fibre
hostile `(9)`.
Normal and optimized transcripts agree and match the stored LF-normalized
output.

This theorem identifies the exact common grammar requested by the binary /
ternary analogy.  It does not produce a common physical LRC atom, a degree-four
Keller exclusion, a faithful finite modular action, or a canonical ternary
orientation.  Those require a lift from the tetrahedron back to the
Bass--Serre cover (or an equivalent `Gamma(3)` cocycle), plus the rooted block
and physical coupling already isolated by THM-3141.
