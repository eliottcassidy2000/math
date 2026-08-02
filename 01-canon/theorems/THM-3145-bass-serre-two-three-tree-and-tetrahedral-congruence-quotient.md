---
id: THM-3145
title: "Bass-Serre two-three tree and tetrahedral congruence quotient"
status: "CANDIDATE + VERIFIED-EXACT; AWAITING INDEPENDENT HOSTILE AUDIT"
source: root/codex-thm3088-push-2026-08-02
depends_on:
  - THM-3141-quartic-v4-modular-congruence-shadow-and-gamma3-sidecar-boundary
  - THM-2056-kelvin-polar-farey-defect-certificate
related:
  - THM-2950-three-conjugate-pair-v-four-torsor-and-quartic-resolvent-frame
  - THM-2951-fifth-compound-reconstruction-and-v-four-phase-scalarization-boundary
script: 04-computation/modular_bass_serre_two_three_tetrahedral_quotient_thm3145.py
output: 05-knowledge/results/modular_bass_serre_two_three_tetrahedral_quotient_thm3145.out
hash_basis: LF-normalized bytes
script_sha256: 3bbd4eb91dfe58948275cca7083a052f095010757f5b70333b342908d5cdfdd3
output_sha256: 78893a91ae111bca1257cb408a6248e29be488952014c9e154029ad1d9a3a1d7
---

# THM-3145 -- Bass--Serre two-three tree and tetrahedral congruence quotient

**CANDIDATE + VERIFIED-EXACT; AWAITING INDEPENDENT HOSTILE AUDIT.**

The binary and ternary tower grammars have a literal common object.  It is
not a homogeneous tree and not the finite four-point torsor: it is the
Bass--Serre `(2,3)`-biregular tree of

```text
G = C2 * C3 = <s,r | s^2=r^3=1>.                         (1)
```

The rooted binary tree and the ternary triple are two different views of
this one tree.  THM-3141's four-point `V4` torsor is its level-three finite
quotient after the missing cycle coordinate has been forgotten.

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
This is the binary fraction-tree grammar: binary branching is the *forward*
view of a rooted trivalent tree, not a second unrelated edge type.

Conversely, a chosen generator `r` cyclically orders the three edges incident
to every `H3` vertex.  This local cyclic star is the ternary triple grammar.
Changing `r` to `r^-1` reverses the cyclic orientation, so the ternary order
is a sidecar and not part of the unpointed abstract tree.

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

Both the subdivided graph and `K4` have first Betti number

```text
12-(6+4)+1 = 6-4+1 = 3.                                 (6)
```

Accordingly `N` is the free rank-three cycle group of this covering.  The
four tetrahedral triangles span a three-dimensional cycle space.  The first
six-edge cycle before suppression is the relation `(SR)^3`, exactly the
`Gamma(3)` coordinate that THM-3141 shows is invisible on the four points.

## 4. Partial-cube and gate boundaries

Every tree is a partial cube: fix a root and map a vertex to the set of edge
cuts on its root path; symmetric difference has size equal to tree distance.
Thus the universal Bass--Serre tree has honest binary cut coordinates.
The suppressed quotient `K4` has triangles and is not bipartite, hence cannot
be a partial cube.  The congruence quotient has glued the cut coordinates;
its ternary triangle is a quotient cycle, not a new hypercube direction.

The lost cycle is already load-bearing for THM-2056's sufficient gate.  In
the one-tail chart put

```text
D(a,b)=max(13|b|,|a-12b|),
d=(43,2),                  T^3 d=(49,2).                  (7)
```

The two primitive slopes are identical in `P1(F3)`, but

```text
91D(d)=2366 > 1853=||d||^2,
91D(T^3d)=2366 <= 2405=||T^3d||^2.                       (8)
```

Thus even the truth value of this fixed sufficient certificate is not a
function of the tetrahedral quotient point.  This says nothing about physical
unsafety of `d`; it proves that a `Gamma(3)` lift is indispensable before the
certificate can be transported.

## 5. Evidence and scope

The companion pins THM-3141 and uses finite permutation arithmetic only.  It
reconstructs `A4`, all `6+4` cosets and `12` incidence edges, the complete
tetrahedral edge list, the regular `V4` action, the cyclic `C3` stabilizer,
four triangles, both Betti-number invoices, and the gate-fibre hostile `(8)`.
Normal and optimized transcripts agree and match the stored LF-normalized
output.

This theorem identifies the exact common grammar requested by the binary /
ternary analogy.  It does not produce a common physical LRC atom, a degree-four
Keller exclusion, a faithful finite modular action, or a canonical ternary
orientation.  Those require a lift from the tetrahedron back to the
Bass--Serre cover (or an equivalent `Gamma(3)` cocycle), plus the rooted block
and physical coupling already isolated by THM-3141.
