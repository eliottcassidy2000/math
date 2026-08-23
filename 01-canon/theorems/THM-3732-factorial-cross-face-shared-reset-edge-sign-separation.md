---
id: THM-3732
title: "Factorial cross-face shared-reset-edge sign separation"
status: >
  FINITE-EXACT + INDEPENDENTLY AUDITED.  In the two named bank-I2 charts
  F12/F13, the smallest shared hostile state has one common reset-directed
  edge, but the strict-positive named-row fibres on that edge are disjoint
  and have no equal positive increment in the inherited raw 22-row gauge.
  Its endpoint has no common next reset edge and the positive fibre sizes
  differ.  The full 239-state common graph has 568 edges and six sinks.
  This excludes the inherited fixed-fibre named-row/common-edge Cech torsor
  model only; it proves no FC(3), HFC(3), or invertibility result.
source: root + lrc14-cover-defect-bridge / 2026-08-22
audit: >
  PASS after narrowing an overbroad torsor claim.  An independent
  reconstruction uses the older THM-3238 coefficient source, a fresh
  partition/upset implementation, and direct L1 reset descent.  It reproduces
  both local fibres, all raw signs, the empty equal-increment relation, all
  239 states, 568 edges, the outdegree histogram, six sinks, and the separate
  HFC boundary-current control.  Normal and optimized streams byte-match all
  frozen transcripts.
depends_on:
  - THM-3238-complete-physical-product-gamma-bank-unique-reset-stitch
  - THM-3286-three-face-availability-helly-defect-and-binary-origin-width
related:
  - THM-3466-factorial-face-stokes-and-keller-boundary-current
script: 04-computation/fc_cross_face_shared_reset_edge_thm3732.py
output: 05-knowledge/results/fc_cross_face_shared_reset_edge_thm3732.out
script_sha256: 12be0cc06ab78fb7eb1d73012436bfce6814c82412ed5721ba50b2f011e6b40d
output_sha256: 3ebae437db2a7b5f90765e0950c766af8f906cf863f3f7bc0b0f3257cad078f3
boundary_control_script: 04-computation/hfc_boundary_primitive_cocycle_control_thm3732.py
boundary_control_output: 05-knowledge/results/hfc_boundary_primitive_cocycle_control_thm3732.out
boundary_control_script_sha256: 09c9fcd704d4fc4035d85680ace0f620ee7b42cc1398e84792e728ab00b76d14
boundary_control_output_sha256: 9ea9186cbd6aba1c64814a60e2b584e5fab1231bd84c4f9d31d5c74add9be5c7
independent_audit_script: 04-computation/fc_cross_face_shared_reset_edge_independent_audit_thm3732.py
independent_audit_output: 05-knowledge/results/fc_cross_face_shared_reset_edge_independent_audit_thm3732.out
independent_audit_script_sha256: 3bee90374c6665166e1ba833f11b605f94b580b92bec79e644be8b1f46033133
independent_audit_output_sha256: 3567968730a68d1b3a30df339a03993739aff7fa5e929ed7667a578d1ae1940f
hash_basis: raw LF bytes
---

# THM-3732 -- a common reset edge has opposite named-row support

**FINITE-EXACT + INDEPENDENTLY AUDITED.**  This is a complete finite graph
and exact-integer response statement for two inherited charts.  It is not a
claim about arbitrary factorial banks.

## 1. The two named bank-I2 charts

The exact profiles are

```text
F12:
  P multiplicities ((1,4),(2,3),(3,2),(4,1),(5,1))
  reset Q12=(1,2,2,3,4,5)

F13:
  P multiplicities ((1,4),(2,3),(3,2),(4,2),(5,2),
                    (6,1),(7,1),(8,1))
  reset Q13=(1,3,3,4,5,6,7,8).                        (1)
```

They are two faces of the same named product-Gamma bank `I2`, not
relabelled copies.

For face `f`, state `n`, and reset-directed one-pole neighbour `n'`,
define

```text
W_f(n->n')={i in {1,...,22}:r_i^f(n')-r_i^f(n)>0}.     (2)
```

The universe and every response in (2) are exact integers.

## 2. The smallest shared hostile

THM-3286 identifies the smallest `F12/F13` common-row hostile as

```text
n0=(3,4,5),                    n1=(1,3,4,5).            (3)
```

The reset-directed neighbours of `n0` are

```text
N12(n0)={(1,3,4,5),(2,3,4,5)},

N13(n0)={(1,3,4,5),(3,3,4,5),
         (3,4,5,6),(3,4,5,7),(3,4,5,8)}.              (4)
```

Thus `e:n0->n1` is the unique common outgoing edge.  On that literal edge,

```text
W_12(e)={2,5,8,9,11,12,14,16,18,22},
W_13(e)={3,4,6,7,10,13,17,19,20,21}.                  (5)
```

Each set in (5) is the full availability set of its face at `n0`.  They
are disjoint and their union has 20 rows.  The omitted rows `1,15` have
negative increments in both charts.  Hence

```text
W_12(e) intersect W_13(e)=empty.                       (6)
```

In the inherited raw 22-row normalization, even the cross-row
equal-increment relation is empty:

```text
{(i,j):Delta_12,i(e)=Delta_13,j(e)>0}=empty.            (7)
```

Equation (7) is gauge-scoped; arbitrary row rescalings are not being
classified.

## 3. Chronological break at the endpoint

At `n1`,

```text
N12(n1)={(1,2,3,4,5)},

N13(n1)={(1,3,3,4,5),(1,3,4,5,6),
         (1,3,4,5,7),(1,3,4,5,8)},                    (8)
```

so

```text
N12(n1) intersect N13(n1)=empty.                       (9)
```

The positive availability cardinalities change as follows:

```text
             at n0       at n1
F12            10           9
F13            10          10.                         (10)
```

Thus a freely chosen bijection between the two ten-point positive fibres at
`n0` would neither be the inherited named-row identity nor raw-response
matching, and it cannot extend as one fixed positive fibre over `n1`.

## 4. Complete common reset graph

Exhausting the 239-state shared physical multiplicity box gives exactly 568
common reset-directed edges and

```text
outdegree:state count
 0:6, 1:41, 2:85, 3:75, 4:28, 5:4.                   (11)
```

The six sinks are

```text
(1,3,4,5),
(1,2,3,4,5), (1,3,3,4,5),
(1,2,2,3,4,5), (1,2,3,3,4,5),
(1,2,2,3,3,4,5).                                      (12)
```

The complete coordinate rule is:

```text
root 1: move c_1 one step toward 1 whenever c_1!=1;
root 2: 3->2;
root 3: 0->1;
root 4: 0->1;
root 5: 0->1.                                          (13)
```

There are no other common moves.  Formula (13) explains the mechanism:
the two resets demand opposite root-2 directions at multiplicity one, while
root 3 loses its common direction immediately after insertion.  The edge
`e` lands at the smallest sink in (12).

The graph was reconstructed independently from its multiplicity boxes and
L1 descent, while the responses were reconstructed from the older THM-3238
coefficient source.  The stable graph and edge-record digests are

```text
8cffe16bf5b2301cf077948da8588edf61edbe837a9015a916b3bea3a963c0e3
b94614818d071bfc341057bafa220e71c887680639c896d06c393389cebc27ab.
```

## 5. Exact categorical consequence

Consider the inherited candidate whose source consists of

```text
(face,state,named row,witness neighbour),
```

whose cross-face transition is identity on the raw state and named row with
a common neighbour, and whose preserved predicate is reset-directed motion
plus strict positive response.

This candidate does not form a fixed-fibre Cech torsor cocycle across
`F12/F13`:

- its identity-row certificate relation is empty on the unique common edge
  by (6);
- raw-gauge equal positive increments supply no alternative relation by (7);
- positive fibre cardinalities disagree at the endpoint by (10);
- chronological continuation by a common physical edge stops by (9).

This conclusion does **not** exclude an arbitrary bijection between separate
finite sets, an imposed regular group action, a variable groupoid/torsor, or
an enriched relation carrying history or magnitude.  None of those
structures is supplied by the inherited bank.

## 6. HFC boundary primitive is a different positive control

Under the HFC-null-candidate premise of THM-3466, local primitives of
`g^2 d(conjugate g)` on labelled oriented boundary edges differ by additive
constants.  Their cyclic sum is the total boundary current, hence zero.  An
orientation and basepoint therefore trivialize this genuine additive
cocycle, although THM-3466's slit-annulus hostile shows that the primitive
can still fail to separate sheets.

The exact rational-triangle control has vertices

```text
(-1,-1), (2,-1), (-1,2)
```

and edge increments

```text
-3i, -3+3i, 3,
```

whose sum is zero; its global-basepoint offsets are `0,-3i,-3`.  This
checks one boundary-current closure moment.  It is not a full HFC example.

No maintained map sends a product-Gamma state/row/face to a labelled
boundary point, primitive constant, or sheet.  Therefore no reset/HFC
cocycle identification is typed.

## 7. Scope and reproduction

The theorem proves no `FC(3)`, `HFC(3)`, `SFC(3)`, GMC, positivity, or
polynomial-invertibility conclusion.

```bash
python3 -B 04-computation/fc_cross_face_shared_reset_edge_thm3732.py
python3 -B -O 04-computation/fc_cross_face_shared_reset_edge_thm3732.py
python3 -B 04-computation/hfc_boundary_primitive_cocycle_control_thm3732.py
python3 -B -O 04-computation/hfc_boundary_primitive_cocycle_control_thm3732.py
python3 -B 04-computation/fc_cross_face_shared_reset_edge_independent_audit_thm3732.py
python3 -B -O 04-computation/fc_cross_face_shared_reset_edge_independent_audit_thm3732.py
```

Each normal/optimized pair must agree byte for byte with its frozen
transcript.
