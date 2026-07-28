---
id: THM-2709
title: "Complete linear target-plane classification on polynomial graph Keller slices"
status: >
  PROVED + VERIFIED-EXACT + INDEPENDENTLY HOSTILE-AUDITED.  Every affine-linear
  target two-plane in the
  THM-2696 quotient is classified on every polynomial graph.  Planes
  containing d are exactly THM-2705.  Of the remaining planes, those
  containing A have one explicit shifted-cubic triangular graph, while a
  plane containing neither A nor d has no nonzero-constant-Jacobian polynomial
  graph.  Hence every surviving map is a polynomial automorphism.  Nonlinear
  target projections, arbitrary polynomial source surfaces, JC(2), and DC(2)
  remain open.
source: root/reflection-quotient-pluecker-graph-2026-07-28; root-long-frontiers independent referee 2026-07-28
depends_on:
  - THM-2696-reflection-completed-s4-relative-different-and-coordinate-invariant-jacobian-gate
  - THM-2702-polynomial-graph-coordinate-projection-keller-classification
  - THM-2705-linear-target-planes-containing-d-polynomial-graph-keller-classification
related:
  - THM-2699-affine-plane-linear-projection-keller-slice-classification
  - THM-2655-quartic-keller-resolvent-v4-quasietale-torsor-and-kummer-class-group-gate
script: 04-computation/jacobian_s4_polynomial_graph_all_linear_planes_thm2709.py
output: 05-knowledge/results/jacobian_s4_polynomial_graph_all_linear_planes_thm2709.out
script_sha256: d851c6671576baac86e44eeae57d61e1aaec28eb02a77409931bba50b6718cb0
output_sha256: 961fab46f659be3db5467257a149b3bbdff8fc1c9d30c2d8d829a060ba65ad53
hash_basis: LF-normalized bytes
---

# THM-2709 -- every linear target plane is classified on every polynomial graph

**PROVED + VERIFIED-EXACT + INDEPENDENTLY HOSTILE-AUDITED.**

THM-2702 classifies the three coordinate target pairs, and THM-2705 classifies
every target plane containing `d`.  There is a short Pluecker completion: the
remaining target planes either contain `A`, where one translated version of
the THM-2702 cubic survives, or contain neither `A` nor `d`, where the Keller
constant contradicts the last coefficient of the graph equation.

This is the complete affine-Grassmannian answer for polynomial graphs in the
fixed THM-2696 quotient coordinates.  It is not the planar Jacobian conjecture:
the source is still a graph in a special threefold quotient and the targets
are still linear functions of `(A,B,d)`.

## 1. Pluecker form of an arbitrary target plane

Work over `C` on

```text
z=f(x,y),                 f in C[x,y],                    (1)
A=x^2-2y,                 B=y^2-2xf,       d=f.           (2)
```

Let `(U,V)` be any linearly independent pair of linear forms in `(A,B,d)`.
Its `2 by 3` coefficient matrix has Pluecker coordinates

```text
P=u_A v_B-u_B v_A,
Q=u_B v_d-u_d v_B,
R=u_d v_A-u_A v_d.                                      (3)
```

Then bilinearity of the determinant gives

```text
Jac(U,V)=P Jac(A,B)+Q Jac(B,d)+R Jac(d,A).               (4)
```

The nonzero projective triple `[P:Q:R]` depends only on the oriented target
plane, up to the determinant of a target basis change.

The geometric flags are exact:

```text
P=0       iff the target plane contains d,
Q=0       iff the target plane contains A.               (5)
```

Indeed, `P` (respectively `Q`) is the determinant of the two row restrictions
to the `(A,B)` (respectively `(B,d)`) coordinate plane.  Its vanishing produces
a nonzero row combination supported on `d` (respectively `A`), and the
converse is immediate.

Fix `kappa!=0`.  The case `P=0` is exactly THM-2705.  It remains to assume
`P!=0`.

## 2. The normalized graph equation

The three coordinate determinants from THM-2702 are

```text
Jac(A,B)=-4(xf_x+x^2f_y+f-xy),
Jac(B,d)=-2(yf_x+ff_y),
Jac(d,A)=-2(f_x+xf_y).                                  (6)
```

Set

```text
a=Q/(2P),             b=R/(2P),
c=-kappa/(4P).                                         (7)
```

Since `P,kappa!=0`, also `c!=0`.  Equations `(4)` and `(6)` show that
`Jac(U,V)=kappa` is equivalent to

```text
(x+ay+b)f_x+(x^2+af+bx)f_y+f-xy=c.                     (8)
```

This equation remembers the target flag: `a=0` is precisely the branch
`Q=0`, hence the planes containing `A`; `a!=0` means the plane contains
neither `A` nor `d`.

## 3. Planes containing `A`: one shifted cubic

Suppose `a=0`.  Put

```text
t=y-x^2/2.                                              (9)
```

The vector field in `(8)` becomes `(x+b)partial_x` in `C[x,t]`, and

```text
partial_x((x+b)f)=xt+x^3/2+c.                          (10)
```

Therefore

```text
(x+b)f=x^2t/2+x^4/8+cx+H(t).                           (11)
```

Polynomiality is equivalent to divisibility by `x+b`.  Evaluating the
numerator at `x=-b` forces

```text
H(t)=-b^2t/2-b^4/8+bc.                                (12)
```

Thus the unique solution is

```text
f=(x-b)t/2+(x-b)(x^2+b^2)/8+c,
t=y-x^2/2.                                             (13)
```

This includes THM-2702's `(A,B)` cubic at `b=0`.

A target basis with Pluecker triple `(P,0,R)` is

```text
U=P A,                 V=B-2b d.                       (14)
```

In the coordinate `(x,t)`, solution `(13)` gives

```text
U=-2Pt,
V=(t+b^2/2)^2-2c(x+b).                                (15)
```

Since `P,c!=0`, this is triangular.  Its inverse is

```text
t=-U/(2P),
x=((t+b^2/2)^2-V)/(2c)-b,
y=t+x^2/2.                                             (16)
```

In particular, the sole graph in this branch is a polynomial automorphism.

## 4. Planes containing neither `A` nor `d` are empty

Assume `a!=0`, and let `s=deg_y(f)` for nonzero `f`.

- If `s>=3`, the nonlinear term `a f f_y` has unique top `y`-degree
  `2s-1`, with nonzero coefficient `a s f_s(x)^2`.  Every other term in
  `(8)` has degree at most `s+1`, a contradiction.
- If `s=2`, the coefficient of `y^3` is

  ```text
  a(f_2'+2f_2^2).                                      (17)
  ```

  The polynomial Riccati equation `f_2'+2f_2^2=0` has only the zero
  polynomial solution: positive degrees cannot match, and a nonzero constant
  leaves a nonzero square.

Hence `deg_y(f)<=1`.  Write

```text
f=u(x)y+v(x).                                          (18)
```

The coefficients of `y^2,y,1` in `(8)` are

```text
a u'=0,
a u^2+a v'+(x+b)u'-x+u=0,
a u v+(x+b)v'+(x^2+bx)u+v-c=0.                        (19)
```

Thus `u` is constant, and the second equation gives

```text
v=x^2/(2a)-(u+a u^2)x/a+v_0.                          (20)
```

Substitution in the last equation of `(19)` has `x^2` coefficient

```text
3(au+1)/(2a),                                         (21)
```

so `u=-1/a`.  After this substitution the entire remaining constant term is
exactly `-c`.  It cannot vanish because `c!=0`.  This contradiction proves
that the `P Q!=0` branch is empty.

The boundary is sharp at zero Jacobian.  If `c=0`, then for arbitrary
`v_0 in C`

```text
f=-(y-x^2/2)/a+v_0                                    (22)
```

solves `(8)`.  Thus it is the nonzero Keller constant, not mere polynomial
invariance of the graph, that empties this Pluecker cell.
In the corresponding normalized row basis the first target coordinate is

```text
U=A-2a d=-2a v_0,
```

so the equality family has literal target-rank collapse; it is not a hidden
zero-Jacobian analogue of the surviving automorphism family.

## 5. Complete classification

Combining the target flags `(5)` with THM-2705 and Sections 3--4 gives:

| Pluecker branch | intrinsic target flag | polynomial graph solutions |
|---|---|---|
| `P=0` | contains `d` | exactly the THM-2705 triangular/affine families |
| `P!=0,Q=0` | contains `A` but not `d` | the unique shifted cubic `(13)` |
| `P Q!=0` | contains neither `A` nor `d` | none |

Every nonzero-constant-Jacobian solution in the entire target Grassmannian is
therefore an explicit polynomial automorphism.

## 6. Boundary and next obstruction

The classification is complete only inside the following box:

```text
source:  one global polynomial graph z=f(x,y);
target:  two affine-linear combinations of A,B,d;
map:     the restriction of the fixed THM-2696 quotient coordinates.
```

It does not cover a source surface that is not a graph in the chosen chart,
a nonlinear target projection, an arbitrary Keller map, `JC(2)`, or `DC(2)`.
Affine target constants do not matter, but polynomial target terms do.

The proof also identifies the next exact probe.  One must leave either the
graph chart or the linear target Grassmannian.  Varying another linear target
plane cannot produce new behavior: its Pluecker flag already lands in the
three-row table above.

## 7. Exact companion

Run

```text
python 04-computation/jacobian_s4_polynomial_graph_all_linear_planes_thm2709.py
python -O 04-computation/jacobian_s4_polynomial_graph_all_linear_planes_thm2709.py
```

Both modes must byte-match

```text
05-knowledge/results/jacobian_s4_polynomial_graph_all_linear_planes_thm2709.out.
```

The companion uses explicit `require` checks.  It verifies the general row
wedge `(4)`, normalized PDE `(8)`, the shifted cubic and polynomial inverse,
the all-degree top-degree/Riccati inequalities, every coefficient in `(19)`,
the final `-c` obstruction, and the sharp zero-Jacobian family `(22)`.  The
finite identities support but do not replace the all-degree proof.

An independent hostile audit first verified promoted THM-2705, then rederived
the Pluecker expansion and both intrinsic target flags.  It separately checked
the signs and scalings in `(8)`, polynomial divisibility and the inverse in the
shifted-cubic branch, the all-degree and Riccati gates, every coefficient in
`(19)`, the terminal `-c` obstruction, the sharp `c=0` family, exhaustiveness,
and scope.  Normal and optimized execution byte-matched the canonical recorded output
SHA-256
`961fab46f659be3db5467257a149b3bbdff8fc1c9d30c2d8d829a060ba65ad53`;
the script SHA-256 is
`d851c6671576baac86e44eeae57d61e1aaec28eb02a77409931bba50b6718cb0`.
The stored transcript has ten lines; the script has zero Python `assert` nodes
and twelve explicit `require` calls.

A second independent referee used the projective normal `(a,b,1)` and row basis
`(A-a d,B-b d)`, providing a separately normalized derivation of the same
three Pluecker cells.  It checked the wedge signs and scalings, the `r>=3`
top-degree gate, the quadratic Riccati boundary, all `r<=1` coefficients, the
zero-Jacobian rank collapse, the Euler inverse, the shifted hostile
specialization, and the theorem's scope.  Its exact companions are

```text
04-computation/jacobian_s4_polynomial_graph_all_linear_target_planes_referee_20260728.py
05-knowledge/results/jacobian_s4_polynomial_graph_all_linear_target_planes_referee_20260728.out
```

with SHA-256 values

```text
840fb2b3ba5f48d897d0867db484cd565d04154edaf8507a03850e51923f692a
76b06fb9357268255969794010da74785822a6182f6ac64661006ba60f224faf
```

Normal and optimized executions of both the canonical companion and the
independent referee byte-match their frozen transcripts.  THM-2705 is now a
proved dependency.

QED.
