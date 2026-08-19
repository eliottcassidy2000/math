---
id: THM-3553
title: "Fixed three-dimensional Keller map polynomial graph-section no-go"
status: >
  PROVED + VERIFIED-EXACT + INDEPENDENTLY DERIVATION-AUDITED.  For the fixed
  three-dimensional Keller map of THM-1300, no polynomial source graph
  z=h(x,y) makes the restricted target pair (F1,F2) have constant planar
  Jacobian.  Constant graphs retain a universal -89 coefficient.  For every
  positive degree d, the top tangential minor is
  -3 S partial_x(S), S=x^3 y^2 h_d, and is nonzero.  Hence this entire fixed
  graph chart cannot descend the known collision to a planar Keller pair or
  be carried into a polynomial target graph in the displayed coordinates.
  Other source/target coordinates and nongraph coordinate hypersurfaces
  remain open.
source: root-2026-08-18-planar-jacobian-counterexample-hostiles
depends_on:
  - THM-1300-jacobian-counterexample-dixmier-A3-explicit
related:
  - THM-3543-quotient-stabilization-collision-versus-ramification-boundary
  - THM-3544-planar-keller-target-pencil-total-degree-six-floor
  - THM-3546-invariant-graph-keller-descent-criterion
  - THM-3550-prime-degree-exclusion-and-pencil-height-eight-floor
  - THM-3554-punctured-kummer-collision-surface-normal-form
script: 04-computation/jc_fixed_threedimensional_graph_section_no_go_thm3553.py
output: 05-knowledge/results/jc_fixed_threedimensional_graph_section_no_go_thm3553.out
script_sha256: e958c30bdfe681a62149b625ba310ba77354918384c1667d1f45cc81e54296f3
output_sha256: f37a82c1cb18a5532d1b39f4061a8919bbc26c25753c570ef9fb4a15d72785af
semantic_sha256: 8a21ebd99f5ba06d83feb6de7561f36bde3a4053b262f1a7abaaab5206c12976
hash_basis: LF-normalized bytes
---

# THM-3553 -- fixed three-dimensional Keller map polynomial graph-section no-go

**PROVED + VERIFIED-EXACT + INDEPENDENTLY DERIVATION-AUDITED.**

Consider the fixed map from THM-1300.  With `u=1+xy`, write

```text
F1=u^3 z+y^2u(4+3xy),
F2=y+3xu^2z+3xy^2(4+3xy),
F3=2x-3x^2y-x^3z.                                      (1)
```

Its three-dimensional Jacobian determinant is `-2`, and the three points

```text
p0=(0,0,-1/4),   p+=(1,-3/2,13/2),   p-=(-1,3/2,13/2)  (2)
```

all map to `(-1/4,0,0)`.

For every polynomial `h in C[x,y]`, restrict `(1)` to the graph `z=h(x,y)`:

```text
P_h=F1(x,y,h),        Q_h=F2(x,y,h).                     (3)
```

Then

```text
Jac(P_h,Q_h) is never constant.                          (4)
```

Consequently no graph in this fixed source chart, including every graph
through any collision pair in `(2)`, yields a planar Keller pair using the
fixed target coordinates `(F1,F2)`.

More strongly, there are no polynomials `h,g in C[x,y]` satisfying

```text
F3(x,y,h(x,y))=g(P_h(x,y),Q_h(x,y)).                    (4a)
```

Indeed, `(4a)` is scheme-theoretic containment of the source graph in a
polynomial target graph.  THM-3546 would then force `Jac(P_h,Q_h)` to be a
nonzero constant, contradicting `(4)`.

## 1. Exact tangential-minor equation

Set `t=xy`.  The chain rule gives

```text
T_h:=Jac(P_h,Q_h)=M0(x,y,h)+Mx(x,y,h)h_x+My(x,y,h)h_y,  (5)
```

where

```text
M0 = -9x u^4 h^2
     -3y u^2(15t^2+25t+7)h
     -y^3(54t^3+189t^2+222t+89),

Mx = -u^2(3h x^2u^2+9t^3+12t^2-t-1),

My = -3u^2(h u^2+y^2(3t+4)).                            (6)
```

Equivalently, `M0,Mx,My` are the three relevant `2 x 2` minors of the
ambient Jacobian, evaluated on the graph.  The exact companion verifies
`(5)--(6)` both by direct differentiation of `(3)` and by the cofactor
formula.

## 2. Every positive graph degree fails at its top layer

Let `d>=1`, and let `h_d` be the nonzero degree-`d` homogeneous part of `h`.
The unique top pieces of `(3)` are

```text
(P_h)_(d+6)=x^3y^3h_d,
(Q_h)_(d+5)=3x^3y^2h_d.                                 (7)
```

Put

```text
S=x^3y^2h_d.                                             (8)
```

Then the homogeneous part of `T_h` of degree `2d+9` is

```text
(T_h)_(2d+9)=Jac(yS,3S)=-3S partial_x(S).               (9)
```

Every monomial of `S` contains `x` to exponent at least three.  In
characteristic zero, `partial_x(S)` is therefore nonzero whenever `h_d` is
nonzero.  Since the polynomial ring is a domain, the right side of `(9)` is
nonzero.  Thus `T_h` cannot be constant for any positive-degree `h`.

If `h` is constant, `(5)--(6)` instead give the universal coefficient

```text
[y^3]T_h=-89.                                            (10)
```

So constant graphs fail as well.  This proves `(4)`.

The top layer also gives an independent constant-graph check.  For `h=c`,

```text
(P_h)_6=x^2y^3(cx+3y),       (Q_h)_5=3x^2y^2(cx+3y),
```

whose Jacobian is `-9x^3y^4(cx+3y)(cx+2y)`, nonzero even at `c=0`.

## 3. Collision interpolation is not the missing ingredient

For an affine graph, the value constraints through any of

```text
(p0,p+),       (p0,p-),       (p+,p-)                   (11)
```

leave one affine parameter.  Direct substitution preserves `(10)` for all
three families.  For degrees two and three, the interpolation spaces have
dimensions four and eight; equation `(9)` eliminates every exact-degree
member.  A specialization with vanishing leading form drops into an already
eliminated lower degree.

These bounded searches are not rejected merely by the inherited degree
floor.  For exact graph degree `d>=1`, equations `(7)` have different
degrees, so every nonzero target-pencil member has degree at least `d+5`.
The quadratic and cubic cells therefore meet THM-3544 with floors seven and
eight before the tangential obstruction kills them.

The only constant graph through a collision pair is `h=13/2`, through
`(p+,p-)`.  Its restricted degrees are `(6,5)`, so it violates THM-3544's
degree-six floor for a nonautomorphic planar Keller pair; `(10)` rejects it
independently.

## 4. Scope and the surviving descent problem

THM-3546 gives a sufficient descent mechanism when a three-dimensional
Keller map carries a polynomial coordinate graph scheme-theoretically to a
target coordinate graph.  The present theorem closes the most literal
realization for `(1)`: source graphs over `(x,y)` and tangential target
coordinates `(F1,F2)`.

It does **not** close:

1. a graph after changing the source polynomial coordinates;
2. a three-dimensional target coordinate change that mixes in `F3`;
3. a polynomially straightened hypersurface that is not a graph in this
   chart;
4. another three-dimensional Keller map with a collision.

Those distinctions are essential.  The obstruction `(9)` is chart-specific,
not a theorem that graph descent can never work.

## 5. Exact verification

Reproduce with

```bash
python3 04-computation/jc_fixed_threedimensional_graph_section_no_go_thm3553.py
python3 -O 04-computation/jc_fixed_threedimensional_graph_section_no_go_thm3553.py
```

The ordinary and optimized transcripts agree byte-for-byte.  The companion
checks the ambient determinant and three-point collision, the compact PDE
against direct differentiation, all three affine collision families, all
six quadratic/cubic pair cells, the target-pencil floors, and the general
leading-form identity `(9)`.  An independent derivation obtained `(9)`
directly from the top restricted coordinates and confirmed that no bounded
coefficient search is being extrapolated into the all-degree claim.

**QED.**
