---
id: THM-2451
title: "Separated-ruling hostile to the two-jet plane-map cone prediction"
status: >
  PROVED + VERIFIED-EXACT + INDEPENDENTLY HOSTILE-AUDITED. For
  arbitrary constant vectors a,b,n in
  C^3, the separated-ruling family A=a+xn, B=b+yn, C=0 makes every
  z-graded Jacobian bracket D_5,...,D_0 vanish. Nevertheless, for
  w=A cross B one has the exact identity
  [w_x,w_y,w]=[a,b,n]^2. Thus whenever a,b,n are independent the
  plane map w is dominant, refuting THM-2446's open prediction (P2)
  even under the stronger hypothesis D_5=D_4=D_3=0. The hostile lies
  in the line-direction stratum of A; a repaired prediction restricted
  to a nondegenerate conic direction map remains open. This gives no
  Keller map, wildness certificate, or planar Jacobian consequence.
source: codex-2026-07-26-twojet-plane-map-cone-hostile
depends_on:
  - THM-2446-twojet-zgraded-jacobian-decomposition-and-cone-system
related:
  - THM-1310-conic-pair-fibers-and-design-equations
script: 04-computation/jacobian_twojet_plane_map_cone_hostile_thm2451.py
output: 05-knowledge/results/jacobian_twojet_plane_map_cone_hostile_thm2451.out
script_sha256: f9a2d9af608887d45d1fb28b86715ed816d8e21a0379d5454884bf3daee88202
output_sha256: 38a93824fec5b5123bf361ea4a18b28c2374c58035bd3da62c6542df3934072a
hash_basis: working-tree bytes (LF)
---

# THM-2451 -- the first three two-jet brackets do not cone the plane map

**PROVED + VERIFIED-EXACT + INDEPENDENTLY HOSTILE-AUDITED.**

THM-2446 opens the `z`-quadratic architecture

```text
F=A(x,y)z^2+B(x,y)z+C(x,y)
```

and records the prediction

```text
D_5=D_4=0 plus the A-part of D_3
    => [w_x,w_y,w]=0,                 w=A cross B.       (P2)
```

The prediction is false.  The failure is not a coefficient accident but
an entire separated-ruling family.

## 1. Counterfamily

Fix constant vectors `a,b,n in C^3` and put

```text
A(x,y)=a+x n,

B(x,y)=b+y n,

C(x,y)=0.                                             (1)
```

Then

```text
A_x=n,       A_y=0,       B_x=0,       B_y=n.         (2)
```

Every bracket in THM-2446's `D_5,D_4,D_3` has a zero column or two
copies of `n`.  Explicitly,

```text
D_5=2[n,0,A]=0,

D_4=[n,0,B]+2[n,n,A]+2[0,0,A]=0,

D_3=[n,n,B]+[0,0,B]
    +2([n,0,A]+[0,0,A]+[0,n,A])=0.                  (3)
```

Since `C=0`, the same inspection gives

```text
D_2=D_1=D_0=0.                                      (4)
```

Thus (1) satisfies not only every natural reading of the proposed
“`A`-part of `D_3`”, but the complete stronger condition
`D_5=D_4=D_3=0`.

## 2. The plane map is dominant

The cross-product plane map is

```text
w=A cross B
  =a cross b+y(a cross n)+x(n cross b),             (5)

w_x=n cross b,                  w_y=a cross n.       (6)
```

The two variable terms in (5) repeat a column in
`[w_x,w_y,w]`, so

```text
[w_x,w_y,w]
 =[n cross b,a cross n,a cross b].                  (7)
```

Using

```text
(n cross b) cross (a cross n)=[a,b,n] n,
```

equation (7) becomes

```text
[w_x,w_y,w]=[a,b,n]^2.                              (8)
```

If `a,b,n` are independent, the right side is nonzero.  Therefore
`w:C^2->C^3` is dominant onto the affine plane chart of its projective
image rather than being constrained to a curve.

The smallest integral witness is

```text
a=(1,0,0),       b=(0,1,0),       n=(0,0,1),

A=(1,0,x),       B=(0,1,y),

w=(-x,-y,1),     [w_x,w_y,w]=1.                    (9)
```

## 3. Mechanism and repaired frontier

The failed implication forgets the **transverse pairing of the two
rulings**.  The top Jacobian brackets see only repeated first-jet
directions: `A` moves in the `x` direction by `n`, while `B` moves in
the `y` direction by the same `n`.  Their cross product converts those
two individually rank-one motions into the independent tangent vectors
`n cross b` and `a cross n`.  The lost mixed coordinate is exactly the
nonzero volume `[a,b,n]`.

This hostile belongs to THM-2446's line-direction stratum:
`[A(x,y)]` traces a projective line and is independent of `y`.  Hence
the strongest immediate repair is:

```text
Does a nondegenerate conic direction map [A],
plus the appropriate D_3 slice, force [w_x,w_y,w]=0?             (10)
```

Equation (10) remains **OPEN**.  Any replacement must explicitly exclude
the separated-ruling kernel (1), rather than merely add more copies of
the already vanishing brackets.

## 4. Scope

The counterfamily has `det J(F)=0`; it is not a Keller map.  It refutes
only the proposed intermediate cone implication.  It neither constructs
nor excludes a wild `z`-quadratic Keller map, and it has no consequence
for planar `JC(2)` or `DC(2)`.

## 5. Exact companion

Run

```text
python3 04-computation/jacobian_twojet_plane_map_cone_hostile_thm2451.py
python3 -O 04-computation/jacobian_twojet_plane_map_cone_hostile_thm2451.py
```

The dependency-free symbolic companion expands all six brackets for
generic `a,b,n`, proves (8), checks the integral witness (9), and verifies
the full Jacobian determinant is zero.  Normal and optimized transcripts
must match

```text
05-knowledge/results/jacobian_twojet_plane_map_cone_hostile_thm2451.out
```

byte-for-byte.

An independent hostile audit reconstructed every bracket directly from
THM-2446 and checked that the example refutes all ordinary readings of
“the `A`-part of `D_3`”: the `C`-free part, the terms containing `A`,
the third-slot term, the complete `D_3` after `C=0`, and the leading
`E_3` expression with `gamma=0`.  It also confirmed that THM-2446 states
no hidden conic-nondegeneracy premise; that condition appears only as a
later hunt-design suggestion.  QED.
