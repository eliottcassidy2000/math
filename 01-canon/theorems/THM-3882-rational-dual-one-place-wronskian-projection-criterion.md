---
id: THM-3882
title: "Rational-dual one-place Wronskian projection criterion"
status: >
  PROVED + VERIFIED-EXACT + INDEPENDENTLY HOSTILE-AUDITED.  For
  an everywhere-immersed rational plane curve, a line section of its dual
  has divisor twice the base fibre of the corresponding point projection
  plus that projection's ramification divisor.  One-point support occurs
  exactly when the projection has degree one and its base fibre is a single
  point of multiplicity d-1.  Immersion therefore forbids one-place lines
  for every primal degree d>=3.  Consequently no irreducible sextic with
  exactly six A2 cusps and four A1 nodes has an affine chart with A1
  normalization.  This closes an equisingularity architecture, not JC(2).
source: root / post-THM-3879 dual-Wronskian reframe, 2026-08-23
audit: >
  INDEPENDENT HOSTILE AUDIT PASS (jc_quartic_c3_construct, 2026-08-23).  The
  audit rederived the homogeneous tangent-map basepoint criterion, the exact
  `2D+Ram` divisor identity, the Riemann--Hurwitz/local-index equivalence,
  and the immersed one-address boundary.  It independently checked the
  Pluecker/reflexivity inference for the full `6A2+4A1` packet and the
  THM-3879 `(3,3)` node specialization.  The proof identifies the line pullback
  with the Wronskian of a point projection, factors its base divisor with
  multiplicity two, invokes Riemann--Hurwitz only after recording the exact
  projection degree, and checks the local immersion boundary.  The exact
  companion verifies the Wronskian factorization, every degree ledger through
  d=40, the Pluecker packet, and the THM-3879 node/projection specialization.
  Normal and optimized runs byte-match the frozen transcript.
depends_on:
  - THM-3879-rational-torus-sextic-c3-packet-one-place-tradeoff
related:
  - THM-3841-cuspidal-jelonek-three-puncture-obstruction
  - THM-3851-tricuspidal-quartic-rank-two-two-place-tradeoff
script: 04-computation/jc2_rational_dual_one_place_wronskian_projection_criterion_thm3882.py
output: 05-knowledge/results/jc2_rational_dual_one_place_wronskian_projection_criterion_thm3882.out
script_sha256: 94374be74db4c08049d18a06316f53a342b20340bdd3dcdee12d1810bb060854
output_sha256: 8c9fa0ede45e359ea08a928506ba227662ce32d8825c9b9f1cc3d392f0d6438b
semantic_sha256: 6452278a5afd71f244134a6871e247195c6af4f80f641c1b05c7f6b6dd6b8fc2
hash_basis: raw LF bytes
---

# THM-3882 -- a one-place dual line is a degree-one projection

**PROVED + VERIFIED-EXACT + INDEPENDENTLY HOSTILE-AUDITED.**
Work over an algebraically closed field `k` of characteristic zero.  Let

```text
nu:P1 -> C subset P2                                           (1)
```

be the normalization of a nonlinear rational plane curve of degree `d>=2`.
Assume that `nu` is an immersion at every normalization point.  Let
`nu_dual:P1 -> C_dual` be the tangent-line map.

For a point `P in P2`, write `H_P` for the line in the dual plane whose
points are the lines through `P`.  Choose two independent linear forms
`L_0,L_1` vanishing at `P`.  Their pullbacks are degree-`d` binary forms.
Factor their complete common base divisor:

```text
L_0 o nu = g x,              L_1 o nu = g y,
D_P=div(g),                  m=deg D_P,
gcd(x,y)=1,                  e=d-m.                            (2)
```

Then `[x:y]` is the resolved projection from `P`, a morphism

```text
phi_P:P1 -> P1                     of degree e,                (3)
```

and the exact divisor identity is

```text
(nu_dual)^* H_P = 2D_P + Ram(phi_P).                           (4)
```

In particular,

```text
supp((nu_dual)^*H_P)={p}
iff
D_P=(d-1)p and deg(phi_P)=1.                                  (5)
```

Because `nu` is immersive, a base fibre supported at one point has
coefficient one there.  Thus `(5)` can occur only for `d=2`.  For every
`d>=3`, every line on the dual curve meets its normalization in at least two
support points.

## 1. A dual line is a projection Wronskian

Move `P` to `[0:0:1]` and write the parametrization as `[X:Y:Z]`.  The line
`H_P` asks whether the tangent to `C` passes through `P`.  In a local
parameter `t`, its pullback is therefore the Wronskian

```text
W_P=X Y'-Y X'.                                                 (6)
```

With `X=gx` and `Y=gy`, the derivative of `g` cancels exactly:

```text
W(gx,gy)=g^2 W(x,y).                                          (7)
```

The zero divisor of the coprime Wronskian `W(x,y)` is precisely the
ramification divisor of `[x:y]`.  Formula `(7)` proves `(4)` locally, hence
globally.  It also explains the otherwise easy-to-miss factor two on the
fibre over `P`.

Immersion is the exact hypothesis which prevents a further cancellation.
The three tangent coordinates have degree `2d-2`; at a nonimmersed
normalization address they acquire a common factor, and the reduced dual
map divides it out.  No such common tangent factor exists here.

The degree ledger is consequently

```text
deg(2D_P)+deg Ram(phi_P)
 =2m+(2e-2)
 =2d-2
 =deg (nu_dual)^*H_P.                                        (8)
```

## 2. Riemann--Hurwitz makes the criterion exact

Suppose the left side of `(5)` is supported at one point `p`.  Both effective
summands in `(4)` must then be supported at `p`.  If `e>=2`, Riemann--Hurwitz
gives

```text
deg Ram(phi_P)=2e-2,                                          (9)
```

whereas the ramification coefficient at one point is at most `e-1`.  Since
`2e-2>e-1`, a degree-`e` self-map of `P1` cannot have all its ramification at
one point.  Hence `e=1`, the ramification divisor vanishes, and `(2)` gives
`m=d-1`.  The support assumption then says `D_P=(d-1)p`.

Conversely, if `D_P=(d-1)p`, then `e=1`, so `(4)` is simply

```text
(nu_dual)^*H_P=2(d-1)p.                                      (10)
```

This proves `(5)`.  Finally, when the fibre of `nu` over `P` has only the
address `p`, immersion says that at least one of the two local coordinates
through `P` has a nonzero first derivative.  Thus the common vanishing order
in `(2)` is one.  Equality `d-1=1` forces `d=2`.  The conic/tangent-line case
is the sharp positive boundary.

## 3. The entire six-cusp/four-node sextic packet is two-place

Let `Gamma` be an irreducible plane sextic having exactly six ordinary `A2`
cusps and four ordinary `A1` nodes and no other singularities.  The classical
Pluecker formulas give

```text
genus(Gamma) = (5*4)/2 - 6 - 4 = 0,
deg(Gamma_dual) = 6*5 - 3*6 - 2*4 = 4,
inflection weight = 3*6*4 - 8*6 - 6*4 = 0.                  (11)
```

Thus `Gamma_dual` is a rational quartic, and the zero inflection weight says
that its normalization map is everywhere immersive.  Biduality identifies
`Gamma` with the dual of that immersed quartic.  Applying `(5)` with `d=4`
shows that every line section of `Gamma` has at least two normalization
support points.

Therefore no member of the whole equisingularity class

```text
degree 6;  singular packet 6A2+4A1                            (12)
```

has an affine chart whose normalization is `A1`.  This is independent of a
chosen torus equation, contact conic, or coordinates.  A one-place
counterexample design cannot preserve the exact THM-3879 singular packet.

## 4. Why the THM-3879 best line has type (3,3)

For the trinodal quartic in THM-3879, the dual line `C=0` corresponds to the
primal node `P=[0:0:1]`.  Its two projection coordinates have the factorization

```text
X=ST(S^2-2T^2),                  Y=ST(S^2-T^2).               (13)
```

Hence

```text
D_P=[S=0]+[T=0],              e=2.                            (14)
```

The residual degree-two projection is ramified once at each of those two
addresses.  Formula `(4)` gives

```text
(nu_dual)^*(C=0)=3[S=0]+3[T=0],                               (15)
```

exactly the pullback `2S^3T^3` found in THM-3879.  The line's two-place
optimality is therefore not an accidental sparse coefficient pattern: its
two node branches each carry one base unit twice and one ramification unit.

## 5. Counterexample-search consequence and scope

This theorem constructs no Keller map and proves no case of `JC(2)`.  It
does remove a full geometric architecture from the search: varying the
contact conic or torus presentation inside the `6A2+4A1` sextic packet can
never repair the one-place defect.

There are now only three honest ways around this obstruction:

1. make the primal normalization nonimmersive, so a common tangent factor is
   cancelled from `(7)`;
2. change the dual/singularity packet, necessarily paying a different
   Pluecker and Cardano divisor ledger; or
3. abandon the dual construction and build the one-place branch directly.

Any such move must still preserve the connected codimension-one-unramified
`C3` packet of THM-3879 while avoiding its globally monogenic cubic order.

Run

```text
python3 04-computation/jc2_rational_dual_one_place_wronskian_projection_criterion_thm3882.py
python3 -O 04-computation/jc2_rational_dual_one_place_wronskian_projection_criterion_thm3882.py
```

and compare both streams byte-for-byte with the frozen output.
