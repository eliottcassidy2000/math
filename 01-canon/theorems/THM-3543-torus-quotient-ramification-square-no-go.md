---
id: THM-3543
title: "Torus-quotient ramification square and planar Keller compression no-go"
status: >
  PROVED + VERIFIED-EXACT + INDEPENDENTLY HOSTILE-AUDITED WITH SCOPE REPAIRS.
  The categorical
  two-dimensional torus quotient of the fixed THM-1300 Keller map is an
  explicit noninjective polynomial map whose Jacobian is the square
  2(2-3v-t)^2.  Every polynomial pair of target invariants factors through
  this quotient, so its pullback Jacobian is either zero or remains divisible
  by that square; polynomial reparametrization cannot turn the quotient into
  a planar Keller map.  This excludes invariant polynomial compression of
  this fixed map, not non-invariant surfaces or JC(2).
source: kps-s183
depends_on:
  - THM-1300-jacobian-counterexample-dixmier-A3-explicit
related:
  - THM-1305-jacobian-equivariant-anatomy-rigidity-decode
  - THM-1345-jc2-equivariant-category-poisson-reframing-dc1-shadow
companion: 04-computation/jacobian_torus_quotient_ramification_square_kps_s183.py
output: 05-knowledge/results/jacobian_torus_quotient_ramification_square_kps_s183.out
script_sha256: 9066e3319ba74c31408cbd505e1ae0f849de00d0357582cafb4ccf4aeaf13e51
output_sha256: 3d9388a73d7bd3308098cba3f2754be9b356f13c095c1f07fd595e03c65345c6
hash_basis: LF-normalized bytes
---

# THM-3543 -- the torus quotient keeps the collision and pays a ramification square

**PROVED + VERIFIED-EXACT + INDEPENDENTLY HOSTILE-AUDITED WITH SCOPE
REPAIRS.**  The most direct
two-dimensional shadow of the fixed three-dimensional Keller
counterexample is indeed noninjective, but it is not Keller.  Its exact
Jacobian is a square supported on the orbit-collapse divisor.  More strongly,
that square divides the Jacobian of every invariant polynomial compression.

The field has characteristic zero throughout.

## 1. Both categorical quotients are polynomial planes

Write the source coordinates as `(x,y,z)` and give them weights

```text
                         (1,-1,-2).                    (1)
```

An invariant monomial satisfies `a-b-2c=0`, hence

```text
x^a y^b z^c=(xy)^b(x^2z)^c.
```

Therefore the source invariant ring is exactly

```text
k[x,y,z]^(G_m)=k[v,t],       v=xy,       t=x^2z.       (2)
```

Write the target coordinates as `(alpha,beta,gamma)`.  Their weights under
THM-1300's equivariance are `(-2,-1,1)`.  An invariant monomial satisfies
`-2i-j+k=0`, so it is

```text
(alpha gamma^2)^i(beta gamma)^j.
```

Thus the target invariant ring is also a polynomial ring,

```text
k[alpha,beta,gamma]^(G_m)=k[T,V],
T=alpha gamma^2,             V=beta gamma.             (3)
```

These are the categorical quotient coordinates; no localization or choice of
gauge is hidden in `(2)` or `(3)`.

## 2. Exact quotient map

For the fixed map `F=(F_1,F_2,F_3)` of THM-1300 put `u=1+v` and

```text
A(v,t)=t u^3+v^2u(4+3v),
B(v,t)=v+3tu^2+3v^2(4+3v),
C(v,t)=2-3v-t.                                      (4)
```

Direct expansion gives the three polynomial identities

```text
x^2F_1=A,                 xF_2=B,                 F_3=xC.  (5)
```

Although `(5)` is suggested by the chart `x != 0`, it is an identity in the
whole polynomial ring.  Consequently the induced quotient map
`G:A^2_(v,t)->A^2_(V,T)` is

```text
G(v,t)=(BC,AC^2).                                      (6)
```

Its Jacobian is

```text
det d(V,T)/d(v,t)=2C^2=2(2-3v-t)^2.                   (7)
```

The companion verifies `(4)`--`(7)` by exact symbolic arithmetic under both
normal and optimized Python execution.

## 3. The triple collision meets a collapsed quotient line

The fixed branch point and the two torus-orbit points of THM-1300 have source
quotients

```text
q_0=(0,0),                 q_1=(-3/2,13/2),            (8)
```

respectively.  The two orbit points already coincide at `q_1`.  Exact
substitution gives

```text
G(q_0)=G(q_1)=(0,0),
det JG(q_0)=8,             det JG(q_1)=0.              (9)
```

In fact the whole line `C=0`, equivalently `t=2-3v`, maps to `(0,0)` under
`G`.  Along this line,

```text
B=4v+6,                    A=(v+1)(v+2),               (10)
```

so `q_1` is the unique point there at which the unfactored covariant `B`
also vanishes, while `A(q_1)=-1/4` remains a unit.  Thus the quotient does
not merely retain the distinguished fixed/orbit collision: it contracts an
entire unstable divisor.  Every polynomial postcomposition still contracts
that line, giving a geometric proof of the no-go below independently of the
Jacobian factorization.

The ramification is already present in the tangent cone.  In local coordinates

```text
r=v+3/2,                    c=C,                        (11)
```

the quadratic initial map at `q_1` is

```text
(V_2,T_2)=(4cr-3c^2/4,-c^2/4),                        (12)
```

whose Jacobian is `-2c^2`, exactly the full local Jacobian after the
orientation-reversing change `(v,t)->(r,c)`.  Thus the factor in `(7)` is not
a high-degree artifact: the quotient turns the orbit doubling into a
quadratic branch at first order.

## 4. Invariant-compression no-go

Let `H=(h_1,h_2)` be any pair of target-invariant polynomials.  By `(3)` there
are polynomials `tilde h_i(V,T)` with

```text
h_i=tilde h_i(beta gamma,alpha gamma^2).               (13)
```

The pair pulled back along `F`, viewed on the source quotient, is therefore
`tilde H o G`.  The chain rule and `(7)` give

```text
det J(tilde H o G)
 = (det J tilde H)(G(v,t)) * 2C(v,t)^2.                (14)
```

If `tilde h_1,tilde h_2` are algebraically dependent, the left side is zero.
If they are algebraically independent, the first factor in `(14)` is nonzero
in characteristic zero, but the nonconstant square `C^2` still divides the
Jacobian.  In neither case is the Jacobian a nonzero constant.

In particular:

1. no pair of invariant target observables compresses this map to a planar
   Keller map;
2. postcomposition by a polynomial automorphism of the target quotient only
   multiplies `(7)` by a nonzero constant; and
3. polynomial automorphisms of the source quotient merely move the divisor
   `C=0` and likewise cannot remove it.

This is the precise obstruction: categorical compression retains the
noninjectivity only by forgetting the character coordinate, and forgetting
that coordinate converts the orbit doubling into ramification.

## 5. Boundary and next counterexample lane

This theorem does **not** prove planar JC, and it does not exclude:

- a non-invariant polynomial surface or graph in `A^3` on which `F` induces a
  planar map;
- a semiconjugate surface that is not the categorical quotient;
- a branched affine modification carrying an additional root/owner sidecar;
- a different three-dimensional counterexample with a transverse invariant
  surface; or
- stable Jacobian/Dixmier constructions that change dimension and map.

Polynomiality is load-bearing.  The rational target recombination

```text
(R,S)=(T/V^2,V)
```

pulls back to `(A/B^2,BC)` and has exact Jacobian `-2/B^2`; it removes the
`C^2` factor at the price of a pole.  Likewise, the theorem concerns the full
torus invariant rings, not the much larger invariant rings of the single
involution `lambda=-1`, and it says nothing about arbitrary semi-invariants.

It does show what any descent of this particular collision must repair.  On
the quotient, a would-be coordinate change needs the impossible polynomial
Jacobian factor `C^(-2)`.  A viable planar-counterexample search should
therefore keep a transverse character or normal multiplier until the final
step and test whether an invariant graph can absorb the square before it is
discarded.  That is a concrete search program, not a claimed construction.

The independent hostile audit reconstructed both invariant rings, the exact
quotient and chain-rule factorization, found the whole-line contraction and
rational escape above, and forced the target-invariant/full-torus/polynomial
scope stated here.
