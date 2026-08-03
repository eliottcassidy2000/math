# Degree seven reopens the cubic wall, then retunes to an odd quartic escape

**Status:** VERIFIED-EXACT WILDCARD ARTIFACT.  This is not a canonical theorem,
does not reserve an ID, and does not change the planar-JC ledger.

## Inheritance pass

- Closest proved mechanism: THM-3237 creates the first degree-nine infinity
  face, and THM-3257 restricts to its primitive wall and uses a dormant
  degree-eight coefficient to expose a cubic reciprocal Newton edge.
- Corrected near miss: simply adding `r*x^7` while retaining THM-3257's old
  `s_*` does **not** tune the cubic carry.  It reopens the preceding `x^97`
  coefficient.
- Canonical hostile: THM-3068 shows that local reciprocal inertia and exact
  Laurent ledgers do not supply a polynomial affine source or the required
  reciprocal cofactor unit.
- Least-used relevant sidecar: THM-3064's pointed branchwise
  inverse-different unit.  THM-3066 proves that a symmetric product quotient
  cannot replace it.

The live concept board was: reciprocal resultant jet; strict-transform wall;
local inertia sign; residual critical locus; branchwise Keller cofactor.

## Exact result

Use THM-3237's response pair and write

```text
V=a^2*x^16*(1+r1/x+r2/x^2+r3/x^3+r4/x^4+...),
B=1+r*x^7+s*x^8+t*x^9.
```

Retain THM-3257's normalized invariants

```text
Kappa=4*r2-3*r1^2,
Rho=r1^4+16*r1*r3-16*r2^2-8*r1,
```

and define

```text
Lambda=4*r2-r1^2,

Sigma=
 128-128*r3+256*r2*r4-256*r3^2-64*r1^2*r4
 +128*r1*r2*r3-64*r2^3-32*r1^3*r3
 +48*r1^2*r2^2-12*r1^4*r2+r1^6.
```

At the primitive wall `t=a=2/Gamma`, the top two coefficients still vanish,
but the next one is

```text
[x^97]K=-2*a^11*(a*Kappa+4*r1*s-8*r).                    (1)
```

Consequently, at THM-3257's old value

```text
s_old=-a*Kappa/(4*r1),
```

one has

```text
[x^97]K=16*a^11*r.                                       (2)
```

Equation (2) is the minimal refutation of the direct fixed-`s` probe.  Degree
seven is transverse to, rather than tangent to, the old cubic wall.

The repaired strict transform retunes

```text
s=s(r)=(8*r-a*Kappa)/(4*r1).                             (3)
```

Along (3), the next coefficient is exactly

```text
[x^96]K=-a^12*(Rho+8*r*Lambda/a)/r1.                     (4)
```

If `Lambda!=0`, there is one cancelling value

```text
r_*=-a*Rho/(8*Lambda),
s_*=(8*r_*-a*Kappa)/(4*r1).                              (5)
```

At (5), the next coefficient is

```text
[x^95]K=-3*a^12*Sigma/(8*Lambda).                        (6)
```

Both accessory fields have nonzero `Lambda` and `Sigma`, so (5)--(6) are
defined and nonzero under every characteristic-zero embedding.  Since the
passport boundary is monic of degree 44, the saturated residual has degree
`51`.

## Exact accessory data

| passport | `Norm(Lambda)` | `Norm(Sigma)` | `Norm(r_*)` | `Norm(s_*)` | residual digest |
|---|---:|---:|---:|---:|---|
| `(4,1,1,1)` | `12448544/30625` | `-157881139404422797328384/28722900390625` | `2796986463061/2489708800000` | `-15118672709/6224272000` | `3c2ac8e2cbd13ee2568f2565bd780e2b0b8bb11facad583cc4b626e1f2265626` |
| `(3,2,1,1)` | `-42471424/275625` | `64559323245631953174528/28722900390625` | `8958651810434375/82564448256` | `-1836115611739609375/14090999169024` | `2808b49be3b3ac590efdd39e1015da59e620d45b11bd103271e21e82bdf6d499` |

The same good accessory reductions give:

| passport | `(p,u)` | `(t_*,s_*,r_*)` | `deg H` | `gcd(B,g)` | `gcd(H,g)` | `gcd(H,H')` | `[x^51]H` |
|---|---:|---:|---:|---:|---:|---:|---:|
| `(4,1,1,1)` | `(113,85)` | `(85,103,91)` | `51` | `0` | `0` | `0` | `41` |
| `(3,2,1,1)` | `(101,64)` | `(89,87,55)` | `51` | `0` | `0` | `0` | `9` |

Thus each characteristic-zero wall residual is squarefree and disjoint from
the inherited owner boundary.  The constant leading `y`-coefficient in the
gradient equation leaves 51 reduced transverse, hence Morse, critical points.

## Odd quartic local inertia

Fix `(r,s)=(r_*,s_*)`, put

```text
delta=t-a,              epsilon=Gamma*t-2=Gamma*delta,
w=1/x,                  Hhat=w^55*H(1/w).
```

The transverse leading derivative is still `-16*a^11`, while (6) is the
`w^4` coefficient.  Therefore

```text
Hhat=
 -16*a^11*delta
 -3*a^12*Sigma*w^4/(8*Lambda)
 +O(delta^2,delta*w,delta*w^2,delta*w^3,w^5),

w^4~-128*Lambda*delta/(3*a*Sigma)
   =-64*Lambda*epsilon/(3*Sigma).                         (7)
```

Equation (7) has four Puiseux branches and one local 4-cycle.  Its sign is
odd.  This is a new exact odd-inertia object, but it acts on escaping roots of
a **critical resultant**, not on normalization sheets of a generically finite
polynomial map.

## Connection contract and stopping boundary

```text
source:      THM-3257 reciprocal critical-resultant jet;
map:         add r*x^7, restrict to t=a, retune s=s(r), then divide/recompute;
target:      the first nonzero strict-transform coefficient and Newton edge;
preserved:   vanishing order, escaping intersection multiplicity, local inertia;
destroyed:   inverse-sheet labels, physical Jacobian, affine/Keller realization;
sidecar:     a genuine second polynomial coordinate and branchwise primitive-
             element cofactors, including their normalized inverse-different units;
hostiles:    fixed-s obstruction (2), THM-3066 product-gauge blindness, and
             THM-3068's disconnected Laurent realization;
cheapest test: compute the pointed cofactor ratio only after a lawful polynomial
               inverse cover is supplied.
```

The result does not decrement the planar-JC frontier.  Fifty-one critical
points remain, so the displayed scalar polynomial has no constant-Jacobian
mate in the current chart.  The odd 4-cycle cannot be retyped as Jelonek
inertia without the missing inverse-cover and cofactor-unit data.

## Reproduction and hashes

```text
python3 04-computation/jc_degree7_retuned_quartic_infinity_wall_wildcard_20260803.py
python3 -O 04-computation/jc_degree7_retuned_quartic_infinity_wall_wildcard_20260803.py
```

Both executions byte-match the frozen output.

```text
script_sha256 = aabdb06854ef9045d31c4fdbe76e2933282f44b49230ce574a77591af9a6df56
output_sha256 = 2ad12c392ff4c3228a7a242cffa58449a734339e102628d3e73386883a9bd0b2
```

The companion pins 14 exact script/output artifacts from THM-3237, THM-3257,
THM-3057, THM-3059, THM-3064, THM-3066, and THM-3068.  The five quartic
companions replay in normal and optimized modes; the new computation has no
assertion node, floating literal, randomness, or discovery cache.
