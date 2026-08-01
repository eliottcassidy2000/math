---
id: THM-3038
title: "Split monogenic order cross-resultant conductor and affine-owner boundary"
status: >
  PROVED CANDIDATE + VERIFIED-EXACT; AWAITING INDEPENDENT HOSTILE AUDIT.
  For a split monogenic quartic order R[T]/((T-a)g), with R a normal
  local domain, g generically separable, and d=g(a) regular, the order is the fibre product of the
  singleton and complementary cubic orders over R/(d).  Its singleton
  conductor slice, normalization quotient, and matching-resolvent order
  defect are the same cyclic module R/(d).  Hence the singleton idempotent
  splits exactly when d is a unit.  This calculation applies literally to
  the monogenic order; a larger graph order requires its own module-membership
  test.  Even a split finite-order idempotent does not prove affine ownership:
  ownership is regularity of every original source coordinate along the
  boundary trace.  No quartic Keller, G1, JC(2), DC(2), or LRC exclusion follows.
source: codex-split-order-conductor-2026-08-01
depends_on:
  - THM-2598-quartic-v4-resolvent-torsor-and-universal-cusp-boundary
  - THM-3037-cusp-braid-s4-lift-dichotomy-and-common-sheet-owner-boundary
related:
  - THM-2455-quartic-swallowtail-scaffold-and-endpoint-corrections
  - THM-2570-jelonek-cusp-cylinder-normalization-and-conductor
  - THM-2633-derangement-character-obstruction-and-d4-keller-exclusion
  - THM-2974-discriminant-cover-integral-order-smith-and-owner-boundary
script: 04-computation/split_monogenic_order_cross_resultant_conductor_thm3038.py
output: 05-knowledge/results/split_monogenic_order_cross_resultant_conductor_thm3038.out
script_sha256: 6042adbe39bf5b33bd4eeaa5f59836b669460210d883f92dbcde018eed4ee5a4
output_sha256: b8231fae9dd680c7e8c95a37f84757195585f76d75a8fdc97a1315a074575b7c
hash_basis: LF-normalized bytes
---

# THM-3038 -- the cross-resultant measures sheet gluing, not affine ownership

**PROVED CANDIDATE + VERIFIED-EXACT; AWAITING INDEPENDENT HOSTILE AUDIT.**

## 1. Inheritance and statement

[THM-3037, the cusp-braid `S4` lift
dichotomy](THM-3037-cusp-braid-s4-lift-dichotomy-and-common-sheet-owner-boundary.md)
leaves a unique common quartic sheet in the local split branch, but warns that
the sheet belongs first to the finite normalization.  [THM-2974, the
integral-order Smith
boundary](THM-2974-discriminant-cover-integral-order-smith-and-owner-boundary.md)
shows independently that an isomorphism of normalized cover algebras need not
identify their graph orders or an affine owner.  The present theorem computes
the exact missing order defect when the quartic is monogenically split, and
then separates it from affine ownership.

Let `R` be a normal local domain, and suppose `g` is squarefree over
`Frac(R)` (equivalently, the finite cubic algebra below is generically
reduced).  Let

```text
g(Y)=Y^3+pY^2+qY+r in R[Y],       a in R,
d=g(a) a nonzero divisor,
C=R[Y]/(g),                       A=R[Y]/((Y-a)g).       (1)
```

Let `epsilon:C -> R/dR` be evaluation at `Y=a`.  Evaluation at the singleton
root and reduction modulo `g` identify the quartic order with

```text
A = R x_(R/dR) C.                                      (2)
```

If `Ctilde` is the normalization of `C` in its total fraction algebra, then

```text
Atilde=R x Ctilde.                                     (3)
```

For the singleton idempotent `e=(1,0)`, one has

```text
A intersect Re = conductor_(Atilde/A) intersect Re = dRe. (4)
```

Consequently

```text
e in A  iff  e is in the conductor  iff  d is a unit,   (5)
e=g(Y)/g(a) when these conditions hold.                 (6)
```

The cross-resultant is therefore an exact gluing coordinate:

```text
(Res_Y(Y-a,g))=(g(a))=(d).                              (7)
```

Over a DVR, the quotient and index formulas are

```text
(R x C)/A = R/(d),
length(Atilde/A)=length(Ctilde/C)+v(d),                 (8)
Disc((Y-a)g)=Disc(g)d^2.                                (9)
```

There is a second, canonical occurrence of the same cyclic defect.  The
standard matching resolvent of the split quartic embeds into the complementary
cubic order, and its cokernel is again `R/(d)`.  Thus raw
quartic--resolvent discriminant equality duplicates the gluing tax; it does
not cancel it.

Finally, suppose a polynomial Keller map has a Zariski-main open immersion

```text
j:A^n -> Xbar                                             (10)
```

into its finite normalization, and let `x_1,...,x_n` denote the original
source coordinates as rational functions on `Xbar`.  Then

```text
j(A^n)=intersection_i Dom(x_i).                          (11)
```

A singleton boundary trace is therefore an actual affine-source owner exactly
when all source coordinates are regular there.  At a height-one trace `eta`,
this is the finite valuation test

```text
v_eta(x_i) >= 0 for every i.                             (12)
```

Equations (5) and (12) are independent gates.

## 2. The fibre product and conductor

Consider the map

```text
Phi:A -> R x C,       h |-> (h(a),h mod g).              (13)
```

Its image satisfies the fibre-product congruence by construction.  If a
polynomial representative lies in the kernel, write it as `sg`.  Evaluation
at `a` gives `s(a)d=0`; regularity of `d` forces `s(a)=0`, so `s` is divisible
by `Y-a`.  Hence the kernel before quotienting is exactly `((Y-a)g)`.

Conversely, let `(b,h mod g)` be compatible.  Then

```text
b-h(a)=kd                                                   (14)
```

for some `k in R`, and `h+kg` maps to the required pair.  This proves (2).
Generic separability gives the reduced total fraction algebra, and normality
of `R` plus its product decomposition gives (3).  The fibre-product and
cross-resultant statements themselves remain valid without generic
separability; that hypothesis is needed precisely for the normalization and
index language.

The fibre-product congruence says `(b,0)` lies in `A` exactly when `b in dR`,
which proves the first equality in (4).  Every element of `dRe` multiplies
`Atilde` into `A`.  Conversely the conductor is contained in `A`, proving the
second equality.  Formula (5) and the explicit idempotent (6) follow.

This last equivalence survives for any finite `R`-subalgebra

```text
A_graph subset R x B                                      (15)
```

that contains diagonal `R`: if `e in A_graph`, then
`e(R x B)=Re subset A_graph`, while the conductor is always contained in the
order.  What does **not** survive is the formula `A_graph intersect Re=dRe`.
A larger graph order may contain `e` through generators absent from the
monogenic suborder.  Thus:

- `d` a unit is always sufficient for the full graph order, because the
  monogenic suborder already contains `e`;
- `d` a nonunit proves absence only from the monogenic suborder;
- the actual full graph order requires direct module membership, by Smith
  form over a DVR or an appropriate local Groebner/module calculation in
  higher dimension.

This is the first truth-surface boundary of the theorem.

## 3. The same defect inside the matching resolvent

Let `y,z,w` be the roots of `g`, so the quartic roots are `a,y,z,w`.  In the
matching containing `{a,y}` and `{z,w}`, the standard resolvent coordinate is

```text
omega_y=ay+zw=y^2+(a+p)y+q.                              (16)
```

Therefore the standard integral matching-resolvent order maps canonically to
the complementary cubic order:

```text
O_res=R[W]/(R_f(W)) -> C,       W |-> omega_y.            (17)
```

In the power bases `(1,W,W^2)` and `(1,y,y^2)`, the matrix of (17) is

```text
[ 1   q       q^2-(2a+p)r ]
[ 0  a+p          pq-r    ]
[ 0   1           a^2+q   ].                             (18)
```

Its determinant is `g(a)=d`.  More sharply, its lower `2 x 2` block is
elementary-equivalent over `R` to

```text
diag(1,-d).                                               (19)
```

The literal unit entry in (18) makes this an exact Smith calculation over any
local domain, not merely an equality of determinants.  Hence

```text
C/O_res = R/(d).                                         (20)
```

Equations (8) and (20) exhibit one cyclic defect module in two places:

```text
quartic singleton gluing  <---- R/(d) ---->  resolver-coordinate loss.
                                                               (21)
```

In particular,

```text
O_res=C  iff  d is a unit  iff  e splits in A.            (22)
```

If `R` is a DVR and

```text
i_comp=length(Ctilde/C),
i_res =length(Ctilde/O_res),
i_4   =length(Atilde/A),                                  (23)
```

then

```text
i_res=i_comp+v(d)=i_4.                                   (24)
```

The equality of raw quartic and resolvent discriminants is the familiar
polynomial identity from [THM-2598](THM-2598-quartic-v4-resolvent-torsor-and-universal-cusp-boundary.md).
Formula (24) explains its integral meaning in the split chart: both graph
orders carry the same extra index.  Comparing those two raw discriminants or
indices cannot prove that `d` is a unit.

## 4. The affine-owner test

Let

```text
A^n --j--> Xbar --Fbar--> A^n                            (25)
```

be the Zariski-main factorization of a polynomial Keller map.  The original
source coordinates define rational functions `x_i` on `Xbar`.  Their common
domain `U=intersection Dom(x_i)` contains `j(A^n)`, and on `U` they define a
morphism `psi:U -> A^n`.  The maps `j psi` and `id_U` agree on the dense open
`j(A^n)`.  Since `Xbar` is separated, they agree on all of `U`; hence every
point of `U` equals `j(psi(point))`.  This proves (11).

For a height-one prime `eta` on the normal space `Xbar`, membership
`x_i in O_(Xbar,eta)` is equivalent to (12).  At a higher-codimension cusp
point one must test the local ring itself, equivalently all height-one
valuations through that point.  This condition is not detected by the finite
order conductor.

The distinction is structural:

```text
finite-order section       asks whether e is integral in the graph order;
affine-source owner        asks whether every original x_i is regular there.
                                                               (26)
```

The second can fail after the first succeeds.

## 5. Sharp hostiles

### 5.1 Raw discriminant equality with a glued singleton

Over `R=C[[t]]`, put

```text
g=T^3-T+t,             f=Tg=T^4-T^2+tT,             d=t. (27)
```

The complementary cubic is etale because

```text
Disc(g)=4-27t^2                                             (28)
```

is a unit.  Nevertheless the singleton idempotent is `g/t`, has a genuine
`t` denominator, and is absent from the monogenic quartic order.  Its
conductor slice is `tRe`.  The standard resolvent is

```text
W^3+W^2-t^2,                                               (29)
```

and both raw discriminants equal

```text
t^2(4-27t^2).                                              (30)
```

Thus normalized complementary sheets and raw quartic--resolvent equality do
not split the singleton.  The missing quantity is precisely the common module
`R/(t)` from (21).

### 5.2 Unit conductor with opposite affine-owner choices

In the completed local cusp model, take

```text
(T-1)(T^3-3uT-2v),          Delta: v^2=u^3.              (31)
```

Here

```text
g(1)=1-3u-2v                                                    (32)
```

is a unit, so the singleton idempotent already belongs to the order and its
conductor.  On `u=s^2,v=s^3`, the cubic becomes

```text
(T-2s)(T+s)^2.                                             (33)
```

Start from the same finite normalization and remove the cubic ramification
locus.  One may either retain the singleton trace above `Delta`, or also
delete that trace while retaining the simple cubic branch `T=2s`.  Both opens
retain all four sheets away from `Delta` and are etale over their images.  The
finite order and conductor are identical, but the singleton owner bit is
opposite.  This is a local Zariski-main hostile, not a polynomial Keller map.

### 5.3 A genuine source-coordinate pole

[THM-2570, the Jelonek cusp-cylinder
normalization](THM-2570-jelonek-cusp-cylinder-normalization-and-conductor.md)
has normalized survivor coordinate

```text
x=2c/theta^2.                                             (34)
```

Its valuation at the conductor prime `theta=0` is `-2`, so (12) excludes
that trace from the affine source.  The finite sections `theta=2,-1` instead
have regular source coordinates.  This is the concrete model for the owner
gate left open by THM-3037.

## 6. Transfer contract and scope

| source | target / map | preserved | destroyed or undecided | required sidecar |
|---|---|---|---|---|
| split quartic monogenic order | fibre product (2) | exact singleton gluing | larger graph generators | full-order membership |
| complementary cubic | matching resolver (16) | normalized algebra and explicit order inclusion | cyclic defect `R/(d)` | unit cross-resultant |
| finite normalization section | affine source via (11) | rational inverse coordinates | boundary ownership | source-coordinate regularity |
| local split `S3` | connected global degree four | global `S4` monodromy by THM-3037 | global owner | boundary incidence |

The `C2*C3`/Farey interpretation is now exact but limited: the local braid
chooses the split `H^1(C2*C3,V4)` class, while (21) measures whether that
finite sheet is separated integrally.  Neither datum supplies the affine
owner valuation (12).  No LRC role, comb, current, or carrier is canonically
identified by this comparison.

No quartic Keller map, forbidden `G1` degree, JC(2), DC(2), SFC(4), or LRC(14)
row is excluded here.  The next physical quartic test has two ordered stages:

1. compute singleton membership in the **full** completed graph order;
2. if it exists, evaluate all original source-coordinate valuations along
   its boundary trace.

## 7. Exact companion

The companion verifies, with explicit runtime checks in ordinary and
optimized Python:

- the cross-resultant and split discriminant identities;
- the coefficient-lattice determinant of the fibre product;
- the singleton-idempotent numerator relation;
- the standard quartic/resolvent discriminant identity;
- the complementary-cubic resolver formula, determinant, and cyclic Smith
  reduction;
- the `C[[t]]` order-mismatch hostile and its actual denominator;
- the cusp-owner factorization and simple/double branches; and
- THM-2570's exact source-coordinate pole order.

The ordinary transcript, optimized transcript, and stored output agree
byte-for-byte after LF normalization.
