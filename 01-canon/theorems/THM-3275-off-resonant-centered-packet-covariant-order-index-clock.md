---
id: THM-3275
title: "Off-resonant centered-packet covariant order-index clock"
status: >
  PROVED + VERIFIED-EXACT.  In the 3-not-dividing-m lane of the actual
  THM-3201 fixed-plus-escaping-C3 graph quartic, put s=pi^(3m) in K and let
  eta=sN(z) be the normalized THM-3271 packet value on the moving cubic
  field.  The fixed/cubic packet cross-resultant is a unit.  Either eta is
  in K, equivalently the THM-3273 collision covariant C is zero, or, writing
  h=w(eta-Tr(eta)/3), one has h>0, 3 does not divide h, and the exact laws
  v_K(s^4C)=h, v_K(disc(mu_eta))=v_K(disc(P_tilde))=2h, and
  length_O_K(O_L/O_K[eta])=h-1.  Thus the first nonbase Laurent character,
  the normalized collision-covariant order, and the moving packet-order
  defect are the same integer in three coordinate systems.  At h=1 the
  normalized covariant is nonunit but the packet order is already maximal;
  an exact h=2 graph-anatomy model has index length one.  This measures the
  order/different correction but not its normalized unit, the true
  chain-rule cofactor, a C3 exclusion, or JC(2).
source: jc-packet-order-clock-2026-08-02
depends_on:
  - THM-3274-graph-quartic-centered-packet-fixed-decoder-and-cofactor-independence
  - THM-3273-quartic-centered-norm-packet-collision-covariant-and-discriminant-factor
  - THM-3201-c3-local-resolvent-splitting-and-matching-newton-gate
related:
  - THM-3272-tame-pure-c3-centered-norm-packet-integral-fixed-projector
  - THM-3064-pointed-cubic-norm-keller-decoder-and-inverse-different-boundary
  - THM-3038-split-monogenic-order-cross-resultant-conductor-and-affine-owner-boundary
script: 04-computation/jc_off_resonant_packet_covariant_conductor_clock_thm3275.py
output: 05-knowledge/results/jc_off_resonant_packet_covariant_conductor_clock_thm3275.out
script_sha256: 3cd3287c5fe93d1ea74fdf5c3fff4a4d76bc3768112a65d45468a7a4aa575b32
output_sha256: 900487b852f0ac2b5e7c4e5840d47177645f3a4a79101a26fe06c83122fa952e
hash_basis: LF-normalized bytes
---

# THM-3275 -- the off-resonant packet covariant is an exact order clock

**PROVED + VERIFIED-EXACT.**

## 1. Local normalization and the moving packet order

Let `(K,v)` be a complete discretely valued field of residue characteristic
different from two and three, with

```text
v(t)=1,                                                   (1)
```

and assume `K` contains a primitive cube root `zeta`.  Let `L/K` be a tame,
totally ramified cyclic cubic extension.  Normalize its value `w` by

```text
w(pi)=1,                 w restricted to K=3v,
sigma(pi)=zeta pi,       pi^3=epsilon t,                 (2)
```

where `epsilon` is a unit.

Now retain the actual fixed-plus-escaping graph-quartic anatomy of
[THM-3201](THM-3201-c3-local-resolvent-splitting-and-matching-newton-gate.md)
in the lane

```text
3 does not divide m.                                    (3)
```

The completed quartic algebra is `K x L`.  Write `z_f in K` for the fixed
depressed root and `z in L` for one moving root.  With the packet map

```text
N(X)=-(20X^3+18pX+27q)/27                               (4)
```

from [THM-3271](THM-3271-universal-quartic-centered-norm-packet-and-local-singleton-projector.md),
put

```text
s=pi^(3m)=(epsilon t)^m in K,
eta_f=sN(z_f) in O_K,             eta=sN(z) in O_L.     (5)
```

THM-3274 gives the residual packet factorization

```text
bar(P_tilde)(W)=(W-A^3)(W-7A^3/27)^3,
P_tilde(W)=s^4P_N(s^(-1)W).                             (6)
```

In particular

```text
eta_f-sigma^i(eta) is a unit for i=0,1,2.               (7)
```

Thus the quartic packet discriminant has positive valuation only because
the three moving packet values collide modulo `t`; the fixed/cubic
cross-resultant contributes no valuation.

## 2. The first nonbase character

For any integral `eta in L`, define its base projection and nonbase depth by

```text
c=(1/3)Tr_(L/K)(eta) in O_K,
h=w(eta-c),                                              (8)
```

provided `eta` is not in `K`.  In the Kummer basis, uniquely write

```text
eta=c+a_1 pi+a_2 pi^2,                  a_1,a_2 in K.   (9)
```

The two possible nonzero term values in `(9)` are respectively congruent to
one and two modulo three, so they cannot tie.  Therefore

```text
h>0,                         3 does not divide h.       (10)
```

Moreover every conjugate difference has exactly this value:

```text
w(eta-sigma(eta))
 =w(eta-sigma^2(eta))
 =w(sigma(eta)-sigma^2(eta))=h.                         (11)
```

Indeed the unique lowest term from `(9)` is multiplied by one of
`1-zeta`, `1-zeta^2`, or `zeta-zeta^2`, all units.  This is the intrinsic
meaning of `h`: it is the first Laurent scale at which the normalized packet
remembers a nontrivial `C3` character.

If `eta` is not in `K`, it generates the prime-degree extension `L/K`.  Let

```text
mu_eta(W)=Norm_(L/K)(W-eta).                            (12)
```

Multiplying the squared differences from `(11)` gives

```text
w(disc(mu_eta))=6h,
v(disc(mu_eta))=2h.                                     (13)
```

This is the discriminant of the **moving cubic minimal polynomial**, not yet
the discriminant of the four-valued scalar packet.

## 3. Six pair quotients recover the same integer

For two depressed quartic roots, THM-3273 gives

```text
N(u)-N(v)=(u-v)Q_p(u,v),
product_(i<j) Q_p(z_i,z_j)=-(2^6/3^18)C(p,q,r).         (14)
```

The off-resonant THM-3201 initials give the following complete valuation
ledger:

```text
pair type       count   w(root difference)  w(packet difference)  w(Q)
fixed/moving      3             -m                   -3m            -2m
moving/moving     3             -m                  h-3m           h-2m.
                                                                    (15)
```

The constants in `(14)` are units.  Summing the last column of `(15)` yields

```text
w(C)=3h-12m.                                             (16)
```

Because `w(s)=3m` and `s^4C` lies in `K`, equations `(2)` and `(16)` give
the exact base-field identity

```text
v(s^4C)=h.                                               (17)
```

There is also a sharp collapsed alternative.  If `eta in K`, its three
moving packet conjugates agree, so `(14)` gives `C=0`.  Conversely, if
`C=0`, some pair quotient vanishes.  Equation `(7)` rules out a fixed/moving
pair, hence two moving packet values agree.  A nontrivial element of `C3`
then fixes `eta`, so `eta in K`.  Therefore

```text
C=0  iff  eta in K;                                     (18)
```

otherwise `(17)` holds with finite `h`.

This strengthens THM-3274's statement that `bar(s^4C)=0`: the order of that
zero is not noise.  It is exactly the first surviving nonbase packet scale.

## 4. Cubic versus quartic packet discriminants

Over the completion, the normalized quartic packet factors as

```text
P_tilde(W)=(W-eta_f)mu_eta(W).                          (19)
```

The discriminant product law is

```text
disc(P_tilde)
 =disc(mu_eta) mu_eta(eta_f)^2.                         (20)
```

By `(7)`, `mu_eta(eta_f)` is a unit.  Consequently `(13)` gives

```text
v(disc(P_tilde))=v(disc(mu_eta))=2h.                    (21)
```

Thus the reduction in `(6)` has zero discriminant, but the characteristic
polynomial itself is separable whenever `C!=0`; its positive discriminant
order is exactly `2h`.  Equation `(21)` does not confuse the inherited
moving triple collision in the special fibre with an exact collision in
characteristic zero.

## 5. Exact order index and different correction

Since `eta` is integral and generates `L`, its monogenic order has
discriminant `(13)`.  Tame total ramification of degree three gives

```text
v(disc(O_L/O_K))=3-1=2.                                (22)
```

For the full-rank inclusion `O_K[eta] subset O_L`, discriminants satisfy

```text
disc(O_K[eta])
 =disc(O_L) [O_L:O_K[eta]]^2.                          (23)
```

Taking ideal valuations in `(13)`, `(22)`, and `(23)` proves

```text
length_(O_K)(O_L/O_K[eta])=h-1.                         (24)
```

Equivalently, the derivative of the packet minimal polynomial has

```text
w(mu_eta'(eta))=2h,                                     (25)
```

whereas a generator of the maximal different has `w`-value two.  Passing
from the packet order to the maximal cubic order therefore requires a
different correction of `w`-value

```text
2(h-1).                                                  (26)
```

This is the precise sense in which `(17)` is a conductor-scale clock.  It
determines the index length and different **valuation**, not the conductor
ideal's full unit data or the normalized inverse-different cofactor.

## 6. The maximal boundary `h=1`

The boundary `h=1` is subtle and useful.  Let `K=k((t))`, where `k` has
characteristic zero and contains `zeta`, take `pi^3=t`, and let `B in k*`.
Consider the depressed local graph-anatomy model

```text
f_1(X)=(X+3B/4)((X-B/4)^3-t^(-1)).                     (27)
```

Before depression this has fixed root zero and moving roots
`B+sigma^i(pi^(-1))`, so it is the `m=1` off-resonant anatomy.  Its normalized
fixed packet value is one.  On a moving root,

```text
eta_1=(7-15Bpi+3B^2pi^2-2B^3t)/27,
h=1.                                                     (28)
```

Exact elimination gives

```text
disc(mu_eta1)
 =-B^6t^2(B^3t+125)^2/19683,                            (29)

t^4C
 =27B^3t(B^3t+125)(8B^6t^2+475B^3t+8000)/64.           (30)
```

Thus the respective valuations are two and one, as `(13)` and `(17)`
predict.  Equation `(24)` gives index length zero:

```text
O_K[eta_1]=O_L.                                         (31)
```

This is a sharp hostile to the tempting implication

```text
normalized C is a nonunit  =>  packet order is nonmaximal.                (32)
```

The normalized covariant must vanish off resonance; a **simple** zero is
exactly the maximal-order boundary.

## 7. A nonmaximal `h=2` hostile

In the same field, take

```text
f_2(X)=X(X^3-3BX-(t^(-1)+B^3t)).                       (33)
```

Its moving root `pi^(-1)+Bpi` has the two opposite `C3` characters.  The
normalized fixed and moving packet values are

```text
eta_f=1+B^3t^2,
eta_2=(7+7B^3t^2-6Bpi^2-6B^2tpi)/27,
h=2.                                                     (34)
```

This time exact elimination gives

```text
disc(mu_eta2)
 =-64B^6t^4(B^3t^2-1)^2/19683,                         (35)

t^4C
 =-27B^3t^2(1000B^6t^4+1757B^3t^2+1000).              (36)
```

Hence

```text
v(t^4C)=2,          v(disc(mu_eta2))=4,
length(O_L/O_K[eta_2])=1.                              (37)
```

Both `(27)` and `(33)` retain the same residual packet
`(W-1)(W-7/27)^3` and a unit fixed/cubic cross-resultant.  They prove that
the integer `h`, invisible in the reduced factorization alone, genuinely
distinguishes maximal from nonmaximal moving packet orders.

These are exact local graph-anatomy models, not polynomial Keller maps or
global `S4` realizations.

## 8. Frontier consequence

The off-resonant information ledger is now

```text
first nonbase packet character depth h
   = normalized collision-covariant order v(s^4C)
   = 1 + packet-order index length
   = one half of the packet discriminant order.         (38)
```

This exposes two concrete next gates.

1. If actual graph incidence forces `v(s^4C)=1`, then the centered packet is
   already a maximal integral primitive for the moving cubic field.  The
   remaining obstruction is purely the normalized inverse-different unit
   isolated by THM-3274/3064.
2. If `v(s^4C)>1`, the positive index length in `(24)` is an additional
   conductor correction that a true graph cofactor must traverse.  Its
   valuation is now known exactly from `C`, so only its unit and affine-owner
   incidence remain unmeasured.

Neither case by itself excludes the `C3` branch.  In particular, `(17)`
cannot repair the norm-one cofactor hostile of THM-3274: that twist preserves
all valuations, including `h` and `(26)`.

## 9. Exact reproduction

Run

```bash
python3 04-computation/jc_off_resonant_packet_covariant_conductor_clock_thm3275.py
python3 -O 04-computation/jc_off_resonant_packet_covariant_conductor_clock_thm3275.py
```

Both modes must reproduce
`05-knowledge/results/jc_off_resonant_packet_covariant_conductor_clock_thm3275.out`
byte for byte.  The companion checks thirteen abstract tame-`C3`
discriminant/index cases, one hundred pair-valuation ledgers, both exact
quartic packet factorizations, the cubic/quartic discriminant separation,
unit cross-resultants, both normalized collision covariants, and the `h=1`
and `h=2` index boundaries.
