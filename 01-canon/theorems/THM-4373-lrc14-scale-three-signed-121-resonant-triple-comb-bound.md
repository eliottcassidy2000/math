---
id: THM-4373
title: "LRC(14) scale-three signed-(1,2,1) resonant triple-comb bound"
status: >
  PROVED + VERIFIED-EXACT. For every primitive triple of distinct positive
  odd speeds prime to three that satisfies a signed relation with absolute
  coefficient pattern (1,2,1), the scale-three failure comb has Haar measure
  at most 6/77. Equality occurs only for (1,5,11), up to permutation. An exact
  endpoint-pair parametrization and a period-three quadrature identity reduce
  the infinite family to a uniform pq>=81 estimate and thirteen explicit
  smaller endpoint pairs. This proves the resonant sector only; the universal
  nonresonant 6/77 bound and LRC(14) remain open.
source: root/adaptive_clock / LRC14 continuation session, 2026-09-03
depends_on:
  - THM-2060-crt-tail-coset-saturation
  - THM-2064-multitail-sheet-capacity-and-dyadic-seam
  - LRCUpTo13
related:
  - THM-4004-lrc14-three-detuned-divisor-comb-profile
  - THM-4150-safe-set-haar-measure-universal-odd-tail-lrc14-transfer
script: 04-computation/lrc14_scale_three_triple_owner_probe_root_20260903.py
output: 05-knowledge/results/lrc14_scale_three_triple_owner_probe_root_20260903.out
script_sha256: 72cd5142eb23570f94ccb909ec4270be05bb2c7c4d48c76adc5d470ff9a07e22
output_sha256: 637a7bafc4baa457d4c8cec959e6c97e06e5066096c62c9ae6ec767873c68226
hash_basis: raw LF bytes
---

# THM-4373 -- the signed-(1,2,1) scale-three sector has maximum `6/77`

**PROVED + VERIFIED-EXACT.**  For a positive integer `w` prime to three and
`j in Z/3Z`, put

```text
D_(w,j)={x in R/Z: ||w(x+j/3)||<1/14}.                 (1)
```

For distinct positive integers `a,b,c`, all odd and prime to three, define
the scale-three failure comb

```text
F_(a,b,c)=disjoint_union_(pi in S_3)
  D_(a,pi(0)) intersect D_(b,pi(1)) intersect D_(c,pi(2)), (2)
```

up to the measure-zero endpoints.  Equivalently, `F_(a,b,c)` is the set of
physical phases at which the three tails strictly spoil the three different
scale-three lifts.

Assume

```text
gcd(a,b,c)=1                                             (3)
```

and, after a permutation and choices of signs,

```text
+/-a +/-2b +/-c=0.                                      (4)
```

Then

```text
mu(F_(a,b,c))<=6/77.                                    (5)
```

Equality holds if and only if

```text
{a,b,c}={1,5,11}.                                       (6)
```

For a nonprimitive signed-`(1,2,1)` triple, common dilation preserves Haar
measure: multiplication by the common factor is Haar-preserving on the
circle, and because that factor is prime to three it merely permutes the
three sheet labels.  Therefore equality occurs exactly when its primitive
reduction is `(1,5,11)`, subject to the standing odd and prime-to-three
hypotheses.

This proves the sharp resonant sector of the triple-comb conjecture.  It does
not prove `(5)` for triples without a signed-`(1,2,1)` relation.

**Body-specific corollary.**  Let `C` be a finite positive core whose safe set

```text
G_C={y in R/Z: ||cy||>=1/14 for every c in C}
```

is nonempty.  If

```text
mu(G_C)>=6/77,                                           (6a)
```

then `3C union {a,b,c}` is `1/14`-lonely for every distinct odd unit triple
in the signed-`(1,2,1)` sector.  Indeed, multiplication by three identifies
the physical safe set for `3C` with a Haar preimage of `G_C`.  The strict
failure comb is a proper open subset of the circle (`x=0` is not in it),
while this body-safe
preimage is nonempty and compact.  Strict inequality in `(6a)` excludes
containment by measure; at equality, a nonempty compact subset of a proper
open subset of the connected circle cannot have the same measure as that
open set.  Hence some body-safe phase lies outside the failure comb, and one
of its three lifts is safe for all three tails.  For a ten-speed core `C`,
the stated nonemptiness is automatic from cited `LRCUpTo13`.

## 1. Endpoint parametrization of every resonant triple

Designate the speed carrying coefficient two in `(4)` by `b`, and sort the
two coefficient-one speeds as `p<q`.  Positivity leaves exactly two forms:

```text
2b=p+q                         (mean branch),
2b=q-p                         (difference branch).     (7)
```

Because `p,q,b` are odd units modulo three, these branches have exact
congruence descriptions:

```text
mean branch       iff q=p  mod 12,       b=(p+q)/2,
difference branch iff q=-p mod 12,       b=(q-p)/2.     (8)
```

Indeed, `b` is odd in the mean branch precisely when `p,q` agree modulo four,
and it is prime to three precisely when they agree modulo three.  In the
difference branch the corresponding conditions are disagreement modulo four
and disagreement modulo three.  The Chinese remainder theorem gives `(8)`.
Conversely, `(8)` makes `b` an odd positive unit and gives `(7)`.  The two
orientations cannot coexist under the unit hypotheses; solving two distinct
relations forces one of the three speeds to be divisible by three.

Since `p,q` are odd,

```text
gcd(p,b,q)=gcd(p,q).                                    (9)
```

Thus primitive resonant triples are parametrized exactly by

```text
p<q,  p and q odd,  3 does not divide pq,
gcd(p,q)=1,  and q=+/-p mod 12,                         (10)
```

with `b` recovered from `(8)`.

## 2. The middle tail is forced by the two endpoint errors

Pass from the physical phase `x` to the body phase

```text
y=3x mod 1.                                             (11)
```

The failure property is invariant under `x -> x+1/3`, so multiplication by
three identifies its quotient with the corresponding set of `y` and
preserves normalized Haar measure.

Put

```text
r=3/14.                                                 (12)
```

When an endpoint tail `w in {p,q}` is eligible, let `n_w` be the unique
nearest integer to `wy`, and write

```text
e_w=wy-n_w,                         |e_w|<r.             (13)
```

In the mean branch, the speed relation is `2b=p+q`; in the difference branch
it is `q=p+2b`.  If all three tails are eligible, their errors obey,
respectively,

```text
2e_b=e_p+e_q,                      2n_b=n_p+n_q,
2e_b=e_q-e_p,                      2n_b=n_q-n_p.         (14)
```

The integer identities are forced rather than assumed: before imposing them,
the absolute value of the error side is strictly below

```text
4r=6/7<1,
```

so the corresponding integer difference must vanish.  Conversely, if the
two endpoint inequalities `(13)` hold and the displayed half-sum or
half-difference of `n_p,n_q` is integral, then `(14)` defines the middle
nearest integer and

```text
|e_b|<=max(|e_p|,|e_q|)<r.                              (15)
```

Thus endpoint eligibility plus one parity condition is equivalent to
eligibility of all three tails.

The parity condition has a useful determinant form.  Put

```text
K=q n_p-p n_q.                                          (16)
```

Since `p,q` are odd, the half-sum or half-difference in `(14)` is integral if
and only if `K` is even.  Write

```text
K=2k,                         k in Z.                   (17)
```

A direct reduction of the exact owner colours

```text
o_w(y)=-w^(-1)n_w mod 3                                (18)
```

using `(8)` and `(14)` gives

```text
{o_p(y),o_b(y),o_q(y)}=Z/3Z
       iff                         3 does not divide k. (19)
```

For example, in the mean branch all three speeds have the same nonzero
residue modulo three; the middle nearest integer is the average of the
endpoint integers, so the three colours are distinct exactly when the two
endpoint colours differ.  This is `K!=0 mod 3`, equivalently `k!=0 mod 3`.
The difference branch gives the same criterion after replacing the average
by the half-difference.

## 3. One determinant, one component, and its exact length

By `(10)`, `gcd(p,q)=1`.  For every integer `k`, the equation

```text
q n_p-p n_q=2k                                         (20)
```

has integer solutions, and all its solutions form one orbit under

```text
(n_p,n_q) -> (n_p+p,n_q+q).                            (21)
```

This orbit is exactly translation of `y` by one.  Hence each `k` contributes
at most one component on the `y`-circle.

The endpoint danger intervals centered at `n_p/p` and `n_q/q` have radii
`r/p` and `r/q`.  Their center separation is

```text
|n_p/p-n_q/q|=2|k|/(pq).                               (22)
```

Since `p<q`, the exact intersection length is

```text
lambda_k=max(0, min(2r/q,
              r/p+r/q-2|k|/(pq))).                     (23)
```

Define

```text
A=3(q-p)/28,                      B=3(q+p)/28,
f(t)=2/(pq) ((B-t)_+-(A-t)_+),       t>=0.              (24)
```

Then `(23)` is exactly

```text
lambda_k=f(|k|).                                        (25)
```

The class `k=0` fails the owner condition `(19)`.  Positive and negative
determinants have equal lengths, so `(19)--(25)` prove the exact master
formula

```text
mu(F_(p,b,q))
 =2 sum_(k>=1, 3 does not divide k) f(k).               (26)
```

This formula depends only on the coefficient-one endpoints `p,q`; the mean
and difference branches have the same measure when they share endpoints.

## 4. Period-three quadrature

For `t>=0`, define

```text
E(t)=sum_(k>=1, 3 does not divide k) (t-k)_+ - t^2/3.   (27)
```

Writing `rho=t mod 3` with `0<=rho<3`, exact summation of the two retained
residue classes gives

```text
E(t)= -rho^2/3                         for 0<=rho<=1,
      rho-1-rho^2/3                    for 1<=rho<=2,
      2rho-3-rho^2/3                   for 2<=rho<3.     (28)
```

The endpoint values agree.  In particular

```text
E(t+3)=E(t),                         -1/3<=E(t)<=0.      (29)
```

Also, from `(24)`,

```text
integral_0^infinity f(t) dt=9/196.                      (30)
```

Substituting `(27)` into `(26)` yields the exact quadrature identity

```text
mu(F_(p,b,q))
 =3/49 + 4/(pq) (E(B)-E(A)).                            (31)
```

Therefore

```text
mu(F_(p,b,q))<=3/49+4/(3pq).                            (32)
```

If `pq>=81`, then

```text
3/49+4/(3pq)<6/77,                                     (33)
```

because `6/77-3/49=9/539` and `27*81>4*539`.

## 5. The thirteen small endpoint pairs

It remains to consider `pq<81` under `(10)`.  Necessarily `p<9`, so direct
use of the congruences `q=+/-p mod 12` leaves exactly

```text
(1,11), (1,13), (1,23), (1,25), (1,35), (1,37),
(1,47), (1,49), (1,59), (1,61), (1,71), (1,73),
(5,7).                                                  (34)
```

Formula `(26)` gives:

| `(p,q)` | `b` | `mu(F)` |
|---:|---:|---:|
| `(1,11)` | `5` | `6/77` |
| `(1,13)` | `7` | `6/91` |
| `(1,23)` | `11` | `12/161` |
| `(1,25)` | `13` | `12/175` |
| `(1,35)` | `17` | `12/245` |
| `(1,37)` | `19` | `2/37` |
| `(1,47)` | `23` | `22/329` |
| `(1,49)` | `25` | `24/343` |
| `(1,59)` | `29` | `24/413` |
| `(1,61)` | `31` | `24/427` |
| `(1,71)` | `35` | `30/497` |
| `(1,73)` | `37` | `30/511` |
| `(5,7)` | `1` | `8/245` |

The unique maximum in this table is `6/77` at `(p,b,q)=(1,5,11)`.
Together with `(33)`, this proves `(5)--(6)`.

## 6. Exact audit and scope

The exact companion retains two independent implementations of `(2)`: direct
rational interval intersection and a complete wall-cell calculation on one
`1/3`-period.  It also:

- verifies `(26)` and `(31)` against all `1,023` primitive resonant triples
  through speed `199`;
- checks the bound and unique equality case on all of them;
- enumerates the thirteen cases `(34)` and their exact measures;
- verifies permutation and common-dilation invariance on the independent
  audit pool; and
- checks all `47,499` primitive odd unit triples through `199`, whose unique
  finite maximum remains `(1,5,11)` and whose nonresonant maximum remains
  `(1,11,43)` with measure `12/301`.

The finite nonresonant census is evidence only.  The universal inequality
`mu(F_(a,b,c))<=6/77` for nonresonant triples remains **OPEN**, so this theorem
does not give an all-tail scale-three transfer or prove LRC(14).
**QED.**
