---
id: THM-4105
title: "Primitive reciprocal phase descent and quantitative arrival"
status: >
  PROVED + VERIFIED-EXACT + INDEPENDENTLY AUDITED. For a primitive integer row, the reciprocal
  exponent equations on an abelian group are equivalent to membership in its
  one-parameter power orbit. On the circle, the maximum reciprocal
  commutator defect controls distance to the physical LRC orbit within an
  explicit Bezout tariff, and conversely. An exact exponent factorization
  retains the signed lift determinant of a binding pair-sum straddle; its
  physical owner sign is intrinsic only below the antipodal boundary. These
  are physical-entry/near-arrival certificates, not LRC(14).
source: codex-lrc14-abc-exponent-reciprocity-20260825
depends_on: []
related:
  - THM-668-pair-sum-ruler-witness-structure
  - THM-4009-euclidean-covering-transference-short-relation-compression
  - THM-4095-exact-arithmetic-field-transport-gapped-pair-margins-and-order-tournament-blindness
  - THM-4106-lrc-pair-owner-reciprocity-and-valuation-tree-decoder
  - THM-4108-abc-conditional-reciprocal-power-radicals-and-lrc-gauge-obstruction
script: 04-computation/reciprocal_phase_arrival_straddle_thm4105.py
output: 05-knowledge/results/reciprocal_phase_arrival_straddle_thm4105.out
script_sha256: 2d05399a87db8b8ff6006622f225bc87e27d542ced4016f96d3b37a69054fd1c
output_sha256: dfaea1aff7b61e2d5389e60196c198c98a5f23ecfb38c2edad0fe559c8517825
hash_basis: raw LF bytes
audit: >
  PASS after adding the antipodal owner-sign firewall and distinguishing exact
  physicality from nonzero-defect near-arrival. The independent audit checked
  every quantifier, both constants, the straddle algebra, hostiles, hashes,
  and normal/-O/frozen transcript identity.
---

# THM-4105 -- primitive reciprocal phase descent and quantitative arrival

**PROVED + VERIFIED-EXACT + INDEPENDENTLY AUDITED.** Reciprocal exponentiation becomes an exact LRC
coordinate only after the objects and gauge are declared. The raw comparison
`a^b ? b^a` is a scalar order; the equations
`z_i^(v_j)=z_j^(v_i)` instead test whether a whole phase packet has one
physical time.

## 1. Exact reciprocal-exponent descent

Let `G` be an abelian group, written multiplicatively, and let

```text
v=(v_1,...,v_r) in Z^r,              gcd(v_1,...,v_r)=1.       (1)
```

For `z=(z_1,...,z_r) in G^r`, the following are equivalent:

```text
(i)  z_i^(v_j)=z_j^(v_i) for every i,j;
(ii) there is a unique x in G such that z_i=x^(v_i) for every i. (2)
```

Choose a Bezout word `c in Z^r` with

```text
sum_i c_i v_i=1
```

and put `x=product_i z_i^(c_i)`. Condition `(i)` gives, for every `j`,

```text
x^(v_j)=product_i z_i^(c_i v_j)
         =product_i z_j^(c_i v_i)=z_j.                       (3)
```

If `x` and `y` both realize the packet, then

```text
x/y=product_i (x/y)^(c_i v_i)=1,                            (4)
```

so the realization is unique. The converse is immediate.

For `G=R/Z`, written additively, `(2)` says exactly

```text
v_j theta_i-v_i theta_j=0 mod 1 for every i,j

iff

theta_i=v_i t mod 1 for one unique t in R/Z.                 (5)
```

Thus the reciprocal commutators characterize the physical one-parameter LRC
orbit of a primitive speed row. They do not assert that the realized packet
is lonely.

### Primitive boundary

Primitivity is load-bearing. For speeds `(2,4)` and phases `(0,1/2)`,

```text
4*0-2*(1/2)=0 mod 1,                                       (6)
```

but no `t` satisfies `(2t,4t)=(0,1/2) mod 1`. Dividing the row by its gcd
changes the commutator to the primitive `(1,2)` equation and exposes the
defect. Raw cross-power compatibility on a nonprimitive row is therefore not
a physical-entry certificate.

## 2. Quantitative arrival tariff

Let `v=(v_1,...,v_r)` be a primitive row of positive integers. For
`theta in (R/Z)^r`, define

```text
A_v(theta)=inf_(t in R/Z) max_i ||theta_i-v_i t||,
D_v(theta)=max_(i,j) ||v_j theta_i-v_i theta_j||,
B(v)=min {sum_i |c_i| : c in Z^r, sum_i c_i v_i=1},
V=max_i v_i.                                               (7)
```

Then

```text
boxed: D_v(theta)/(2V) <= A_v(theta) <= B(v)D_v(theta).    (8)
```

For the upper bound, take a minimizing Bezout word and set
`t=sum_i c_i theta_i`. In `R/Z`,

```text
theta_j-v_j t
 =sum_i c_i(v_i theta_j-v_j theta_i),                     (9)
```

so the circle triangle inequality gives the right side of `(8)`. Conversely,
if `e_i=theta_i-v_i t`, then

```text
||v_j theta_i-v_i theta_j||
 <=v_j||e_i||+v_i||e_j||
 <=2V max_k ||e_k||.                                      (10)
```

Maximizing over pairs and infimizing over `t` gives the left side. At defect
zero, `(8)` recovers the exact descent theorem.

The tariff `B(v)` is the price of choosing one global time from pairwise
defects. It may be large. A small commutator defect is therefore a controlled
arrival statement, not by itself a uniform LRC margin.

## 3. Exact exponentiation lift of a binding straddle

Let `u,w` be distinct positive speeds with `uw>1`, let `alpha,beta` be
integer lifts, and put

```text
Q=u+w,
P=alpha+beta,
D=u beta-w alpha>0,
t=P/Q,
M=D/Q.                                                     (11)
```

Then

```text
u t-alpha=M,                  w t-beta=-M.                (12)
```

When `0<D<Q/2`, this is strict local pair-sum straddle geometry. At
`D=Q/2`, the two errors are antipodal and their physical owner sign is not
canonical, although the chosen integer lifts and the algebra below remain
valid. In either case it has the exact reciprocal-power factorization

```text
boxed:
(u^beta/w^alpha)^Q=(u^w/w^u)^P (uw)^D.                   (13)
```

Indeed,

```text
beta Q=wP+D,                  alpha Q=uP-D,               (14)
```

which is precisely the exponent ledger on the two sides of `(13)`. In
logarithmic form,

```text
M=[log(u^beta/w^alpha)-t log(u^w/w^u)]/log(uw).           (15)
```

Reversing the oriented owner pair--that is, swapping the labelled pairs
`(u,alpha)` and `(w,beta)`--inverts the reciprocal ratios and negates `D`;
the ordinary LRC margin retains only `|D|/Q`. Hence `(13)` is a lossless
multiplicative lift of the local owner sign for `0<D<Q/2` once a strict
binding straddle has been located. At `D=Q/2` it records only the selected
lift orientation: for `(u,w,t)=(1,3,1/2)`, lifts `(0,2)` and `(1,1)` give
opposite signed determinants for the same physical antipodal packet.

For AP13, the canonical active straddle

```text
(u,w,alpha,beta)=(1,13,0,1)
```

has `P=D=1` and `t=M=1/14`.

## 4. The scalar comparator is a hostile control

For positive reals,

```text
C(a,b)=a^b/b^a,
phi(a)=log(a)/a,
sign log C(a,b)=sign(phi(a)-phi(b)).                     (16)
```

Thus raw `a^b ? b^a` is a scalar-potential preorder. Its edge values satisfy
the exact weighted zero-curl law

```text
C(a,b)^c C(b,c)^a C(c,a)^b=1.                           (17)
```

Any cyclic exponent tournament must therefore introduce a contextual
normalization or another coordinate and must ledger the information changed
by that operation. THM-4107 identifies one such pairwise gcd gauge.

## 5. Transfer boundary and generated obligations

The theorem changes the physical-entry layer of the LRC proof architecture:

```text
phase packet
   -- reciprocal commutators --> exact orbit membership / arrival defect
   -- binding straddle --------> signed determinant factorization
   -- loneliness test ---------> still requires all runner clearances.      (18)
```

In particular:

- exact physicality requires vanishing reciprocal commutators; nonzero defect
  gives only the quantified distance-to-orbit bound `(8)`;
- the signed factorization `(13)` does not locate the binding pair and does
  not synchronize the other eleven speeds;
- radical support of `(13)` is not an arrival decoder; THM-4108 gives the
  precise ABC boundary;
- a useful next computation is the minimum-`l1` Bezout tariff on the 165
  valuation orbits and 17 live `11+2` types, followed by an intersection with
  the actual owner cells rather than their scalar summaries.

No statement here proves LRC(14), ABC, or an IUT-to-ABC implication.

## 6. Exact audit

Reproduce with

```text
python -B 04-computation/reciprocal_phase_arrival_straddle_thm4105.py
python -B -O 04-computation/reciprocal_phase_arrival_straddle_thm4105.py
```

The normal and optimized streams equal the frozen transcript. The audit uses
only integer and `Fraction` arithmetic. It checks `354,466` finite-group
states, all `3,010` compatible primitive phase packets, `361,622`
quantitative pointwise gates, the nonprimitive hostile, `3,049` straddle
factorizations, and the AP13 owner-sign control.
