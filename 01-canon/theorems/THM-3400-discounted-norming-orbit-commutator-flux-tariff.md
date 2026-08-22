---
id: THM-3400
title: "Discounted norming-orbit commutator-flux tariff"
status: >
  PROVED + VERIFIED-EXACT + INDEPENDENTLY AUDITED.  Let
  c>=2 and let uniformly bounded all-iterate defects
  E_n=cV^*Q^{*n}V-T^{*n} come from an isometry V and a contraction Q.  If
  kappa=||T||>c, every norming vector pays discounted commutator leakage at
  least kappa(kappa-c).  Consequently leakage at most epsilon gives the
  optimal robust bound (c+sqrt(c^2+4epsilon))/2.  Global commutation is only
  the zero-leakage special case.  The tariff and robust bound are attained
  for every epsilon by the exact family T=kappa J_2,Q=J_2,V=I.  A
  subcritical noncommuting family has strictly negative leakage, while a
  power-bounded nonnormal involution has leakage kappa^2+1 and a scalar
  hostile shows why an infinite uniform tail remains necessary.  This
  strictly strengthens THM-3390 but supplies no Crouzeix completion.
source: root-2608-crouzeix-puzzle-2026-08-15
audit: algebraic recurrence, discounted-tail, tariff, equality-family, noncommuting, periodic, and finite-prefix exact controls; independent immutable-file proof, replay, hash, semantic-digest, AST, routing, dependency, and hostile audit complete
depends_on:
  - THM-3390-all-iterate-commuting-completion-norm-bound
related:
  - THM-2222-scalar-transfer-parity-tower-and-four-checkpoint-survivor-reduction
  - THM-3315-tournament-cut-switching-centered-coronal-walk-compiler
script: 04-computation/discounted_norming_orbit_commutator_flux_tariff_thm3400.py
output: 05-knowledge/results/discounted_norming_orbit_commutator_flux_tariff_thm3400.out
script_sha256: 778784e967ee26d5973dec01f9ed87c048da99589183103aec3b655a53203a82
output_sha256: 14b74082e9dc653f15b5e7d08a9bbb3c8518d06ea97eab26cd800144e195a2db
semantic_sha256: 6cfccf40aefd4df5b0dc41c8a6c17a16ae39e36adb1581121cfae0349b79f6dc
hash_basis: LF-normalized bytes
---

# THM-3400 -- discounted norming-orbit commutator-flux tariff

**PROVED + VERIFIED-EXACT + INDEPENDENTLY AUDITED.**

## 1. Inheritance and the changed coordinate

[THM-3390](THM-3390-all-iterate-commuting-completion-norm-bound.md) proves
that a uniformly bounded all-iterate completion controls `||T||` when every
defect commutes with `T`.  Its proof never uses the full operator equation.
It sees only one discounted real commutator current along a norming orbit.

The closest hostile is still `Q=0,T=MJ_2`: the defects are bounded because
`J_2^2=0`, but noncommutation permits `M>c`.  The corrected near miss is
therefore not "commutation can be omitted."  It is that commutation can be
replaced by the exact scalar leakage it killed.  The least-used sidecar is
the norming vector itself.

| field | exact connection |
|---|---|
| source | the all-iterate positive completion of THM-3390 |
| target | one discounted real commutator current on a norming orbit |
| map | contract `[T^n,E_n]` and `[T^n,E_(n+1)]T` against one norming vector, then discount by `||T||^(-n)` |
| preserved | exactly the noncommuting term in the singular-vector telescope |
| destroyed | operator-level commutation, behavior off the norming orbit, and each individual time after summation |
| required sidecars | every iterate and a uniform defect bound |
| cheapest decisive tests | `T=kappa J_2,Q=J_2` for equality and `T=MJ_2,Q=0` for positive leakage |

## 2. The quantitative theorem

Let `H` be a finite-dimensional complex Hilbert space, let `T in L(H)`, and
let `c>=2`.  Suppose there are a Hilbert space `K`, an isometry

```text
V:H -> K,
```

and a contraction `Q in L(K)` such that

```text
E_n=cV^*Q^{*n}V-T^{*n},                 n>=1,          (1)
sup_(n>=1)||E_n||=B<infinity.                            (2)
```

Put `kappa=||T||`.  When `kappa>1`, choose any unit norming vector `x`, so

```text
T^*T x=kappa^2 x.                                      (3)
```

With `[A,B]=AB-BA`, define

```text
C_n(x)=Re <[T^n,E_n]x,x>,
D_n(x)=Re <[T^n,E_(n+1)]T x,x>,
L_n(x)=kappa C_n(x)-D_n(x),                            (4)

Lambda_x=sum_(n>=1) kappa^(-n)L_n(x).                 (5)
```

The series in `(5)` converges absolutely.  If `kappa>c`, then every norming
vector satisfies the sharp **leakage tariff**

```text
Lambda_x>=kappa(kappa-c).                              (6)
```

Consequently, for `epsilon>=0`, if there is a norming vector with

```text
Lambda_x<=epsilon,                                     (7)
```

then

```text
||T||<=[c+sqrt(c^2+4epsilon)]/2.                       (8)
```

When `||T||<=1`, `(8)` is automatic and no series is needed.  The proof and
conclusion remain valid on an arbitrary Hilbert space whenever `T` attains
its norm.

Two useful transparent corollaries are

```text
L_n(x)<=0 for every n  ==>  ||T||<=c,
E_nT=TE_n for every n  ==>  Lambda_x=0 and ||T||<=c.   (9)
```

Thus THM-3390 is the zero-leakage special case.  The first implication in
`(9)` is strictly more general: it asks for neither equality nor global
operator commutation.

## 3. Proof

Use the convention that the inner product is linear in its first variable.

### 3.1 Bounded defects type the infinite current

Put

```text
P_n=cV^*Q^{*n}V.                                      (10)
```

Since `Q` is a contraction and `V` is an isometry,

```text
||P_n||<=c,
||T^n||=||T^{*n}||=||P_n-E_n||<=c+B.                  (11)
```

Therefore

```text
|C_n(x)|<=2(c+B)B,
|D_n(x)|<=2kappa(c+B)B,
|L_n(x)|<=4kappa(c+B)B.                               (12)
```

For `kappa>1`, equation `(12)` proves absolute convergence of `(5)`.  This
is why a finite prefix or merely pointwise finite defects cannot replace
the uniform all-iterate hypothesis.

### 3.2 Exact commutator ledger

Define

```text
m_n=Re <E_nT^n x,x>,
y_n=(P_(n+1)T-kappa P_n)x,
a=kappa^2-kappa.                                      (13)
```

No commutation is assumed.  Moving `E_n` and `E_(n+1)` through `T^n` while
retaining their commutators gives

```text
r_n:=Re <T^n(kappa E_n-E_(n+1)T)x,x>
    =kappa m_n-m_(n+1)+L_n(x).                        (14)
```

On the other hand, `(3)`, `(10)`, and `(13)` give

```text
kappa E_n-E_(n+1)T
 =-P_(n+1)T+kappa P_n+(kappa^2-kappa)T^{*n}
```

after application to `x`.  Hence

```text
r_n=a||T^{*n}x||^2-Re <T^n y_n,x>
   =a||T^{*n}x-y_n/(2a)||^2-||y_n||^2/(4a)
   >=-||y_n||^2/(4a).                                 (15)
```

Combining `(14)` and `(15)`, dividing by `kappa`, and iterating through `N`
gives

```text
m_1>=kappa^(-N)m_(N+1)
     -sum_(n=1)^N kappa^(-n)||y_n||^2/(4a)
     -sum_(n=1)^N kappa^(-n)L_n(x).                   (16)
```

Equation `(11)` and the uniform defect bound make `m_(N+1)` bounded, so the
first term tends to zero.

### 3.3 The tariff

Set

```text
w=Q^*VT x-kappa Vx.                                   (17)
```

Then

```text
y_n=cV^*Q^{*n}w,                    ||y_n||<=c||w||.   (18)
```

Letting `N` tend to infinity in `(16)` and summing the geometric series
yields

```text
m_1>=-c^2||w||^2/[4kappa(kappa-1)^2]-Lambda_x.        (19)
```

Independently, contraction of `Q`, `||Tx||=kappa`, and
`P_1=E_1+T^*` give

```text
||w||^2
 <=2kappa^2-(2kappa/c)(m_1+kappa^2).                  (20)
```

Substituting `(19)` into `(20)` gives the exact quantitative remainder

```text
[1-c/(2(kappa-1)^2)]||w||^2
 <=(2kappa/c)[Lambda_x-kappa(kappa-c)].                (21)
```

If `kappa>c>=2`, the left coefficient is positive because

```text
(kappa-1)^2>(c-1)^2>=c/2.                             (22)
```

The left side of `(21)` is then nonnegative, forcing `(6)`.  Under `(7)`,

```text
kappa(kappa-c)<=epsilon,
```

whose positive quadratic root is exactly `(8)`.

## 4. Optimality for every leakage budget

Let

```text
J=[0 1]
  [0 0],             H=K=C^2,       V=I,
Q=J,                 T=kappa J,     kappa>=c.          (23)
```

Take `x=e_2`.  Since `J^2=0`,

```text
E_1=(c-kappa)J^*,                 E_n=0 for n>=2,
w=0,
C_1=kappa(kappa-c),              D_1=0,
L_1=kappa^2(kappa-c),
Lambda_x=kappa(kappa-c).                              (24)
```

Thus `(6)` is equality for every `kappa>=c`.  Given any `epsilon>=0`, choose

```text
kappa=[c+sqrt(c^2+4epsilon)]/2.                        (25)
```

Then `Lambda_x=epsilon` and equality holds in `(8)`.  Neither the tariff nor
the robust norm bound can be improved within the stated hypotheses.

## 5. Strict weakening and sharp failure boundaries

### 5.1 Favorable leakage need not be commutation

In `(23)`, instead take `1<s<c` and `T=sJ`.  Then

```text
E_1=(c-s)J^*,                   [T,E_1]!=0,
Lambda_x=s(s-c)<0.                                    (26)
```

The one-sided condition in `(9)` holds strictly even though the only
nonzero defect does not commute with `T`.  This proves that the new gate is
genuinely weaker, not a restatement of THM-3390.

### 5.2 A bounded periodic hostile pays visible leakage

Let `kappa>c`, take `Q=0,V=I`, and put

```text
T=[0       kappa]
  [1/kappa 0].                                         (27)
```

Then `T^2=I`, so `E_n=-T^{*n}` is uniformly bounded although `T` is
nonnormal and `||T||=kappa`.  For `x=e_2`, even leakage terms vanish and

```text
L_(2j+1)=kappa(kappa^2-kappa^(-2)).                   (28)
```

Therefore

```text
Lambda_x
 =sum_(j>=0)kappa^(-(2j+1))L_(2j+1)
 =kappa^2+1.                                          (29)
```

Uniform all-iterate defects alone do not imply the norm bound.  The missing
coordinate is exactly the positive observer leakage.

### 5.3 Omitting leakage or the infinite uniform tail

For `M>c`, `Q=0,V=I,T=MJ` gives bounded defects

```text
E_1=-MJ^*,                 E_n=0 for n>=2,
Lambda_x=M^2>0,            ||T||=M.                   (30)
```

This is the minimal-dimensional hostile to deleting the leakage hypothesis.

For every finite `N`, the scalar choice `Q=0,V=1,T=M>c` has zero
commutator leakage and a finite cap `M^N` on `E_1,...,E_N`, but
`E_n=-M^n` is unbounded on the infinite tail.  Thus no arbitrarily long
finite prefix with an unspecified finite cap can replace `(2)`.

## 6. Scope

The theorem is a self-contained operator inequality.  It does not construct
an all-iterate completion for a given operator, prove scalar or complete
Crouzeix, improve the sharp mass `c=2`, or control a non-norm-attaining
operator without an approximation argument.  Lorist--Schwenninger
[arXiv:2608.03841v1](https://arxiv.org/abs/2608.03841) is a very recent cited
preprint motivating THM-3390; no claim from that preprint is used here.

The creative gain is diagnostic: any putative norm excess must now carry a
quantified scalar incompatibility.  This is the operator analogue of a
closure tariff, not a transfer into LRC, FC, tournament, or JC objects.

## 7. Exact companion

Run

```bash
python3 04-computation/discounted_norming_orbit_commutator_flux_tariff_thm3400.py
python3 -O 04-computation/discounted_norming_orbit_commutator_flux_tariff_thm3400.py
```

The standard-library companion verifies the recurrence ledger, the tariff
and robust equality family on a rational grid, the strictly favorable
noncommuting family, the square-zero omission hostile, the exact infinite
periodic sum `(29)`, and the finite-prefix boundary.  It uses exact rational
arithmetic only.  The computation is a hostile referee; the theorem follows
from `(11)--(22)`.
