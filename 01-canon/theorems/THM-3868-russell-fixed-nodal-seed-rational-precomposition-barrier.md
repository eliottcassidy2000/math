---
id: THM-3868
title: "Russell fixed nodal seed rational precomposition barrier"
status: >
  PROVED + VERIFIED-EXACT + INDEPENDENTLY HOSTILE-AUDITED.  For the fixed
  minimal nodal seed of THM-3860, no rational controls Z,S in Frac(B) can make
  the composite both regular on the Russell surface and Keller.  A monic cubic
  recovers S integrally, then Z polynomially; closure of the Russell Poisson
  bracket and B*=k* give the contradiction.  This covers genuinely mixed
  Z_s!=0 controls, but not a different seed or nonrational controls.
source: root / mixed_normal_pole_scout / planar JC mixed-normal lane, 2026-08-23
audit: >
  INDEPENDENT HOSTILE PROOF/REPLAY AUDIT PASSED.  The independent companion's
  101 checks rederive the monic cubic and its generic irreducibility,
  discriminant/cusp and derivative/different, normal integral recovery,
  Poisson closure and the constant-unit trap.  An independent convolution
  verifies formal rows through order eight, and a Laurent-UFD audit proves the
  genuinely mixed hostile has a height-one pole avoiding r,z,D,G,s.  Normal
  and optimized executions of both companions byte-match frozen LF outputs.
depends_on:
  - THM-3785-linear-higher-pole-russell-pseudoplane-maximal-observable
  - THM-3846-formal-arm-darboux-lift-and-algebraization-gate
  - THM-3860-russell-higher-normal-rational-lifts-and-vertical-pole-barrier
related:
  - THM-3801-cubic-etale-normalization-nonmonogenic-and-companion-sheet-gate
  - THM-3843-russell-arm-birational-immersion-and-forced-self-identification
  - THM-3849-russell-arm-conductor-polynomial-and-residual-contact-graph
  - THM-3862-russell-finite-completion-nonmonogenic-branch-contract
script: 04-computation/jc2_russell_fixed_nodal_seed_rational_precomposition_barrier_thm3868.py
output: 05-knowledge/results/jc2_russell_fixed_nodal_seed_rational_precomposition_barrier_thm3868.out
script_sha256: be15ae2ac76a118bb009dadf2ce898b9d780ec5572e2b94156c9ecb736f337b1
output_sha256: caa9ef246b0e4be6a88ee07b56eee1a76330f12c0bb60dbf0da7b665986137d2
semantic_sha256: 53ae645f196938bccdb3400f6f54082e8802e53a9e2e32e1f73040eb5b2feb37
independent_script: 04-computation/jc2_russell_fixed_nodal_seed_rational_precomposition_barrier_independent_audit_thm3868.py
independent_output: 05-knowledge/results/jc2_russell_fixed_nodal_seed_rational_precomposition_barrier_independent_audit_thm3868.out
independent_script_sha256: 144ac480ae96038bae4fe909ddf93d76e5baf7d3b6022ff18cda61016952beb7
independent_output_sha256: 7ea79598e49a557516cb415914e0088a48794ed6c994a1cc4879aaf2ab9d8fd8
independent_semantic_sha256: 64a1358ae53b5f8b00f234161db9a3812eb8055f40716aa449f3dafa07009821
hash_basis: raw LF bytes
---

# THM-3868 -- the fixed nodal seed admits no rational regular Keller precomposition

**PROVED + VERIFIED-EXACT + INDEPENDENTLY HOSTILE-AUDITED.**  Let `k` be
algebraically closed of characteristic zero, fix `c in k*`, and retain the
smooth Russell surface and its completed arm

```text
B=k[r,z,e]/(r^2e-z^3+c^3r),                 Y=Spec(B),          (1)
Bhat=k[s][[z]],                             {z,s}=1,
s=e/[3(c^3+er)],
e=3c^3s+9s^2z^3,                  r=z^3/(c^3+3sz^3).          (2)
```

Fix the minimal nodal seed in independent controls `(Z,S)`:

```text
A_0=9c^6S^2-Z/(3c^3),
C_0=27c^9S^3-3c^3S-(3/2)SZ.                                (3)
```

There are **no** `Z,S in Frac(B)` such that

```text
A=A_0(Z,S) in B,             C=C_0(Z,S) in B,               (4)
{A,C}=1.                                                          (5)
```

Thus every rational precomposition of this fixed seed is excluded,
including genuinely mixed controls with `Z_s!=0`.  This is stronger than
THM-3860's vertical-subclass barrier.  It does not exclude a different
higher-normal seed, controls outside `Frac(B)`, or arbitrary pairs in `B^2`
not factoring through `(3)`.  No JC(2) conclusion is claimed.

## 1. Monic recovery and the local pole mechanism

The first equation of `(3)` gives

```text
Z=27c^9S^2-3c^3A.                                            (6)
```

Substituting into the second gives

```text
C=-(27/2)c^9S^3+((9/2)c^3A-3c^3)S,                          (7)
```

or the monic recovery cubic

```text
S^3+[(2-3A)/(9c^6)]S+2C/(27c^9)=0.                         (8)
```

If `(4)` holds, `(8)` makes the rational function `S` integral over `B`.
THM-3785 proves that `B` is smooth, hence normal, so `S in B`; equation
`(6)` then puts `Z in B`.  Therefore

```text
A_0(Z,S),C_0(Z,S) in B and Z,S in Frac(B)  ==>  Z,S in B.   (9)
```

This packages the hostile valuation calculation.  If at a DVR
`ord(S)=n<0`, regularity of `A` could only occur with `ord(Z)=2n` and
leading relation `Z=27c^9S^2`.  The two lowest terms of `C` then have
coefficient

```text
27c^9-(3/2)(27c^9)=-(27/2)c^9 != 0,                         (10)
```

so `C` still has order `3n<0`.  Poles in the two controls cannot cancel in
both nodal coordinates.

## 2. The cubic derivative is the seed different

For independent `(A,C)`, put

```text
f(T)=T^3+pT+q,
p=(2-3A)/(9c^6),                  q=2C/(27c^9).              (11)
```

Equation `(8)` is `f(S)=0`, and

```text
Disc(f)=-4[(2-3A)^3+27C^2]/(729c^18).                       (12)
```

The generic cubic is irreducible.  Indeed `k(A,C)=k(p,q)`, and at the
`q`-infinity valuation a hypothetical rational root has integral order `m`.
The orders of its three terms are `3m,m,-1`; for every integer `m` their
minimum occurs uniquely, which is impossible in a zero sum.  Thus the
finite degree-three seed cover is branched over the cusp

```text
(3A-2)^3=27C^2.                                              (13)
```

Together `(6),(8)` and generic irreducibility identify

```text
k[S,Z]=k[A,C,S]=k[A,C][T]/(f).                             (D1)
```

Its power-basis derivative, hence its monogenic different generator, is
exactly the seed Jacobian factor:

```text
f'(S)=3S^2+(2-3A)/(9c^6)
     =(Z+2c^3)/(9c^9)
     =[2/(9c^6)](1+Z/(2c^3)).                               (14)
```

This is a literal bridge to THM-3862's nonmonogenic finite-completion
contract: the density factor here is the cubic derivative/different, not
merely an analogous discriminant statistic.

## 3. Poisson closure and the unit contradiction

Direct differentiation in `(3)` gives

```text
J_(Z,S)(A_0,C_0)=1+Z/(2c^3).                                (15)
```

The bracket of `(2)` restricts to `B`; on its generators,

```text
{z,r}=-3r^2,                 {z,e}=3c^3+6er,
{r,e}=9z^2.                                                   (16)
```

Hence `{B,B} subset B`.  Under `(4)`, recovery `(9)` puts `Z,S` in `B`,
so the chain rule and `(5),(15)` give an identity in `B`:

```text
1=(1+Z/(2c^3)){Z,S}.                                        (17)
```

The first factor is therefore a unit.  THM-3785 proves `B*=k*`, so `Z` is
constant.  Then `{Z,S}=0`, contradicting `(17)`.  This proves `(4),(5)`.

Equivalently, with `w=1/(2c^3)` and

```text
T=(1+wZ)S,                                                   (18)
```

one has `{Z,T}=(1+wZ){Z,S}`.  The density equation alone is an ordinary
rational Keller equation; the hard sidecar is divisibility of `T` by
`1+wZ`.  Integral recovery restores that sidecar globally, where constant
units make it impossible.

## 4. Formal mixed equations through order six

The contradiction is global, not formal.  Write

```text
Z=z-(w/2)z^2+uz^3+vz^4+pz^5+qz^6+tz^7+...,
S=s+hz^2+kz^3+lz^4+mz^5+nz^6+...,                           (19)
```

where every displayed coefficient is a function of `s` and primes below
mean `d/ds`.  The equation

```text
(1+wZ)J_(z,s)(Z,S)=1                                        (20)
```

has the following exact rows:

```text
E_2=3u+h'-(3/2)w^2,                                         (21)

E_3=4v+k'+4wu+(1/2)w^3,                                    (22)

E_4=5p+l'+5wv-(5/2)w^2u
    +(3u-(3/2)w^2)h'-2hu',                                 (23)

E_5=6q+m'-3w^2v+6wp+3wu^2
    +(3u-(3/2)w^2)k'-(2wh+3k)u'
    +((1/2)w^3+4wu+4v)h'-2hv',                             (24)

E_6=7t+n'-(7/2)w^2p+7wq+7wuv
    +(3u-(3/2)w^2)l'-(2wh+3k)v'
    +((1/2)w^3+4wu+4v)k'
    +(w^2h-3wk-4l)u'
    +(-(5/2)w^2u+5wv+5p)h'-2hp'.                           (25)
```

Each row is set equal to zero.  The first says

```text
u=w^2/2-h'/3,                                                (26)
```

so arbitrary polynomial `h` gives polynomial `u`.  More generally, if

```text
Z=sum_(i>=1)a_i(s)z^i, a_1=1,
S=s+sum_(j>=2)b_j(s)z^j,                                   (27)
```

the only new terms at order `z^N` are

```text
(N+1)a_(N+1)+b_N'.                                         (28)
```

All others come from earlier rows.  One may choose `b_N` freely and solve
for `a_(N+1)` in characteristic zero.  Polynomial coefficients stay
polynomial.  Thus no finite formal-density obstruction exists, and the jet
`(u,h)` alone contains no global denominator data.

## 5. Exact genuinely mixed rational hostile

Let `eta in k*` and define

```text
D=1-eta s z^2,                x=z/D^2,             y=sD^3.  (29)
```

Directly, `J_(z,s)(x,y)=1`.  Compose this rational Hamiltonian control with
THM-3860's Mobius density control:

```text
H=D^2+(w/2)z,                 G=D^2+(3w/2)z,
Z=z/H,                        S=sH^3/(DG).                    (30)
```

Equivalently `Z=x/(1+wx/2)` and
`S=y(1+wx/2)^3/(1+3wx/2)`, so

```text
(1+wZ)J_(z,s)(Z,S)=1                                           (31)
```

exactly.  The first jets are

```text
Z=z-(w/2)z^2+(w^2/4+2eta s)z^3
    -(w^3/8+2eta ws)z^4
    +(w^4/16+(3/2)eta w^2s+3eta^2s^2)z^5+...,

S=s+(3w^2s/4-3eta s^2)z^2-w^3s z^3
    +(3w^4s/2+(3/4)eta w^2s^2+3eta^2s^3)z^4+... .          (32)
```

Thus `Z_s` first occurs at order three and `(21)` holds with

```text
u=w^2/4+2eta s,                 h=3w^2s/4-3eta s^2.          (33)
```

For a named pole specialize `w=1/(2c^3)` and use

```text
Y_(rz)=Spec B_(rz)=Spec k[z^(+-1),r^(+-1)],
s=1/(3r)-c^3/(3z^3).                                       (34)
```

Clearing Laurent units from `H` gives

```text
N_eta=4c^3(3rz+eta c^3r-eta z^3)^2+9r^2z^3.                (35)
```

This is a nonzero nonunit Laurent polynomial.  (Its `r^0z^6` term is
nonzero for `eta!=0`, while it also has terms involving `r`.)  Hence the
Laurent UFD has a height-one factor of `H`.  At such a factor,
`H`, the functions `D,G,s,z` are generically units because

```text
H-D^2=(w/2)z,                 G-H=wz,                        (36)
```

and `H|_(s=0)=1+wz/2`.  Here `s` is a Laurent unit times the irreducible
element `z^3-c^3r`, so the
last specialization really proves that no prime factor of `H` is an `s`
factor.  If the chosen factor has multiplicity `mu>0` in `H`, then

```text
ord(Z)=-mu,                   ord(S)=3mu,
ord(A_0(Z,S))=-mu.                                            (37)
```

At `eta=0`, `(35)` reduces to

```text
N_0=9r^2z^2(z+4c^3),                                        (38)
```

the vertical divisor of THM-3860.  The genuinely mixed coordinate moves
the divisor; it does not cancel it.

## 6. Boundary and next test

The candidate mechanism is

```text
regular nodal composite
 -> monic recovery in normal B
 -> controls in B
 -> seed different is a unit
 -> constant-unit contradiction.                             (39)
```

It closes this fixed rational-precomposition strategy.  The cheapest live
test is therefore an alternative higher-normal seed with the same arm and
first-normal packet: eliminate its hidden control and test whether recovery
is still monic and whether its Jacobian factor is still a power-basis
derivative.  The first failure of integral recovery is the first genuine
escape from `(39)`.  A nonrational control requires a separate normalization
and valuation audit.  JC(2) remains open.

## 7. Reproduction

```bash
python3 -B 04-computation/jc2_russell_fixed_nodal_seed_rational_precomposition_barrier_thm3868.py
python3 -B -O 04-computation/jc2_russell_fixed_nodal_seed_rational_precomposition_barrier_thm3868.py
python3 -B 04-computation/jc2_russell_fixed_nodal_seed_rational_precomposition_barrier_independent_audit_thm3868.py
python3 -B -O 04-computation/jc2_russell_fixed_nodal_seed_rational_precomposition_barrier_independent_audit_thm3868.py
```

Each primary run matches
`05-knowledge/results/jc2_russell_fixed_nodal_seed_rational_precomposition_barrier_thm3868.out`;
each independent run matches
`05-knowledge/results/jc2_russell_fixed_nodal_seed_rational_precomposition_barrier_independent_audit_thm3868.out`.
The scripts perform 50 and 101 exact checks without random or floating-point
input.

**QED.**
