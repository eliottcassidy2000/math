---
id: THM-2570
title: "Deep-augmentation Parseval energy and nonlinear holotopy obstruction"
status: >
  PROVED CANDIDATE + VERIFIED-EXACT; INDEPENDENT HOSTILE AUDIT PENDING.
  A lawful diagonal-zero C_13^3 table of mass rho has deep target-zero
  anchor energy at least rho^2/(12*13^6), sharply.  The six displayed
  THM-2562 duty replicas form a rank-one Gram packet of energy at least
  d_g^2 rho^2/(2*13^6) for every gain.  One adaptively selected deep colour,
  common to all gains, has anchor at least rho/(12*13^3).  No prescribed
  colour has a positive uniform floor.  This is coefficient-level
  partial-bare energy, not a physical common-carrier observable; no row is
  excluded and LRC(14) remains open.
source: root-holotopy-2026-07-27-deep-energy
depends_on:
  - THM-2562-canonical-duty-commutator-line-rank-and-anchor-rigidity
  - THM-2567-deep-coloured-duty-replica-cycle-and-augmentation-cancellation
related:
  - THM-2383-polarized-complete-subcube-gram-tomography
  - THM-2563-paired-dipole-deep-target-corner-and-partial-bare-boundary
script: 04-computation/lrc14_deep_augmentation_energy_thm2570.py
output: 05-knowledge/results/lrc14_deep_augmentation_energy_thm2570.out
script_sha256: 40809dbb1f03bf48541ad32183ab7150e7d691d974d15939a7f47790b271ef20
output_sha256: 7c52420949b10bd70338ac004aa5b7f3ae8fcb03937aa4ff0a175c1168e45cd5
hash_basis: LF-normalized bytes
---

# THM-2570 -- the cancelled coloured cycle has positive Gram energy

**PROVED CANDIDATE + VERIFIED-EXACT; INDEPENDENT HOSTILE AUDIT
PENDING.**

THM-2567 finds a genuine linear cancellation: every nonzero deep colour has
a rigid duty face, but the sum over all thirteen colours is zero.  Squaring
before forgetting colour changes the answer.  The resulting Parseval class
has a sharp universal floor, and all gains see the same quantitatively
selected deep colour.

This is a nonlinear holotopy obstruction only in the finite algebraic sense.
It is not yet a physical Gram measurement on the LRC carrier.

## 1. The deep anchor profile

Let `p=13`, and let

```text
H:F_p^3 -> Q_(>=0),                 H not identically zero,

H(t,s,t)=0                          for every s,t.             (1)
```

Use THM-2567's normalized transform and physical target reindexing

```text
B(m,b,h)
 =p^(-3) sum_(r,s,t) H(r,s,t) zeta^(mr+bs+ht),

A_m(q_1,q_2)=B(m,q_1,q_2-m).                                (2)
```

Put

```text
rho=sum_(r,s,t)H(r,s,t)>0,

G(u)=sum_(s,t)H(u+t,s,t),

a_m=A_m(0,0)=p^(-3)sum_u G(u)zeta^(mu).                     (3)
```

Then

```text
G(0)=0,                    sum_u G(u)=rho.                    (4)
```

THM-2567 proves qualitatively that every `a_m`, `m!=0`, is nonzero.  The
present theorem computes the sharp aggregate size.

## 2. Exact Parseval identity and sharp floor

The nonzero-colour anchor energy is exactly

```text
E_deep:=sum_(m!=0)|a_m|^2

 =p^(-6)(p sum_u |G(u)|^2-rho^2)                             (5)

 =rho^2/(12p^6)
   +p^(-5)sum_(u!=0)|G(u)-rho/12|^2.                         (6)
```

Consequently

```text
E_deep >=rho^2/(12*13^6),                                   (7)
```

with equality if and only if

```text
G(1)=...=G(12)=rho/12.                                      (8)
```

### Proof

Parseval for the unnormalized Fourier transform gives

```text
sum_m |sum_u G(u)zeta^(mu)|^2=p sum_u|G(u)|^2.              (9)
```

The zero colour is `rho`, so deleting it and dividing by `p^6` proves
(5).  Substitute `G(0)=0` and complete the square around the mean of the
twelve remaining entries to obtain (6).  This proves (7) and its equality
case.

Equivalently, THM-2567's augmentation law gives

```text
sum_(m!=0)a_m=-a_0=-rho/p^3,                                (10)
```

and Cauchy on twelve colours gives (7) directly.  Thus the aggregate floor
needs only the diagonal zero and nonzero mass; rationality and nonnegativity
are needed for THM-2567's stronger *every-colour* nonvanishing.

The hostile

```text
H(r,s,t)=1_(s!=0)1_(r!=t)                                  (11)
```

has `rho=1872`, `G(u)=156` for `u!=0`, and

```text
a_m=-12/169                         for every m!=0,

E_deep=1728/28561.                                          (12)
```

It attains equality in (7).  The stronger three-zero-plane hostile from
THM-2567 also has uniform off-zero `G` and attains equality.

## 3. Rank-one duty-replica Gram packet

For a gain `g!=0`, let `K_g` and `d_g` be the unnormalized canonical duty
commutator and anchor jump from THM-2562.  Restrict its output to the six
uncontaminated target coordinates and twelve nonzero deep colours:

```text
R_g(j,m):=(K_g A_m)(jg),          j=1,...,6, m!=0.           (13)
```

THM-2562 gives the exact outer-product factorization

```text
R_g(j,m)=-d_g a_m.                                         (14)
```

Hence `R_g` has rank one and

```text
||R_g||_F^2
 =6d_g^2 E_deep
 >=d_g^2 rho^2/(2*13^6).                                   (15)
```

Its row Gram matrix is

```text
d_g^2 E_deep J_6,                                          (16)
```

with sole nonzero eigenvalue `6d_g^2 E_deep`.  Equality in (15) is attained
by (11).  Sharpness refers to these six displayed coordinates; the full
commutator can have additional energy.  For example, on (11) and the gain
`g=(1,0)`, its full squared norm is `25/24` times (15)'s displayed energy.

The gain classes are

```text
d_g=229692  on 24 axes,
    =440232  on 12 antidiagonal gains,
    =440244  on the other 132 gains.                         (17)
```

Therefore the stack over all gains is still rank one in the replica/deep
factorization, while its displayed energy obeys

```text
sum_(g!=0)||R_g||_F^2
 >=14587701710688/4826809 * rho^2.                           (18)
```

For the normalized duty `nu=kappa n` of THM-2562, multiply (15) and (18)
by `kappa^2`.

## 4. One quantitative colour works for every gain

Averaging (7) over the twelve nonzero colours gives one `m_*!=0` with

```text
|a_(m_*)|>=rho/(12*13^3)=rho/26364.                         (19)
```

The same `m_*` is independent of `g`.  Equations (14) and (19) therefore
give, simultaneously for every gain and every `j=1,...,6`,

```text
|R_g(j,m_*)|>=d_g rho/26364.                                (20)
```

This adaptive common-colour bound is sharp at (11).  It is stronger than
selecting a fresh deep colour for each gain, but it does not prescribe the
colour in advance.

## 5. Application to the paired-dipole corner

For THM-2563, let `rho_0` be its source weight and `M=sum R_(r,s,t)` its
positive table mass.  Its capacity bounds give

```text
M>=63rho_0                 when exactly one graft role is the guard,

M>=81rho_0                 when both roles are ordinary.                    (21)
```

Applying (7) with `rho=M` yields

```text
E_deep >=1323/19307236 rho_0^2       in the guard case,

E_deep >=2187/19307236 rho_0^2       in the ordinary case.                  (22)
```

The displayed per-gain energy is six times `d_g^2` times these floors.
The constants in (22) inherit the coarse `63/81` capacity bounds and are not
claimed sharp for the full THM-2563 geometry.

## 6. No prescribed-colour amplitude floor

Qualitative every-colour survival does not admit a uniform quantitative
upgrade.  Put `theta=2pi/13` and choose positive rationals

```text
delta_n ->(2cos theta)^(-1).
```

For fixed mass `rho`, set

```text
c_n=rho/(12+2delta_n),

G_n(1)=G_n(-1)=c_n(1+delta_n),

G_n(u)=c_n                    for u!=0,+-1,

G_n(0)=0.                                                   (23)
```

Then

```text
a_1(G_n)
 =c_n p^(-3)(-1+2delta_n cos theta) ->0.                    (24)
```

No term is zero: otherwise the rational profile would contradict the
`Phi_13` argument of THM-2567 (equivalently, `2cos theta` would be rational).
These profiles are realized even with all three THM-2563 zero loci.  For
`u!=0,-1`, place mass `G_n(u)` at `(r,s,t)=(u+1,1,1)`; place the `u=-1`
mass at `(1,1,2)`.  Thus `r=0`, `s=0`, and `r=t` all remain empty.

So (19) is the correct uniform conclusion: one adaptive colour is large,
while a prescribed colour can be arbitrarily small.

## 7. What the quadratic class preserves and loses

Linear deep augmentation kills the THM-2567 cycle, but the Hermitian Gram
packet (16) cannot cancel.  This is the exact nonlinear survivor.  It also
shows a limit: stacking all gains raises the energy but not the rank, so
another gain does not manufacture a deep selector.

THM-2383 explains what would make the Gram class observable: a lawful
complex quadrature against a spanning oriented reference on the **same**
carrier.  No such reference is supplied here.  Squaring derived
partial-bare coefficients is not itself a Boolean overlap, Abel current, or
Hall edge.  The missing sidecar remains the common deep-refined
Reynolds/endpoint carrier and its physical phase reference.

## 8. Exact companion and stopping boundary

Run

```bash
python3 04-computation/lrc14_deep_augmentation_energy_thm2570.py
python3 -O 04-computation/lrc14_deep_augmentation_energy_thm2570.py
```

Both executions must reproduce

```text
05-knowledge/results/lrc14_deep_augmentation_energy_thm2570.out
```

byte-for-byte.  The dependency-free exact companion exhausts all `531,440`
nonzero ternary off-diagonal profiles and all `4,095` Boolean profiles,
checks every nontrivial anchor on the latter, replays all `2,028`
diagonal-free singleton cells, verifies both equality hostiles, all `168`
gain classes, the global factor in (18), the two THM-2563 floors, and the
strict triangle geometry behind the prescribed-colour infimum.  It executes
`1,647,551` explicit integer/Fraction assertions under both normal and
optimized Python.

This theorem upgrades an exact linear cancellation to a sharp positive
coefficient-level energy and a common quantitative colour.  It does not
construct the common carrier, a physical Gram measurement, a canonical
selector/current, a scalar-row exclusion, or LRC(14). **QED.**
