---
id: THM-2525
title: "Unit-guard collision floor and word--owner cross-cospan collapse"
status: >
  PROVED + VERIFIED-EXACT + INDEPENDENTLY AUDITED.  A nonzero Boolean response supported in the
  unit guard C_H={||Hx||>1/7}, with 13 not dividing H, has depth-one K_13
  collision drift at least 3/13 of its mass.  The constant is sharp, and
  equality means that almost every occupied root fibre contains exactly ten
  occupied predecessors.  Every translated chi_7 bank then has a signed
  scalar of modulus at least 845/364 of the event mass; at the guard-aligned
  slope the sharp constant is 1014/305.  Applied to every THM-2349 live word
  event, this upgrades qualitative first-collision nonvanishing to uniform
  relative floors.  However, the apparently ordered cross-cospan from the
  word-bearing owner event to its bare owner is pointwise identical to the self-cospan:
  the delayed word is already constant on every depth-one root fibre.  Its
  cross profile and translated chi_7 bank are therefore antipodally even and
  contain no source/arrival orientation.  More generally, unshifted
  target-zero scalar-comb events are reflection-even.  A genuinely odd
  Boolean repair must retain a root-visible pre-collision or chiral target
  sidecar.  No row exclusion or LRC(14) proof is claimed.
source: codex-2026-07-27-unit-guard-cross-cospan
depends_on:
  - THM-2349-first-depth-one-delayed-shallow-restart
  - THM-2522-intrinsic-collision-depth-toothpick-descent-and-late-owner-decoupling
  - THM-2524-translated-chi7-hamilton-polarization-inversion
related:
  - THM-2365-lawful-target-coshift-and-h-drift-dichotomy
  - THM-2478-delayed-owner-handoff-graft-and-deep-sheet-rebase-boundary
  - THM-2518-perron-inverse-branch-owner-word-cospan-recovery
script: 04-computation/lrc14_unit_guard_cross_cospan_thm2525.py
output: 05-knowledge/results/lrc14_unit_guard_cross_cospan_thm2525.out
script_sha256: 975e97e86fe79849a875a8edda147cf0f6e60bf5621fa62b66fba24f9d419f24
output_sha256: f5ba3ff8b9d2507773818eab9c40d72b4217aa4312c2f8707e4ffd5988efe0c1
hash_basis: working-tree bytes (LF)
---

# THM-2525 -- the live word--owner cospan is quantitatively active but unoriented

**PROVED + VERIFIED-EXACT + INDEPENDENTLY AUDITED.**

THM-2522 proves that the unit guard forbids a live word event from being a
complete depth-one replica.  The same geometry has a quantitative form.  A
thirteen-point root fibre places at most ten points in the guard, so every
occupied fibre pays at least three missing predecessors.  This yields an
exact relative floor for the first collision:

```text
nonzero Boolean F supported in C_H
  -> D_13(F;1,1) >= (3/13) measure(F) > 0.                    (1)
```

The most obvious attempt to turn this collision into an ordered Boolean
object is to distinguish the terminal-word refinement `F=E W` from its bare
owner `E`.  That attempt fails for an exact reason.  At the depth-one
collision, every delayed word factor is already common to the two
predecessors.  Boolean idempotence then makes the ordered cross product
literally equal to the old self product, before integration or Fourier
transform:

```text
F(x)E(x+u/13)=F(x)F(x+u/13).                                  (2)
```

Thus the semantic leg labels do not add an orientation coordinate.  The
positive gain in (1) is real, but it lives in the same even collision profile
already inverted by THM-2524.

## 1. The thirteen-point unit-guard invoice

Put

```text
C_H={x in R/Z: ||H x||>1/7},                 13 does not divide H.  (3)
```

For every `x`, multiplication by `H` permutes the depth-one root phases, so

```text
#{u in F_13:x+u/13 in C_H}
 =#{u in F_13:Hx+u/13 in C_1}.                           (4)
```

The set `C_1` is one open circular arc of length `5/7`.  Eleven points of a
thirteen-grid have circular span at least `10/13`, and

```text
10/13>5/7.                                                   (5)
```

Consequently

```text
#{u:x+u/13 in C_H}<=10.                                     (6)
```

Away from the finite endpoint phases the count is exactly nine or ten,
because `13(5/7)=65/7` lies strictly between nine and ten.  Only the upper
bound in (6) is needed below.

Let `F` be a nonzero rational Boolean step function with

```text
support(F) subset C_H                                             (7)
```

up to null endpoints.  Define its depth-one collision profile

```text
B_u=integral_T F(x)F(x+u/13)dx,               u in F_13,      (8)

rho=integral_T F=measure(support(F)).                            (9)
```

Since `F` is Boolean, `B_0=rho`.  Fubini and (6) give

```text
sum_u B_u
 =integral_T F(x) sum_u F(x+u/13)dx
 <=10 rho.                                                   (10)
```

Therefore THM-2519/2522's depth-one drift obeys

```text
D_0
 :=B_0-(1/13)sum_u B_u
 >=(3/13)rho
 >0.                                                         (11)
```

This is a relative floor.  THM-2349 proves only row-dependent positivity of
`rho`, so (11) is not an absolute uniform mass bound over the `165` rows.

There is an exact equality boundary.  Put

```text
m_F(x)=sum_u F(x+u/13).                                      (12)
```

This is constant on each root orbit and at most ten on every occupied
orbit.  Subtracting the right side of (11) from the left gives

```text
D_0-(3/13)rho
 =(1/13)integral_T F(x)(10-m_F(x))dx.                        (13)
```

Hence

```text
D_0=(3/13)rho
 iff m_F(x)=10 for almost every x with F(x)=1.               (14)
```

Equivalently, almost every occupied root fibre has exactly ten occupied
predecessors.  If every occupied fibre has at most nine, the stronger floor
`D_0>=(4/13)rho` follows by the same proof.

The constant `3/13` is sharp even inside the shallow unit-guard geometry.
Take `H=1` and

```text
F=1_(D_13 intersection C_1),
D_13={x:||13x||<1/14}.                                      (15)
```

On every active fibre `z=13x in D_1`, precisely the roots

```text
r=2,3,...,11
```

lie in `C_1`; the roots `0,1,12` do not.  Thus

```text
rho=10/91,
D_0=30/1183=(3/13)rho.                                      (16)
```

The strict/open endpoint convention changes only null sets.

## 2. Every rational root colour survives

The vector `B=(B_u)` is rational.  Equation (11) makes it nonconstant.  If
one nontrivial normalized root transform vanished,

```text
Bhat(a)=(1/13)sum_u B_u zeta_13^(-au)=0,
                                                       a!=0, (17)
```

then the rational polynomial `sum_u B_u X^u` would be divisible by
`Phi_13(X)=1+X+...+X^12`.  Its thirteen coefficients would all be equal,
contrary to (11).  Hence

```text
Bhat(a)!=0                         for every a in F_13^*.     (18)
```

In fact these coefficients are positive: the self-correlation identity of
THM-2519 writes each one as the squared norm of the corresponding root
character.  Their sum is exactly `D_0`, so (11) also gives

```text
sum_(a!=0) Bhat(a)>=(3/13)rho.                               (19)
```

No comparable lower bound for each individual algebraic conjugate follows
from (19) alone.

## 3. The word-bearing/bare-owner cross-cospan collapses

The collapse is an abstract Boolean identity.  Let `E,W:T->{0,1}` be
measurable and suppose

```text
W(x+u/13)=W(x)                  for every u in F_13 a.e.      (20)
```

Put

```text
F=E W.                                                         (21)
```

For the apparently ordered profiles

```text
C^(F,E)_u=integral_T F(x)E(x+u/13)dx,
C^(E,F)_u=integral_T E(x)F(x+u/13)dx                          (22)
```

and the self profile `B` of (8), Boolean idempotence and (20) give

```text
F(x)F(x+u/13)
 =E(x)W(x)E(x+u/13)W(x+u/13)
 =F(x)E(x+u/13).                                             (23)
```

Thus

```text
C^(F,E)_u=B_u.                                                (24)
```

Every real self-correlation is antipodally even:

```text
B_(-u)=B_u.                                                   (25)
```

Reversing the two legs in a cross-correlation reverses `u`, so (24)--(25)
also give

```text
C^(E,F)_u=C^(F,E)_u=B_u.                                     (26)
```

The two nominally ordered cospans are identical.  This conclusion does not
require `E` or `W` to be reflection-even; root invariance and Boolean
idempotence suffice.

The same statement can be expressed on THM-2524's predecessor module.  Put

```text
p^F_r(z)=F((z+r)/13)-(1/13)sum_j F((z+j)/13),
p^E_r(z)=E((z+r)/13)-(1/13)sum_j E((z+j)/13).                 (27)
```

Then the centred cross profiles satisfy

```text
b^(F,E)_u
 =(1/13)integral_T sum_r p^F_r(z)p^E_(r+u)(z)dz
 =C^(F,E)_u-(1/13)sum_v C^(F,E)_v,                           (28)
```

and similarly for `b^(F,F)`.  Equations (24)--(26) imply

```text
b^(F,E)=b^(E,F)=b^(F,F),
(b^(F,E))^-=0.                                                (29)
```

For every `tau!=0`, THM-2524 therefore gives not a new ordered bank but the
exact old one:

```text
R_tau^(F,E)
 =13 A_tau b^(F,E)
 =13 A_tau b^(F,F)
 =R_tau^(F,F).                                                (30)
```

It is lossless on the nonconstant even profile and has all twelve nonzero
root modes by (18), but its odd/source-order component is identically zero.

## 4. Uniform application to all `165` live rows

On every first-depth-one row, THM-2349 supplies

```text
E=1_(E_j),
W(x)=1_(Q_(j,sigma))(13^k x),
F=E W,
k>=2,
rho=measure(E_(j,sigma,k))>0.                                (31)
```

It also proves

```text
support(F) subset E_j subset A_0 subset C_H,
13 does not divide H.                                        (32)
```

The word is root-constant because

```text
13^k(x+u/13)=13^k x+13^(k-1)u=13^k x        mod 1.           (33)
```

Sections 1--3 apply verbatim.  Thus every live row has the quantitative
first-collision invoice

```text
D_0>=(3/13)measure(E_(j,sigma,k))>0,                          (34)
```

all twelve rational collision colours, and a lossless translated `chi_7`
bank.  But marking one leg as `F` and the other as `E` does not turn that
bank into an emitted owner--word edge.  On every contributing depth-one
pair, (33) puts the same terminal word on both predecessors; equation (23)
is the exact algebraic record of that semantic loss.

## 5. A sharp signed scalar survives the collapse

Although the word/bare label creates no new profile, the unit guard forces a
quantitative entry in the signed translated bank.  First, an elementary
bound works for every slope.  Put

```text
N_tau=A_tau^5-39A_tau^3+299A_tau.                            (35)
```

For `tau=1`, the zeroth row is

```text
(N_1)_(0,*)
 =(0,85,15,-15,-5,-85,5,5,-85,-5,-15,15,85),

||(N_1)_(0,*)||_1=420.                                      (36)
```

Every other nonzero slope only permutes this row.  THM-2524's inverse says

```text
D_0=b_0=(1/4225)(N_tau R_tau)_0.                             (37)
```

Consequently (11), (36), and (37) give, for every `tau!=0`,

```text
max_t |R_tau(t)|
 >=(4225/420)D_0
 >=(845/364)rho.                                             (38)
```

The guard geometry permits a stronger sharp estimate after aligning the
slope with `H`.  Split the circle into the thirteen inverse branches of
`z=13x` and put

```text
e_r(z)=F((z+r)/13),
n(z)=sum_r e_r(z),
c_t(z)=sum_r e_r(z)e_(r+t)(z).                               (39)
```

Then

```text
B_t=(1/13)integral_T c_t(z)dz,
rho=(1/13)integral_T n(z)dz,
R_tau=integral_T A_tau c(z)dz.                               (40)
```

The last equality uses `R_tau=13A_tau(B-average(B))` and the fact that
`A_tau` kills constants.

Reindex each root fibre by the unit `H`.  The complement of `C_H` is a
closed arc of length `2/7`, so the zero set of every nonempty bit profile
`e(z)` contains at least three consecutive roots in this `H`-gauge.  Choose

```text
tau H=1                         mod 13                         (41)
```

and reorder the translation coordinate by the same unit.  This turns the
last identity in (40) into the `A_1` bank without changing its maximum
absolute entry.

Here is the exact finite lemma.  Let `e in {0,1}^13` be nonzero and suppose
its zero set contains a cyclic run of length three.  Put

```text
c_t=sum_r e_r e_(r+t),
r=A_1 c,
m=sum_r e_r.                                                  (42)
```

Then

```text
9r_0+368r_1-340r_4-198r_5 >=234m.                            (43)
```

This is an exact finite statement.  The companion exhausts all `5,434`
admissible nonzero masks, using integer arithmetic; the minimum slack is
zero and there are `78` equality masks.  Since the coefficient `l1` norm in
(43) is

```text
9+368+340+198=915,                                           (44)
```

integrating (43) in (40) yields

```text
915 max_t |R_tau(t)|
 >=234 integral_T n(z)dz
 =3042 rho.
```

Therefore the guard-aligned slope satisfies the stronger bound

```text
max_t |R_tau(t)| >=(1014/305)rho.                             (45)
```

The constant in (45) is sharp for the stated Boolean unit-guard class.  An
exact primal certificate uses the three masks

```text
0x25b  (mass 6),
0x2db  (mass 7),
0x3ff  (mass 10)                                             (46)
```

with probabilities

```text
52/429,                 39/429,                 338/429.     (47)
```

All three zero sets contain the common run `{10,11,12}`.  As a normalized
finite mixture it has

```text
rho=305/429,
max_t |R_1(t)|=1014/429,
```

so equality holds in (45).  For a physical realization take `H=1` and a
small rational base interval inside `D_1`.  Its guard complement is exactly
the cyclic run `{12,0,1}`.  Rotate the three masks so their common zero run is
`{12,0,1}`, partition a rational subset of that interval of measure `delta`
in the proportions (47), and use the masks on its inverse branches.  Cyclic
rotation does not change an autocorrelation or its bank.  The resulting
circle event has

```text
rho=delta 305/429,
max_t |R_1(t)|=delta 1014/429,
```

and is a rational Boolean subset of `D_13 intersection C_1`.  Thus the sharp
certificate is not merely an abstract fractional mask.

Equations (38) and (45) are quantitative signed-scalar survivors.  The
selected `R_tau(t)` can have either sign, is a linear combination of Boolean
intersections, and remains the self-bank by (30).  Neither bound makes it a
positive Boolean atom or an ordered source/arrival edge.

## 6. Reflection no-go and the first legal escape

There is a second, slightly broader no-go with deliberately narrow scope.
Let `mathcal B_0` be the target-zero Boolean algebra generated by the
unshifted centred scalar combs

```text
1_(||n x||<alpha),
1_(||n x||>alpha),                                            (48)
```

using finite Boolean operations and integer dilations.  Every member is
reflection-even:

```text
f(-x)=f(x).                                                   (49)
```

Consequently, for any two real `f,g in mathcal B_0`, the unweighted root
cross profile satisfies

```text
K_u=integral_T f(x)g(x+u/13)dx,
K_(-u)=K_u.                                                   (50)
```

Indeed, substitute `x -> -x`.  This applies to the unshifted target-zero
owners, word factors, guard factors, and danger-comb factors used in
THM-2349.  It does **not** apply after a lawful nonzero target co-shift, to
an absolute inverse-branch sheet with its carry retained, or to a signed
endpoint current.

Equations (23) and (50) identify the cheapest legal escapes.  A Boolean odd
cross-cospan must retain at least one of:

1. a pre-collision factor which is not constant under `x -> x+u/13`;
2. a reflection-chiral sidecar, such as one fixed nonzero THM-2365 target
   co-shift rather than its `+/-` quotient; or
3. an absolute ordered ancestry sheet/carry before THM-2518's simultaneous
   sheet sum.

A signed endpoint-current pairing can also be odd, but it is not a Boolean
repair.  These are necessary escape coordinates, not assertions that their
odd components are nonzero on every live row.  The cheapest next decisive
test is therefore a target-coshifted word/owner cross bank with the target
label retained under reflection; another unshifted word/bare-owner pairing
cannot work.

## 7. Exact companion and stopping boundary

Run

```bash
python3 04-computation/lrc14_unit_guard_cross_cospan_thm2525.py
python3 -O 04-computation/lrc14_unit_guard_cross_cospan_thm2525.py
```

Both executions reproduce

```text
05-knowledge/results/lrc14_unit_guard_cross_cospan_thm2525.out
```

byte-for-byte.  The dependency-free companion verifies the exact unit-guard
phase-cell counts, all `8,099` nonempty thirteen-bit fibres allowed by the
ten-point carrier bound, the collision floor and equality boundary, the
cross/self collapse and reversal symmetry, the sharp ten-of-thirteen
profile, and the rational cyclotomic all-colour criterion.  It also derives
the inverse-numerator row and its `l1` norm, exhausts all `5,434`
three-zero-run profiles for (43), and verifies both exact sides of the sharp
primal/dual certificate (45)--(47).

The theorem improves the live collision from qualitative support to a sharp
relative energy invoice and identifies the literal Boolean reason the
obvious ordered construction fails.  It does not produce a chiral target
sidecar, preserve an absolute branch sheet, transfer the signed bank into a
scalar-cover row, exclude any of the `165` rows, or prove LRC(14). **QED.**

An independent audit rederived the ten-point guard invoice, Fubini floor and
equality boundary, word/bare pointwise collapse, cyclotomic all-colour step,
inverse-numerator row and normalization, `H tau=1` gauge transport, exact
dual inequality, and matching rational shallow-carrier realization.  It also
reproduced normal, optimized, and stored companion output byte-for-byte and
confirmed the recorded hashes.
