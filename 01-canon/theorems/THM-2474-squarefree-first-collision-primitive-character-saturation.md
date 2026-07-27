---
id: THM-2474
title: "Squarefree first-collision primitive-character saturation"
status: >
  PROVED + VERIFIED-EXACT + INDEPENDENTLY HOSTILE-AUDITED. Let N be
  squarefree and let two nonnegative rational step packets remain
  disjoint after every codimension-one Perron quotient P_(N/p), but
  overlap with mass I after P_N. Their N-root correlation is supported
  on units. Every primitive character is nonzero; their signed sum is
  mu(N)I, their squared energy is at least I^2/phi(N), and some has
  modulus at least I/phi(N). Every prescribed unit residue has a common
  original Fourier landing in each endpoint-Prony block. This upgrades
  every interior THM-2313 Fano corner from some 91-unit class to all 72
  classes and removes THM-2318's unselected-root-residue debt. The
  multiplier-four three-prime control has the explicit grade-three,
  root-four atom 13^3*17 and an incident m=1 deepest edge, an incidence
  already covered by THM-2327 on every positive canonical shallow-owner
  word stratum in the 150-row strict bank. For
  nonsquarefree N, singleton-versus-unit root packets satisfy every
  stated predecessor zero but all primitive transforms vanish because
  mu(N)=0, so squarefreeness is sharp for these hypotheses. Pair
  incidence outside THM-2327's positive canonical shallow-owner word
  strata, target gain, base phase, and forcing a one-shot corner on all
  rows remain open; no scalar row is excluded.
source: codex-2026-07-27-squarefree-collision-saturation
depends_on:
  - THM-2313-biprime-pareto-collision-frontier-and-91-unit-current-shell
  - THM-2318-one-shot-three-prime-mobius-amplifier
  - THM-2323-primitive-fixed-colour-cross-correlation-and-same-gauge-word-alignment
  - THM-2398-prime-cyclic-rational-restoration-dichotomy
  - THM-2449-coprime-owner-anova-and-delta-replica-boundary
  - THM-2471-owner-first-collision-weighted-root-service-and-temporal-atom-boundary
related:
  - THM-2293-quadratic-root-energy-raises-the-ancestry-shell
  - THM-2299-rooted-current-service-energy-and-base-phase-no-go
  - THM-2306-owner-normalized-disjoint-support-and-first-collision-shell
  - THM-2319-crt-unit-bispectrum-needle-and-mixed-polarization-no-go
  - THM-2326-vertexwise-septimally-primitive-c3-degree
  - THM-2327-two-colour-marked-unit-c3-triangle
script: 04-computation/lrc14_squarefree_collision_saturation_thm2474.py
output: 05-knowledge/results/lrc14_squarefree_collision_saturation_thm2474.out
script_sha256: bea849bc574326a1b71be8137eceb955a45e98cb87fc78c2b6ee67add9500c40
output_sha256: bfc0882e004d0df0ea8f55a2ee9012d33865dfcf9d263ee7c77d7d5c6ef048f5
hash_basis: working-tree bytes (LF)
---

# THM-2474 -- a full squarefree collision corner fires every unit colour

**PROVED + VERIFIED-EXACT + INDEPENDENTLY HOSTILE-AUDITED.**

THM-2313 and THM-2318 discovered the correct collision geometry: add one
Perron direction for each prime whose divisibility must be removed, and
take a fully minimal corner of the resulting collision upper set.  Möbius
inclusion--exclusion then isolates a nonzero aggregate over frequencies
coprime to the product.  Their landing argument selected *some* unit
residue progression.

At a squarefree corner the physical root packet says more.  Every
codimension-one predecessor zero erases one coordinate hyperplane of the
root difference table.  The correlation is consequently supported on the
unit group itself.  Rational Galois conjugacy and a Ramanujan sum then
force not one but **every** primitive character.

```text
codimension-one predecessor zeros
  -> correlation supported on (Z/NZ)^*
  -> primitive Galois orbit is all-or-all
  -> its Ramanujan sum is mu(N)I != 0
  -> every primitive colour survives.                        (1)
```

This is the recursive toothpick law behind the one-prime `13` tooth, the
`7 x 13` Fano sheet, and the three-prime Kakeya needle.

## 1. Squarefree root-correlation theorem

Let

```text
N=product_(i=1)^d p_i>1                             (2)
```

be a product of distinct primes.  Let `U,V:T->R_(>=0)` be nonzero
rational-valued step functions with rational breakpoints.  For every
prime `p|N`, assume the codimension-one predecessor is disjoint:

```text
integral_T (P_(N/p)U)(P_(N/p)V)=0.                  (3)
```

Assume the full quotient collides:

```text
I=integral_T(P_NU)(P_NV)>0.                         (4)
```

On the canonical oriented `Z/NZ` root chart put

```text
A(y,r)=U((y+r)/N),

F(y,r)=V((y+r)/N),                                 (5)
```

and define the unnormalized difference correlation

```text
C(s)=integral_T sum_(r mod N) A(y,r+s)F(y,r)dy.    (6)
```

Then

```text
C(s)>=0,

C(s)=0                         if gcd(s,N)>1,

sum_s C(s)=N^2 I.                                  (7)
```

For `k mod N`, use normalized amplitudes

```text
alpha_k(y)=(1/N)sum_r A(y,r)zeta_N^(-kr),

phi_k(y)=(1/N)sum_r F(y,r)zeta_N^(-kr),

J(k)=integral_T alpha_k(y)conjugate(phi_k(y))dy.   (8)
```

Every primitive colour survives:

```text
J(k)!=0                         for gcd(k,N)=1.      (9)
```

The primitive sector has the exact signed ledger and quantitative floors

```text
sum_(gcd(k,N)=1) J(k)=mu(N)I,                       (10)

sum_(gcd(k,N)=1)|J(k)|^2>=I^2/phi(N),

max_(gcd(k,N)=1)|J(k)|>=I/phi(N).                  (11)
```

Here `mu(N)=(-1)^d`; in particular its magnitude is one.

## 2. Why predecessor zeros erase every nonunit shift

Fix `p|N`.  Since the integrand in (3) is nonnegative, (3) implies

```text
(P_(N/p)U)(z)(P_(N/p)V)(z)=0                      (12)
```

for almost every `z`.  For `t mod p`, evaluate at `z=(y+t)/p`.  Apart
from the positive factor `p/N`, the two predecessor values are the block
sums

```text
sum_(a mod N/p) U((y+t+pa)/N),

sum_(b mod N/p) V((y+t+pb)/N).                     (13)
```

Their product is zero.  Nonnegativity makes every cross term in that
product zero.

If `p|s`, then `r` and `r+s` have the same residue `t mod p`.  The term

```text
A(y,r+s)F(y,r)                                     (14)
```

is one of the vanished cross terms in (13).  Hence `C(s)=0`.  Doing this
for every `p|gcd(s,N)` proves the support assertion in (7).

Summing (6) over all shifts gives

```text
sum_s C(s)
 =integral_T (sum_r A(y,r))(sum_r F(y,r))dy
 =N^2 integral_T(P_NU)(P_NV)
 =N^2I,                                             (15)
```

which completes (7).

The proof reveals the geometric object: the positive shifts form a set
of **Kakeya needles avoiding every coordinate hyperplane** in the CRT
root cube.  Marginal spectra cannot certify this; the simultaneous
crosshair zero is load-bearing.

## 3. Galois saturation and the Ramanujan ledger

Finite inversion gives

```text
J(k)=(1/N^2)sum_s C(s)zeta_N^(-ks).                 (16)
```

Every `C(s)` is rational because the packets and all breakpoints are
rational.  The primitive values in (16) are therefore one Galois orbit:

```text
sigma_u(J(k))=J(uk),                gcd(u,N)=1.      (17)
```

Thus either all primitive `J(k)` vanish or none do.

Sum (16) over primitive `k`.  The inner sum is the Ramanujan sum

```text
c_N(s)=sum_(gcd(k,N)=1)zeta_N^(-ks).                (18)
```

By (7), only unit `s` occur; for a unit,

```text
c_N(s)=mu(N).                                       (19)
```

Equations (7), (16), and (19) give (10).  Squarefreeness makes
`mu(N)=+-1`, so the sum cannot vanish.  The all-or-all law (17) proves
(9).  Cauchy--Schwarz and the triangle inequality on the `phi(N)` values
prove (11). QED.

For `d=1`, (10)--(11) are exactly THM-2471's `-I`, `I^2/12`, and
`I/12` invoices.  The parity of the number of prime directions explains
the alternating sign of successive collision derivatives.

## 4. Every prescribed unit class lands syndetically

Physicalizing (5) gives the original functions `U,V`.  The root-gauge
identity is

```text
J(k)=sum_(h in Z)
 U_hat(k+Nh)conjugate(V_hat(k+Nh)).                 (20)
```

The series is absolutely convergent by Cauchy--Schwarz.  Equation (9)
therefore makes the endpoint-product sequence on every prescribed unit
class `q=k+Nh` nonzero.

Let `J_U,J_V` be the numbers of nonzero jumps of `U,V`, and put

```text
L=J_UJ_V.                                           (21)
```

After multiplying by `(2*pi*q)^2`, the product in (20) is an exponential
sum on at most `L` endpoint-difference nodes.  Restricting to one fixed
progression `q=k+Nh` preserves that bound.  A nonzero `L`-node exponential
sum cannot vanish at `L` consecutive indices.  Hence for every unit `k`
and every integer `H`, some

```text
H<=h<=H+L-1                                        (22)
```

satisfies

```text
U_hat(k+Nh)V_hat(k+Nh)!=0.                          (23)
```

Taking the least positive representative `1<=k<=N-1` and `H=0` gives

```text
1<=q<=NL-1,

q congruent k mod N.                               (24)
```

If `U=P_DU_0` and `V=P_DV_0`, equation (23) is equivalently

```text
(U_0)_hat(Dq)(V_0)_hat(Dq)!=0.                     (25)
```

Thus the theorem selects the CRT colour *before* endpoint-Prony; it does
not merely discover a favourable class after aggregation.

## 5. Collision Newton frontiers and spectral rank

Let distinct primes index a commuting Perron staircase

```text
I_(r_1,...,r_d)
 =integral_T
   (product_i P_(p_i)^r_i U_0)
   (product_i P_(p_i)^r_i V_0).                    (26)
```

The positive set is upper.  At any Pareto-minimal point, let `S` be the
coordinates whose entries are positive.  Pull back one step in precisely
those coordinates.  Every immediate predecessor in `S` is zero and the
full `product_(i in S)p_i` quotient is positive.  Sections 1--4 apply with

```text
N_S=product_(i in S)p_i.                            (27)
```

The support size `|S|` is therefore an exact **collision spectral rank**:

```text
|S|=1: one-prime tooth;

|S|=2: Fano/CRT sheet;

|S|=3: Kakeya needle;

general |S|: every one of phi(N_S) full-rank colours.       (28)
```

Axis-born overlap and full interior collision are now mathematically
distinct, rather than two visual descriptions of the same aggregate.

## 6. Upgrades to the existing LRC collision frontiers

### 6.1 THM-2313: every interior Fano corner has all 72 colours

At an interior Pareto-minimal THM-2313 cell `(s,t)`, put

```text
U=P_13^(s-1)P_7^(t-1)A,

V=P_13^(s-1)P_7^(t-1)B.                            (29)
```

The west and south predecessor zeros are exactly (3) for `N=91`, and the
corner mass is (4).  Therefore all

```text
phi(91)=72                                           (30)
```

primitive CRT colours survive.  For every prescribed unit class `k mod
91`, (24) improves THM-2313 equation (21) from “some unit class” to that
exact class.

At its grade-three hostile corner `(s,t)=(3,7)`, THM-2313 needs
`n congruent 9 mod 13` to match the shallow character `4`.  Choose any
unit CRT class with that thirteenth residue.  Theorem (24) now lands it.
The extra factor `7^6` and the lack of pair incidence remain; colour
alignment alone does not erase those coordinates.

### 6.2 THM-2318: the one-shot cube has all primitive colours

At a fully minimal THM-2318 cell `(s,1,1)`, with auxiliary prime `ell`,
put

```text
U=P_13^(s-1)A,

V=P_13^(s-1)B,

N=91ell.                                           (31)
```

The three vanishing faces are exactly (3) for `p=13,7,ell`, and the top
cell is (4).  Hence every one of

```text
phi(91ell)=72(ell-1)                               (32)
```

primitive colours survives.  The prescribed THM-2293 root character may
now be fixed modulo `13` and completed to a unit class modulo `91ell` by
CRT before applying (24).  This removes item 1, “Root residue,” from
THM-2318 Section 5.  It does not identify a target gain or force such a
corner for every live row.

### 6.3 The multiplier-four control already has the desired edge

On THM-2318's exact carrier take `ell=65537` and its minimal corner
`(3,1,1)`.  The normalized source/current interval packets have centers

```text
F: {1,15}/16,             half-width epsilon,

C: {7,9}/16,              half-width 169epsilon,

epsilon=10^(-12).                                  (33)
```

The unit class

```text
n=17,

17 congruent 4 mod 13                              (34)
```

already fires at `h=0`.  Indeed `13^2*17=2873 congruent 1 mod 8`; the
two-center cosine factors are nonzero, and both width sines have strictly
positive arguments below `pi`.  Undoing the owner normalization `c=13`
gives the common word-current/source atom

```text
X=13^3*17=37349,

nu_13(X)=3,

X/13^3 congruent 4 mod 13.                          (35)
```

For the deepest speed

```text
c_3=742586=2*13^5,
```

put

```text
Y=X+c_3=13^3*355.                                  (36)
```

The source coefficient at `Y` is also nonzero: after `P_13`
normalization its frequency is

```text
Y/13=59995 congruent 3 mod 8,                       (37)
```

so the center cosine survives, while the tiny-width sine is nonzero.
Both `X/13^3=17` and `Y/13^3=355` are congruent to `4 mod 13`.
Thus the canonical control contains an incident deepest edge

```text
Y-X=1*c_3,

gcd(1,91)=1,                                       (38)
```

inside the required grade/root class.  This is a positive control, not a
new incidence theorem: THM-2327 already gives the stronger word-marked
`91`-unit incident triangle on every positive canonical shallow-owner word
stratum in the `150`-row strict bank.  Incidence remains a debt outside
that scope, including the repeated-first and alternative resonance branches.

## 7. Squarefreeness is sharp

Let `N` be nonsquarefree.  Realize the following constant root packets by
the rational circle step functions

```text
U=1_[0,1/N),

V=sum_(gcd(r,N)=1)1_[r/N,(r+1)/N).                 (39)
```

On the oriented `N`-root chart they are

```text
A(r)=1_(r=0),

F(r)=1_(gcd(r,N)=1).                                (40)
```

For every prime `p|N`, the `A` packet occupies residue zero modulo `p`,
while `F` occupies only nonzero residues.  Hence every codimension-one
predecessor product in (3) is zero.  The full quotient has positive
service, and its difference correlation is the unit mask.

For every primitive `k`, however,

```text
N^2J(k)
 =sum_(gcd(s,N)=1)zeta_N^(-ks)
 =c_N(k)
 =mu(N)
=0.                                                (41)
```

Thus every primitive colour vanishes despite all stated predecessor zeros
and positive top collision.  Additional prime-power face data could repair
this, but the codimension-one theorem as stated cannot.  This is a minimal
Boolean hostile, not a failure caused by irrational weights or seams.

## 8. Scope and stopping ledger

The theorem proves

```text
fully minimal squarefree collision corner
  -> unit-supported physical root correlation
  -> every primitive character
  -> every prescribed unit Fourier progression
  -> bounded common landing in that progression.                 (42)
```

It closes:

- THM-2313's unselected unit class at every interior corner;
- THM-2318's unselected root residue, conditional on its one-shot corner;
- the corresponding one-, two-, and three-prime spectral-rank dictionary.

It does not:

- force a one-shot corner for every one of the `165` rows;
- remove a pre-existing factor of seven in the owner;
- add pair incidence beyond THM-2327's strict-row scope;
- retain target gain or component base phase;
- intertwine THM-2471's owner-collision root with the THM-2365 deep endpoint;
- exclude a scalar row or prove LRC(14).

## 9. Exact companion

Run

```text
python3 04-computation/lrc14_squarefree_collision_saturation_thm2474.py
python3 -O 04-computation/lrc14_squarefree_collision_saturation_thm2474.py
```

The dependency-free exact companion:

- builds cyclotomic polynomials and reduces every primitive transform
  exactly for squarefree moduli `6,10,15,30,42,91`;
- verifies the Ramanujan signed and energy ledgers using only integer and
  `Fraction` arithmetic;
- checks the crosshair/unit-support implication on explicit CRT root cubes;
- exhausts nonsquarefree hostiles at `4,8,9,12`, where every primitive
  transform vanishes;
- replays the corrected THM-2313 `(5,4)` two-shift Fano control;
- checks the four exact THM-2318 corner inequalities; and
- verifies the `n=17`, `X=37349`, `Y=X+c_3`, grade/root, cosine, width,
  and unit-incidence arithmetic in (33)--(38).

Both transcripts must reproduce

```text
05-knowledge/results/lrc14_squarefree_collision_saturation_thm2474.out
```

byte-for-byte.

An independent hostile audit rederived the unit-support lemma, Galois orbit,
Ramanujan normalization, prescribed-class Prony landing, recursive frontier
specialization, and nonsquarefree Möbius-zero hostile.  It also separated
the explicit multiplier-four incident edge from any unsupported global
pair-incidence claim.

QED.
