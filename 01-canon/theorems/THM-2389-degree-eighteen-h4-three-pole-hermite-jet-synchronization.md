---
id: THM-2389
title: "Degree-eighteen H4 three-pole Hermite-jet synchronization"
status: >
  PROVED CANDIDATE + VERIFIED-EXACT; INDEPENDENT AUDIT PENDING.
  A putative degree-eighteen H4 survivor has a rational normalization
  whose three smooth points above infinity may be placed at
  0,1,infinity. In that coordinate y=A_3/[t(t-1)] and
  v=C_6/[t(t-1)]^2. The three v/y^2 slopes are exactly the three roots
  of one fixed separable cubic. Their missing linear jets impose six
  Hermite conditions and leave only the affine freedom
  C_6=C_*+kappa[t(t-1)]^2. Orders two through five at each pole
  reconstruct B,C,D,W through explicit nonzero pivots. If the three
  reconstructions agree, the full degree-eighteen residual is
  lambda[t(t-1)]^6; one sixth-order pole equation is therefore
  equivalent to the complete spectral identity. This is a lossless
  sparse coefficient reduction, not an H4 exclusion or a proof of
  JC(2) or DC(2).
source: codex-2026-07-26-h4-three-pole-jets
depends_on:
  - THM-2332-degree-eighteen-genus-zero-square-class-and-dessin-trap
  - THM-2386-degree-eighteen-h4-common-root-elimination
related:
  - THM-2373-degree-eighteen-rational-charged-section-atlas
  - THM-2386-degree-eighteen-h4-common-root-elimination
  - THM-2387-degree-eighteen-h4-elliptic-three-isogeny-atlas
script: 04-computation/jc2_degree18_h4_three_pole_jets_thm2389.py
output: 05-knowledge/results/jc2_degree18_h4_three_pole_jets_thm2389.out
script_sha256: fe462d935d743a81606eb7edfdca27bea5f9b0e50f8e9378b6667a6e81b04694
output_sha256: 225723bd7b8c6ad5e0f19832eb8e2958798654e4c4966f43cbd7f3e2273d4051
secondary_script: 04-computation/jc2_degree18_h4_three_pole_jet_synchronization_thm2389.py
secondary_output: 05-knowledge/results/jc2_degree18_h4_three_pole_jet_synchronization_thm2389.out
secondary_script_sha256: abc8a0f9ce1095230e25f24d1bbb1974a85868c87a2244f57914ed6d51a363ad
secondary_output_sha256: 821fc5e5c24512f6bb7edf79cbe9283d3326b007af63139b799fc0d71310766c
hash_basis: working-tree bytes (LF)
---

# THM-2389 -- synchronize the three infinity-pole jets

**PROVED CANDIDATE + VERIFIED-EXACT; INDEPENDENT AUDIT PENDING.**

THM-2332 reduces every remaining degree-eighteen Keller trajectory to
the four-simple square-class identity

```text
4P(y)^3+49Q(y)^2=H_4(y)S_4(y)^2,                 (1)

H_4 squarefree,
```

for the fixed structured quartic and sextic

```text
P
 =245y^4+1890 B y^2-24300 B^2+122472 D,

Q
 =539y^6+11340 B y^4+183708 C y^3
  +(72900 B^2-367416 D)y^2
  +(2361960 BC+2480058 W)y.                       (2)
```

The normalization is rational, the map to the `y`-line has degree
three, and the three points above infinity are smooth and unramified.
The raw coefficient identity in a rational parametrization has nineteen
coefficients. This theorem replaces it by three short local jets and
one scalar lock.

Throughout, use boldface names in prose for the spectral parameters
`(B,C,D,W)` when they might be confused with the denominator polynomial.

## 1. The three-pole normalization

Choose a coordinate `t` on the rational normalization that sends the
three points above `y=infinity` to

```text
t=0, 1, infinity.
```

Put

```text
B_2(t)=t(t-1).
```

The degree-three map has a simple pole at each of these points, so

```text
y=A_3(t)/B_2(t),                                  (3)

A_3=a_0+a_1t+a_2t^2+a_3t^3,

a_0 a_3 A_3(1) !=0.
```

The depressed spectral coordinate `v` has pole order exactly two at
each point above infinity. Hence

```text
v=C_6(t)/B_2(t)^2,                                (4)

deg C_6=6.
```

Define the regular infinity coordinates

```text
x=B_2/A_3=1/y,

r=C_6/A_3^2=v/y^2.                                (5)
```

With

```text
alpha=16/964467,
beta =64/703096443,
```

the depressed cubic of THM-2332 becomes

```text
Phi(r,x)
 =r^3
  +alpha[
     245+1890 B x^2+(-24300 B^2+122472 D)x^4
   ]r
  +beta[
     539+11340 B x^2+183708 C x^3
     +(72900 B^2-367416 D)x^4
     +(2361960 BC+2480058 W)x^5
   ]=0.                                           (6)
```

There is no term of order one or six in `x`. This sparse hole is the
source of both the Hermite reduction and its sharp last obstruction.

## 2. The fixed slope cubic

At `x=0`, equation (6) is

```text
f(s)
 =s^3+(80/19683)s+704/14348907=0.                 (7)
```

Its discriminant is

```text
Disc(f)=-94208/282429536481
       =-2^12*23/3^24 !=0.                        (8)
```

It is exactly the depressed version of THM-2332's infinity polynomial

```text
L(a)=1127-138915a+1607445a^2-26040609a^3:

f(s)=-L(s+5/243)/26040609.                        (9)
```

Thus the three values

```text
s_0   =C_6(0)/A_3(0)^2,

s_1   =C_6(1)/A_3(1)^2,

s_inf =[t^6]C_6/a_3^2                            (10)
```

are a permutation of the three distinct roots of `f`. They cannot
repeat: the three points over infinity are distinct smooth points of
the normalization, and (7) gives one point for each simple slope.

Since

```text
Phi_r(s,0)=f'(s)!=0,

Phi_x(s,0)=0,
```

the implicit expansion at every infinity point starts as

```text
r=s+O(x^2).                                       (11)
```

The absent linear term is load-bearing.

## 3. Six Hermite conditions leave one coefficient

Write

```text
C_6=sum_(j=0)^6 c_j t^j.
```

Equation (11) at `0`, `1`, and `infinity` is exactly

```text
c_0=s_0 a_0^2,

c_1=2s_0 a_0a_1,

c_6=s_inf a_3^2,

c_5=2s_inf a_3a_2,                                (12)

C_6(1)=s_1 A_3(1)^2,

C_6'(1)=2s_1 A_3(1)A_3'(1).
```

These are six independent affine conditions on the seven coefficients
of `C_6`. Here is an explicit solution. Put

```text
Sigma
 =s_1 A_3(1)^2-(c_0+c_1+c_5+c_6),

Tau
 =2s_1 A_3(1)A_3'(1)-(c_1+5c_5+6c_6).            (13)
```

Then, for one free scalar `kappa`,

```text
c_2=3Sigma-Tau+kappa,

c_3=Tau-2Sigma-2kappa,

c_4=kappa.                                        (14)
```

Conversely, (12)--(14) satisfy all six conditions.

The freedom has an intrinsic form. The difference of any two solutions
has double zeros at `0` and `1`, vanishing coefficients of `t^6,t^5`,
and therefore degree at most four. Hence

```text
C_6=C_*+kappa B_2^2.                              (15)
```

The kernel is exactly one-dimensional. It is not an arbitrary
coefficient-count heuristic.

## 4. Orders two through five reconstruct the spectral parameters

At each pole `i in {0,1,infinity}`, use `x=B_2/A_3` itself as local
parameter and write

```text
r
 =s_i+k_i x^2+l_i x^3+m_i x^4+n_i x^5+o_i x^6
  +O(x^7).                                        (16)
```

All six jet coefficients are rational functions of the five
three-pole variables

```text
a_0,a_1,a_2,a_3,kappa
```

and of the chosen slope ordering. No branch root in the `y`-line is
introduced.

For brevity write

```text
f_i'=3s_i^2+80/19683,

g_2(s)
 =1890alpha s+11340beta
 =160(243s+8)/1240029,

h_3=183708beta=256/15309,

g_4(s)
 =122472alpha s-367416beta
 =128(243s-4)/15309,

h_5=2480058beta=128/567.                          (17)
```

Every pivot in (17) is nonzero on a root of `f`:

```text
f_i' !=0                                          by (8),

g_2(s)=0  => s=-8/243 and f(s)=-64/531441,

g_4(s)=0  => s= 4/243 and f(s)= 64/531441,

h_3 h_5!=0.                                       (18)
```

Comparing the coefficients of `x^2,x^3,x^4,x^5` in (6) therefore
reconstructs, successively,

```text
B_i
 =-f_i' k_i/g_2(s_i),                             (19)

C_i
 =-f_i' l_i/h_3.                                  (20)
```

If the three values in (19) agree, call the common value `B`. If the
three values in (20) agree, call it `C`. Then

```text
D_i
 =-1/g_4(s_i) [
      f_i'm_i+3s_i k_i^2
      +alpha(1890Bk_i-24300B^2s_i)
      +beta(72900B^2)
   ],                                             (21)
```

and, after the three `D_i` agree with common value `D`,

```text
W_i
 =-1/h_5 [
      f_i'n_i+6s_i k_i l_i
      +1890alpha B l_i
      +2361960beta BC
   ].                                             (22)
```

Thus one full spectral tuple `(B,C,D,W)` exists through order five if
and only if

```text
B_0=B_1=B_inf,

C_0=C_1=C_inf,

D_0=D_1=D_inf,

W_0=W_1=W_inf.                                    (23)
```

This is eight explicit synchronization equations in the five
three-pole variables, one system for each of the six slope orderings.

## 5. One sixth-order equation is the full global identity

Assume (23), and use its common tuple in (6). The coefficient of `x^6`
at pole `i` is

```text
Xi_i
 =f_i'o_i
  +3s_i(2k_i m_i+l_i^2)+k_i^3
  +alpha[
     1890B m_i+(-24300B^2+122472D)k_i
   ].                                             (24)
```

There is no contribution from `Q` at this order because the structured
sextic (2) has zero constant term.

Now homogenize (6):

```text
N(t)
 =C_6^3
  +alpha P^h(A_3,B_2)C_6
  +beta Q^h(A_3,B_2),                             (25)

Phi=N/A_3^6.
```

The degree of `N` is at most eighteen. Equations (7), (11), and
(19)--(23) say that `Phi=O(x^6)` at all three poles.

At `t=0` and `t=1`, the numerator `A_3` is a unit and `x` is a
uniformizer. Therefore

```text
t^6(t-1)^6 divides N.                             (26)
```

At infinity, put `q=1/t`. Since

```text
A_3^6 is asymptotic to a_3^6 q^(-18),

x is asymptotic to q/a_3,

Phi=O(q^6),
```

we get

```text
deg N<=12.                                        (27)
```

The divisor in (26) already has degree twelve. Hence

```text
N=lambda B_2^6                                   (28)
```

for one scalar `lambda`. Dividing by `A_3^6` gives

```text
Phi=lambda x^6.
```

Consequently

```text
Xi_0=Xi_1=Xi_inf=lambda.                          (29)
```

It is enough to impose one last equation:

```text
Xi_0=0.                                           (30)
```

Equations (23) and (30) are **equivalent** to the complete homogenized
identity (25), after the slope and Hermite conditions. This proves the
lossless sparse reduction.

The final lock is sharp. If a constant term were allowed in `Q(y)`,
then its homogenization would add an arbitrary multiple of `B_2^6`.
All three jets through order five would remain unchanged. Thus no
criterion stopping at order five can remove `lambda`; (30) detects
exactly the coefficient hole special to (2).

## 6. Scope and non-overlap

This theorem changes the live `H_4` coefficient problem from

```text
nineteen global coefficients in A_3,C_6,B,C,D,W
```

to

```text
six slope orderings;

five variables a_0,a_1,a_2,a_3,kappa;

eight pole-synchronization equations;

one sixth-order scalar lock.                       (31)
```

It is compatible with the root-free parameter charts of THM-2373 but
does not spend or replace their weighted scaling action. The coordinate
`t` normalizes the three infinity points of the spectral curve, whereas
THM-2373 normalizes one adjacent pair among `(B,C,D,W)`.

THM-2386 now proves

```text
Res_y(P,Q)!=0
```

for every genuine `H_4` survivor. Thus a current elimination consumer
of (23)--(30) may saturate by this resultant from the outset. The
pole-jet proof itself never divides by `Res_y(P,Q)` and remains an exact
identity on the larger formal coefficient locus; it does not reprove
THM-2386 or empty the coprime lane.

THM-2387 targets the branch quartic, its Cardano divisor, and an
elliptic three-isogeny atlas. This theorem constructs none of those
objects: it retains the degree-six spectral numerator and synchronizes
its infinity jets. The two reductions are complementary.

The remaining `S_3` action permuting `0,1,infinity` merely permutes the
six slope orderings and the three copies of (19)--(24). No ordering is
canonically distinguished.

The theorem is an exact necessary-and-sufficient parametrization of the
spectral identity for a putative rational `H_4` normalization. It does
not prove that the nine equations in (31) are empty, restore the Keller
one-form or flux sidecars, close degree eighteen, or prove `JC(2)` or
`DC(2)`.

## 7. Coordinate-free binary-form restatement

The normalization `0,1,infinity` is computationally economical, but the
mechanism does not depend on it. Put

```text
w=(1701/2)v.
```

Then the depressed equation is

```text
w^3+12P(y)w+56Q(y)=0.                              (32)
```

On an arbitrary coordinate `(X:Z)` of the rational normalization write

```text
y=n/d,                   w=r/d^2,                  (33)
```

with coprime binary cubics `n,d`, squarefree `d`, and a binary sextic
`r`. The three roots of `d` are the three labelled infinity poles.
Their rescaled slopes

```text
alpha_i=(1701/2)s_i
```

are the distinct roots of

```text
alpha^3+2940alpha+30184=0,

Disc=-126247730112.                                (34)
```

Near a pole use the intrinsic germs

```text
x=d/n=1/y,                 a=r/n^2=w/y^2.
```

Equation (32), after multiplication by `x^6`, has no linear `x` term.
Thus

```text
a(0)=alpha_i,                 a'(0)=0.              (35)
```

For a fixed pole labelling, the six confluent conditions (35) put `r`
on one affine line. Indeed, the difference of two solutions has a
double zero at every root of squarefree `d`, hence is exactly

```text
r-r_tilde=tau d^2.                                 (36)
```

Likewise, after the order-two through order-five reconstructions agree,
the homogenized degree-eighteen residual has order at least six at all
three roots of `d`. It is therefore exactly

```text
lambda d^6.                                        (37)
```

The single sixth-order lock kills `lambda`. Equations (36)--(37) are
the coordinate-free content of the two kernel statements
`span(B_2^2)` and `span(B_2^6)` above. They also identify precisely what
the `S_3` quotient forgets: the bijection between the three poles and
the three roots of (34), not any intrinsic tournament orientation.

## 8. Exact companions

The dependency-free exact companion:

- verifies (7)--(9) and the two hostile pivot values in (18);
- checks the explicit Hermite family (12)--(15) and proves its
  homogeneous constraint matrix has rank six with kernel `B_2^2`;
- expands (6) through order six and verifies every reconstruction
  coefficient (19)--(24) on exact rational controls;
- proves by exact row reduction that the eighteen local vanishing
  conditions on a degree-at-most-eighteen numerator have rank eighteen
  and kernel exactly `B_2^6`; and
- checks a split finite-field control at the first good prime `59`,
  where the three slope roots are `5,11,43`.

Run

```bash
python3 04-computation/jc2_degree18_h4_three_pole_jets_thm2389.py
python3 -O 04-computation/jc2_degree18_h4_three_pole_jets_thm2389.py
```

Both transcripts must match

```text
05-knowledge/results/jc2_degree18_h4_three_pole_jets_thm2389.out
```

byte-for-byte after LF normalization. Every executable check raises
explicitly under optimized Python. The declared hashes have been
reproduced; independent audit remains pending before status promotion.

The secondary coordinate-free companion independently verifies
(32)--(37), derives the triangular reconstruction formulas after the
integer rescaling, checks `span(d^2)` and `span(d^6)` on exact confluent
packets, and reconstructs formal branches in
`Q[alpha]/(alpha^3+2940alpha+30184)`. Run

```bash
python3 04-computation/jc2_degree18_h4_three_pole_jet_synchronization_thm2389.py
python3 -O 04-computation/jc2_degree18_h4_three_pole_jet_synchronization_thm2389.py
```

Both transcripts must byte-match

```text
05-knowledge/results/jc2_degree18_h4_three_pole_jet_synchronization_thm2389.out
```

after LF normalization.
