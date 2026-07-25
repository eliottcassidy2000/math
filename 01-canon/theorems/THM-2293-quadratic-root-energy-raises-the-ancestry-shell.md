---
id: THM-2293
title: "Quadratic root energy raises the ancestry shell"
status: >
  PROVED + VERIFIED-EXACT. On every one of the 150 strict first-depth-one
  scalar profiles, summed quadratic proper-root energy retains one shallow
  ancestry label and has a signed whole-class covariance with the named
  deepest-successor comb. It therefore lands an exact terminal Fourier mode
  m*c_3 with 0<|m|<=2533 and gcd(m,91)=1. A Perron-weighted refinement
  retains actual Fourier atoms A,A' of one original exclusive-owner set:
  both have exact thirteen-adic valuation b, have the same nonzero root
  residue after division by 13^b, and satisfy A-A'=m*c_3 for a
  gcd(m,91)=1 multiplier with 0<|m|<=578982. On the 120 interior rows the
  small terminal mode either annotates THM-2284's existing c_1/c_3 pivot or
  extends its anchored relation flag by one unit augmented row. Positive
  Jackson smoothing gives the actual source pair a row-dependent finite word
  address, without a uniform height. Separately, at explicit uniform
  bandwidth it produces a possibly different finite-polynomial pair on the
  same shallow label and a bounded equality of height at most
  11685923626556462293960132370114. That equality may be the tautological
  c_3-only channel, so it does not yet force a nonzero new relation. On all
  120 interior rows the fixed root character selected by THM-2276's c_1
  carry also lands some c_3 terminal difference, but with no uniform
  multiplier bound from the current label data. No profile is excluded and
  LRC(14) remains open.
source: codex-2026-07-25-quadratic-root-shell-lift
depends_on:
  - THM-2273-shallow-owner-flow-and-deep-successor-gap-spread
  - THM-2276-shallow-owner-residue-aligned-crossing
  - THM-2278-two-shallow-proper-root-spectrum-and-gap-ancestry-activation
  - THM-2284-thirteen-adic-anchored-rank-three-plucker-lift
related:
  - THM-2201-cyclic-root-fibre-hasse-jet-transition-carrier
  - THM-2269-marked-expiration-root-spectrum-and-branch-state-no-go
  - THM-2282-thirteen-adic-saturation-and-unit-anchored-minor
  - THM-2283-mixed-rank-two-safe-torus-floor-and-scalar-rank-three-harvest
script: 04-computation/lrc14_quadratic_root_energy_shell_lift_thm2293.py
output: 05-knowledge/results/lrc14_quadratic_root_energy_shell_lift_thm2293.out
script_sha256: d5aae265938fa4c0bdae522a202861ece81ba4297d50cea6dcacd062f8df88d5
output_sha256: e83b01be617ad9e6b442c26ea1121781b92eb8d6bfe647ebca510639e0e9304e
hash_basis: working-tree bytes (LF)
---

# THM-2293 -- quadratic root energy raises the ancestry shell

**PROVED + VERIFIED-EXACT.**

The nontrivial root characters in THM-2278 retain shallow ancestry, but a
single such character lives one thirteen-adic shell below any ordinary
post-expiration target. The lawful repair is quadratic:

```text
same nonzero root character times its conjugate
  -> the root character cancels
  -> two lifts in the same residue are differenced
  -> the difference rises by one thirteen-adic shell.          (1)
```

The deepest-successor safe-gap condition then forces this quadratic energy
to have a signed Fourier interaction with the deepest-successor danger comb.
The interaction can be landed at a bounded multiplier. There are two useful
forms:

```text
unweighted image masks:
  a small terminal mode |m|<=2533;

Perron-weighted source masks:
  two actual original-source Fourier atoms A,A'
  with A-A'=m c_3 and |m|<=578982.                    (2)
```

The first has the smaller coefficient. The second retains the actual integer
lifts of an original exclusive-owner set.

## 1. Strict-profile setup

Use THM-2273 and THM-2278's strict scalar notation

```text
T(x)=13x mod 1,

c_1=13u_1,       c_2=13^b u_2,       c_3=13^c u_3,

2<=b<c,                         5<=c<=19,             (3)
```

where `u_1,u_2,u_3` are thirteen-units. Let `E_1,E_2` be the two
shallow exclusive-owner pieces and put

```text
G_j=T^b(E_j),             Y_j=T(G_j),       Y=Y_1 union Y_2.
                                                               (4)
```

The middle-owner separation and deepest-successor exclusion are

```text
G_1 subset D_(u_2)^c,     G_2 subset D_(u_2),

Y subset D_s^c,

s=c_3/13^(b+1)=13^d u_3,       d=c-b-1>=0.           (5)
```

THM-2273 gives

```text
measure(Y)>=Y_0:=5696989/76962600.                   (6)
```

For `y in R/Z`, label its thirteen roots by

```text
x_r(y)=(y+r)/13,                  r in F_13,          (7)
```

and use THM-2278's binary masks and root transforms

```text
m_(j,y)(r)=1_(G_j)(x_r(y)),

n_j(y)=sum_r m_(j,y)(r),

M_(j,kappa)(y)
 =sum_r m_(j,y)(r) zeta^(-kappa r),

zeta=exp(2*pi*i/13),       1<=kappa<=12.             (8)
```

On its active fibres the first mask is nonempty and proper, while the
second has occupancy one or two:

```text
n_1 in {1,...,12},         n_2 in {1,2}.             (9)
```

## 2. A conditional-covariance shell lemma

The operation used twice below is elementary and exact.

Let

```text
s=13^d u,         13 does not divide u,
Q=13^(d+1),                                           (10)
```

and define the translation average

```text
A_Q F(x)=(1/Q)sum_(r=0)^(Q-1)F(x+r/Q).               (11)
```

It is the orthogonal Fourier projection onto the modes divisible by `Q`.
As `r` runs modulo `Q`,

```text
s(x+r/Q)=sx+ur/13 mod 1.                             (12)
```

Every residue modulo thirteen occurs `13^d` times. A centered arc of
length `1/7`, which lies strictly between `1/13` and `2/13`, contains one
or two points of every translated thirteen-grid. Consequently

```text
A_Q 1_(D_s)(x) in {1/13,2/13}                       (13)
```

away from a null endpoint set.

Let `W` be a nonnegative `L^2` function satisfying

```text
W 1_(D_s)=0.                                         (14)
```

Orthogonal projection and (13) give

```text
sum_(Q does not divide n)
  W_hat(n) overline((1_(D_s))_hat(n))

 =-<A_Q W,A_Q 1_(D_s)>
 <=-(1/13)integral W.                                (15)
```

The Fourier transform of the comb is

```text
(1_(D_s))_hat(n)=0,                     s does not divide n,

(1_(D_s))_hat(ms)
 =sin(pi*m/7)/(pi*m),                  m!=0.         (16)
```

Moreover,

```text
Q does not divide ms       iff       13 does not divide m. (17)
```

Thus (15) is a negative signed aggregate over exactly the multipliers

```text
m!=0,       7 does not divide m,       13 does not divide m. (18)
```

This is a whole-class statement. No individual Fourier coefficient has yet
been selected.

We also need a quantitative tail. If

```text
0<=W<=C,                 mu=integral W,              (19)
```

then Bessel and `W^2<=CW` give

```text
sum_m |W_hat(ms)|^2<=C mu.                           (20)
```

The elementary inscribed-dodecagon estimate gives

```text
pi>31/10.                                             (21)
```

For an exact rational verification, `sqrt(3)<1733/1000`, and hence

```text
36(2-sqrt(3))>36(2-1733/1000)
              =9612/1000
              >961/100.
```

The unit-circle dodecagon has half-perimeter
`6sqrt(2-sqrt(3))`, proving (21). Therefore

```text
sum_(|m|>=L+1)|(1_(D_s))_hat(ms)|^2
 <200/(961L).                                        (22)
```

Equations (20)--(22), followed by Cauchy--Schwarz, make every multiplier
bound below coefficient-independent.

## 3. Small terminal mode from the unweighted root word

Sum the twelve nontrivial character energies separately on the two shallow
labels:

```text
W_j(y)=sum_(kappa=1)^12 |M_(j,kappa)(y)|^2.          (23)
```

Finite Parseval on `F_13` gives the exact occupancy formula

```text
W_j=13n_j-n_j^2.                                     (24)
```

It follows from (9) that

```text
0<=W_1<=42,             0<=W_2<=22.                 (25)
```

At every point of `Y`, at least one active proper mask contributes at least
twelve. Hence, writing `mu_j=integral W_j`,

```text
mu_1+mu_2>=12 measure(Y)>=12Y_0.                    (26)
```

Choose `j in {1,2}` so that

```text
mu_j/C_j
 >=(mu_1+mu_2)/(42+22)
 >=12Y_0/64
 =5696989/410467200,                                (27)

C_1=42,                    C_2=22.
```

This weighted choice, rather than an unweighted pigeonhole, retains the
shallow label without paying the worse pointwise cap on both labels.

By (5), `W_j 1_(D_s)=0`. Apply (15). Suppose every term in (18) with

```text
0<|m|<=2533
```

vanished. The next two absolute multipliers cannot contribute:

```text
7 divides 2534,              13 divides 2535.        (28)
```

Thus the whole aggregate in (15) would lie in `|m|>=2536`. Equations
(20)--(22), with `L=2535`, would give an absolute value strictly smaller
than `mu_j/13`, because the exact rational comparison is

```text
mu_j/C_j
 >=5696989/410467200

 >33800/(961*2535)
  =40/2883,                                           (29)

5696989/410467200-40/2883
 =1910429/394458979200
 >0.
```

This contradicts (15). Therefore there are

```text
j in {1,2},        kappa in {1,...,12},

m in Z,            0<|m|<=2533,       gcd(m,91)=1   (30)
```

such that

```text
Fourier(W_j,ms)!=0,

Fourier(|M_(j,kappa)|^2,ms)!=0.                      (31)
```

The second assertion follows because `W_j` is the sum over `kappa`.
Pulling (31) back through `T^(b+1)` sends the frequency `ms` to

```text
13^(b+1)ms=m c_3.                                   (32)
```

Thus `m c_3` is an exact Fourier atom of a quadratic root-energy function
that retains one shallow ancestry label and one nonzero root character.
The whole signed covariance (15), not merely the existence in (31), is the
certificate.

### 3A. A fixed-character variant

The summation over characters in (23) is needed only for a uniform
multiplier bound. It is not needed for existence.

Suppose one specified label `j` has positive mass, and fix any specified

```text
kappa in {1,...,12}.                                (F1)
```

The proper-root lemma in THM-2278 applies on every active fibre, so

```text
W_(j,kappa):=|M_(j,kappa)|^2

is positive on Y_j,        integral W_(j,kappa)>0.   (F2)
```

Moreover its support is contained in `Y_j subset D_s^c`. Applying (15)
directly to `W_(j,kappa)` gives a strictly negative whole-class sum.
Therefore some, not a priori bounded,

```text
m in Z,        gcd(m,91)=1                           (F3)
```

satisfies

```text
Fourier(|M_(j,kappa)|^2,ms)!=0.                      (F4)
```

Pullback through `T^(b+1)` again makes this a terminal quadratic mode
`m c_3`. A quantitative bound for `m` would require a quantitative lower
bound for this particular label's mass; the available theorem gives only a
floor for the sum of the two labels.

On every one of the `120` interior profiles, THM-2276 proves that the
`c_1`-private label `E_1` is open and nonempty, hence has positive mass. It
also selects a primitive two-coordinate crossing with normalized carry

```text
kappa_carry=A(p)/13 mod 13 in F_13^*.                (F5)
```

Taking `j=1` and `kappa=kappa_carry` in (F1)--(F4) glues that prescribed
normalized carry character to some terminal `c_3` difference. This is
character-value alignment on the same `c_1`-private label. It does not
reconstruct the discarded intermediate root digits, and it has no uniform
multiplier height until the `c_1`-private mass itself has a uniform floor.

## 4. Perron weighting recovers two actual source atoms

The preceding mode belongs to the quadratic image carrier. To recover
actual Fourier atoms of an original exclusive-owner set, keep the Perron
weights instead of replacing positive image density by a binary support.

Put

```text
f_j=1_(E_j),             e_j=measure(E_j),

h_j=P^b f_j,                                             (33)
```

where

```text
P f(y)=(1/13)sum_(r=0)^12 f((y+r)/13).                (34)
```

Then `0<=h_j<=1`, `integral h_j=e_j`, and

```text
support(h_1) subset G_1 subset D_(u_2)^c,

support(h_2) subset G_2 subset D_(u_2).              (35)
```

Define the weighted root transforms, their zero-character sums, and their
total nonconstant energies by

```text
N_(j,kappa)(y)=sum_r h_j(x_r(y))zeta^(-kappa r),

S_j(y)=sum_r h_j(x_r(y)),

V_j(y)=sum_(kappa=1)^12 |N_(j,kappa)(y)|^2.          (36)
```

Finite Parseval gives

```text
V_j
 =13sum_r h_j(x_r)^2-S_j^2.                          (37)
```

The root support sizes in (35) are at most twelve and two, respectively.
Cauchy on each fibre therefore gives

```text
V_1>=S_1^2/12,

V_2>=(11/2)S_2^2.                                   (38)
```

Both `S_j` are supported on `Y_j subset D_s^c`, a set of measure at most
`6/7`, and

```text
integral S_j=13e_j.                                  (39)
```

A second Cauchy inequality yields

```text
integral V_1>=(1183/72)e_1^2,

integral V_2>=(13013/12)e_2^2.                       (40)
```

The quadratic in (37) is convex in each root weight on `[0,1]`. Its maximum
is therefore attained at a binary vertex. The same occupancy table as
(24) gives

```text
0<=V_1<=42,                 0<=V_2<=22.              (41)
```

THM-2273 supplies

```text
e_1+e_2>=L_0:=5696989/367580070.                    (42)
```

The sharp quadratic allocation between the coefficients in (40) is

```text
(1183/72)e_1^2+(13013/12)e_2^2

 >=(13013/804)(e_1+e_2)^2.
```

Consequently, if `nu_j=integral V_j`,

```text
nu_1+nu_2
 >=nu_0
 :=(13013/804)L_0^2
 =227189785662847/58436012221844400.                 (43)
```

Choose a label with

```text
nu_j/C_j>=nu_0/64.                                   (44)
```

Again `V_j1_(D_s)=0`, so (15) applies. If all the unit multipliers through
`578982` vanished, the tail estimate with `L=578982` would contradict
(15). The exact comparison is

```text
nu_0/64
 >33800/(961*578982),                                (45)

nu_0/64-33800/(961*578982)

 =296921577147599
  /346814897688821607884467200
 >0.
```

Therefore there are

```text
j in {1,2},        kappa in {1,...,12},

m in Z,            0<|m|<=578982,      gcd(m,91)=1  (46)
```

such that

```text
Fourier(|N_(j,kappa)|^2,ms)!=0.                      (47)
```

This coefficient is a genuine autocorrelation of original-source Fourier
atoms. Indeed, for every `L^2` function `h`,

```text
N_kappa(y)
 =13 exp(2*pi*i*kappa*y/13)
   sum_(r in Z) h_hat(kappa+13r)exp(2*pi*i*r*y).     (48)
```

Thus a nonzero coefficient at `n=ms` in (47) contains at least one nonzero
summand

```text
h_j_hat(kappa+13(r+n))
  overline(h_j_hat(kappa+13r)).                     (49)
```

The Perron identity is

```text
h_j_hat(q)=f_j_hat(13^b q).                         (50)
```

Put

```text
A =13^b(kappa+13(r+n)),

A'=13^b(kappa+13r).                                 (51)
```

Equations (47)--(51) prove

```text
f_j_hat(A)!=0,                 f_j_hat(A')!=0,

nu_13(A)=nu_13(A')=b,

A/13^b congruent A'/13^b congruent kappa mod 13,

A-A'
 =13^(b+1)n
 =13^(b+1)ms
 =m c_3.                                             (52)
```

This is the promised actual-lift toothpick pair. The two frequencies agree
at the common shallow root residue and then remain congruent through the
entire stalk from depth `b` to depth `c`; their first terminal split is the
bounded thirteen-unit multiple `m u_3`.

## 5. Finite Jackson realization and its exact stopping boundary

There are two distinct finite-word conclusions. First, each particular
actual source pair from Section 4 acquires a finite coordinate address by
coefficientwise Jackson convergence, but the bandwidth depends on the row
and the selected atoms. Second, one explicit uniform bandwidth preserves
the nonzero covariance head and produces a bounded pair in the smoothed
polynomial carrier. The latter pair need not be the former source pair.

Write the scalar row as

```text
w_*=(H,q_1,...,q_5,c_1,c_2,c_3) in Z^9.             (J1)
```

Up to null endpoints, each `f_j=1_(E_j)` is a product of nine factors

```text
f_j(t)=product_(i=1)^9 chi_(j,i)(w_i t),             (J2)
```

where every `chi_(j,i)` is the indicator of one centered interval or its
complement.  The guard coordinate uses radius `1/7`; every other coordinate
uses radius `1/14`.

For `N>=2`, use the normalized squared-Fejer kernel

```text
P_N(t)=sum_(k=0)^(N-1) exp(2*pi*i*k*t),

J_N(t)=|P_N(t)|^4 / integral |P_N|^4.                (J3)
```

It is nonnegative, has integral one, and has degree

```text
H_N=2N-2.                                             (J4)
```

For each base factor put `q_(j,i)=J_N*chi_(j,i)` and define

```text
Q_(j,N)(t)=product_(i=1)^9 q_(j,i)(w_i t).           (J5)
```

The elementary Jackson first-moment estimate gives

```text
0<=q_(j,i)<=1,

||q_(j,i)-chi_(j,i)||_1<3/(2N).                      (J6)
```

For completeness, `integral |P_N|^4=N(2N^2+1)/3`; the bounds

```text
|P_N(t)|<=min(N,1/(2||t||))
```

and a split at `||t||=1/(2N)` give

```text
integral ||t||J_N(t)<3/(4N).
```

Translation changes an interval indicator by at most `2||t||` in `L^1`,
which proves (J6). The product telescope and the fact that both sides of
(J5) lie in `[0,1]` then give

```text
||Q_(j,N)-f_j||_2^2
 <=||Q_(j,N)-f_j||_1
 <27/(2N).                                            (J7)
```

Every nonzero Fourier coefficient of a base interval of length `1/7` or
`2/7` vanishes at multiples of seven. The squared-Fejer multiplier is
strictly positive throughout its support. Consequently every Fourier
frequency of `Q_(j,N)` has a word representation

```text
A=a.w_*,

a in Xi_(H_N)^9,

Xi_H={0} union {k in Z:0<|k|<=H and 7 does not divide k}. (J8)
```

Fix the actual atoms `A,A'` selected in Section 4. Equation (J7) implies
coefficientwise convergence

```text
Fourier(Q_(j,N),q) -> f_j_hat(q)       for every fixed q. (J8a)
```

Thus, for all sufficiently large row- and atom-dependent `N`, both
`Fourier(Q_(j,N),A)` and `Fourier(Q_(j,N),A')` are nonzero. Equation (J8)
then gives finite words `a,a'` with

```text
a.w_*=A,          a'.w_*=A',          (a-a').w_*=m c_3. (J8b)
```

This gives the actual Section 4 edge a finite coordinate address. It gives
no uniform bandwidth or height, because neither the absolute sizes of
`A,A'` nor their nonzero coefficient magnitudes have a uniform lower bound.

We now prove the separate uniform statement. Put

```text
h_(j,N)=P^b Q_(j,N)                                  (J9)
```

and form `N_(j,kappa,N)` and `V_(j,N)` from (36), replacing `h_j` by
`h_(j,N)`. The Perron operator is an `L^2` contraction. For any thirteen
root weights in `[0,1]`, their total nontrivial character energy is at most
`42`. Finite Parseval and Cauchy--Schwarz therefore give the pointwise and
integrated comparison

```text
|V_(j,N)-V_j|
 <=2 sqrt(42)
   sqrt(sum_(kappa=1)^12
        |N_(j,kappa,N)-N_(j,kappa)|^2),

||V_(j,N)-V_j||_2
 <=26 sqrt(42)||h_(j,N)-h_j||_2

 <78 sqrt(63/N).                                      (J10)
```

The factor `26` is exact in this ledger: root Parseval contributes `13`,
and integration over the thirteen inverse branches contributes a second
factor `13`.

We next freeze a quantitative nonzero head before smoothing. Let

```text
L=578982,

T_L=200/(961L)=100/278200851,

r_0=nu_0/64
    =227189785662847/3739904782198041600.             (J11)
```

Choose the same label `j` as in (44), put `C=C_j`, `r=nu_j/C`, and let
`H_j` denote its covariance sum over the unit multipliers
`0<|m|<=L`. Equations (15), (20), and (22) imply

```text
|H_j|
 >C(r/13-sqrt(rT_L))

 >=C(r-169T_L)/26

 >=g_0
 :=(22/26)(r_0-169T_L)

 =296921577147599
  /409872151814061900227097600
 >0.                                                   (J12)
```

The middle inequality uses
`sqrt(r)+13sqrt(T_L)<=2sqrt(r)`, valid because
`r>=r_0>169T_L`; the last uses `C in {22,42}`.

Take the certificate-minimal integer for this coarse ledger,

```text
N_0=2921480906639115573490032947784,

H_0=2N_0-2
    =5842961813278231146980065895566.                 (J13)
```

Exact rational arithmetic gives

```text
g_0^2/4-78^2*63/N_0

 =1512520856980458045994136291
  /81799118868347936026371657539647828685067698505281323958615849827440899472752640000
 >0.                                                   (J14)
```

Hence (J10) is less than `g_0/2`. Comparing the finite covariance heads by
Parseval gives, explicitly,

```text
|H_(j,N_0)-H_j|
 <=||V_(j,N_0)-V_j||_2
   (sum_(0<|m|<=L,gcd(m,91)=1)
      |(1_(D_s))_hat(ms)|^2)^(1/2)

 <=||V_(j,N_0)-V_j||_2 ||1_(D_s)||_2

 <g_0/2.                                                (J15)
```

Here Parseval gives `||1_(D_s)||_2=1/sqrt(7)<1`. Thus the
`N_0`-smoothed head is still nonzero. Therefore, for the same shallow
label—but not necessarily the same character or multiplier as the actual
pair in Section 4—some

```text
kappa_N in {1,...,12},

m in Z,       0<|m|<=578982,       gcd(m,91)=1       (J16)
```

satisfies

```text
Fourier(|N_(j,kappa_N,N_0)|^2,ms)!=0.                (J17)
```

Expand (J17) as in (48)--(51), now using the finite polynomial
`Q_(j,N_0)`. It supplies two nonzero Fourier atoms `B,B'` of that polynomial
with

```text
nu_13(B)=nu_13(B')=b,

B/13^b congruent B'/13^b congruent kappa_N mod 13,

B-B'=m c_3.                                          (J18)
```

Choose nonzero convolution summands for those atoms. By (J8), there are

```text
a,a' in Xi_(H_0)^9
```

such that

```text
a.w_*=B,          a'.w_*=B',          (a-a').w_*=m c_3. (J19)
```

Equivalently,

```text
r:=a-a'-m e_(c_3)
```

lies in the scalar relation lattice and obeys

```text
||a-a'||_infinity
 <=2H_0
 =11685923626556462293960131791132,

||r||_infinity
 <=2H_0+578982
 =11685923626556462293960132370114.                  (J20)
```

This is a genuine uniform bounded-word realization in the smoothed shell
carrier. It is not an identification with the actual source pair from
Section 4, and it is not yet a genuine new-relation theorem. The first
failed implication is exact:

```text
a-a'=m e_(c_3)                                       (J21)
```

is compatible with every conclusion above. In that channel `r=0`; the two
words differ only in the already named deepest coordinate. Even when
`r!=0`, the coefficient of `c_3` in `r` may vanish. The same-root residue
does not rule out (J21), because `c_3` is already deeper than the residue
level `b`.

Thus the next useful operation is a connected or cumulant-style projection
which annihilates the one-coordinate `c_3` autocorrelation while retaining
the negative ancestry/deepest covariance. Without that extra sidecar,
finite smoothing proves the bounded equality (J19), not a nonzero
cross-coordinate relation.

There is an exact binary-relation reframe. For one label and residue, take
as vertices the nonzero source Fourier atoms

```text
Gamma_(j,kappa)={
 A:
 f_j_hat(A)!=0,
 nu_13(A)=b,
 A/13^b congruent kappa mod 13
}.                                                     (J22)
```

Join two vertices when their difference is `m c_3` with `gcd(m,91)=1`.
Section 4 proves that one selected graph has an edge. The qualitative part
of Section 5 gives that exact edge a row-dependent finite word address.

For the uniform statement, define instead

```text
Gamma^(N)_(j,kappa)={
 B:
 Fourier(Q_(j,N),B)!=0,
 nu_13(B)=b,
 B/13^b congruent kappa mod 13
}.                                                    (J23)
```

Section 5 proves that some `Gamma^(N_0)_(j,kappa_N)` on the same label has
an explicitly bounded edge. It may be a different residue and edge from
the actual graph in Section 4. Neither statement points the actual edge at
any previously marked atom. In particular, a later theorem which merely
produces one bounded common atom supplies a marked vertex, not an incident
shell edge. The exact composition target is therefore a marked-degree or
connectivity lemma for the actual `Gamma_(j,kappa)`.

This carrier is an undirected graph, not naturally a tournament: reversing
an edge replaces `m` by `-m`, and no target predicate distinguishes the two
orientations. A legitimate orientation would need an extra phase, time, or
owner-transfer sidecar.

## 6. The anchored Plucker interface on the 120 interior rows

Now assume the interior scope

```text
3<=b<=c-2.                                           (53)
```

THM-2284 supplies three scalar relation rows

```text
p,r,s_* in Z^9
```

and a three-column set

```text
J={c_1,k,l}                                          (54)
```

whose minor is a thirteen-unit and satisfies

```text
0<|Delta_3|<=54991358114.                            (55)
```

Let

```text
a_0=m e_(c_3)                                        (56)
```

be the coefficient row annotated by the small terminal atom in Section 3.
It is not an integer relation:

```text
a_0.w_*=m c_3!=0.                                   (57)
```

If `c_3 in J`, the existing relation pivot is already anchored at both
`c_1` and the analytically selected terminal column `c_3`.

If `c_3 notin J`, append `a_0` to the three relation rows and append `c_3`
to the three pivot columns. The resulting four-by-four determinant is

```text
Delta_4=+/-m Delta_3,

Delta_4!=0 mod 13,

0<|Delta_4|
 <=2533*54991358114
 =139293110102762.                                   (58)
```

Thus every interior strict row has an anchored augmented flag joining
`c_1` to the root-energy-selected deepest column. This is an augmented
relation/spectral flag, not four-dimensional integer relation rank.

## 7. Filtered-module and Bockstein-shaped interpretation

For the scalar coefficient row `w_*`, define

```text
partial(a)=a.w_*,

F^t={a in Z^9:13^t divides partial(a)}.              (59)
```

Then

```text
Z^9=F^0 superset F^1 superset F^2 superset ...,

Lambda_*=ker(partial)=intersection_(t>=0)F^t.        (60)
```

The THM-2284 relation rows are cycles in every `F^t`. The terminal analytic
row in (56) has the exact grade

```text
a_0 in F^c minus F^(c+1),                            (61)
```

because `m` and `u_3` are thirteen-units.

At the frequency level, a nontrivial root atom at the common shallow stage
has exact valuation `b`. Two atoms with the same root character have the
form

```text
13^b(kappa+13r),       13^b(kappa+13r').
```

Their difference lies one step deeper, in `13^(b+1)Z`. Equations (46)--(52)
show that the deepest-comb covariance continues this lift all the way to
the exact grade `c` and lands it on `m c_3`.

This is **Bockstein-shaped**, not a claim that a new cohomology theory has
been constructed:

```text
linear nontrivial character:
  remembers ancestry but lies in the shallow graded piece;

character times conjugate:
  kills the character, differences two same-residue lifts,
  and enters the next filtration step;

deepest-comb covariance:
  selects the exact terminal graded piece.           (62)
```

The relation lattice acts as a gauge on coefficient representatives without
changing `partial(a)`. Therefore a Plucker address can normalize a
representative, but it cannot raise the scalar frequency to another shell.
The quadratic operation is load-bearing.

## 8. Exact loss ledger and non-consequences

The maps proved above are:

```text
source:
  THM-2278's two labelled proper-root words on the common-time
  deepest-gap carrier, or the Perron-weighted indicators of E_1,E_2;

target:
  a signed whole-class covariance with D_s and a terminal mode m c_3;

map:
  sum character squares, project away the 13^(d+1)-periodic sector,
  and use the deepest comb's exact Fourier support;

preserved:
  one shallow label, one nonzero root character, common time b+1,
  exact depths b and c, the deepest blocker label, and a bounded
  multiplier prime to both seven and thirteen;

additionally preserved in Section 4:
  two actual integer Fourier lifts of one original E_j;

additionally preserved after Section 5:
  a row-dependent finite word address for the exact Section 4 source edge,
  and a separate uniform pair of bounded nine-coordinate words in the
  smoothed carrier, with every nonzero word digit prime to seven;

destroyed:
  the complex phase under squaring, the individual deepest-gap index,
  a pointwise first-switch time, the absolute sizes of the unsmoothed
  atoms A,A', and hard support under finite smoothing;

still unresolved after finite smoothing:
  whether the word difference uses any coordinate other than c_3, and
  whether its relation-lattice remainder is nonzero, and whether the
  unpointed actual shell edge is incident to any separately marked atom;

needed next sidecar:
  a connected/cumulant channel killing the c_3-only autocorrelation, or
  a marked-degree/atom-graph connectivity theorem, or a proof that the
  terminal relative class cannot coexist with the anchored relation flag
  in the same integral filtration.                                     (63)
```

The linear stopping boundary is exact. A source Fourier mode that can pair
with a target pulled back through `T^k` must be divisible by `13^k`.
Nontrivial root characters at depth `b` have exact valuation `b`, so they
are orthogonal to every ordinary target at `k>=b+1`. Relation shifts do not
change the scalar frequency. Consequently THM-2269/2278's linear root
energy cannot be glued directly to THM-2291/2288/2289-style delayed
handoffs; the quadratic shell lift or a different signed whole-class
operation is necessary.

Neither (31), (52), nor (58) gives positive safe mass, contradicts the
cover, forces a handoff at expiration, or excludes a valuation profile.
The theorem applies to the `150` strict profiles only. The `15`
repeated-first profiles are outside the two-shallow input. LRC(14) remains
open.

## 9. Exact reproduction

Run

```bash
python3 04-computation/lrc14_quadratic_root_energy_shell_lift_thm2293.py
python3 -O 04-computation/lrc14_quadratic_root_energy_shell_lift_thm2293.py
```

The companion checks the full `150=120+15+15` strict-profile split, every
occupancy inequality, the translated-root count, the dodecagon lower bound
for `pi`, both rational energy floors, both multiplier-tail margins, the
excluded `2534/2535` boundary, all valuation identities, the augmented
determinant bound, and the finite-smoothing head gap, bandwidth, word height,
and first tautological stopping boundary. Normal and optimized runs produce
the stored transcript byte for byte. QED.
