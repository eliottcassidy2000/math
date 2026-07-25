---
id: THM-2276
title: "Shallow pair crossing, ancestry phase cone, and the 696-multiplier hard bank"
status: >
  RESERVED / COMPLETE PROOF CANDIDATE UNDER AUDIT. Candidate: every one of
  the 120 interior first-depth-one scalar profiles has a primitive
  two-coordinate guard/blocker or unit/blocker crossing of height at most
  9841. Its carry is exactly plus or minus m*c_1, where 1<=m<=757 and m
  is a thirteen-unit. The strict c_1-private locus is nonempty and open.
  Its multiplicity pushforward has the exact Fourier law
  g_hat(h)=13*(1_E)_hat(13h), while its image support has marked energy in
  the normalized carry's nonzero residue class. For m=1,2,3, the owner
  phase cone forces nonzero exact atoms in the original ancestry, descended
  multiplicity, and image support. The boundary is sharp: each of the 696
  possible thirteen-unit multipliers 4<=m<=757 has an exact two-interval
  local c_1-private carrier with the same pair relation and every-residue
  root energy but zero ancestry and descended atoms at the relation
  frequency. The local witness does not assert a global scalar cover. No
  profile is excluded, and this file is not a proved dependency until
  independently audited.
source: codex-2026-07-25-shallow-owner-phase-cone
depends_on:
  - THM-2266-depth-one-deep-pair-centered-signed-dual-and-relation-atlas
  - THM-2268-two-shell-private-owner-trident-and-raw-carry-cocycle-no-go
  - THM-2269-marked-expiration-root-spectrum-and-branch-state-no-go
related:
  - THM-2255-valuation-separated-pair-cap-and-exclusive-owner-mass
  - THM-2263-thirteen-adic-gap-pair-spectrum-and-profile-sharp-owner-floor
  - THM-2274-mixed-scalar-relative-rank-harvest-and-adaptive-pair-crossing
  - THM-2275-mixed-scalar-relation-and-guard-blocker-crossing
  - THM-2273-shallow-owner-flow-and-deep-successor-gap-spread
  - THM-2279-shallow-blocker-anchored-relation-minor
  - THM-2282-thirteen-adic-saturation-and-unit-anchored-minor
script: 04-computation/lrc14_shallow_owner_phase_cone_thm2276.py
output: 05-knowledge/results/lrc14_shallow_owner_phase_cone_thm2276.out
script_sha256: 751e24f7a159a7b87c2dbc1ee96b44863eeee2784ad657fd703924c652e4eedd
output_sha256: eec21897d98afebc02076063b3646a1f050598417005eb7ec8fa8f627c7990af
hash_basis: working-tree bytes (LF)
---

# THM-2276 -- shallow pair crossing and the ancestry phase-cone boundary

**RESERVED / COMPLETE PROOF CANDIDATE UNDER AUDIT.**

Use

```text
D_a={x in R/Z:||ax||<1/14},
C_H={x in R/Z:||Hx||>1/7},
T(x)=13x mod 1.                                      (1)
```

Assume a live scalar row in the counterexample branch of THM-2198 and
THM-2266. In particular,

```text
H,q_1,...,q_5 are positive thirteen-units;
H is odd;
q_1,...,q_5 are pairwise distinct;
c_1,c_2,c_3 are distinct positive multiples of thirteen;

C_H subset union_(i=1)^5 D_(q_i)
             union union_(j=1)^3 D_(c_j)             (2)
```

almost everywhere. Consider one of THM-2266's `120` interior
first-depth-one profiles:

```text
c_1=13u_1,       c_2=13^b u_2,       c_3=13^c u_3,

3<=b<=c-2,       5<=c<=19,
13 does not divide u_1u_2u_3.                        (3)
```

For a scalar relation

```text
z=(alpha_H,alpha_1,...,alpha_5,beta_1,beta_2,beta_3)
  in Z^9,

alpha_H H+sum_i alpha_i q_i+sum_j beta_j c_j=0,      (4)
```

fix the sign convention

```text
A(z)=alpha_H H+sum_i alpha_iq_i
    =-sum_j beta_jc_j.                               (5)
```

Call `A(z)` its guard/unit-versus-blocker carry and call `z` a crossing
when `A(z)!=0`.

The candidate proves that there is a primitive two-coordinate relation `p`
and an integer `m` such that

```text
||p||_infinity<=9841,
A(p)=plus or minus m c_1=plus or minus 13m u_1,
1<=m<=757,             13 does not divide m.         (6)
```

In particular,

```text
p is a crossing,
nu_13(A(p))=1,
kappa:=A(p)/13=plus or minus m u_1,
kappa=-beta_1(p)u_1.                                 (7)
```

The last equality is exact, not merely a congruence.

There is also an exact ancestry statement. Let `E_1` be the full strict
`c_1`-private locus and define

```text
g(y)=sum_(r=0)^12 1_(E_1)((y+r)/13),

F={y:g(y)>0}=T(E_1),       B=T(F),       f=1_F.      (8)
```

Then `E_1` is nonempty and open, `F subset D_(u_1)`, and `F,B` have
positive measure. The multiplicity pushforward `g` and the nonlinear
support indicator `f` retain complementary information:

```text
g_hat(h)=13(1_(E_1))_hat(13h),                 h in Z, (9)

sum_(q congruent kappa mod 13)|f_hat(q)|^2
 >=4 measure(B)/28561
 >0.                                                   (10)
```

For the three small multipliers the phase cone upgrades both statements to
exact atoms. Put `K=A(p)=13kappa`. If `m in {1,2,3}`, then

```text
Re (1_(E_1))_hat(K)
 >=cos(pi*m/7)measure(E_1)>0,

Re g_hat(kappa)
 >=13cos(pi*m/7)measure(E_1)>0,

Re f_hat(kappa)
 >=cos(pi*m/7)measure(F)>0.                          (11)
```

The `3/4` boundary is exact using the retained local data. For every

```text
m in H={4<=m<=757:13 does not divide m},             (12)
```

and every strict profile `(1,b,c)`, there is an explicit local scalar row,
the primitive pair relation

```text
13q_1-mc_1=0,              A=m c_1,                 (13)
```

and an open `c_1`-private set `E` for which `T:E->F` is one-to-one, the
pushforward multiplicity equals `1_F`, all nonzero root residues of `F`
carry positive energy, yet

```text
(1_E)_hat(13m)=g_hat(m)=(1_F)_hat(m)=0.              (14)
```

The exact hard bank has

```text
|H|=696.                                             (15)
```

## 1. The shallow pair is already the desired crossing

THM-2266 proves that one of the six reduced pairs

```text
(H,u_1),             (u_1,q_i), 1<=i<=5             (16)
```

has coprime product at most `757`. There are two cases.

For `(H,u_1)`, put

```text
g_0=gcd(H,u_1),       a=H/g_0,       d=u_1/g_0.      (17)
```

Then

```text
gcd(a,d)=1,         ad<=757,        13 does not divide ad,

13d H-a c_1=0.                                      (18)
```

Take `p_H=13d`, `p_(c_1)=-a`, and all other coefficients zero. Its
carry, with exactly the convention (5), is

```text
A(p)=13dH=a c_1=13a u_1.                            (19)
```

Thus (6)--(7) hold with the positive sign and `m=a`.

For `(u_1,q_i)`, put

```text
g_0=gcd(u_1,q_i),     a=u_1/g_0,       d=q_i/g_0.    (20)
```

Then

```text
gcd(a,d)=1,         ad<=757,        13 does not divide ad,

d c_1-13a q_i=0.                                    (21)
```

Take `p_i=-13a`, `p_(c_1)=d`, and all other coefficients zero. Now

```text
A(p)=-13a q_i=-d c_1=-13d u_1,                      (22)
```

so (6)--(7) hold with the negative sign and `m=d`.

The two vectors are primitive because

```text
gcd(13d,a)=1             in (18),
gcd(d,13a)=1             in (21).                   (23)
```

Their respective heights are

```text
max(13d,a),             max(d,13a),                 (24)
```

and each is at most `13*757=9841`. This proves (6)--(7).

The height-`462` mixed crossing from THM-2275 and a `{-1,0,1}`
perturbation are unnecessary for the stated residue alignment. Adding that
crossing would weaken the height to `10303` without preserving its
coordinatewise seven-unit Fourier provenance. The direct pair is the
strictly stronger object here.

## 2. Private support versus ancestry multiplicity

Define

```text
E_1={
 x:
 ||Hx||>1/7,
 ||c_1x||<1/14,
 ||q_ix||>1/14 for every i,
 ||c_2x||>1/14,
 ||c_3x||>1/14
}.                                                   (25)
```

Every condition in (25) is strict, so `E_1` is open. THM-2268 supplies at
least ten strict `c_1`-private torsion points, hence

```text
E_1 is nonempty and measure(E_1)>0.                  (26)
```

Because `c_1=13u_1`, owner transport gives

```text
x in E_1
  => ||13u_1x||<1/14
  => T(x) in D_(u_1).
```

Therefore `F=T(E_1) subset D_(u_1)`. The circle covering `T` is open and
nonsingular, so `F` and `B=T(F)` have positive measure.

The function `g` in (8) is the ancestry multiplicity, not merely an image
bit. A branchwise change of variables gives, for every integer `h`,

```text
g_hat(h)
 =sum_(r=0)^12 integral_0^1
    1_(E_1)((y+r)/13)exp(-2*pi*i*h*y)dy

 =13 integral_(R/Z)1_(E_1)(x)exp(-2*pi*i*13h*x)dx

 =13(1_(E_1))_hat(13h).                             (27)
```

The branch phase `exp(2*pi*i*h*r)` is one. This proves (9). Passing from
`g` to `f=1_{g>0}` discards multiplicity and is nonlinear; neither object
may silently replace the other.

THM-2269 applies to the marked predecessor support `F`. A unit danger comb
meets any thirteen-root fibre in at most two sheets. Every nonzero
character of each occupied one/two-sheet mask is nonzero, and its
quantitative global form is

```text
sum_(q congruent k mod 13)|f_hat(q)|^2
 >=4 measure(B)/28561,             1<=k<=12.         (28)
```

Both `m` and `u_1` are thirteen-units, so

```text
k=kappa mod 13=plus or minus m u_1 mod 13
```

is nonzero. Substitution in (28) proves (10). This aligns the selected
pair relation with a positive marked residue class, but not yet with a
particular exact atom.

## 3. The low-frequency ancestry phase cone

Let `v` be a positive integer and let `X subset D_v` be measurable with
positive measure. For every `x in X`, choose

```text
theta(x)=vx mod 1 in (-1/14,1/14).                  (29)
```

For either sign and every positive integer `m`,

```text
Re exp(-2*pi*i*(plus or minus mv)x)
 =cos(2*pi*m*theta(x)).                              (30)
```

When `1<=m<=3`, the whole phase interval lies strictly inside the right
half-plane:

```text
|2*pi*m*theta(x)|<pi*m/7<=3pi/7<pi/2.               (31)
```

Cosine is even and decreasing on `[0,pi]`, so

```text
Re (1_X)_hat(plus or minus mv)
 >=cos(pi*m/7)measure(X)>0.                         (32)
```

Take `(X,v)=(E_1,c_1)`. Since `K=plus or minus m c_1`, this proves the
first line of (11). Equation (9), with `h=kappa`, gives the second line.
Taking `(X,v)=(F,u_1)` gives the third line.

Thus the exact relation-to-ancestry map is

```text
pair carry K
  -> phase-cone atom (1_(E_1))_hat(K)
  -> multiplicity atom g_hat(K/13).                 (33)
```

The support atom `f_hat(K/13)` is separately positive in the same small
multiplier range. This is stronger than the residue-class statement (10).

## 4. Sharp normalized cancellation at every hard multiplier

Fix `m>=4`, put `a=1/(4m)`, and choose

```text
0<epsilon<min(1/(4m), 1/14-1/(4m)).                 (34)
```

The second upper bound is positive exactly when `m>=4`. In the centered
representative interval for the circle, put

```text
I_+=(a-epsilon,a+epsilon),
I_-=(-a-epsilon,-a+epsilon),
F_(m,epsilon)=I_+ union I_-.                        (35)
```

The intervals are disjoint, their closures lie in `D_1`, and their equal
Fourier amplitudes at frequency `m` have opposite quarter-turn phases:

```text
(1_(F_(m,epsilon)))_hat(m)

 =sin(2*pi*m*epsilon)/(pi*m)
    [exp(-2*pi*i*m*a)+exp(2*pi*i*m*a)]

 =sin(2*pi*m*epsilon)/(pi*m)
    [exp(-pi*i/2)+exp(pi*i/2)]

 =0.                                                 (36)
```

The coefficient at `-m` is its complex conjugate and also vanishes.

Use the centered real lift in (35) and divide both intervals by thirteen.
After reducing modulo one, call the resulting set `E_(m,epsilon)`. Then

```text
E_(m,epsilon) subset D_13,
T:E_(m,epsilon)->F_(m,epsilon) is one-to-one,

g_E(y):=sum_(r=0)^12 1_E((y+r)/13)
       =1_(F_(m,epsilon))(y).                        (37)
```

The pushforward identity (9) now makes (36) equivalent to

```text
(1_E)_hat(13m)
 =g_E_hat(m)/13
 =(1_F)_hat(m)/13
 =0.                                                 (38)
```

Moreover `F subset D_1`, so every occupied thirteen-root fibre has at most
two marked sheets. Its image under `T` has positive measure, and the
one/two-sheet character proof of THM-2269 gives positive Fourier energy in
every nonzero residue class. If `13` does not divide `m`, the energetic
class `m mod 13` nevertheless contains the zero atom in (36).

This already proves the abstract stopping boundary. The next section
shows that the same carrier can satisfy every *local* labelled eligibility
condition, not only owner support.

## 5. An exact local private-flow witness with the pair relation

Fix `m` in the hard bank (12), and fix any strict profile

```text
2<=b<c,             5<=c<=19.                       (39)
```

Set

```text
u_1=1,              u_2=u_3=2m,

c_1=13,             c_2=2m*13^b,       c_3=2m*13^c,

q_(k+1)=m+13k,                         0<=k<=4,

H=m                  if m is odd,
H=m+13               if m is even.                  (40)
```

All six guard/unit speeds are thirteen-units, `H` is odd, the five `q_i`
are pairwise distinct, and the blockers have the required strict profile.
The scalar hypotheses require the `q_i` to be pairwise distinct; they do
not require `H` to be distinct from every `q_i`.
Choose `r in {0,...,12}` so that

```text
mr congruent 6 mod 13.                              (41)
```

At the two normalized carrier centers put

```text
y_sigma=sigma/(4m),
x_sigma=(r+y_sigma)/13 mod 1,        sigma in {-1,+1}. (42)
```

The shallow blocker owns both points:

```text
||c_1x_sigma||=||y_sigma||=1/(4m)<1/14.             (43)
```

For `q=q_(k+1)`, equation (41) gives

```text
qx_sigma congruent
 [6+sigma*q/(4m)]/13 mod 1.                         (44)
```

Since `0<=k<=4` and `m>=4`,

```text
1/4<=q/(4m)<=7/2,

5/2<=6+sigma*q/(4m)<=19/2.                          (45)
```

For every numerator `N` in the interval (45),

```text
min(N,13-N)>=5/2>13/7.
```

The norm in (44) is therefore strictly larger than `1/7`, and in
particular larger than `1/14`. The chosen `H` equals `q_1` when `m` is odd
and `q_2` when `m` is even, so (45) also proves

```text
||Hx_sigma||>1/7.                                   (46)
```

Finally, for `j=2,3`,

```text
c_jx_sigma
 congruent 2m*13^(lambda_j-1)y_sigma
 =sigma*13^(lambda_j-1)/2
 congruent 1/2 mod 1,                               (47)
```

because every power of thirteen is odd. Thus both deep blockers are
strictly safe.

All inequalities (43)--(47) are strict. By continuity, one may choose
`epsilon>0` small enough that the two intervals

```text
y in (y_sigma-epsilon,y_sigma+epsilon)
```

retain every labelled condition. For an explicit rational choice, take
half the minimum of the disjointness/owner-support bounds in (34) and all
center margins divided by the corresponding slope in `y`. The companion
constructs this number exactly for all `696*150=104400` cases.

Let `F` be these two `y` intervals and lift them through the fixed root

```text
x=(r+y)/13.                                         (48)
```

Their lift `E` is an open, positive-measure, strict `c_1`-private set:
the guard is safe, all five unit combs and both deep blockers are safe, and
only `c_1` is dangerous. The map `T:E->F` is one-to-one, so its ancestry
multiplicity is exactly `g=1_F`.

The first unit speed is `q_1=m`, and the primitive scalar relation

```text
13q_1-mc_1=0                                       (49)
```

has height `max(13,m)<=757`, blocker coefficient `-m`, and carry

```text
A=13q_1=mc_1=13m.                                  (50)
```

Equations (36), (38), and (50) give the simultaneous exact failure

```text
(1_E)_hat(A)=0,
g_hat(A/13)=0,
(1_F)_hat(A/13)=0,                                  (51)
```

while the marked root spectrum of `F` has positive energy in the aligned
nonzero residue `m mod 13`.

This is an honest local source-eligibility, ancestry-branch, and pair-
relation witness. It is **not** asserted to satisfy the global cover (2).
It proves that the missing high-multiplier input is genuinely global or
phase-localized; merely retaining all local owner labels does not supply it.

## 6. Count and fixed-carrier no-go

There are `754` integers from `4` through `757`; exactly

```text
floor(757/13)=58
```

are divisible by thirteen. Hence

```text
|H|=754-58=696.                                     (52)
```

The `3/4` boundary is geometric. At `m=4`, the cancellation centers
`plus or minus 1/16` lie inside `D_1`. At `m=3`, the required centers
`plus or minus 1/12` lie outside `D_1`, and the positive half-plane
certificate (32) forbids cancellation by any positive-measure subset of
`D_1`.

There is also a fixed-carrier no-go. For `F=D_1`,

```text
(1_F)_hat(n)=sin(pi*n/7)/(pi*n),        n!=0.        (53)
```

It vanishes at every nonzero multiple of seven. Multiplication by seven
permutes the thirteen residue classes, so each residue class contains
infinitely many vanishing exact atoms of this one carrier, even though
each nonzero residue class has positive root energy. Residue energy can
never imply support at an arbitrary member of the residue class.

## 7. Connection and loss ledger

The candidate connection is

```text
source:
  THM-2266's bounded shallow reduced pair and THM-2268's strict
  c_1-private points;

map:
  lift the pair to its primitive scalar crossing, push the private-locus
  indicator through T with multiplicity, and evaluate the owner phase cone
  at the original and normalized carries;

preserved:
  the c_1 label, exact scalar sign convention, primitive coefficients,
  height 9841, valuation one, multiplier m, owner unit u_1, ancestry
  multiplicity, image support, and marked root residue;

positive output:
  exact original-ancestry, descended-multiplicity, and image-support
  Fourier atoms when m=1,2,3;

lost:
  THM-2275's seven-unit convolution provenance, current blocker service,
  sibling gluing, and any uniform width of the private intervals;

hard boundary:
  696 possible multipliers, each with an exact local private-flow
  cancellation witness.                                             (54)
```

Even in the low branch, a nonzero marked-ancestry coefficient is not by
itself a contradiction to the scalar cover. The pair relation came from a
covariance bound, not as a nonzero common convolution summand of the
guard/unit and blocker safe approximants. A closure theorem must still
couple this ancestry atom to current service or to a signed cover identity
at the same exact frequency.

The local witnesses in Section 5 do not arise from, and do not refute, a
global counterexample. They freeze the exact first missing implication.
No valuation profile is excluded here, and LRC(14) remains open.

## 8. Exact verification

Run

```bash
python3 04-computation/lrc14_shallow_owner_phase_cone_thm2276.py
python3 -O 04-computation/lrc14_shallow_owner_phase_cone_thm2276.py
```

The primary companion uses only integer and rational arithmetic. It checks
the complete thirteen-unit reduced-pair atlas, both primitive lift types,
the height `9841`, the exact `699=3+696` multiplier split, rational
two-interval cancellation at every hard multiplier, all one/two-sheet
root-character obstructions, and all `104400` local strict-profile
witnesses including their scalar labels, pair relation, center margins,
and an explicit admissible rational half-width. The normal and optimized
transcripts are required to be byte-identical.

Nothing may depend on this file until its proof, dependency scope, and
companion have been independently audited and the reserved status has been
promoted.
