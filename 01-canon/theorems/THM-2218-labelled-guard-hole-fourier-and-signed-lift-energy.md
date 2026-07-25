---
id: THM-2218
title: "Labelled guard-hole Fourier transform and signed lift energy"
status: >
  PROVED + VERIFIED-EXACT. A scalar thirteen-lift guard-hole vector is an
  exact sum of phasewise cyclic correlations. Its full lift-label data has
  equivalent cyclotomic Fourier, integral group-algebra, and confluent
  Hasse descriptions; the scalar family sum retains only the constant
  mode. The theorem proves the root-digit laws for the depth-three
  profiles, a signed top-k energy inequality with equality cases, an exact
  common-lift regret identity, and a family-knapsack bound. A complete
  integer stress audit exhibits 22 locally negative fibres whose global
  deficit is repaired by regret 3,668, and verifies the mixed hostile
  family-knapsack margin 16,361/13. These identities do not supply the
  missing uniform all-depth correlation estimate and do not prove LRC(14).
source: klein-2026-07-24-labelled-guard-hole-fourier-energy
depends_on:
  - THM-2204-scalar-depth-223-thirteen-lift-capacity-law
related:
  - THM-2201-cyclic-root-fibre-hasse-jet-transition-carrier
  - THM-2205-scalar-depth-113-exact-lift-capacity-exclusion
  - THM-2207-scalar-depth-123-labelled-guard-hole-exclusion
  - THM-2216-residual-capacity-hinge-gram-law
script: 04-computation/lrc14_labelled_guard_hole_fourier_energy_thm2218.py
output: 05-knowledge/results/lrc14_labelled_guard_hole_fourier_energy_thm2218.out
script_sha256: 753523f776970d3c1c608a062b132914e40cb7744b7e96bb7c092651b18cf2f0
output_sha256: 9799cd4ae6a75842ecf76746179555e5a789d8bd41eafb8ecd74e3d74ab8e039
hash_basis: working-tree bytes (LF)
---

# THM-2218 -- labelled guard-hole Fourier transform and signed lift energy

Sections 1--7 are analytic and `PROVED`. Section 8 is a separate
`VERIFIED-EXACT` finite stress audit: it tests the usefulness and sharp
boundary of the proved inequalities, but it is not another exclusion
proof. THM-2205 and THM-2207 close `(1,1,3)` and `(1,2,3)` independently
by complete exact audits.

We use the notation and reversed guard trivialization of THM-2204.

## 1. Conventions

Put

```text
p=13,                  G=F_13,
Q=13^m,                N=13Q.
```

The guard is normalized to one. A primitive image phase is represented by

```text
r in {1,...,Q-1},       p does not divide r,
y_r=r/Q.
```

All sheet expressions below are in `G`. Thus `r^(-1)` means the inverse of
the residue of `r` modulo `p`, and every exponent of a group-algebra
variable is reduced modulo `p`.

For a finite set `P` of primitive image phases, retain THM-2204's notation

```text
B(r) = guard-safe sheets above r,
F(r) = G\B(r),                                      (A1)

I_a(r) = Z intersection
         (ar/Q-p/14, ar/Q+p/14),                    (A2)

delta = a^(-1) mod p,                               (A3)

A_(a+lQ)(r)
 = {delta(d+lr):d in I_a(r)},       l in G.         (A4)
```

Here `a` is a unit modulo `Q`. The set `I_a(r)` has one or two integers.
There are no multiplicities in (A4), because `delta` is a unit and the two
possible integers are consecutive.

Define

```text
E_a(P)=sum_(r in P)|I_a(r)|,                        (A5)

b_(a,l)(P)
 =sum_(r in P)|A_(a+lQ)(r) intersection F(r)|,     (A6)

C_(a,l)(P)=E_a(P)-b_(a,l)(P).                      (A7)
```

Thus `C_(a,l)` is the guard-safe capacity and `b_(a,l)` is its labelled
guard-hole loss. The value `E_a(P)` is independent of `l`.

For Fourier transforms choose

```text
zeta=exp(2*pi*i/p),
hat(f)(h)=sum_(s in G) f(s) zeta^(-hs),             (A8)
hat(b)(h)=sum_(l in G) b_l zeta^(-hl).              (A9)
```

These transforms are unnormalized.

## 2. Exact cyclic-correlation and Fourier formula

For each phase define the endpoint and hole functions

```text
e_(a,r)(s)=1_{s in {delta*d:d in I_a(r)}},
f_r(s)=1_{s in F(r)}.                               (A10)
```

Use the cyclic correlation

```text
(e star f)(x)=sum_(s in G)e(s)f(s+x).               (A11)
```

Then (A4) gives the exact phasewise recursion

```text
b_(a,l)(P)
 =sum_(r in P)(e_(a,r) star f_r)(delta*r*l).        (A12)
```

This is a recursion in the lift label: every phase contributes one cyclic
correlation, sampled after the phase-dependent automorphism
`l -> delta*r*l`.

### Fourier identity

For every `h in G`,

```text
hat(b_a)(h)
 =sum_(r in P)
    hat(f_r)(h*a*r^(-1))
    sum_(d in I_a(r)) zeta^(h*d*r^(-1)).            (A13)
```

**Proof.** If `c=e star f`, then directly from (A8) and (A11),

```text
hat(c)(j)=hat(f)(j)hat(e)(-j).                      (A14)
```

If `v` is nonzero, the transform of `l -> c(vl)` at frequency `h` is
`hat(c)(hv^(-1))`. In (A12), `v=delta*r`, so

```text
j=h(delta*r)^(-1)=h*a*r^(-1).
```

Finally,

```text
hat(e_(a,r))(-j)
 =sum_(d in I_a(r))zeta^(j*delta*d)
 =sum_(d in I_a(r))zeta^(h*d*r^(-1)).
```

Summing over `r` proves (A13). QED.

At `h=0`, (A13) becomes

```text
hat(b_a)(0)
 =sum_(r in P)|F(r)||I_a(r)|
 =sum_(l in G)b_(a,l)(P).                           (A15)
```

Since `C_(a,l)=E_a-b_(a,l)`,

```text
hat(C_a)(0)=pE_a(P)-hat(b_a)(0),                    (A16)
hat(C_a)(h)=-hat(b_a)(h),             h nonzero.    (A17)
```

Equation (A16) is THM-2204's family-sum law. Equations (A13) and (A17)
identify its missing coordinates: the twelve nonconstant Fourier modes.
Parseval gives the exact two-sided energy

```text
sum_(l in G)|b_(a,l)-bar(b_a)|^2
 =1/p sum_(h nonzero)|hat(b_a)(h)|^2.               (A18)
```

The phases of the twelve nonconstant modes in (A13), not merely their
magnitudes, are load-bearing.

### Power-spectrum obstruction

This loss has a sharp abstract witness. Let `x in R^13` be nonconstant and
let `tau_j x` be its cyclic translate by `j`. Compare two packets of
thirteen local vectors:

```text
aligned packet:       x,x,...,x;
rotating packet:      tau_0 x,tau_1 x,...,tau_12 x. (A18a)
```

Every local vector in both packets has the same mean and the same Fourier
magnitudes. Nevertheless, for every `1<=k<=12`,

```text
H_k(sum aligned packet)=13H_k(x),
H_k(sum rotating packet)=k sum_l x_l,

13H_k(x)>k sum_l x_l.                               (A18b)
```

Here `H_k` is the top-`k` support function defined formally in (A42a)
below. Indeed, the rotating aggregate is the constant vector
`(sum_l x_l)1`. The average of the `k`-subset sums of `x` is
`k(sum_l x_l)/13`, so their maximum is at least that average. Equality
would make all `k`-subset sums equal; exchanging one chosen and one
unchosen coordinate then forces all coordinates of `x` equal, contrary to
the hypothesis. Thus no recursion which keeps only local zero modes and
power spectra can determine or sharply reconstruct the aggregate top-five
statistic. Relative Fourier phase, or an equivalent rooted sidecar, is
necessary.

## 3. Integral group algebra and the confluent Hasse chart

Let

```text
mathcal B_(a,P)(X)
 =sum_(l in G)b_(a,l)(P)X^l
 in Z[X]/(X^p-1).                                   (A19)
```

The coefficient formula (A6) has the exact integral form

```text
mathcal B_(a,P)(X)
 =sum_(r in P)
   sum_(d in I_a(r))
   sum_(s in F(r))
      X^(r^(-1)(a*s-d)).                            (A20)
```

**Proof.** A triple `(r,d,s)` contributes to the coefficient of `X^l`
precisely when

```text
s=delta(d+lr).
```

Multiplying by `a=delta^(-1)` and solving for `l` gives

```text
l=r^(-1)(a*s-d) mod p.                              (A21)
```

This is exactly (A20). QED.

Evaluation at `X=zeta^(-h)` recovers (A13). Thus (A20) is an exact
integer carrier; no reduction modulo `p` has occurred.

Now reduce (A19) modulo `p` and put

```text
X=1+epsilon,
F_p[C_p]=F_p[epsilon]/(epsilon^p).                  (A22)
```

Write

```text
mathcal B_(a,P)(1+epsilon)
 =sum_(j=0)^(p-1) J_j(a;P) epsilon^j.               (A23)
```

Then

```text
J_j(a;P)
 =sum_(l in G)b_(a,l) binom(l,j)

 =sum_(r in P)
   sum_(d in I_a(r))
   sum_(s in F(r))
      binom(r^(-1)(a*s-d),j)              mod p.    (A24)
```

The Pascal matrix `(binom(l,j))_(0<=l,j<p)` is triangular with diagonal
one, so the full jet `(J_0,...,J_12)` reconstructs the vector
`(b_(a,l) mod p)_l`. Deck translation is triangular unipotent exactly as
in THM-2201.

There is an important scope boundary:

> The Hasse jet (A24) reconstructs the accumulated correlation vector only
> modulo thirteen. The integers `b_(a,l)` can be much larger than thirteen.
> Order inequalities therefore require the integral polynomial (A20), the
> cyclotomic Fourier values with their phases, or an equivalent integer/
> 13-adic lift. The reduced Hasse jet alone does not determine signs or
> order statistics.

This is the precise Fourier/Hasse relation: Fourier separates the thirteen
semisimple characters over characteristic zero; Hasse jets separate their
confluent collision at `X=1` in characteristic thirteen.

## 4. Two-blocker inclusion-exclusion

Let `Omega` be the primitive phase universe under consideration and put

```text
P_(u,v)=Omega\(D_u union D_v).                       (A25)
```

For any subset `S` of phases, let `mathcal B_(a;S)` denote (A20) with
`P=S`. Linearity and

```text
1_(P_(u,v))
 =1_Omega-1_(Omega intersection D_u)
          -1_(Omega intersection D_v)
          +1_(Omega intersection D_u intersection D_v)
```

give

```text
mathcal B_(a;P_(u,v))
 =mathcal B_(a;Omega)
  -mathcal B_(a;Omega intersection D_u)
  -mathcal B_(a;Omega intersection D_v)
  +mathcal B_(a;Omega intersection D_u intersection D_v).  (A26)
```

Consequently (A26) holds coordinatewise for `b_(a,l)`, modewise for
`hat(b_a)(h)`, and jetwise for `J_j(a;P)`.

No distinctness is required. If `u=v`, the last three terms reduce to one
subtraction, as they should. Open/closed endpoint choices make no
difference on the power-of-thirteen torsion layers used by THM-2204,
because those layers are coprime to seven and fourteen.

## 5. Root-digit recursion for the depth-three profiles

Assume

```text
Q=pQ_0,                 Q_0=13^(m-1),
```

with `m>=2`. Every primitive phase has a unique form

```text
r=s+kQ_0,               k in G,
s primitive mod Q_0.                                  (A27)
```

Because `Q_0=0 mod p`,

```text
r=s mod p,              r^(-1)=s^(-1) mod p.        (A28)
```

Define the local labelled kernel

```text
K_(a;s,k)(X)
 =sum_(d in I_a(s+kQ_0))
  sum_(z in F(s+kQ_0))
       X^(s^(-1)(a*z-d)).                            (A29)
```

If

```text
P_s={k in G:s+kQ_0 in P},                            (A30)
```

then (A20) disintegrates exactly as

```text
mathcal B_(a,P)(X)
 =sum_s sum_(k in P_s) K_(a;s,k)(X).                (A31)
```

The two blocker types have different digit laws. Let `beta` be the
coefficient of a blocker on the image-phase circle modulo `Q`.
For the lower phase define the corresponding integer window

```text
I_beta^(0)(s)
 =Z intersection
  (beta*s/Q_0-p/14,beta*s/Q_0+p/14).                (A31a)
```

* If `p` does not divide `beta`, then above the lower phase `s/Q_0` its
  blocked digit set is

  ```text
  L_beta(s)
   ={-beta^(-1)d mod p:d in I_beta^(0)(s)},         (A32)
  ```

  a singleton or a two-digit chord.

* If `beta=p gamma`, then every digit is blocked when
  `s/Q_0 in D_gamma`, and no digit is blocked otherwise.              (A33)

To verify (A32), write the thirteen phase roots as
`(s/Q_0+k)/p` and repeat THM-2198's integer-window calculation:
`d=pn-beta*k`, so `k=-beta^(-1)d`. Equation (A33) is the
all-or-nothing divisible-blocker law.

After the unique depth-three blocker disappears, an actual blocker of
valuation `lambda` has image coefficient `13^(lambda-1)u`. Therefore:

```text
(1,1,3):  P_s is G minus the union of two singleton/chords;
(1,2,3):  P_s is empty when the divisible bit is on,
           and otherwise G minus one singleton/chord;
(2,2,3):  P_s is either all of G or empty.           (A34)
```

This explains structurally why `(2,2,3)` scalarizes and why the two
shallower-depth profiles retain a labelled digit correlation.

## 6. Signed top-k lift energy

The ordinary Parseval variance is unnecessarily pessimistic when one lift
is a large negative outlier. The correct elementary refinement separates
the two signs.

### Lemma (signed top-k bound)

Let `c_0,...,c_12` be real numbers, put

```text
S=sum_l c_l,                 mu=S/13,
x_l=c_l-mu,

V_+=sum_l max(x_l,0)^2,
V_-=sum_l max(-x_l,0)^2.                           (A35)
```

For every `1<=k<=12`,

```text
sum of the k largest c_l
 <=k*mu+min(sqrt(k V_+),sqrt((13-k)V_-)).           (A36)
```

The statement is unchanged by ties: the left side is the maximum over all
`k`-subsets, and every maximizing subset obeys the same inequality.

**Proof.** Let `L` be any `k`-subset. Since `sum_l x_l=0`,

```text
sum_(l in L)x_l
 <=sum_(l in L)max(x_l,0)
 <=sqrt(k V_+),                                      (A37)

sum_(l in L)x_l
 =-sum_(l notin L)x_l
 <=sum_(l notin L)max(-x_l,0)
 <=sqrt((13-k)V_-).                                  (A38)
```

The last steps are Cauchy--Schwarz. Take the smaller bound and maximize
over `L`. QED.

### Equality and boundary cases

If every `c_l` is equal, then `V_+=V_-=0` and (A36) is equality for every
`L`.

Assume the vector is nonconstant.

* Equality in the positive bound in (A37) requires exactly `k` positive
  deviations, all equal, with all positive energy supported on `L`; no
  negative deviation may lie in `L`.
* Equality in the negative bound in (A38) requires exactly `13-k`
  negative deviations, all of equal magnitude, with all negative energy
  supported on the complement of `L`.
* If the two square-root bounds tie and (A36) is equality, the vector is
  two-level: `k` equal positive deviations on `L` and `13-k` equal
  negative deviations off `L`, with the magnitudes forced by zero sum.

If only one square-root term is the minimum, equality in (A36) needs only
the equality conditions for that active side. Ties cause no issue in the
inequality. In every nonconstant equality case above, however, the
maximizing `L` is unique: it is the set of exactly `k` positive deviations
in the positive case, or the complement of the exactly `13-k` negative
deviations in the negative case. Only the constant case has multiple sharp
maximizing subsets.

### Exact integer certificate for five lifts

For integer capacities `C_l`, define

```text
S=sum_l C_l,
z_l=13C_l-S,
Z_+=sum_l max(z_l,0)^2,
Z_-=sum_l max(-z_l,0)^2.                            (A39)
```

Then (A36) at `k=5` is

```text
Top5(C)
 <=(5S+sqrt(min(5Z_+,8Z_-)))/13.                   (A40)
```

Let `R` be an integer residual sheet mass and put

```text
D=13R-5S.                                           (A41)
```

The fully integral sufficient certificate

```text
D>0,
D^2>min(5Z_+,8Z_-)                                  (A42)
```

implies `Top5(C)<R`. Equality in the second line of (A42) gives only the
nonstrict conclusion `Top5(C)<=R`; it is not an uncovered-point
certificate. If `D<=0`, (A40) supplies no strict conclusion even though
the actual top five may still miss `R`.

For a THM-2204 family, apply this with `C_l=E_a-b_(a,l)`. The centered
vectors of `C` and `-b` coincide, so the exact Fourier phases in (A13), or
the integer coefficients in (A20), determine `V_+` and `V_-`. Their sum is
the Parseval energy, but the split into signs is strictly stronger for
order statistics.

### Lemma (common-lift regret)

For a vector `x in R^13`, write

```text
H_k(x)=max_(L subset G, |L|=k) sum_(l in L)x_l.     (A42a)
```

If a capacity vector is decomposed over lower fibres as

```text
C=sum_s C^(s),
```

define its common-lift regret by

```text
Reg_k((C^(s))_s)
 =sum_s H_k(C^(s))-H_k(sum_s C^(s)).                (A42b)
```

Then `Reg_k>=0`. If the residual mass also decomposes as `R=sum_s R_s`,
one has the exact margin identity

```text
R-H_k(C)
 =sum_s (R_s-H_k(C^(s)))+Reg_k((C^(s))_s).          (A42c)
```

Moreover, including all tie cases,

```text
Reg_k=0
```

if and only if there is one `k`-subset `L` which is a maximizer in
(A42a) for every local vector `C^(s)`.

**Proof.** For every fixed `L`,

```text
sum_(l in L)sum_s C_l^(s)
 <=sum_s H_k(C^(s)).
```

Maximizing the left side proves nonnegativity. Subtract and add the local
envelope to obtain (A42c). Equality in the displayed inequality after the
maximum is possible exactly when some common `L` attains every local
maximum; this argument does not select or break ties. QED.

Thus negative local margins are not fatal. A proof may instead extract a
positive turn tax from the incompatibility of the locally winning lift
sets. Fourier phases in (A13) are one coordinate system for this regret.

## 7. Family-knapsack bound for the global top five

In the finite sign-class audit, partition all unit sign classes into their
thirteen-lift base families. To retain the cyclic order, use the raw
oriented lifts

```text
q_l=a+lQ mod N.                                     (A43)
```

Replacing a raw lift by its canonical sign representative does not change
capacity because `D_q=D_(-q)`, but it can reverse the displayed cyclic
order. Fourier and Hasse calculations must use (A43), not sorted sign
representatives.

For each base family `a`, define

```text
T_a(k)=sum of its k largest capacities,             (A44)

U_a(k)=k*mu_a+
       min(sqrt(kV_(a,+)),sqrt((13-k)V_(a,-))),
U_a(0)=0.                                           (A45)
```

Let `Top_K(all lifts)` be the sum of the `K` largest capacities over the
disjoint union of all lift families. Then

```text
Top_K(all lifts)
 =max_(sum_a k_a=K) sum_a T_a(k_a)
 <=max_(sum_a k_a=K) sum_a U_a(k_a),                (A46)
```

where `0<=k_a<=min(13,K)`.

**Proof.** Any `K`-subset chooses a unique number `k_a` from each family.
For fixed occupancies, its best contribution is the sum of the familywise
top-`k_a` values. This proves the equality. Apply (A36) in each family to
obtain the inequality. Ties merely create more maximizing subsets or
occupancy vectors. QED.

The right side is a six-state max-plus recursion when `K=5`:

```text
d_0(0)=0,             d_0(j)=-infinity for j>0,

d_i(j)=max_(0<=k<=j)
       (d_(i-1)(j-k)+U_(a_i)(k)).                   (A47)
```

There is a rational, integer-auditable version. With `S_a,z_(a,l),Z_(a,+),
Z_(a,-)` as in (A39), put

```text
M_a(k)=min(kZ_(a,+),(13-k)Z_(a,-)),

tilde(U)_a(k)
 =(kS_a+ceil(sqrt(M_a(k))))/13.                     (A48)
```

Then `T_a(k)<=U_a(k)<=tilde(U)_a(k)`. Running (A47) on the integer
numerators in (A48) gives a rigorous rational upper bound, with no
floating-point comparison. Equality in (A46) requires both an optimizing
occupancy vector and equality in every active family bound; otherwise the
knapsack bound is strict.

## 8. Finite stress audit, kept separate from the proof

Everything in this section is `VERIFIED-EXACT` finite evidence. It is not
used in Sections 1--7.

### 8.1. The spread-2460 family

For THM-2204's residual associated with the diagonal depth-two pair `(4,4)`
and base `a=1098`, use the oriented lifts `q_l=1098+2197l`. The exact
vectors are

```text
C =
(2460,2456,2454,2450,2456,2452,   0,
 2452,2450,2456,2452,2452,2454),

E_a=3190,

b=E_a-C =
(730,734,736,740,734,738,3190,
 738,740,734,738,738,736).                           (A49)
```

The residual has

```text
1716 image phases:
  1228 guard fibres of size 9,
   488 guard fibres of size 10,

R=15932,
S=sum C=29444,
actual family Top5=12282,
actual family margin=3650.                          (A50)
```

For (A39)--(A42),

```text
z =
(2536,2484,2458,2406,2484,2432,-29444,
 2432,2406,2484,2432,2432,2458),

Z_+=72261760,
Z_-=866949136,
D=59896,

D^2=3587530816,
min(5Z_+,8Z_-)=361308800,
integer certificate gap=3226222016.                 (A51)
```

The signed bound is

```text
Top5(C)<=12786.77881374741,
R-bound=3145.22118625259.                           (A52)
```

For comparison, the ordinary two-sided projection/variance bound is

```text
Top5(C)
 <=5*mean(C)+sqrt((5*8/13)sum_l(C_l-mean(C))^2)
 =15459.81547207350,
```

leaving only `472.18452792650`. The sign split removes the false cost of
the single zero-capacity outlier.

### 8.2. Exact obstruction to a fibrewise top-five proof

Disintegrate the same residual one level further, over the `156` primitive
lower phases modulo `169`. For each lower phase `s`, let

```text
R_s = residual guard-safe sheet mass above s,
C_l(s) = the contribution of that lower fibre to lift l,
m_s=R_s-Top5_l(C_l(s)).                             (A53)
```

The exact hostile values are

```text
number of lower fibres with m_s<0: 22,
min_s m_s=-5,
sum_s m_s=-18.                                      (A54)
```

Two minimum rows are

```text
s=4:
R_s=120,
C(s)=(25,7,25,13,20,25,0,25,20,13,25,7,25),
Top5=125;

s=6:
R_s=120,
C(s)=(25,25,25,20,13,7,0,7,13,20,25,25,25),
Top5=125.                                           (A55)
```

Thus no proof obtained by summing nonnegative per-lower-fibre top-five
margins can handle even the already-closed `(2,2,3)` hostile row.

The exact repair is the common-lift regret. Put

```text
local envelope
 =sum_s Top5_l(C_l(s))
 =R-sum_s m_s
 =15950,

global common-lift Top5
 =Top5_l(sum_s C_l(s))
 =12282,

regret=15950-12282=3668.                            (A56)
```

Then

```text
global margin
 =sum_s m_s+regret
 =-18+3668
 =3650.                                             (A57)
```

The missing positive quantity is therefore not local capacity. It is the
incompatibility of the locally winning lift labels across lower fibres,
equivalently the phase alignment retained by (A13)/(A20).

### 8.3. Depth-three-profile stress rows

The scratch audit tested (A42) on the seven particularly concentrated base
families

```text
a in {1,1098,799,599,1015,597,1013}.                (A58)
```

For each of `a=1,1098`, it checked all

```text
514605 unordered (1,1,3) blocker pairs,
 79092 ordered   (1,2,3) blocker pairs.             (A59)
```

For each of the other five bases in (A58), it repeated both complete pair
universes. There were zero failures of the exact integer certificate
(A42). This is a seven-family stress test, not a proof over all `1014`
base families.

On the exact mixed-profile hostile residual with blocker labels `(799,46)`,

```text
R=13526,
actual global Top5=11918,
actual margin=1608.                                 (A60)
```

The global family-knapsack bound (A46), evaluated over all `1014` base
families, is

```text
Top5(all lifts)<=159477/13
                =12267.461538461539,

R-bound=16361/13
       =1258.538461538461.                          (A61)
```

The five actual winning families are

```text
(799,599,1015,597,1013),
```

with exact maxima

```text
(2604,2472,2292,2288,2262).                         (A62)
```

The unrounded signed `k=1` bounds equal the first four maxima and give
`2262.276870424941` for the last. Their direct sum therefore leaves
`1607.723129575059`. Equation (A61) instead uses the fully global dynamic
program and the integer ceiling in (A48).

### 8.4. Reproduction status

The analytic proof in Sections 1--7 needs no computation.

The spread row and the underlying THM-2204 capacities are independently
reproduced by the canonical artifact

```text
04-computation/lrc14_scalar_depth223_thirteen_lift_capacity_thm2204.py
05-knowledge/results/lrc14_scalar_depth223_thirteen_lift_capacity_thm2204.out
```

whose THM-2204 frontmatter records

```text
script SHA256:
9c16e1e7a69834f9304c877f1627232f374688aa7f8d2372c260dfdbfa056ac8

output SHA256:
0b89f80f6f82c63bbc6700239f133ff1e4c9ed3b2fe6389bc604f6b8de229bb8.
```

The additional spread, lower-fibre, seven-family exhaustive, and mixed
family-knapsack checks are frozen in the canonical artifacts

```text
04-computation/lrc14_labelled_guard_hole_fourier_energy_thm2218.py
05-knowledge/results/lrc14_labelled_guard_hole_fourier_energy_thm2218.out
```

with hashes

```text
script SHA256:
753523f776970d3c1c608a062b132914e40cb7744b7e96bb7c092651b18cf2f0

output SHA256:
9799cd4ae6a75842ecf76746179555e5a789d8bd41eafb8ecd74e3d74ab8e039.
```

The normal and `python3 -O` transcripts are byte-identical to the stored
output. Reproduce with

```text
python3 04-computation/lrc14_labelled_guard_hole_fourier_energy_thm2218.py
python3 -O 04-computation/lrc14_labelled_guard_hole_fourier_energy_thm2218.py
```

## 9. Scope and consequence

The proved contribution is an exact analytic recursion and an inequality:

```text
source:
  a base unit class a and a two-blocker residual P;

map:
  sum phasewise cyclic correlations and take either the integral
  group polynomial, its cyclotomic Fourier values, or its confluent
  Hasse reduction;

preserved:
  lift labels, Fourier phase, blocker inclusion-exclusion, and the
  common-lift order statistic through the signed-energy/knapsack bound;

lost by the scalar family sum:
  all twelve nonconstant modes;

lost by Fourier magnitudes alone:
  cross-phase alignment;

lost by the reduced Hasse jet alone:
  integer lifts and hence signs/order;

cheapest hostile control:
  (A49)--(A57).                                     (A63)
```

The root-digit law (A34) identifies the exact new state in `(1,1,3)` and
`(1,2,3)`. The signed-energy and family-knapsack lemmas turn that state
into rigorous order-statistic certificates. By themselves they do not
prove either valuation profile empty without a uniform estimate or a
complete exact audit over every residual and base family; those exclusions
are supplied independently by THM-2205 and THM-2207. QED for Sections
1--7.
