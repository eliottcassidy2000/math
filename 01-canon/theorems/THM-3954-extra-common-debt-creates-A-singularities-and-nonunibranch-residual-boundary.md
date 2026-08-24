---
id: THM-3954
title: "Extra common debt creates A-singularities and non-unibranch residual boundary"
status: >
  PROVED + VERIFIED-EXACT + INDEPENDENTLY HOSTILE-AUDITED.
  In the universal nonconstant-ratio A1 internal-split packet, arbitrary
  common multiplier c is now divided into two exact cases. If c is coprime
  to the three finite color forms R,(S+omega^2 R),(S-omega R), THM-3951's
  repeated-incidence forest applies. If c vanishes to order m at a finite
  point of one color, the natural cubic surface has completed local ring
  k[[x,y,t]]/(xy-t^(3m)), hence is already normal there, while its one
  irreducible residual ramification prime has two normalization addresses
  collapsing to that point. It is therefore non-unibranch and cannot be an
  affine-plane boundary. This provisionally closes every c in this A1
  packet, not other internal-split grammars or JC(2).
source: jc-extra-debt-local / post-THM-3951 universal shared-color audit, 2026-08-24
audit: >
  INDEPENDENT HOSTILE AUDIT PASS (jc-zero-debt-lift and jc-cohn3709,
  2026-08-24). Both audits reconstructed the primitive linear-in-P domain
  proof, reduced-discriminant index argument for global normality, residual
  prime and birational-image calculation, exact A_(3m-1) completed local
  ring at every shared color, two smooth normalization addresses of that one
  residual prime, and the exhaustive THM-3951/THM-3920 gcd dichotomy. Normal
  and optimized runs byte-match the frozen 92-gate output, all hashes agree,
  documentation checks pass, and no repair was required.
depends_on:
  - THM-3950-a1-internal-split-denominator-debt-and-equianharmonic-shadow
  - THM-3951-affine-plane-boundary-incidence-forest-and-equianharmonic-survivor-nonentry
  - THM-3920-affine-plane-boundary-unibranch-depressed-cubic-chart-obstruction
script: 04-computation/jc2_extra_debt_flank_normalization_20260824.py
output: 05-knowledge/results/jc2_extra_debt_flank_normalization_20260824.out
script_sha256: 8a6c3ce8b6478abfc8d373ed60e45ce3328578a73ad9176a6d105b942e0f6dae
output_sha256: 80f68ffe8a3021a8d1c28407f80f2ee7b3893bf5f12c8f34b847b55bc975bce6
semantic_sha256: 01393277071248a2c4b22382d343627adef25b29ca676db9d99455f107d01bf3
hash_basis: raw LF bytes
---

# THM-3954 -- extra common debt does not separate the residual boundary

**PROVED + VERIFIED-EXACT + INDEPENDENTLY HOSTILE-AUDITED.** Work over an
algebraically closed field `k` of characteristic zero.
Fix

```text
omega^2+omega+1=0,       delta=omega-omega^2,       delta^2=-3.       (1)
```

The candidate closes the common-multiplier boundary left open by THM-3951.
Its scope is the natural cubic over `A2_(P,t)` attached to the universal
nonconstant-ratio `A1_t` packet of THM-3950. It does not classify a different
root gauge, splitting across several cube-difference rows, a non-`A1`
primary branch, or an arbitrary cubic field.

## 1. Universal packet and statement

Let `R,S,c in k[t]` satisfy

```text
gcd(R,S)=1,             R/S nonconstant,             c!=0.           (2)
```

Put

```text
U=S+omega^2 R,                 V=S-omega R,
r0=RUV,
A0=(S-R)V^2,                   B0=(S+R)U^2,
D0=(1-omega^2)B0-(1-omega)A0,
E0=(B0-A0)A0B0.                                             (3)
```

The exact THM-3950 identities compress these to

```text
D0=3r0+delta S^3,             E0=r0^2(2r0-D0).              (4)
```

For the arbitrary extra common multiplier, set

```text
A=cA0,          B=cB0,          r=cr0,
q=cD0 P+c^3E0.                                                   (5)
```

Consider the natural depressed-cubic surface

```text
X_c=Spec B_c,
B_c=k[P,t,Z]/(F),
F=Z^3-3PZ+q.                                                  (6)
```

The sign change `T=-Z` identifies `(6)` with THM-3951's natural cubic. The
provisional conclusion is:

> For every data set `(2)-(6)`, `X_c` is a normal integral finite-flat cubic
> surface. Its graph and residual ramification loci are two distinct prime
> divisors. There is no dense open `A2` in `X_c` on which the finite map
> `X_c -> A2_(P,t)` is etale. Equivalently, this entire natural-cubic
> realization of the nonconstant-ratio `A1` packet has no same-function-field
> planar Keller chart.

The proof separates `gcd(c,RUV)=1`, already handled by THM-3951, from the
extra-debt case. The latter is the new content below.

## 2. The displayed cubic surface is already normal

The coefficient of `P` in `(6)` is `cD0-3Z`; the remaining term is
`Z^3+c^3E0`. Their gcd in `k[t,Z]` is one. Indeed, the first polynomial is
primitive and linear in `Z`, and substituting `Z=cD0/3` in the second gives

```text
c^3 K,                    K=D0^3/27+E0.                    (7)
```

The polynomial `K` is not identically zero. A nonconstant map
`R/S:P1_t -> P1` hits the three distinct colors

```text
R=0,                         U=0,                         V=0. (8)
```

At most one of these three fibres is supported only at `t=infinity`, so one
has a finite point. At every point of `(8)`, coprimality in `(2)` gives
`S!=0`, and `(4)` gives

```text
r0=E0=0,              D0=delta S^3,
K=(delta S^3)^3/27!=0.                                  (9)
```

Thus `(7)` is nonzero. Gauss's lemma applied to the primitive polynomial,
linear in `P`, proves that `F` is irreducible. Hence `B_c` is a domain and,
because `F` is monic in `Z`, it is finite free of rank three over `k[P,t]`.

Its discriminant is

```text
disc_Z(F)=-27H,                 H=q^2-4P^3,               (10)
```

with exact factorization

```text
H=(c^2r0^2-P)N,
N=4P^2-c^2(D0^2-4r0^2)P+c^4r0^2(D0-2r0)^2.              (11)
```

The quadratic discriminant is

```text
disc_P(N)=c^4(D0-2r0)^3(D0+6r0).                        (12)
```

Modulo squares its last two factors are

```text
(R-S)(R+S)(R+delta S)(3R+delta S).                       (13)
```

THM-3950 proves that `(13)` is not a square in `k(t)`: otherwise one
component of `P1_t` would map nonconstantly to the fixed smooth `j=0`
elliptic double cover, contradicting Riemann--Hurwitz. Therefore `N` is
irreducible over `k(t)`, hence over `k[t,P]` by primitivity. It is distinct
from the graph factor because

```text
N(c^2r0^2,t)=-4delta c^4r0^3S^3!=0                     (14)
```

as a polynomial. Thus `H` is reduced.

Now `B_c` is a hypersurface, hence `S2`. Away from `(10)` its finite map is
etale. At the generic point of either factor in `(11)`, the order
discriminant has valuation one. Passage from this order to the maximal order
changes that valuation by twice the local index, so the index is zero.
The corresponding height-one local rings of `B_c` are DVRs. This proves
`R1`, and Serre's criterion proves

```text
X_c is normal.                                             (15)
```

This global argument is complemented by the explicit completed local ring
at every extra-debt point in Section 4.

## 3. There is one residual prime, not two boundary primes

The relative derivative of `(6)` is

```text
F_Z=3(Z^2-P).                                             (16)
```

On `P=Z^2`, the remaining equation factors as

```text
F|_(P=Z^2)=(Z-cr0)Q_c,
Q_c=-2Z^2+c(D0-2r0)Z+c^2r0(D0-2r0).                     (17)
```

This gives the graph ramification prime

```text
E_Gamma=(P-Z^2,Z-cr0)                                   (18)
```

and the residual ramification prime

```text
E_R=(P-Z^2,Q_c).                                        (19)
```

The latter is **one prime**. Its quadratic discriminant in `Z` is `c^2`
times `(12)` without the square `c^4`; its square class is again `(13)`.
Thus `Q_c` is irreducible in `k[t,Z]`, and

```text
B_c/E_R ~= k[t,Z]/(Q_c)                                  (20)
```

is a domain. Moreover

```text
Q_c(Z)Q_c(-Z)=N(Z^2,t).                                  (21)
```

At the generic point, `(16)` and `F=0` give `q=2Z^3` and `P=Z^2`, hence

```text
Z=q/(2P).                                                 (22)
```

Therefore `E_R` maps birationally to the irreducible residual branch
`N=0`; it is the unique residue-degree-one ramification prime there, not a
pair of conjugate boundary primes. Similarly, `(14)` shows that `(18)` and
`(19)` are distinct. Since `X_c` is already the finite normalization, neither
prime can disappear under a further ambient normalization.

## 4. A shared color gives an exact `A_(3m-1)` point

Let `alpha in A1_t` be a finite zero of one of `R,U,V`, and suppose

```text
m=ord_alpha(c)>=1.                                       (23)
```

Only one of `R,U,V` can vanish there: two would force `R(alpha)=S(alpha)=0`.
Consequently `S(alpha)!=0`, and `(9)` says `D0` and `K` are units in the
completed local ring `k[[tau]]`, where `tau=t-alpha`.

Make the exact polynomial coordinate changes

```text
W=3Z-cD0,
X=P-(W^2+3cD0W+3c^2D0^2)/27.                            (24)
```

Direct expansion, with no discarded higher terms, gives

```text
F=-XW+c^3K.                                              (25)
```

Write `c=tau^m u` with `u` a unit. Rescaling `X` by the unit `u^3K`
turns `(25)` into

```text
completed O_(X_c,x_alpha)
   ~= k[[tau,X,W]]/(XW-tau^(3m)).                        (26)
```

Thus the extra-debt point is the rational double point `A_(3m-1)`. In
particular it is a normal surface point. For a simple common factor
`m=1`, it is exactly an `A2` point.

Equation `(26)` answers the normalization question sharply: ambient
normalization does **not** separate the graph and residual curves. There is
only one point of `X_c` above `(P,t,Z)=(0,alpha,0)`.

## 5. The one residual prime has two branches at that point

In the function field of `(19)`, put

```text
xi=Z/c.                                                    (27)
```

The denominator-free equation is

```text
Q_0=-2xi^2+(D0-2r0)xi+r0(D0-2r0)=0,                     (28)
Q_c(cxi)=c^2Q_0.                                         (29)
```

The element `xi` is integral over `k[t]`, so adjoining it to `(20)` gives a
finite birational partial normalization. At the shared color point,
`r0=0` and `D0` is a unit, hence

```text
Q_0(alpha,xi)=xi(D0(alpha)-2xi).                         (30)
```

It has the two distinct simple roots

```text
xi=0,                         xi=D0(alpha)/2.             (31)
```

The partial normalization is regular at both points, so the full
normalization of the single irreducible curve `E_R` has two distinct points
there. Under its finite map back to `E_R`, both have

```text
Z=cxi=0,                         P=c^2xi^2=0,             (32)
```

and hence both map to the one normal surface point in `(26)`. Therefore

```text
E_R is non-unibranch at x_alpha.                          (33)
```

This is stronger than failure of a clean pairwise-intersection test. The
two addresses in `(31)` are two branches of **one** irreducible ramification
prime, not two primes that could be assigned to different boundary
components.

## 6. Exhaustion of every common multiplier

Suppose the function field of `(6)` admitted plane coordinates `x,y` with

```text
k(x,y)=Frac(B_c),          P,t in k[x,y],
Jac_(x,y)(P,t) in k*.                                      (34)
```

Normalization-form Zariski Main would identify `A2_(x,y)` with an open
`U0 subset X_c`. Because `(34)` is etale, every ramification prime of the
finite map `X_c -> A2_(P,t)`, including `(18)` and `(19)`, must lie in
`X_c minus U0`.

There are now exactly two cases.

1. If `gcd(c,RUV)=1`, THM-3951 gives at least two distinct clean incidences
   between the two boundary primes `E_Gamma` and `E_R`. Its boundary-tree
   forest lemma excludes `U0 ~= A2`.
2. If `gcd(c,RUV)!=1`, there is a finite common zero `alpha` over the
   algebraically closed field. Sections 4--5 make the mandatory boundary
   prime `E_R` non-unibranch at `x_alpha`. THM-3920 says every irreducible
   boundary curve of a normal affine surface containing `A2` is unibranch.
   This again excludes `U0`.

These cases exhaust every nonzero polynomial `c`. The apparent way around
THM-3951--making one clean color collision non-clean--therefore produces a
strictly stronger self-incidence obstruction.

## 7. Minimal extra debt and the explicit hostile

For `R=1,S=Y`, write

```text
U=Y+omega^2,                 V=Y-omega,
g0=UV=Y^2-delta Y-1.                                     (35)
```

The literal two-clean-intersection hypothesis of THM-3951 first fails when
`U|c` or `V|c`; a single linear factor is the minimal asymmetric extra debt.
To corrupt both flanks symmetrically requires `UV|c`, whose minimal member is

```text
c=UV.                                                     (36)
```

At `U=0` and `V=0`, the factors `(A,B)` have vanishing-order vectors
`(1,3)` and `(3,1)`, respectively. Both zeros of `c` are simple, so `(26)`
gives two `A2` surface points. The residual equations in `(30)` are exactly

```text
U=0:  -xi(delta+2xi)=0,
V=0:  -xi(-delta+2xi)=0.                                (37)
```

Thus each `A2` point receives two residual normalization addresses. The
graph prime also passes through the point, but `(33)` already excludes the
candidate before the two-prime forest is used.

No extra multiplier actually survives. What is minimal is only the debt
needed to evade the *literal clean-intersection premise*: one flank factor,
or `UV` for the symmetric attempt. The cost is an `A2` singularity and a
non-unibranch residual boundary.

## 8. Reproduction and open boundary

Run

```bash
python3 04-computation/jc2_extra_debt_flank_normalization_20260824.py
python3 -O 04-computation/jc2_extra_debt_flank_normalization_20260824.py
```

Both runs must byte-match
`05-knowledge/results/jc2_extra_debt_flank_normalization_20260824.out`. The
92 exact gates freeze the universal factorization, four-color discriminant,
natural-cubic coordinate change, three color specializations, residual-prime
norm, explicit `c=UV` hostile, and both collapsed normalization addresses.

This theorem closes only the one-factor, nonconstant-ratio `A1` packet in
the natural `A2_(P,t)` cubic realization. Simultaneous splitting across
several cube-difference factors, genuinely bivariate factors, non-`A1`
primary branches, and the planar Jacobian conjecture remain **OPEN**.
