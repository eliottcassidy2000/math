---
id: THM-2844
title: "Sharp support-cut and signed adjacent-ray orientation boundary"
status: >
  PROVED + VERIFIED-EXACT + INDEPENDENTLY HOSTILE-AUDITED.  Among
  coordinatewise-positive adjacent-difference orthants,
  THM-2830's cubic orientation holds for every pair exactly when one
  support cut separates them.  On the first signed adjacent ray
  V_t=t d_1+d_2 with t>-1, the Pascal quotient is monotone exactly for
  t>=alpha, where alpha is the unique root in (-2/3,-1/2) of
  5t^3+30t^2+57t+24.  Primitive integer and Pascal-positive hostiles show
  both boundaries are sharp.  Below alpha the orientation certificate
  fails, but the first three factorial moments still detect the plane;
  this is not an SFC or GMC counterexample.
source: root/signed-adjacent-ray-orientation-boundary-2026-07-28
depends_on:
  - THM-2824-arbitrary-three-slot-factorial-moment-divisibility-and-atomic-orientation-boundary
  - THM-2830-disjoint-positive-adjacent-cone-factorial-moment-three-detection
related:
  - THM-2830-disjoint-positive-adjacent-cone-factorial-moment-three-detection
  - THM-2841-all-order-adjacent-difference-factorial-tensor-positivity
script: 04-computation/gmc_signed_adjacent_pascal_boundary_thm2844.py
output: 05-knowledge/results/gmc_signed_adjacent_pascal_boundary_thm2844.out
script_sha256: be181a810745af5cd93abfa74c698010328fac94b772a075e0a27bc799ec9029
output_sha256: 2eb5e5455f9ad8ee51c2f2a00da48da02e3445551e0d214d17dcbc8b3fbae2c3
hash_basis: LF-normalized bytes
---

# THM-2844 -- sharp support-cut and signed adjacent-ray boundary

**PROVED + VERIFIED-EXACT + INDEPENDENTLY HOSTILE-AUDITED.**

Let

```text
L(s^n)=n!,                 f_n=s^n/n!,
d_n=f_(n+1)-f_n,                                      (1)

A_n(V)=L(d_nV),           B_n(V)=L(d_nV^2),
R_n(V)=B_n(V)/A_n(V),                                (2)

D(U,V)=2L(V^3)L(UV)-3L(UV^2)L(V^2).                  (3)
```

The theorem identifies both the maximal coordinate orthants on which
`D>=0` universally and the first exact signed extension of THM-2830's
Pascal-quotient certificate.

## 1. A support cut is necessary and sufficient

Let `I` be a finite nonempty subset of `Z_(>=0)` and `J` a finite
nonempty subset of `Z_(>=1)`.  Then

```text
D(U,V)>=0                                                   (4)
```

for every nonzero

```text
U=sum_(i in I)lambda_i d_i,       lambda_i>=0,
V=sum_(j in J)mu_j d_j,           mu_j>=0                 (5)
```

if and only if

```text
max I<min J.                                               (6)
```

For a singleton upper ray `V=d_j`, THM-2830's global quotient
monotonicity gives the sharper sign law

```text
D(d_i,d_j)
 =6A_i A_(j-1)(R_(j-1)-R_i),

sign D(d_i,d_j)=sign(j-i-1).                              (7)
```

Sufficiency in `(6)` is THM-2830 with the cut `b=min J`.
If `(6)` fails, choose `i=max I` and `j=min J`; the singleton
specialization in `(7)` is negative.  In the separated orthant, equality
occurs exactly at an occupied adjacent cut:

```text
U=lambda d_(b-1),                V=mu d_b.                (8)
```

If either adjacent label is absent, the inequality is strict.  Thus a
support cut is not cosmetic; every overlap or reversal contains a
negative singleton witness.  For example,

```text
D(d_2,d_1)=-228.                                         (9)
```

## 2. The maximal first signed adjacent ray

Put

```text
V_t=t d_1+d_2,                    t>-1,
u=t+1.                                                   (10)
```

Then every linear and quadratic Pascal response stays positive:

```text
A_n=(n+1)(t+(n+2)/2)>0,                                 (11)

6B_n=
 n^5+10n^4+47n^3+122n^2+168n+96
 +u(5n^4+38n^3+121n^2+184n+108)
 +6u^2(n^3+6n^2+13n+10)>0.                             (12)
```

Define

```text
Delta_n=A_nB_(n+1)-A_(n+1)B_n.                          (13)
```

The first step is

```text
Delta_0=2P(t),
P(t)=5t^3+30t^2+57t+24.                                (14)
```

For every `n>=1`, with `m=n-1`,

```text
4Delta_n=
 [m^6+17m^5+119m^4+427m^3+788m^2+632m+104]
 +u[6m^5+86m^4+494m^3+1410m^2+1980m+1080]
 +u^2[12m^4+140m^3+608m^2+1152m+792]
 +u^3[8m^3+72m^2+208m+184]
 >0.                                                    (15)
```

Since

```text
P'(t)=15(t+2)^2-3>0                         for t>=-1,
P(-2/3)=-58/27,             P(-1/2)=19/8,              (16)
```

there is a unique

```text
alpha in (-2/3,-1/2),
alpha=approximately -0.5820748606339202,                (17)
```

with `P(alpha)=0`.  Therefore

```text
R_0,R_1,R_2,... is nondecreasing
 iff t>=alpha.                                          (18)
```

At `t=alpha`, only `R_0=R_1`; every later step is strict.
Direct derivative integration also gives

```text
D(d_0,V_t)=6Delta_0=12P(t).                             (19)
```

Thus the same algebraic endpoint is the exact maximal
`D(d_0,V_t)>=0` interval inside this Pascal-positive signed slice.

## 3. Primitive integer and Pascal-positive hostiles

The primitive adjacent integer ray

```text
V_*=3d_2-2d_1                                           (20)
```

has

```text
A_n=(n+1)(3n+2)/2>0,
B_n=(n+2)(3n^4+29n^3+123n^2+253n+208)/2>0,
Delta_0=-116,                      D(d_0,V_*)=-2088.    (21)
```

All later quotient steps are positive.  Indeed, for `n>=1`,

```text
4Delta_n=
 27n^6+351n^5+1863n^4+4901n^3
 +6066n^2+2344n-464>0,                                 (22)
```

as is immediate after shifting `n=m+1`.  The ray `(20)` is
coefficient-`l_1`-height minimal among primitive integer rays in
`Z d_1 direct-sum Z d_2`:
heights below five allow no negative ratio smaller than `-1/2`, while
`alpha<-1/2`; height five first admits `-2/3`.

Even strict positivity of every `A_n` and `B_n` is insufficient for
quotient monotonicity.  The globally minimal primitive example is

```text
V_dagger=2d_3-d_1.                                      (23)
```

It obeys

```text
A_n=(n+1)(n^2+5n+3)/3>0,

B_n=(n+2)
 (n^6+26n^5+273n^4+1519n^3
  +4796n^2+8151n+5814)/18>0,                           (24)

R_0=646,              R_1=1715/3,       Delta_0=-446. (25)
```

For `n=m+1`,

```text
Delta_n=(2/27)(
 m^9+42m^8+765m^7+7920m^6+51198m^5
 +213129m^4+566273m^3+912909m^2
 +791487m+269703)>0.                                   (26)
```

Any signed primitive vector of coefficient height two has `A_0=0`, so
`(23)` is minimal.  This refutes the tempting implication

```text
all A_n,B_n positive  =>  B_n/A_n monotone.            (27)
```

The square map loses the total-positivity coordinate.  Notice also that

```text
D(d_0,V_dagger)=24792>0.                               (28)
```

One local quotient dip need not reverse every longer secant.

## 4. Certificate failure is not moment failure

For `U=d_0`, `V=V_t`, and `u=t+1>0`, THM-2824's first binary
quadratic/cubic division invariant is

```text
I_1=-2(7u^3+43u^2+97u+81)<0.                           (29)
```

Consequently

```text
L(H)=L(H^2)=L(H^3)=0    implies H=0
```

for every complex `H in span_C(d_0,V_t)` throughout the whole
Pascal-positive range `t>-1`, including `alpha>t>-1` where `(19)` is
negative.

The endpoint `alpha` is therefore a sharp boundary for one
orientation/quotient certificate.  It is not a factorial-moment
counterexample, an SFC counterexample, or a failure of Gaussian
moment-six detection.

## 5. Exact companion

The exact companion:

1. symbolically derives the direct `A_n,B_n,Delta_n,D`, and `I_1`
   identities;
2. exhausts `156` singleton sign cells and `3,969` Cartesian support
   pairs;
3. tests `1,112` exact rational parameters around the algebraic endpoint;
4. scans primitive adjacent integer rays through coefficient height `50`,
   with `(-2,3)` the unique first hostile;
5. rederives both primitive hostiles by literal divided-power
   multiplication; and
6. checks the height-two, endpoint, equality, reversed-cut, and
   certificate-survival boundaries.

All truth-bearing gates are explicit exceptions and all arithmetic is
integer, rational, or symbolic polynomial arithmetic.  Reproduce with

```text
python 04-computation/gmc_signed_adjacent_pascal_boundary_thm2844.py
python -O 04-computation/gmc_signed_adjacent_pascal_boundary_thm2844.py
```

Both modes byte-match the stored transcript.

## 6. Independent hostile audit

An independent audit rederived the support-cut iff against THM-2830's
proper transport-ordered interlaced subcone, every symbolic response and
endpoint identity, both primitive hostiles, the equality cases, and the
`I_1` survivor.  It restricted the first height-minimality statement to
the fixed lattice `Z d_1 direct-sum Z d_2`, exactly the universe proved
and scanned.  It independently replayed normal and optimized companions
against the stored transcript and declared LF hashes.

**QED.**
