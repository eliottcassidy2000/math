---
id: THM-2945
title: "Nonnegative complete-intersection norm and repeated-divisor gate"
status: >
  PROOF-COMPLETE CANDIDATE + VERIFIED-EXACT; AWAITING INDEPENDENT
  HOSTILE AUDIT.  A nonzero real polynomial H which is nonnegative on
  [0,infinity) has interior zero set exactly the positive-real root set
  of gcd(H,H').  Strict positivity on the closed ray additionally
  requires H(0)>0.  For a real ternary (2,3,4) complete intersection
  with no real Q=C point, the locally dehomogenized multiplication norm
  is nonnegative; Poisson's formula and polynomial continuity give a
  globally oriented nonnegative resultant.  Positive denominator
  clearings and positive extraneous flag factors preserve its
  positive-ray zero set, so factorial first-window failure reduces
  conditionally to a repeated-divisor root test.  No universal
  terminal-pole divisor formula or new width closure is claimed here.
source: codex-gmc-uniform-width-extension-2026-07-29
depends_on:
  - THM-2843-four-slot-projective-divisibility-and-resolvent-reduction
  - THM-2925-general-width-terminal-pole-cancellation-and-macaulay-degree-law
  - THM-2942-macaulay-extraneous-flag-factor-and-pluecker-mutation
related:
  - THM-2824-arbitrary-three-slot-factorial-moment-divisibility-and-atomic-orientation-boundary
  - THM-2943-width-seven-eight-two-chart-macaulay-resultant-closure
  - THM-2944-width-nine-ten-two-chart-macaulay-resultant-closure
script: 04-computation/gmc_nonnegative_norm_repeated_divisor_gate_thm2945.py
output: 05-knowledge/results/gmc_nonnegative_norm_repeated_divisor_gate_thm2945.out
script_sha256: cad3d52fc9bfd1a848b55ad298aa2588a3e193a5eaca42009f77049009b037e7
output_sha256: 1fa19116bb21d63676b00ba16ef161025f33fdb14e6160d7952fbddcefa05313
thm2843_dependency_sha256: 4832a9e4cda1608473e3a4bcfcb880a7fa3f1c6db47b7e939b4b5183cb9549aa
thm2925_dependency_sha256: 83d70a95f0943992d0e4b7027eede431d4dc968b66655e37b43fd0acfc692e47
hash_basis: LF-normalized bytes
---

# THM-2945 -- nonnegative complete-intersection norm and repeated-divisor gate

**PROOF-COMPLETE CANDIDATE + VERIFIED-EXACT; AWAITING INDEPENDENT
HOSTILE AUDIT.**

## 1. The positive-ray lemma

Let

```text
H in R[t]\{0},                    H(t)>=0 for t>=0,    (1)
S=gcd(H,H'),                     S monic.             (2)
```

Then

```text
{xi>0:H(xi)=0}={xi>0:S(xi)=0}.                       (3)
```

Consequently,

```text
H(t)>0 for every t>0
 iff S has no positive real root,                     (4)
```

and

```text
H(t)>0 for every t>=0
 iff H(0)>0 and S has no positive real root.          (5)
```

To prove `(3)`, write at an interior real zero

```text
H(t)=(t-xi)^m U(t),                  U(xi)!=0.         (6)
```

If `m` were odd, the sign would change across `xi`, contradicting
`(1)`.  Hence `m>=2` is even and `H'(xi)=0`.  This proves the forward
inclusion.  The reverse inclusion is immediate from the definition of
the gcd.

The zero polynomial is a separate identically-degenerate branch; its
gcd with its derivative is not used here.

## 2. The endpoint is not an interior point

Equation `(5)` cannot be weakened to a gcd test alone.  The polynomial

```text
H(t)=t                                                   (7)
```

is nonnegative on `[0,infinity)` and has

```text
gcd(H,H')=1,                                           (8)
```

but vanishes at the endpoint.  More generally, a one-sided
nonnegative polynomial may have an odd-order zero at `0`.  Thus every
factorial application must either verify `H(0)>0` directly or hand the
unshifted support to a separate theorem.

Nonnegativity is also load-bearing.  The sign-indefinite polynomial

```text
H(t)=t-1                                               (9)
```

has a simple positive root and gcd `1`.  The first failed implication
without `(1)` is precisely the even-multiplicity conclusion in `(6)`.

## 3. Positive clearing and positive multiplier transport

Let `N` be a rational function on the ray and suppose

```text
N(t)=P(t)/D(t),               D(t)>0 for t>=0,
N(t)>=0 for t>=0.                                     (10)
```

Let `E in R[t]` satisfy

```text
E(t)>0 for t>=0,                                      (11)
```

and assume

```text
H=E*P                                                  (12)
```

is nonzero.  Then

```text
H(t)>=0,
Z_[0,infinity)(H)=Z_[0,infinity)(N).                  (13)
```

The lemma applies to `H`.  Thus the positive-ray zero problem for `N`
can be decided from `H(0)` and the positive real roots of
`gcd(H,H')`, without first separating the positive multiplier `E`.

The complete factorization of the gcd is not invariant under this
operation.  For example,

```text
N=1,                    E=(t+1)^2                     (14)
```

gives

```text
H=(t+1)^2,              gcd(H,H')=t+1.                (15)
```

The extra repeated root is negative and irrelevant to `(3)`.  Hence a
nontrivial gcd does not obstruct positivity on the ray, and a scan of
negative repeated roots is not by itself an intrinsic formula for the
norm.

Strict positivity in `(10)--(11)` is essential for preservation of
the zero set.  If the proposed clearing multiplier is only
nonnegative, then

```text
N=1,                    E=(t-3)^2                     (16)
```

manufactures a false positive-ray repeated divisor at `t=3`.

## 4. Complete-intersection norms are nonnegative

Let `Q_t,C_t,F_t` be real ternary forms of degrees `2,3,4`, depending
algebraically on a real parameter `t`.  Assume on a dense parameter
open set that

```text
Q_t,C_t form a regular sequence,
V_R(Q_t,C_t)=empty.                                   (17)
```

Choose locally a real linear form `L` missing the finite scheme
`V(Q_t,C_t)`, and set

```text
A_(t,L)
 =R[x0,x1,x2]/(Q_t,C_t,L-1),
g_(t,L)=F_t/L^4.                                      (18)
```

Poisson's product formula gives, with one convention-dependent sign
`epsilon in {+1,-1}`,

```text
epsilon*Res(Q_t,C_t,F_t)
 =Res(Q_t,C_t,L)^4 det_R(m_g:A_(t,L)->A_(t,L)).       (19)
```

The real Artin-local factors of `A_(t,L)` occur in conjugate pairs.
On a conjugate pair, multiplication by `g` contributes

```text
det_C(m_g at P) det_C(m_g at Pbar)
 =|det_C(m_g at P)|^2>=0.                             (20)
```

This remains true for nonreduced local factors.  The other factor on
the right of `(19)` is a fourth power.  Therefore

```text
epsilon*Res(Q_t,C_t,F_t)>=0                           (21)
```

on every such local chart.

No single dehomogenizing `L` is asserted to work globally.  The
admissible `L`-charts cover the dense regular-sequence locus, and
`(19)` agrees on their overlaps because its left side is the same
homogeneous resultant.  Polynomial continuity extends `(21)` across
all isolated exceptional parameters of this one-parameter regular
locus, whether the exception comes from a failed selected chart or a
nonregular pair.  In the factorial positive-definite case below the
quadratic is geometrically irreducible, so nonregularity is exactly
`Q_t|C_t`; there it is enough that `Q_t` not divide `C_t` over `R(t)`.

If divisibility holds generically, the resultant is identically zero
and the nonzero-polynomial hypothesis in Section 1 fails.  That is a
separate common-component branch, not a repeated-root conclusion.

## 5. Factorial first-window typing

Fix a normalized translated four-slot support

```text
(n,n+a,n+b,n+M),             0<a<b<M,                 (22)
```

and eliminate the mean as in THM-2843.  Let

```text
Q_n=L(H^2),             C_n=L(H^3),             F_n=L(H^4) (23)
```

be the resulting ternary forms.

The THM-2925 rational continuation is obtained by dividing the Gamma
moment coefficients by common positive Gamma factors.  Hence its
quadratic is a positive rescaling of the Gram form of the distinct
functions

```text
s^(n+d)/Gamma(n+d+1),                 d in {0,a,b,M},  (23a)
```

in `L^2(R_+,e^(-s)ds)`.  Thus `Q_n` remains positive definite for
every real `n>=0`, not only for integer depths.  THM-2824, as used in
THM-2843, supplies a physical specialization with `Q` not dividing
`C`; therefore divisibility cannot hold over `R(n)`.  The local
complete-intersection argument applies on the generic parameter locus
and extends across its isolated exceptional parameters.  With the
fixed resultant convention, choose `epsilon` once so that

```text
N_(M,a,b)(n)
 =epsilon*Res(Q_n,C_n,F_n)>=0                         (24)
```

on the nonnegative physical ray wherever the rational continuation is
defined.

THM-2925 supplies denominator clearings made from factors `n+r` with
`r>0`; their resultant powers are strictly positive for `n>=0`.
THM-2942 further identifies the two variable chart coordinates

```text
p_J0=-q200*K,                   p_J1=-P_alt,           (25)
```

and, in the finite width bank it audits, the normalized common
extraneous factor

```text
G_extr=gcd(K,P_alt)                                   (26)
```

has positive constant term and nonnegative coefficients.  Thus

```text
G_extr(n)>0                         for n>=0.          (27)
```

After any positive THM-2925 clearing, one may therefore work with a
polynomial of the form

```text
H_(M,a,b)(n)
 =positive clearing * G_extr(n) * N_(M,a,b)(n).       (28)
```

It is `H`, not necessarily the bare resultant, that a common-chart gcd
computation naturally returns.  Equations `(24)`, `(27)` and Section 3
show that `H` is nonnegative and has exactly the same nonnegative-ray
zeros as the complete-intersection norm.

Consequently, provided `H` has been constructed exactly and is
nonzero,

```text
first-window failure at a real n>0
 iff gcd(H,H') has a positive real root at n.          (29)
```

For integer-depth closure it suffices to prove that the gcd has no
positive real root and then check `H(0)>0` separately.  This is a
reduction, not a claim that the repeated divisor has a uniform formula.

## 6. An exact conjugate-pair model

The companion includes a literal six-point model of `(19)--(21)`.
On the conic coordinate put

```text
c_6(u,v)
 =(u^2+v^2)(u^2+4v^2)(u^2+9v^2),                     (30)
f_8(u,v)
 =v^6(u^2+(t-2)v^2).                                 (31)
```

The sextic has the three conjugate pairs

```text
u/v=+-i, +-2i, +-3i,                                  (32)
```

and no root at infinity.  Exact elimination gives

```text
Res(c_6,f_8)
 =(t-3)^2(t-6)^2(t-11)^2,                            (33)
gcd(Res,partial_t Res)
 =(t-3)(t-6)(t-11).                                  (34)
```

Thus every interior zero is visibly a conjugate-pair collision and is
detected by the derivative gcd.

## 7. Finite factorial scout, not proved content

An independent exact scout, whose constructor is not yet a stable
dependency of this theorem, reports the following stronger pattern.
For every normalized support `(0,a,b,M)` of widths `6<=M<=10`, form

```text
H=G_extr*N                                            (35)
```

after the chosen positive clearing and remove scalar content.  Across
all

```text
10+15+21+28+36=110                                   (36)
```

supports, the reported radical is

```text
rad gcd(H,H')
 ~ product_(j=1)^a (2n+2j-1)
   product_(r=a)^M (n+r),                             (37)
```

independent of `b`.  Every displayed root is negative.  The reported
width digests begin

```text
M=6: bde57c31...,       M=7: 5149bd24...,
M=8: ed574f22...,       M=9: 41ac672e...,
M=10:07317515....                                      (38)
```

This is **SCOUT**, not a dependency and not part of the proved
candidate.  Formula `(37)` concerns `H=G_extr*N`, not the bare norm:
the positive extraneous multiplier can contribute negative repeated
roots exactly as in `(14)--(15)`.  A stable constructor, full hashes,
normal/optimized replay, and direct endpoint checks are required before
`(37)` can close these widths canonically.

No all-width extrapolation of `(37)` is made.  In particular, this
theorem does not assert a universal terminal-pole divisor law.

## 8. Exact companion and scope

The exact companion:

- hash-pins the audited THM-2843 and THM-2925 sources;
- verifies the conjugate-pair and nonreduced-local norm signs;
- proves `(33)--(34)` by exact resultant and gcd computation;
- checks `1,296` factor-structured nonnegative ray polynomials,
  including `324` simple-endpoint cases, `720` negative repeated-factor
  cases, and `1,152` cases with an interior zero;
- checks the endpoint hostile `(7)`, sign hostile `(9)`, positive
  clearing hostile `(14)`, and non-strict clearing hostile `(16)`.

The companion does not reproduce the provisional factorial scan
`(35)--(38)`.  The abstract proof is Sections 1--5; the computation
audits its algebra and sharp boundaries.

Run

```text
python 04-computation/gmc_nonnegative_norm_repeated_divisor_gate_thm2945.py
python -O 04-computation/gmc_nonnegative_norm_repeated_divisor_gate_thm2945.py
```

Both modes must byte-match the stored transcript and the LF-normalized
hashes above.

**QED, pending independent hostile audit.**
