---
id: THM-3097
title: "Translated-support Monge compactification and finite first-bad-prefix spectrum"
status: >
  PROVED CANDIDATE + VERIFIED-EXACT; INDEPENDENT AUDIT REQUESTED.  At every
  fixed width t, translating an arbitrary factorial support sufficiently far
  makes its physical first-window resultant uniformly positive, independently
  of every internal gap.  The minimal bad prefix bank is finite.  Since widths
  one through three are universally good, bad t-slot supports in [0,X] are
  O_t(X^(t-4)); in particular the four-slot bad set is finite.  Conditional
  SFC through width b improves the exponent to t-b-1.  This is existential
  and non-effective and does not prove SFC(4).
source: root-gmc-translated-support-compactification-2026-08-01
depends_on:
  - THM-2824-arbitrary-three-slot-factorial-moment-divisibility-and-atomic-orientation-boundary
  - THM-3069-one-normal-remote-terminal-suspension-and-physical-tropical-flag
  - THM-3073-upper-triangular-resultant-norm-and-torsion-character-death-barcode
  - THM-3093-arbitrary-gap-remote-cluster-monge-flag-compactification
  - THM-3096-physical-support-complete-intersection-and-arbitrary-radial-pair-radical
related:
  - THM-3018-factorial-conjecture-as-a-simplex-moment-problem
  - THM-3086-arbitrary-cluster-composition-chambers-and-alternant-clutch-holotopy
script: 04-computation/gmc_translated_support_compactification_thm3097.py
output: 05-knowledge/results/gmc_translated_support_compactification_thm3097.out
script_sha256: ab8920f64ec986a7566d14e0cabe211e8535f141032e5e6d46daf5d33d831da3
output_sha256: e8c8828f5f4586bc80260ef546be8882c8956ba8f8e0a032b6471de5b31808ed
hash_basis: LF-normalized bytes
---

# THM-3097 -- translated-support compactification and bad-prefix spectrum

**PROVED CANDIDATE + VERIFIED-EXACT; INDEPENDENT AUDIT REQUESTED.**

THM-3093 compactifies every internal gap of a remote cluster above a fixed
good child.  There is a second boundary face that it deliberately excludes:
translate the **entire** support, so the child has width zero and the degree
one moment becomes the first normal equation.  That face is even cleaner.
There is no lower-system error at all, and the degree-one response is the
constant row that anchors the Monge staircase.

Combining this translation face with the gap-composition faces gives more
than a remote-tail theorem.  It makes the arithmetic bad-support locus a
finite collection of cylinders over minimal bad prefixes.  The first three
prefix widths are empty, so the bad locus has arithmetic codimension at
least four.

## 1. Physical translated system

Let

```text
L(s^n)=n!,                    f_n=s^n/n!.                 (1)
```

Fix `t>=2`, a translation `N>=0`, and arbitrary distinct integer gaps

```text
0=h_1<h_2<...<h_t.                                     (2)
```

For `x=(x_1,...,x_t)`, normalize the first `t` moment forms by

```text
F_r(x)=L((sum_j x_j f_(N+h_j))^r)/L(f_N^r),
                                      1<=r<=t.           (3)
```

In particular,

```text
F_1=x_1+...+x_t.                                        (4)
```

Choose once and for all a determinant-one linear coordinate change whose
last coordinate is `w=F_1`.  The uneliminated resultant

```text
R_t(N;h)=Res(F_1,...,F_t)                               (5)
```

is therefore the physical mean-zero first-window resultant obtained by
restricting `F_2,...,F_t` to `w=0`.  This fixes its orientation rather than
remembering it only up to sign.

There are constants

```text
N_tr(t)<infinity,                   eta_t>0              (6)
```

depending only on `t`, such that for every `N>=N_tr(t)` and every gap vector
`(2)`,

```text
R_t(N;h)>0.                                             (7)
```

More precisely, after the normalization in Section 3 its normal resultant
is at least `eta_t`, uniformly over arbitrary internal gap diameter.  For
`t=1`, set `R_1=Res(x_1)=1`; there is no obstruction.

## 2. The degree-one anchored Monge flag

For an integer `h>=0`, put

```text
V_r(h)=(rN+1)_(rh)/(N+1)_h^r.                           (8)
```

Thus `V_1(h)=1`.  Across an adjacent gap interval define

```text
phi_j(r)=r^(-1) log[V_r(h_(j+1))/V_r(h_j)]

 =integral_(h_j)^(h_(j+1))
  [psi(r(N+x)+1)-psi(N+x+1)] dx.                       (9)
```

Its derivative is

```text
partial_r phi_j(r)=integral_(h_j)^(h_(j+1))
 (N+x)psi'(r(N+x)+1) dx>0.                             (10)
```

At the endpoint `N=h_j=0` the integrand vanishes only at one point; every
positive-length integer edge remains strict.  Hence for `1<=a<b`,

```text
[V_a(h_(j+1))/V_a(h_j)]^b
 <[V_b(h_(j+1))/V_b(h_j)]^a.                           (11)
```

Set `r_i=i` and choose dual potentials

```text
lambda_1=1,
lambda_(j+1)/lambda_j
 =[V_j(h_j)/V_j(h_(j+1))]^(1/j),

a_i=V_i(h_i)lambda_i^i.                                (12)
```

Since `V_1=1`,

```text
lambda_1=lambda_2=1>lambda_3>...>lambda_t.              (13)
```

The pure normalized root weights are at most one and expose exactly the
staircase

```text
{(i,i),(i,i+1):i<t} union {(t,t)}.                      (14)
```

This is the finite symbolic relation state.  The equality on the first edge
is structural, not a loss of strict response.

## 3. Exact coefficient contraction and covariance

For a multi-index `alpha` with `|alpha|=r`, the coefficient of `x^alpha` in
`F_r` is its multinomial coefficient times

```text
Q_(r,alpha)
 =(rN+1)_(sum_j alpha_j h_j)
   /product_j (N+1)_(h_j)^alpha_j.                      (15)
```

The function `s -> log(rN+1)_s` is discretely convex and vanishes at zero.
Jensen at the barycentre gives

```text
Q_(r,alpha)^r<=product_j V_r(h_j)^alpha_j.              (16)
```

Define

```text
B_i(x)=a_i^(-1)F_i(lambda_1x_1,...,lambda_tx_t),
Etilde_t(N;h)=Res(B_1,...,B_t).                         (17)
```

Equations `(11)--(16)` bound every coefficient of `B_i` by its multinomial
coefficient.  No lower-child atom or physical outer error exists.

The resultant has total variable degree `t!` and degree `t!/i` in equation
`i`.  Every `lambda_i` occurs once as a variable scale and once in `a_i`, so
its powers cancel exactly:

```text
R_t(N;h)
 =product_(i=1)^t V_i(h_i)^[t!/i] Etilde_t(N;h).        (18)
```

The displayed carrier is positive.  Formula `(18)` is the magnitude and
orientation sidecar omitted by the scale-free Monge flag.

## 4. Compactification of every gap vector

Suppose uniform positivity failed along translations `N_nu->infinity`.
After passing to a subsequence, every adjacent difference

```text
Delta_j=h_(j+1)-h_j                                    (19)
```

is fixed or tends to infinity, while all bounded coefficients in `(17)`
converge.  The fixed edges form consecutive blocks.

For fixed degrees `a<b`, the response across an interval of length `Delta`
satisfies

```text
g_Delta(a)/g_Delta(b)<=c_(a,b)^Delta,       c_(a,b)<1. (20)
```

This includes `a=1`; for example one may take the explicit bound

```text
g_Delta(1)/g_Delta(b)<=(b!)^(-Delta/b).                 (21)
```

Consequently a divergent edge kills every backward column, and Jensen kills
every mixed coefficient using such a column.  Only arbitrary forward terms
remain.

On a fixed block `I=[a,b]`, write `d_j=h_j-h_a`.  The exact common-offset
identity

```text
V_(r,N)(H+d)=V_(r,N)(H)V_(r,N+H)(d)                    (22)
```

rebases the block even when its common origin `H` also diverges.  Since
`N+H->infinity`, its diagonal restriction tends to powers of linear forms
whose line matrix, up to cancelling positive column multipliers, is

```text
L_(i,j)=r_i^d_j/r_i^d_i,                   i,j in I.    (23)
```

Its generalized alternant has the gap-independent floor

```text
det L_I
 =det(r_i^d_j)_(i,j in I)/product_(i in I)r_i^d_i

 >=Vandermonde(r_a,...,r_b)
   /product_(i=a)^b r_i^(i-a)>0.                        (24)
```

Indeed, after removing the ordinary Vandermonde, the remaining Schur
polynomial contains the required ordered monomial with coefficient one.

Iterate THM-3073 from the rightmost block.  Arbitrary forward entries do not
change the resultant, and if `P` is the consecutive-block partition then

```text
lim Etilde_t(N;h)=product_(I in P)(det L_I)^t!>0.        (25)
```

There are only `2^(t-1)` partitions, and `(24)` gives a positive row-only
floor on all of them.  Compactness contradicts the assumed sequence and
proves `(6)--(7)`.  Any block-order reversal has sign raised to `t!`, which
is even for `t>=2`; `t=1` was calibrated separately.

## 5. Support-configuration holotopy

The ordered `t`-slot chamber has two kinds of faces at infinity:

1. the **translation face**, where all slots move right and Sections 1--4
   give the width-zero Monge flag;
2. the **gap-composition faces**, where the first divergent gap leaves a
   fixed lower child and THM-3093 (`q>=2`) or THM-3069 (`q=1`) supplies the
   remote normal block.

After positive carrier normalization the resultant section extends to every
face.  On a gap face its zero divisor is exactly the lower-child resultant:
the normal generalized alternant never vanishes.  Thus the closure of a bad
locus can meet infinity only over a bad child stratum.  This is a genuine
compactified holotopy of support configurations, not only a subsequence
argument.

The `n=0` child boundary is included.  In THM-3093's proof every denominator
is `(n+1)_C`, every response and Jensen inequality remains strict for
`n>=0`, and the same outer contraction applies.  THM-3069 already states
`n>=0`.  No division by the child minimum is used.

## 6. Finite minimal bad-prefix banks

Let

```text
B_t={0<=a_1<...<a_t : R_t(a_1;0,a_2-a_1,...,a_t-a_1)=0}.  (26)
```

Call `c=(c_1,...,c_m)` a **minimal bad prefix** if `c in B_m` and every
proper initial prefix is good.  Let `C_m` be the set of such prefixes.

Each `C_m` is finite.  Otherwise take an infinite sequence of distinct
members.  The translation theorem bounds `c_1`, so pass to a subsequence
with fixed minimum.  Let `j<m` be the maximal fixed prefix length; then
`c_(j+1)->infinity`.  The fixed `j`-prefix is good by minimality.  If
`m-j>=2`, THM-3093 makes the whole support good uniformly over every later
gap; if `m-j=1`, THM-3069 does the same.  Both contradict badness.

Widths one and two are elementary, and THM-2824 proves width three.
Therefore

```text
C_1=C_2=C_3=empty.                                     (27)
```

Every bad `t`-support has a unique first bad prefix, so

```text
B_t subset union_(m=4)^t union_(c in C_m)
 {t-supports extending c}.                              (28)
```

This is an upper cover: an extension of a bad prefix need not itself be bad.
For supports in `[0,X]`, a fixed prefix `c` of length `m` has exactly

```text
binom(X-c_m,t-m)                                        (29)
```

extensions when `c_m<=X`, and zero otherwise.  Hence, for a finite constant
`K_t`,

```text
#B_t(X)
 <=sum_(m=4)^t sum_(c in C_m,c_m<=X) binom(X-c_m,t-m)
 <=K_t(X+1)^(t-4).                                     (30)
```

Since the full chamber has `binom(X+1,t)` supports, bad supports have density
`O_t(X^-4)`.  Equivalently, the bad locus has arithmetic/Minkowski
codimension at least four: the finite first bad prefix is the relation state,
and the free tail coordinates are its carry directions.

More generally, if SFC is known through width `b`, then `C_m` is empty for
`m<=b`, and

```text
#B_t(X)=O_t(X^(t-b-1)).                                 (31)
```

If `b=t-1`, the bad width-`t` set is finite.  Unconditionally,

```text
B_4=C_4 is finite.                                      (32)
```

This does not say that `B_4` is empty.

## 7. Uniform complete-intersection and radial consequences

On every good support, THM-3096 turns `F_1,...,F_t` into a complete
intersection with socle degree

```text
D_t=t(t-1)/2,                         E_t=D_t+1.         (33)
```

For the normalized systems `B` in `(17)`, consider the degree-`E_t` map

```text
T_B: direct_sum_(j=1)^t S_(E_t-j) -> S_(E_t),
     (A_j) |-> sum_j A_j B_j.                           (34)
```

It is surjective.  The coefficient vectors of all translated `B` lie in a
compact multinomial box; Sections 3--4 keep their closure inside the
positive-resultant locus.  Therefore the smallest singular value of `(34)`
has a positive width-only floor.  Minimal-norm right inverses give uniformly
bounded normalized Bezout certificates for every `x_i^E_t`, independently
of every internal gap.  Undoing `(17)` costs only the displayed
`lambda_i`, `a_i`, and carrier scalings; no hidden condition number remains.

For positive supports, THM-3096 consequently gives arbitrary-complex-radial
pair-radical certificates through moment `2t` on all but

```text
O_t(X^(t-4))                                             (35)
```

of the `t`-slot supports in `[1,X]`.  This genuinely handles arbitrary
radial coefficient phases on the good-support population.  It does not
handle neutral-channel cancellation or general charge patterns.

## 8. Sharp boundaries and exact evidence

The thresholds, finite prefix banks, and constants above are existential;
the compactness proof does not enumerate them.  Making `(24)--(25)` and the
resultant Lipschitz constant effective is the next finite-core problem.
Growing width is outside the theorem.

The bad-child boundary is load-bearing.  If a projective child witness kills
not only its first `m` forms but all forms through target width `t`, then
setting every new coefficient to zero makes **every** `t`-slot extension bad.
This literally saturates a cylinder of order `X^(t-m)`.  If later child forms
cut that witness out, zero-tail persistence disappears, but mixed-tail roots
still require a transverse detection-depth sidecar.  Suspension alone cannot
improve the exponent.

As a purely logical sharpness control, declaring an abstract support bad iff
its first four entries are `(0,1,2,3)` obeys the present good-prefix remote
implications and has `Theta(X^(t-4))` bad extensions.  It is not asserted to
be a physical factorial counterexample.

The exact companion verifies:

1. `840` strict response inequalities including `r=1` and `N=0`;
2. `7,602` full-degree multivariate Jensen coefficients;
3. `540` staircase weights, including `lambda_1=lambda_2`;
4. `372` exact composition-face/block-triangular controls;
5. `209` common-offset and resultant-covariance identities;
6. `20` rank-two/rank-three limiting face formulas;
7. `45` exact finite-prefix extension counts; and
8. `15` codimension-four abstract-cylinder controls, plus repeated-gap zeros.

Run

```text
python 04-computation/gmc_translated_support_compactification_thm3097.py
python -O 04-computation/gmc_translated_support_compactification_thm3097.py
```

Both modes must equal the stored transcript after LF normalization.

This theorem does not prove SFC(4), SFC in any higher width, the full
arbitrary-radial GMC(2), NC2 outside the two-charge resultant-good family,
LRC(14), JC(2), or DC(2).

**QED, pending independent audit and status promotion.**
