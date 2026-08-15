---
id: THM-3079
title: "Newton-PF row transform and strict integer-mesh terminal-minus-one cone"
status: >
  PROVED + VERIFIED-EXACT + INDEPENDENTLY HOSTILE-AUDITED.  Every finite
  PF-infinity row filter with positive constant
  coefficient preserves reversed strict total positivity of the reciprocal-
  Gamma kernel on its positive shifted strip.  If P has only negative real
  roots, its falling-factorial coefficient polynomial is PF-infinity: under
  P->(x+r)P it evolves by C->(z+r)C+zC', whose roots strictly interlace.
  Consequently every strict integer shape mesh with terminal prefix -1 has
  the checkerboard sign at every generalized Hankel order whenever the common
  base satisfies alpha_0+u_min+v_min>M, where M is the residual polynomial
  degree.  In particular H_n=(a+m)_n/(a)_n^2 is closed for every integer m>=1
  and every minor with a+u_min+v_min>m.  The shallow strip, nonintegral gaps,
  and terminal prefixes below -1 remain open.
source: codex-gmc-strict-prefix-hostile-2026-08-01
depends_on:
  - THM-3056-product-gamma-reciprocal-hypergeometric-and-hankel-reversal
related:
  - THM-3053-beta-gamma-prefix-transport-and-multiplicative-holotopy-cone
  - THM-3065-reciprocal-beta-gap-contiguous-hankel-wall
script: 04-computation/gmc_laguerre_pf_strict_reciprocal_beta_gap_thm3079.py
output: 05-knowledge/results/gmc_laguerre_pf_strict_reciprocal_beta_gap_thm3079.out
script_sha256: 9971f45a464c7f2fcd6a4c086a278d095f2b9eb8fb14e692497b961ba79527e5
output_sha256: 64b1800ec7fe672ebc268389408ac24694d5436d8a1484fc84af6e25e8fa16b7
independent_script: 04-computation/gmc_newton_pf_thm3079_independent_audit.py
independent_output: 05-knowledge/results/gmc_newton_pf_thm3079_independent_audit.out
independent_script_sha256: 57e2e7333db66f596429f2f428284ec6e3a84d7dcdc8302436a8b73b980585f2
independent_output_sha256: ce94015d6dc7fc140a565ce152b71c46f1f10b1428dac8ad4d7a8f5eaa6cae34
independent_semantic_sha256: 7c413b0b7d5c568f77e795c74150b91337621ab018801bff8450721dfa0cfaeb
audit: >
  An independent Fraction/Leibniz implementation reconstructs the Newton
  coefficients by finite differences, rechecks Toeplitz and generalized
  minors without importing the primary companion, and supplies an exact
  below-strip hostile.  A separate line audit rederived the reciprocal-Gamma
  orientation, Newton interlacing, distinguished Cauchy--Binet term, and tail
  shift.  The positive-strip hypothesis remains load-bearing.
hash_basis: LF-normalized bytes
---

# THM-3079 -- Newton-PF row transforms close a strict integer-mesh cone

**PROVED + VERIFIED-EXACT + INDEPENDENTLY HOSTILE-AUDITED.**

THM-3065 proves that nonpositive dual prefixes force the checkerboard sign at
order two but not at higher order: a zero prefix permits exact order-three and
order-four hostiles.  Its strict-prefix repair survived the finite banks but
had no theorem beyond order two.

A mixed all-order cone is now closed.  Integer reverse-Beta gaps are not
handled by a generic Hadamard-product closure.  Their residual is a polynomial
with negative roots; its Gregory--Newton coefficients form a Pólya-frequency
filter acting on the row coordinate of the mandatory reciprocal-Gamma
baseline.  The single Laguerre edge is only the first instance.

## 1. Statement

Let `m>=1` be an integer and let `a>m`.  Put

```text
H_n=(a+m)_n/(a)_n^2,                         n>=0.       (1)
```

For every order `r>=1` and every pair of strictly increasing nonnegative
integer tuples

```text
u_1<...<u_r,                v_1<...<v_r,                 (2)
```

one has

```text
sign det[H_(u_i+v_j)]_(i,j=1)^r
  =(-1)^(r(r-1)/2).                                      (3)
```

In particular every such minor is nonzero.  A positive geometric factor
`c^n` may be inserted in `(1)` without changing `(3)`.

There is a useful abstract form.  Let `M>=0`, let `a>M`, and let
`c_0,...,c_M` have `c_0>0` and a totally nonnegative lower Toeplitz matrix
`[c_(p-q)]`, with `c_j=0` outside `0<=j<=M`.  Then

```text
F_n=sum_(j=0)^M c_j Gamma(a)/Gamma(a+n-j)              (3a)
```

obeys the same strict generalized-Hankel sign law `(3)`.  Thus a finite
PF-infinity filter is an all-order closure operation for the reversed
reciprocal-Gamma kernel, provided its entire shifted row strip stays positive.
The rest of the theorem supplies the nontrivial Laguerre filter which turns
`(3a)` into `(1)`.

On the two-shape mesh `(a,a+m)`, the exponent inventory of `(1)` is

```text
e=(-2,+1),                    E=(-2,-1).                 (4)
```

Thus `(3)` is an all-order theorem inside the strict-dual-prefix cone which is
not covered by THM-3056's coordinatewise-negative product-Gamma theorem.  The
positive exponent in `(4)` is essential to that distinction.

Section 6 proves the stronger tail form: for arbitrary tuples `(2)`, it is
enough that

```text
a+u_1+v_1>m.                                           (4a)
```

## 2. The Laguerre falling-factorial filter

Write `x^(underline j)=x(x-1)...(x-j+1)`.  The elementary connection formula
between rising and falling factorials is

```text
(x+1)_m=sum_(j=0)^m c_(m,j) x^(underline j),            (5)

c_(m,j)=m! binom(m,j)/j!>0.                             (6)
```

For completeness, `(5)--(6)` follow by induction.  Multiplication by
`x+m+1` sends

```text
x^(underline j)
  -> x^(underline (j+1))+(m+j+1)x^(underline j),        (7)
```

and the coefficients in `(6)` satisfy exactly that recurrence, beginning
with `c_(0,0)=1`.

The coefficient generating polynomial is

```text
C_m(z)=sum_(j=0)^m c_(m,j)z^j=m!L_m(-z),                (8)
```

where `L_m` is the ordinary Laguerre polynomial.  Its `m` roots are positive
and simple, so every root of `C_m` is negative and simple.  Hence

```text
C_m(z)=c_(m,0) product_(ell=1)^m(1+lambda_ell z),
lambda_ell>0.                                           (9)
```

This makes the orientation of the Pólya-frequency claim explicit.  Define the
bi-infinite lower Toeplitz matrix

```text
T_(p,q)=c_(m,p-q),
c_(m,j)=0 outside 0<=j<=m.                              (10)
```

Each factor `(1+lambda_ell z)` in `(9)` gives a nonnegative lower-bidiagonal
Toeplitz path matrix.  Their ordered product is `(10)`, up to the positive
factor `c_(m,0)`.  Therefore every minor of `T` in increasing row and column
order is nonnegative.  This is the concrete PF-infinity sidecar used below;
real-rootedness is not being invoked as a free-standing slogan.

## 3. The reciprocal-Gamma kernel on the shifted strip

Set

```text
r_x=Gamma(a)/Gamma(a+x),                   x>-a.       (11)
```

THM-3056 proves the required reciprocal-Gamma sign law for nonnegative row
nodes.  The identical proof works on the full strip in `(11)`, and we record
the extension because the negative rows are load-bearing here.

Let `x_1<...<x_r` all satisfy `x_i>-a`, and let the `v_j` be as in `(2)`.
Put `V=v_r`.  Multiplication of row `i` by the positive quantity
`Gamma(a+x_i+V)/Gamma(a)` changes column `j` into

```text
(a+x_i+v_j)_(V-v_j).                                  (12)
```

Writing `y_i=a+x_i>0`, these are evaluations of the nested positive-factor
polynomials

```text
Q_v(y)=product_(t=v)^(V-1)(y+t).                       (13)
```

Reverse the columns.  Their degrees are then strictly increasing, and the
factor list `V-1,V-2,...,0` makes their coefficient matrix a product of
nonnegative tail-bidiagonal matrices.  Cauchy--Binet with the positive
generalized Vandermonde matrix `[y_i^d]` is strictly positive: the term using
the actual degree set has coefficient minor one.  Undoing the reversal gives

```text
sign det[r_(x_i+v_j)]=(-1)^(r(r-1)/2).                 (14)
```

Thus the reversed-column kernel is strictly totally positive on every finite
part of the strip `x>-a`.

## 4. Exact PF row transform

The Pochhammer rectangle identity gives

```text
(a+m)_n/(a)_n=(a+n)_m/(a)_m.                          (15)
```

Apply `(5)` with `x=a+n-1`.  Since

```text
(a+n-1)^(underline j) r_n=r_(n-j),                    (16)
```

where `(11)` supplies the meaning when `n<j`, equations `(1),(5),(15)` yield
the exact row transform

```text
H_n=1/(a)_m sum_(j=0)^m c_(m,j) r_(n-j).              (17)
```

This is the place where the discrete integer gap enters.  It turns the
reverse-Beta factor into a finite PF filter; no such finite row transform is
asserted for a nonintegral gap.

Fix `(2)` and let

```text
Z={u_i-j:1<=i<=r, 0<=j<=m}                            (18)
```

with duplicates removed and the remaining integers put in increasing order.
Define

```text
A_(i,z)=c_(m,u_i-z),                                  (19)
K_(z,j)=r_(z+v_(r+1-j)).                              (20)
```

Then `(17)` says literally

```text
[H_(u_i+v_(r+1-j))]=1/(a)_m A K.                      (21)
```

The orientation and strictness are now transparent.

- `A` is a submatrix of the increasing-coordinate Toeplitz matrix `(10)`, so
  every maximal minor of `A` is nonnegative.
- Every `z in Z` satisfies `z>=-m>-a`.  Hence every maximal minor of `K` is
  strictly positive by the reversed-column form of `(14)`.
- Choose in Cauchy--Binet the nodes `z_i=u_i`.  Then
  `[A_(i,u_j)]` is lower triangular: entries above the diagonal vanish and
  every diagonal entry is `c_(m,0)=m!`.  Its determinant is `(m!)^r>0`.

Therefore every term in the Cauchy--Binet expansion of `(21)` is
nonnegative, and the distinguished term is strictly positive.  The reversed
determinant is positive.  Undoing the column reversal proves `(3)`.

Nothing in this last argument used the Laguerre formula after total
nonnegativity of `(10)` and positivity of `c_(m,0)` were known.  Replacing
`m,c_(m,j)` by `M,c_j` proves the abstract PF row-transform criterion `(3a)`.

## 5. Why this repairs one zero-cut failure mode

The single reverse-Beta sequence

```text
B_n=(a+m)_n/(a)_n                                      (22)
```

has dual prefixes `(-1,0)`.  THM-3065 shows that at the integer wall it is a
degree-`m` Gregory--Newton polynomial and its Hankel rank stops at `m+1`.
Multiplying `(22)` coefficientwise by the reciprocal-Gamma sequence
`1/(a)_n` changes the terminal prefix from zero to `-1` and produces `(1)`.
The new Gamma factor does more than perturb the zero wall: equations
`(5)--(21)` identify it as the base kernel on which the same Gregory--Newton
polynomial acts by a PF row transform.  That transform restores strictness at
every order.

This explains one precise role of a strict prefix.  It does **not** prove that
every strict inventory admits such a Laguerre transform.  Section 6 shows that
the larger terminal-minus-one integer-mesh residual is still a negative-root
polynomial and hence has a Newton-PF factorization.  Terminal prefixes below
minus one or nonintegral gaps leave additional reciprocal-Gamma or
nonpolynomial Gamma-ratio factors, for which no finite Toeplitz factorization
is asserted.

That baseline is itself exact and universal.  For an arbitrary increasing
shape mesh and integer inventory

```text
H_n=product_(j=0)^N(alpha_j)_n^(e_j),
E_j=sum_(i=0)^j e_i<0,
S_j=-E_j>=1,                                           (22a)
```

summation by parts gives

```text
H_n=(alpha_N)_n^(-S_N)
    product_(j=0)^(N-1)
      [(alpha_(j+1))_n/(alpha_j)_n]^(S_j).             (22b)
```

Write `S_j=1+T_j`.  The unit flow telescopes all the way across the mesh:

```text
H_n=1/(alpha_0)_n
    * (alpha_N)_n^(-T_N)
    * product_(j=0)^(N-1)
      [(alpha_(j+1))_n/(alpha_j)_n]^(T_j),
T_j>=0.                                                (22c)
```

So integral strictness is exactly **one full reciprocal-Gamma baseline plus
a nonpositive-prefix residual**.  For `(4)`, the residual is one reverse-Beta
edge of integer length `m`; `(5)--(21)` prove that its polynomial action is a
PF filter.  Formula `(22c)` reframes the general open problem: determine when
the residual in `(22c)` acts on the baseline through a PF falling-factorial
filter, or what stronger sidecar replaces that filter.

## 6. Newton-PF closure for an integer mesh with terminal prefix minus one

The Laguerre identity is an instance of a general discrete stability lemma.

### Newton-PF lemma

Let

```text
P(x)=kappa product_(ell=1)^M(x+r_ell),
kappa>0, r_ell>0,                                     (23)

P(x)=sum_(j=0)^M c_j x^(underline j),
C(z)=sum_(j=0)^M c_jz^j.                              (24)
```

Then every `c_j>0` and every zero of `C` is negative and simple.  In
particular `(c_j)` is PF-infinity in the lower-Toeplitz orientation `(10)`.

This has a short induction which exposes the mechanism.  Suppose `C` is the
falling-coefficient polynomial of `P`.  Multiplication by `x+r`, `r>0`, uses

```text
(x+r)x^(underline j)
 =x^(underline (j+1))+(j+r)x^(underline j),           (25)
```

so the new coefficient polynomial is

```text
D(z)=(z+r)C(z)+zC'(z).                                (26)
```

The coefficient recurrence `d_j=c_(j-1)+(j+r)c_j` is strictly positive.
If the roots of `C` are `-rho_i`, `rho_i>0`, then away from those roots

```text
D(z)/C(z)=z+r+sum_i z/(z+rho_i),
d/dz[D/C]=1+sum_i rho_i/(z+rho_i)^2>0.                (27)
```

Across each pole the quotient jumps from `+infinity` on the left to
`-infinity` on the right.  Equation `(27)` therefore gives one simple zero
in every interval between consecutive poles, one to the left of all poles,
and one between the rightmost pole and zero because `(D/C)(0)=r>0`.  These
are all `deg(C)+1` zeros of `D`, and all are negative.  Starting from the
positive constant polynomial proves the lemma.

### Multi-edge strict-prefix theorem

Return to `(22a)` and suppose additionally that

```text
E_N=-1,
d_j=alpha_(j+1)-alpha_j is a positive integer.         (28)
```

Then `S_N=1`, so `T_N=0` in `(22c)`.  Every residual edge is the polynomial

```text
(alpha_(j+1))_n/(alpha_j)_n
 =(alpha_j+n)_(d_j)/(alpha_j)_(d_j).                  (29)
```

Set

```text
M=sum_(j=0)^(N-1) T_j d_j.                            (30)
```

Apart from a positive constant, the residual in `(22c)` is a degree-`M`
polynomial in `x=alpha_0+n-1`; all of its roots are negative because its
linear factors are

```text
x+(alpha_j-alpha_0)+t+1,
0<=t<d_j, each repeated T_j times.                    (31)
```

The Newton-PF lemma makes its falling coefficients a PF-infinity filter.
Apply the abstract criterion `(3a)` with `a=alpha_0`.  This proves the
**global all-minor theorem** whenever `alpha_0>M`: every order and every pair
of strictly increasing nonnegative row/column tuples obeys the checkerboard
sign law.

There is a sharper per-minor tail statement.  Given `(2)`, put
`s=u_1+v_1`, `u'_i=u_i-u_1`, and `v'_j=v_j-v_1`.  Factoring the positive
common value `H_s` from every entry replaces every shape by `alpha_j+s` and
leaves the integral gaps, the `T_j`, the roots `(31)`, and the degree `M`
unchanged.  The abstract criterion therefore applies under the precise strip
condition

```text
alpha_0+u_1+v_1>M.                                   (32)
```

Thus `(28)--(32)` also prove the **eventual-offset corollary** for every
`alpha_0>0`: every order and arbitrary row/column gaps obey the checkerboard
sign throughout the tail `u_1+v_1>M-alpha_0`.  This is distinct from the
global all-minor theorem above, which needs `alpha_0>M`.  The two-shape
inventory `(4)` has `T=(1,0)` and `M=m`, so `(32)` is precisely `(4a)`.

## 7. Boundary and open extension

The conditions `a+u_1+v_1>m` and `(32)` are sufficient positive-strip
conditions, not asserted failure boundaries.  They ensure that every shifted
row in the PF transform lies in `z>-a`, where `(14)` is strictly sign-regular
term by term.  Below that strip, individual shifted Gamma rows can meet zeros
or leave the positive chamber before their PF sum cancels them.  No
continuation across that boundary is claimed without a nonvanishing argument.

The strip condition cannot simply be deleted from the abstract criterion.
Take `a=1/2` and the PF filter `C(z)=1+3z` (`M=1`).  At the first base,

```text
F_0=r_0+3r_(-1)=1+3(a-1)=-1/2,                        (34)
```

although `F_n>0` for `1<=n<=6`.  Thus even the order-one sign fails below the
positive strip; this is a boundary hostile, not merely a proof artifact.

Likewise, the following stronger statement remains **OPEN**:

```text
all cumulative integer prefixes E_j<0
  => every generalized Hankel minor has sign (-1)^(r choose 2).       (33)
```

The theorem proves `(33)` on the terminal-minus-one integer-mesh tail `(32)`.
It does not cover terminal exponent `E_N<=-2`, nonintegral shape gaps, the
finite shallow strip `alpha_0+u_1+v_1<=M`, a generic Hadamard-product closure,
or any new physical GMC/SFC/NC2 consequence.

## 8. Exact evidence

The dependency-free exact companion checks:

- `363` Laguerre/falling-factorial equalities for gaps `0..10`;
- `223,398` finite Toeplitz minors in the orientation `(10)`, with no negative
  values;
- `342` exact row-transform identities `(17)`;
- `36` literal nonconsecutive matrix products and complete Cauchy--Binet
  expansions comprising `6,372` rational terms, including the distinguished
  strict term;
- `6,444` direct generalized Hankel determinants of orders `2..6`, computed
  without calling the row-transform determinant path;
- `30` non-Laguerre generalized-minor controls for three generic PF filters
  presented directly as products of positive bidiagonal factors;
- `225,273` exact Newton-expansion and Toeplitz-minor controls for four
  multi-edge integer meshes, together with `65` direct generalized minors;
  one mesh has `alpha_0<=M` and is tested only after the tail shift `(32)`; and
- `128` additional exact order-`3..6` hostile controls on four strict-prefix
  meshes with two through five shapes, explicitly labelled finite-only.

Reproduce with

```text
python 04-computation/gmc_laguerre_pf_strict_reciprocal_beta_gap_thm3079.py
python -O 04-computation/gmc_laguerre_pf_strict_reciprocal_beta_gap_thm3079.py
```

Both modes byte-match the stored twelve-line transcript.  The LF-normalized
script has `18,568` bytes in `459` lines and the output has `636` bytes in twelve
lines; their hashes are pinned in the frontmatter.

The independent companion does not import the primary script.  It reconstructs
five Newton families, including rational and repeated roots, and checks `130`
exact reconstruction cells, `15,018` Toeplitz minors, and `4,430` generalized
minors.  Three multi-edge meshes add `30` transform identities and `2,658`
generalized minors.  It also certifies the hostile `(34)` and records an
independent proof-direction audit.  Reproduce with

```text
python 04-computation/gmc_newton_pf_thm3079_independent_audit.py
python -O 04-computation/gmc_newton_pf_thm3079_independent_audit.py
```

Both modes byte-match the stored eight-line transcript; all independent hashes
are pinned in the frontmatter.

**QED.**
