---
id: THM-3075
title: "Laguerre PF row transform and strict reciprocal-Beta integer gap"
status: >
  PROVED CANDIDATE + VERIFIED-EXACT; INDEPENDENT IMMUTABLE-FILE AUDIT
  REQUESTED.  If m>=1 is an integer and a>m, then every strict generalized
  Hankel minor of H_n=(a+m)_n/(a)_n^2 has sign
  (-1)^(r choose 2).  This is a genuinely mixed strict-dual-prefix inventory
  (-2,+1), with cumulative prefixes (-2,-1).  The mechanism is an exact
  Laguerre/PF-infinity row transform of the reciprocal-Gamma kernel: the
  falling-factorial coefficients of (x+1)_m have generating polynomial
  m!L_m(-z), hence a totally nonnegative Toeplitz matrix.  Cauchy--Binet and
  one distinguished triangular minor give strictness at every order and for
  arbitrary row/column gaps.  More generally, every finite PF-infinity row
  filter with positive constant coefficient preserves this reversed strict
  total positivity on the matching positive Gamma strip.  The bound a>m is a
  proved positive-strip hypothesis, not claimed sharp; the general strict-
  dual-prefix conjecture remains open.
source: codex-gmc-strict-prefix-hostile-2026-08-01
depends_on:
  - THM-3056-product-gamma-reciprocal-hypergeometric-and-hankel-reversal
related:
  - THM-3053-beta-gamma-prefix-transport-and-multiplicative-holotopy-cone
  - THM-3065-reciprocal-beta-gap-contiguous-hankel-wall
script: 04-computation/gmc_laguerre_pf_strict_reciprocal_beta_gap_thm3075.py
output: 05-knowledge/results/gmc_laguerre_pf_strict_reciprocal_beta_gap_thm3075.out
script_sha256: 18b72d0ba72fe8fecbe2035b0ca18e292d26065a0f2ac701ee5867b03ace3059
output_sha256: 2a7f7906e9bf2802e5a27258777493cbf02ec2c652ccee0f4a0ed6373a848c96
hash_basis: LF-normalized bytes
---

# THM-3075 -- a Laguerre row transform closes one strict reciprocal-Beta cone

**PROVED CANDIDATE + VERIFIED-EXACT; INDEPENDENT AUDIT REQUESTED.**

THM-3065 proves that nonpositive dual prefixes force the checkerboard sign at
order two but not at higher order: a zero prefix permits exact order-three and
order-four hostiles.  Its strict-prefix repair survived the finite banks but
had no theorem beyond order two.

The first mixed all-order family is now closed.  The integer gap is not handled
by a generic Hadamard-product closure.  Instead, its Gregory--Newton polynomial
is a Laguerre Pólya-frequency filter acting on the row coordinate of the
reciprocal-Gamma kernel.

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
every strict inventory admits such a Laguerre transform.  With more shapes,
the mandatory reciprocal-Gamma baseline remains, but the residual positive
polynomial or nonintegral Gamma ratio need not have the Toeplitz factorization
`(9)`.

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

## 6. Boundary and open extension

The hypothesis `a>m` is a sufficient positive-strip condition, not an asserted
failure boundary.  It ensures that every shifted row in `(18)` lies in
`z>-a`, where `(14)` is strictly sign-regular term by term.  Below that strip,
individual shifted Gamma rows can meet zeros or leave the positive chamber
before the sum in `(17)` cancels them.  No continuation in `a` is claimed
without a nonvanishing argument.

Likewise, the following stronger statement remains **OPEN**:

```text
all cumulative integer prefixes E_j<0
  => every generalized Hankel minor has sign (-1)^(r choose 2).       (23)
```

The theorem proves `(23)` only for the two-shape inventory `(4)`, integral
shape gap, and `a>m`.  It does not prove a Hadamard-product closure, a
nonintegral-gap result, an arbitrary-mesh result, or any new physical
GMC/SFC/NC2 consequence.

## 7. Exact evidence

The dependency-free exact companion checks:

- `363` Laguerre/falling-factorial equalities for gaps `0..10`;
- `223,398` finite Toeplitz minors in the orientation `(10)`, with no negative
  values;
- `342` exact row-transform identities `(17)`;
- `36` literal nonconsecutive matrix products and complete Cauchy--Binet
  expansions comprising `6,372` rational terms, including the distinguished
  strict term;
- `6,444` direct generalized Hankel determinants of orders `2..6`, computed
  without calling the row-transform determinant path; and
- `30` non-Laguerre generalized-minor controls for three generic PF filters
  presented directly as products of positive bidiagonal factors; and
- `128` additional exact order-`3..6` hostile controls on four strict-prefix
  meshes with two through five shapes, explicitly labelled finite-only.

Reproduce with

```text
python 04-computation/gmc_laguerre_pf_strict_reciprocal_beta_gap_thm3075.py
python -O 04-computation/gmc_laguerre_pf_strict_reciprocal_beta_gap_thm3075.py
```

Both modes byte-match the stored eleven-line transcript.  The LF-normalized
script has `13,726` bytes in `358` lines and the output has `562` bytes in eleven
lines; their hashes are pinned in the frontmatter.

**QED, pending independent immutable-file audit and status promotion.**
