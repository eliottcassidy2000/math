---
id: THM-2982
title: "First-gap wall-stripped norm-core strict ULC through width thirty-four"
status: PROVED + VERIFIED-EXACT + INDEPENDENTLY HOSTILE-AUDITED
source: codex-gmc-first-gap-ulc-m6-m34-2026-07-30
depends_on:
  - THM-2969-first-gap-wall-stripped-resultant-norm-core-atlas
  - THM-2973-first-gap-wall-stripped-norm-core-continuation-through-thirty-one
  - THM-2978-first-gap-wall-stripped-norm-core-at-thirty-two
related:
  - THM-2927-general-width-flagged-macaulay-leading-coefficient
  - THM-2964-general-pure-factorial-moment-resonance-ladder
  - THM-2942-macaulay-extraneous-flag-factor-and-pluecker-mutation
  - THM-2986-first-gap-transfer-one-checkerboard-sign
  - THM-2988-first-gap-self-curvature-negativity
  - THM-2991-pf-infinity-arbitrarily-delayed-newton-ratio-return
script:
  - 04-computation/gmc_first_gap_wall_stripped_norm_core_strict_ulc_through_thirty_four_thm2982.py
  - 04-computation/gmc_first_gap_structural_top_jet_sidecar_thm2982.py
output:
  - 05-knowledge/results/gmc_first_gap_wall_stripped_norm_core_strict_ulc_through_thirty_four_thm2982.out
  - 05-knowledge/results/gmc_first_gap_structural_top_jet_sidecar_thm2982.out
script_sha256:
  - 645353f04d2143f91b7b213e7223761c65c72a81c0ebb522184643fdd97ca24b
  - ca8d03d4ff3b647797623346eee285374b1595942ae66295a19b050ac8bf4547
output_sha256:
  - 305378242e0427a7eacd3e3292ab183b3cc9174b2dae515c623489924213fc9d
  - b23b5cf27d41fb4a20d841f9bfe65f7d983df1c29e8819e34bbd7ee86bc63fd0
hash_basis: LF-normalized bytes
---

# THM-2982 -- first-gap wall-stripped norm-core strict ULC through width thirty-four

**PROVED + VERIFIED-EXACT + INDEPENDENTLY HOSTILE-AUDITED.**

The finite theorem and the conditional all-width frontier are deliberately
separated.

## Statement

For an integer `6<=M<=34`, let `N_M(n)` be the primitive wall-stripped
pure-resultant core attached to normalized support

```text
{0,1,2,M}.                                                (1)
```

For `M<=26` use PROVED THM-2969, for `27<=M<=31` use PROVED
THM-2973, and for `M=32` use PROVED THM-2978.  At `M=33,34` apply the same
two frozen THM-2969 Macaulay charts and exact flag/Smith/seam division.  Write

```text
N_M(n)=sum_(i=0)^d a_(M,i)n^i.                           (2)
```

Then every coefficient is positive and, for every `1<=i<d`,

```text
a_(M,i)^2 i(d-i)
  > a_(M,i-1)a_(M,i+1)(i+1)(d-i+1).                     (3)
```

Equivalently, `a_(M,i)/binom(d,i)` is strictly log-concave.  Thus the binary
homogenization

```text
sum_(i=0)^d a_(M,i)x^i y^(d-i)                           (4)
```

is in the interior of the two-variable Lorentzian cone.  This strictly
strengthens the PF2 conclusions of THM-2969/2973/2978 on their common finite
range.

The exact global minimum of the ratio in `(3)` follows this ridge.  The index
`i` is the true ascending power of `n`, and `k=d-i` is its distance from the
leading edge.

| `M` | `d` | minimizing `i` | `k` | minimum ratio minus one |
|---:|---:|---:|---:|---:|
| 6 | 121 | 61 | 60 | `8.0622846897e-3` |
| 7 | 144 | 80 | 64 | `6.0817438656e-3` |
| 8 | 164 | 97 | 67 | `4.9140166823e-3` |
| 9 | 184 | 116 | 68 | `4.0300499732e-3` |
| 10 | 205 | 136 | 69 | `3.3227746231e-3` |
| 11 | 226 | 157 | 69 | `2.7746894403e-3` |
| 12 | 244 | 178 | 66 | `2.4222778974e-3` |
| 13 | 268 | 202 | 66 | `2.0347593827e-3` |
| 14 | 288 | 225 | 63 | `1.7921408578e-3` |
| 15 | 308 | 248 | 60 | `1.5894253025e-3` |
| 16 | 329 | 271 | 58 | `1.4191207676e-3` |
| 17 | 351 | 296 | 55 | `1.2714738048e-3` |
| 18 | 369 | 317 | 52 | `1.1724047484e-3` |
| 19 | 392 | 343 | 49 | `1.0610547527e-3` |
| 20 | 412 | 366 | 46 | `9.8351860874e-4` |
| 21 | 431 | 388 | 43 | `9.1701007159e-4` |
| 22 | 453 | 413 | 40 | `8.5025942043e-4` |
| 23 | 475 | 438 | 37 | `7.9233028798e-4` |
| 24 | 493 | 460 | 33 | `7.5169288135e-4` |
| 25 | 516 | 486 | 30 | `7.0269453801e-4` |
| 26 | 536 | 509 | 27 | `6.6758691618e-4` |
| 27 | 556 | 533 | 23 | `6.3370632597e-4` |
| 28 | 577 | 557 | 20 | `6.0255281498e-4` |
| 29 | 599 | 583 | 16 | `5.7262656351e-4` |
| 30 | 617 | 604 | 13 | `5.5101670437e-4` |
| 31 | 639 | 630 | 9 | `5.2559781199e-4` |
| 32 | 660 | 654 | 6 | `5.0413591167e-4` |
| 33 | 680 | 678 | 2 | `4.8437799398e-4` |
| 34 | 701 | 700 | 1 | `4.6598363906e-4` |

Hence `M=34` is the first width in this exact atlas at which the globally
tight Newton three-circuit is attached to the leading boundary.

The two new primitive cores are

| `M` | degree | SHA-256 of primitive coefficient tuple | wall invoice |
|---:|---:|:---|:---|
| 33 | 680 | `e86a6beb0403fa1a7cfa020813073a102d94ae7a17e4ed9c4eb673eeb7665040` | cubic wall `r=25` removed by the prescribed `c_M`; no extra Smith correction |
| 34 | 701 | `080a8e37065be58ccb1d7c6f551f3febc35c960f59d76d3b3335dfff658f0d91` | no quadratic, cubic, quartic, or matrix-sporadic correction |

The exact width-34 minimum ratio is

```text
27655303142730846948837614910195097534441250912778059699455304305433884593670175322013125
----------------------------------------------------------------------------------------- .
27642422226229352775679708685499047535273900449378640304403997574607551772198286144871066
```

In particular the first-gap moment nonvanishing conclusion already proved
through width 32 extends to widths 33 and 34.  No arbitrary-width or
arbitrary-radial conclusion is asserted.

## Proof

### 1. Frozen dependency and wall map

The companion checks the LF hashes of the immutable THM-2969, THM-2973, and
THM-2978 companions before importing the THM-2969 engine.  The inherited
roles are exact:

- THM-2969 supplies the two full Macaulay charts, flag, local Smith and seam
  factors, core reconstruction, positivity, and widths `6..26`;
- THM-2973 supplies widths `27..31` and the load-bearing quartic correction
  `(M,r)=(31,25)`;
- THM-2978 supplies the independently audited width-32 core; and
- this companion reconstructs widths 33 and 34 from the same full charts and
  audits their complete wall invoices.

At each width the worker rebuilds the full primitive core rather than reading
stored coefficient data.  It checks the canonical core digest and exact
positive coefficient tuple.  In particular, the strict inequalities below
are properties of the same cores used in the proved nonvanishing chain, not
of a fitted or truncated surrogate.

### 2. Exact strict-ULC scan and orientation repair

The `python-flint` method `fmpz_poly.coeffs()` returns coefficients in native
constant-first order.  The companion additionally checks

```text
tuple(core.coeffs()) = tuple(core[i] for i=0,...,d),       (5)
```

then evaluates both sides of `(3)` as integers at every internal index.  All
slacks are strictly positive.  The global ratio minimum, its ascending index,
and its distance from the leading edge are compared with the frozen ridge
table; an exact rational enclosure marker and a SHA-256 of the exact ratio
are included in every compact transcript row.

This orientation check is load-bearing.  A preliminary scratch transcript
had reversed the native coefficient tuple while retaining an `ascending`
header.  Reversal preserves the truth and value of the ULC inequalities, but
it replaces `i` by `d-i` and therefore falsifies the reported ridge.  The
repaired direct checks give, for example,

```text
M=17: i=296, k=55;        M=32: i=654, k=6,              (6)
```

not the reversed indices `55` and `6`.

Positive coefficients with no internal zeros and strict log concavity after
binomial normalization are exactly the binary strict-Lorentzian criterion,
which proves `(4)`.

### 3. Independent width-34 Jacobi prediction

For a degree-`d` polynomial `f`, put

```text
u(f)   =[n^(d-1)]f/[n^d]f,
ell2(f)=2[n^(d-2)]f/[n^d]f-u(f)^2.                       (7)
```

Both are additive on products.  If the selected width-34 Macaulay chart is

```text
M(n)=n^Delta(M0+n^(-1)M1+n^(-2)M2+...),                 (8)
```

Jacobi's identity gives

```text
u(det M)=tr(M0^(-1)M1),
ell2(det M)=2tr(M0^(-1)M2)-tr((M0^(-1)M1)^2).            (9)
```

The companion builds only the top three coefficient slices of the three
moment forms, hence three exact `36 x 36` rational matrices.  It subtracts
the logarithmic jets of `q^6 c K` and the explicit `W=B E/H` wall product.
If `v=ell2+u^2=2a_(d-2)/a_d`, the leading circuit ratio is

```text
((d-1)/d) u^2/v.                                         (10)
```

This independent calculation exactly equals the full-core width-34 fraction
in the statement.  It uses neither the other 699 coefficients nor the ULC
scan.  The full scan separately proves that this edge circuit is the global
minimum.

### 4. The pure-resultant Hessian is edge-local

The second trace in `(9)` also has a closed coordinate-free source.  Put the
leading response rows into coordinates in which the three pure forms are

```text
(X_2^2,X_3^3,X_4^4).                                   (11)
```

More generally, write

```text
F_j(t)=X_j^j+tG_j+t^2H_j+O(t^3),
alpha_j=[X_j^j]G_j,       beta_j=[X_j^j]H_j,
g_(jk,a)=[X_j^(j-a)X_k^a]G_j.                          (12)
```

For `{j,k,l}={2,3,4}` and `m_j=24/j`, the exact local formula is

```text
ell2(Res)=sum_j m_j(2beta_j-alpha_j^2)
 -2 sum_(j<k) l sum_(a=1)^min(j,k) a g_(jk,a)g_(kj,a). (13)
```

Indeed, torus weight permits at quadratic order only a base-monomial self
term or a pair of opposite displacements along one Newton-simplex edge.  The
coefficient of the latter is forced by the binary reduction

```text
Res(X^j+uX^(j-a)Y^a,
    Y^k+vY^(k-a)X^a,
    Z^l)
  =(1-uv)^(al).                                         (14)
```

This also gives the exact equality/failure boundary: a cross response
vanishes precisely when at least one of its two opposite edge coefficients
vanishes; no non-edge quadratic pair can contribute.  There are seven
possible edge/transfer pairs for degrees `(2,3,4)`.  In the first-gap
factorial slices, six are active and the `(3,4),a=3` pair vanishes at every
audited width.

More explicitly, let `p_k` be the response-dual columns and put

```text
b_jk=B_j dot p_k,                 c_jk=C_j dot p_k.      (14a)
```

The first subleading factorial form gives

```text
g_(jk,1)=(j-1)(2b_jj b_jk-c_jk)/2,
g_(jk,2)=(j-1)b_jk^2/2,
g_(jk,a)=0 for a>=3.                                    (14b)
```

Thus every active transfer-two opposite product is a product of nonzero
squares, while the transfer-three zero is structural rather than numerical.
For the transfer-one term, set `E=x d/dx` and let

```text
lambda_k(x)=A(x) dot p_k                             (14c)
```

be the cardinal response (`lambda_k(j)=delta_jk` at `j=2,3,4`).  Direct
differentiation rewrites `(14b)` as

```text
g_(jk,1)=-(j-1)/2 (E^2+E)(lambda_k/lambda_j)|_(x=j).    (14d)
```

On `6<=M<=50`, exact arithmetic gives the checkerboard sign
`sign g_(jk,1)=(-1)^(row(j)+row(k))`.  PROVED THM-2986 subsequently proves
that sign for every first-gap width `M>=6` by an exact cardinal-response
argument.  That all-width theorem is not used to prove the finite atlas here.

The durable structural sidecar transforms the actual first three slices to
`(11)` and matches `(13)` to the `36 x 36` Jacobi trace at every
`6<=M<=34`: `29` widths and `58` exact first/second-jet equalities.
The same sidecar also checks the closed factorial slices and the wider finite
curvature census below; its reproduction surface and hashes are frozen in
**Exact evidence**.

Formula `(13)` closes the resultant-Hessian term in the symbolic edge
program.  It does not by itself prove the wall-stripped sign: the factorial
slice coefficients and explicit Smith/seam wall jets must still be
substituted and controlled.

Those slice coefficients also admit a closed all-width reduction.  For

```text
T_j(d_1,...,d_j;n)
 =(jn+1)_(sum d_i)/product_i(n+1)_(d_i),                (15)
P1(s)=s(s+1)/2,              P2(s)=s(s+1)(2s+1)/6,
```

put `S=sum d_i`.  The first two logarithmic coefficients at infinity are

```text
lambda1=P1(S)/j-sum_i P1(d_i),
lambda2=-P2(S)/(2j^2)+sum_i P2(d_i)/2,                  (16)
```

so

```text
T_j=j^S(1+lambda1/n+(lambda1^2/2+lambda2)/n^2+...).    (17)
```

The common clearing denominator

```text
D_j=product_(r=1)^(M-1)(n+r)^(j-1)(n+M)^(j-2)          (18)
```

has the same explicit `P1/P2` log jets.  Taking the signed inclusion
`d_i in {0,1,2,M}` and multiplying by `(18)` therefore gives all three form
slices using only finite sums of polynomials in `M` times powers of
`2^M,3^M,4^M`.

The durable structural sidecar checks all `31` form coefficients at all
`29` widths, hence `899` coefficient records and `2,697` scalar top-slice
equalities.

Finally the response-coordinate determinant simplifies to

```text
-2(3*2^M-3*3^M+4^M-1).                                 (19)
```

Combining `(13)`, `(16)--(19)`, and the finite wall sums reduces the `k=1`
strict-ULC numerator, after clearing a known positive denominator, to an
explicit rational-exponential expression in
`M,2^M,3^M,4^M`.  On each residue class modulo `30`, the floor and quartic-
wall indicators are fixed polynomial data.  This is an exact reduction of
the arbitrary-width edge sign problem, not a proof that its numerator is
positive.

### 5. Resultant-to-moment consequence

For widths 33 and 34, the exact core and wall invoices retain the proved
factorization

```text
R_M = W_M N_M                                             (20)
```

up to a nonzero rational unit.  Every coefficient of `N_M` is positive and
`W_M(n)>0` for `n>=0`, so `R_M(n)!=0`.  The denominator-cleared quadratic,
cubic, and quartic factorial-moment forms therefore have no common
projective zero.  Restoring the eliminated mean gives the same first-four-
moment nonvanishing implication as THM-2969/2973/2978.

## Exact evidence

The intended canonical reproduction is

```text
python 04-computation/gmc_first_gap_wall_stripped_norm_core_strict_ulc_through_thirty_four_thm2982.py --workers 4 --output .scratch/thm2982.normal.out
python -O 04-computation/gmc_first_gap_wall_stripped_norm_core_strict_ulc_through_thirty_four_thm2982.py --workers 4 --output .scratch/thm2982.opt.out
```

The compact atlas record digest is

```text
9a496cb20febf0e88b039afcca56540ff0fe352f11aeaa79b965dd7922237e8e.
```

The compact normal and optimized transcripts are byte-identical (`9,229`
UTF-8/LF bytes).  Their frozen LF hashes are

```text
script  645353f04d2143f91b7b213e7223761c65c72a81c0ebb522184643fdd97ca24b
output  305378242e0427a7eacd3e3292ab183b3cc9174b2dae515c623489924213fc9d
```

The single durable structural sidecar is reproduced by

```text
python 04-computation/gmc_first_gap_structural_top_jet_sidecar_thm2982.py --output .scratch/thm2982-structural.normal.out
python -O 04-computation/gmc_first_gap_structural_top_jet_sidecar_thm2982.py --output .scratch/thm2982-structural.opt.out
```

Its normal and optimized transcripts are LF-identical (`8,869` bytes), with
record digest
`c3b116e42666d32eb6712091413c550db6d0ed920cb9d513b439ff1a0cbc2b5e`
and frozen LF hashes

```text
script  ca8d03d4ff3b647797623346eee285374b1595942ae66295a19b050ac8bf4547
output  b23b5cf27d41fb4a20d841f9bfe65f7d983df1c29e8819e34bbd7ee86bc63fd0
```

The independent hostile audit replayed both companions normally and under
`-O` against the stored outputs.  It rederived all `11,909` internal strict
ULC circuits and checked coefficient orientation directly.  A separate
Macaulay reconstruction of widths `33,34` reproduced both primitive-core
digests, the `679+700` new strict circuits, the width-34 edge fraction, and
the exact wall census: width 33 has only the prescribed cubic root `r=25`,
width 34 has no quadratic, cubic, or quartic root, and neither core retains
an integer factor `(n+r)` for `1<=r<=M`.  The structural formulas, finite
`M<=50` scope, hashes, and conditional boundary were also checked
independently.

## Conditional all-width frontier -- not a proved dependency

The sharp local sign target is now a curvature-versus-wall invoice.  In the
pure-resultant formula `(13)`, an exact finite census on `6<=M<=50` finds all
`135` self curvatures negative, all `270` active opposite-edge products
positive (and hence entering `(13)` negatively), and exactly `45` structural
zeros, always `(3,4),a=3`.  Thus the pure resultant `ell2` is negative at all
45 tested widths.  After division by `W`, the only positive contribution to
the core `ell2` is the explicit wall square-sum in Section 3.

The transfer-two positivity and transfer-three zero have the algebraic proof
`(14b)`.  PROVED THM-2986 proves the uniform checkerboard sign for transfer
one, and PROVED THM-2988 proves the uniform negativity of all three self
curvatures.  Together with `(13)--(14b)`, those later theorems show that the
pure-resultant `ell2` is negative for every first-gap width `M>=6`.  They do
not control the positive square-sum introduced by wall division.

The `M=6..50` sign census in the durable structural transcript remains an
independent finite check of those later all-width theorems.  Because this
companion proves no all-width sign statement itself, its frozen transcript
correctly prints `all_width_sign=NOT_ASSERTED`.

The same top-three-slice/Jacobi calculation is exact as an algebraic
calculation at any chosen width.  Outside the fully reconstructed finite
atlas, however, any edge-sign conclusion is **CONDITIONAL** on all of the
following:

1. continuation of the THM-2942 chart/resultant degree and flag identity;
2. continuation of the baseline local Smith and seam invoices;
3. only the THM-2964 quartic ladder `M=10t+1,r=8t+1` and the known matrix
   sporadic `(12,5)` adding correction walls; and
4. absence of a new matrix-only or seam correction.

Even under those hypotheses, edge positivity is not the whole ULC theorem.
An all-width result still needs a **no-return theorem**: after the active
Newton circuit reaches `k=1` at width 34, every deeper circuit must stay no
tighter than the edge circuit.  In each residue class modulo 30, clearing the
positive generalized-Vandermonde denominators turns the edge trace into a
rational-exponential expression in `M,2^M,3^M,4^M`.  Dominant-base separation
plus a finite prefix is a concrete next certificate; it has not yet been
proved.  PROVED THM-2991 shows why this sidecar must be family-specific:
PF-infinity, Hurwitz stability, strict ULC, and an arbitrarily long improving
leading prefix can coexist with a later Newton-ratio return.

The strict ULC atlas therefore proves no arbitrary-width first-gap theorem,
no arbitrary-radial GMC(2), and no NC2 closure.

**QED.**
