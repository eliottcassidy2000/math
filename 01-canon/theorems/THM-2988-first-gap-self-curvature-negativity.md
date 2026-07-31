---
id: THM-2988
title: "First-gap self-curvature negativity"
status: PROVED / VERIFIED-EXACT / INDEPENDENTLY HOSTILE-AUDITED
source: codex-gmc-self-curvature-2026-07-30
depends_on:
  - THM-2927-general-width-flagged-macaulay-leading-coefficient
related:
  - THM-2986-first-gap-transfer-one-checkerboard-sign
  - THM-2982-first-gap-wall-stripped-norm-core-strict-ulc-through-thirty-four
script: 04-computation/gmc_first_gap_self_curvature_negativity_thm2988.py
output: 05-knowledge/results/gmc_first_gap_self_curvature_negativity_thm2988.out
script_sha256: 6fc6ccef77e95525961aff6c2d0d69acd9cbe79fcbdd0073fa4337b954d9924b
output_sha256: d20fd90b98a95922fb8fc92d9d506caec13f4bcb01dca05672f6f6e2092ca83e
hash_basis: LF-normalized bytes
---

# THM-2988 -- first-gap self-curvature negativity

**PROVED / VERIFIED-EXACT / INDEPENDENTLY HOSTILE-AUDITED.**

This theorem closes the three self terms in the local pure-resultant
log-Hessian.  It does not control the wall square-sum, prove edge strict ULC
after wall division, prove interior no-return, handle arbitrary radial
coefficients, or prove GMC(2).

## Statement

For a first-gap width M>=6, put

    A_M(x)=(1-x^M, x-x^M, x^2-x^M)

and let R_M be the response matrix with rows A_M(2), A_M(3), A_M(4).
Transform the first three factorial slices of the degree-j form, for
j in {2,3,4}, into response coordinates in which its leading form is X_j^j:

    F_j=X_j^j+n^(-1)G_j+n^(-2)H_j+O(n^(-3)).

Define

    alpha_j=[X_j^j]G_j,
    beta_j =[X_j^j]H_j,
    kappa_j=2beta_j-alpha_j^2.                           (1)

Then, for every integer M>=6,

    kappa_2<0,  kappa_3<0,  kappa_4<0.                  (2)

All three inequalities are strict.  There is no equality width in the
first-gap range.

## Proof

### 1. Sparse adjugate formula

Write U=2^M and V=3^M, so 4^M=U^2.  Directly,

    det R_M=-2(U^2+3U-3V-1).                            (3)

This determinant is nonzero for M>=6.  Let v_j be the j-th column of
adj(R_M).  If F_j^(r) denotes the r-th factorial slice before the response
change, homogeneity gives

    alpha_j=F_j^(1)(v_j)/(det R_M)^j,
    beta_j =F_j^(2)(v_j)/(det R_M)^j.                   (4)

Consequently,

    kappa_j=C_j/(det R_M)^(2j),
    C_j=2F_j^(2)(v_j)(det R_M)^j-F_j^(1)(v_j)^2.        (5)

The denominator in (5) is positive.  It remains to prove -C_j>0.

The slices in (4) are fully explicit.  For directions d_1,...,d_j drawn
from {0,1,2,M}, put S=sum d_i and

    T_j=(jn+1)_(S)/product_i (n+1)_(d_i).

With

    P1(s)=s(s+1)/2,
    P2(s)=s(s+1)(2s+1)/6,

the first two logarithmic coefficients of T_j are

    lambda1=P1(S)/j-sum_i P1(d_i),
    lambda2=-P2(S)/(2j^2)+sum_i P2(d_i)/2.              (6)

Thus

    T_j=j^S(1+lambda1/n+(lambda1^2/2+lambda2)/n^2+...). (7)

Signed inclusion over d_i in {0,1,2,M}, followed by the common clearing
denominator, constructs F_j^(1),F_j^(2) exactly.  Equations (3)--(7)
therefore construct each C_j as a finite polynomial in M,U,V without a
symbolic resultant or a dense transformed determinant.

### 2. Grouped exponential form

Collect the negative numerator as

    -C_j=sum_(a,b) q_(j,a,b)(M)(2^a 3^b)^M.             (8)

The exact grouped representations are frozen by their SHA-256 digests in the
table below.  The number of expanded monomials and the number of distinct
exponential bases are also shown, so no hidden truncation occurs.

| j | expanded terms | groups | grouped-representation digest |
|:---:|---:|---:|:---|
| 2 | 110 | 25 | 9d82d202148eb5755f8fca6d0a33b234f1bcd41a76ec329450e2bda9c7b5321c |
| 3 | 210 | 49 | ca4e661821f3e31ad6876be528b768fc4d77ffe20ee310125aee3a0777d9e49c |
| 4 | 384 | 81 | b6e21bdae6499779cdedc4af15612d2ed33120fb69a8371d7042bde3d7b6249f |

For each j the largest base B and runner-up base B' have B'/B=3/4.
The complete dominant and lower-term invoices are

| j | dominant base B | dominant coefficient p_j(M) | B' | negative dominant tariff | absolute lower tariff | tail N |
|:---:|---:|:---|---:|---:|---:|---:|
| 2 | 256 | 8(2M^3-3M^2+M+144)/3 | 192 | 8 | 1483244/3 | 54 |
| 3 | 4096 | 4(32M^3+16M+1467)/3 | 3072 | 0 | 1011680732 | 74 |
| 4 | 65536 | 32M(8M^2-8M-1) | 49152 | 288 | 201611683520 | 87 |

Every dominant polynomial has degree 3, while every lower coefficient has
degree at most 4.  If L_j is the dominant leading coefficient and A_j the
negative dominant tariff, then for M>=N

    p_j(M)/M^3 >= L_j-A_j/M=:f_j.                       (9)

At the displayed tail starts,

    f_2=140/27,  f_3=128/3,  f_4=7328/29.              (10)

If T_j is the absolute lower tariff, (8)--(10) give

    (-C_j)/(B^M M^3)
      >= f_j-T_j M(3/4)^M.                              (11)

The exact positive margins after division by f_j are

| j | 1-(T_j/f_j)N(3/4)^N |
|:---:|:---|
| 2 | 62767441708390834036983714788449 / 811296384146066816957890051440640 |
| 3 | 16801722156309940837727838130657657862796799 / 5708990770823839524233143877797980545530986496 |
| 4 | 172505339160395148182314587285988197943367611485671213 / 2741730303580379285656730228261100023817900460781600768 |

All are positive.  Moreover ((M+1)/M)(3/4)<1 from every displayed N,
so the lower majorant in (11) decreases thereafter, while the dominant floor
in (9) increases.  Hence -C_j>0 for every M>=N.

The exact companion evaluates (8) without rounding on each remaining width
6<=M<N.  The prefix lengths are respectively 48, 68, and 81; the exact value
lists have digests

    j=2: 419708bf803054e3fce4a437bee828347013c57d9694327acc1e08e5548b108c
    j=3: 214f438938c26dc2ce7377c8fab427a364fc010f2570bfc47dca4d0c1fa97f7f
    j=4: 3436040bea2b49ccd44c401a6583d8bdc306c3cd9527cb6a7831baa6accdd4f5

Every value is strictly positive.  Thus -C_j>0 for every M>=6.  The positive
denominator in (5) proves (2).

QED.

## Pure-resultant consequence and remaining boundary

For the pure-power resultant, the exact local second logarithmic jet is

    ell2(Res)
      =sum_j (24/j)kappa_j
       -2 sum_(j<k) l sum_a a g_(jk,a)g_(kj,a),          (12)

where l is the remaining degree.  THM-2988 makes every self term in (12)
strictly negative.  The separate THM-2986 candidate proves all transfer-one
opposite products positive for every M>=6; the transfer-two products are
direct products of squares, and the only transfer-three pair vanishes
structurally.  Once both candidates are independently promoted, (12) gives

    ell2(Res)<0 for every first-gap width M>=6.           (13)

Statement (13) is a composition consequence, not a dependency used in this
proof.  It is still not the wall-stripped ULC edge: division by the explicit
wall product adds a positive square-sum to ell2, and strict ULC also charges
the degree-normalized first log jet.  Interior no-return remains separate.

## Exact evidence

Reproduce with

    python 04-computation/gmc_first_gap_self_curvature_negativity_thm2988.py --output .scratch/thm2988.normal.out
    python -O 04-computation/gmc_first_gap_self_curvature_negativity_thm2988.py --output .scratch/thm2988.opt.out

The normal and optimized transcripts are LF-identical (2,542 bytes).  The
companion uses optimization-stable explicit guards.  Besides constructing
and certifying all grouped expressions, it independently compares (5) with
the directly transformed factorial slices at every M=6..50: 135 exact
equalities, with digest

    ed2dc44be0eea5701584da1c025a64c3766e4f098135c9dcc1acf9a388e06fbd.

The three-record digest is

    48b18cb73843dc63b7c48d4929e2290e3641713b37d127dbfcd659abc56e9f0e.

Frozen LF hashes are

    script  6fc6ccef77e95525961aff6c2d0d69acd9cbe79fcbdd0073fa4337b954d9924b
    output  d20fd90b98a95922fb8fc92d9d506caec13f4bcb01dca05672f6f6e2092ca83e

An independent hostile audit rederived the factorial-slice and monic-clearing
normalizations, the adjugate denominator and its sign, and the exhaustive
grouping by the distinct bases `2^a 3^b`.  It independently transformed the
actual factorial forms at eight widths through `M=100`, obtaining 24 further
exact negative curvature equalities, and checked the tariffs, `3/4` tail
monotonicity, all finite-prefix counts, and exact margins.  Normal, optimized,
and stored transcripts and their hashes matched.  No load-bearing defect
remained.
