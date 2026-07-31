---
id: THM-2986
title: "First-gap transfer-one checkerboard sign"
status: PROVED / VERIFIED-EXACT / INDEPENDENTLY HOSTILE-AUDITED
source: codex-gmc-transfer-one-checkerboard-2026-07-30
depends_on:
  - THM-2927-general-width-flagged-macaulay-leading-coefficient
related:
  - THM-2982-first-gap-wall-stripped-norm-core-strict-ulc-through-thirty-four
  - THM-2942-macaulay-extraneous-flag-factor-and-pluecker-mutation
script: 04-computation/gmc_first_gap_transfer_one_checkerboard_sign_thm2986.py
output: 05-knowledge/results/gmc_first_gap_transfer_one_checkerboard_sign_thm2986.out
script_sha256: f2e333bf4b6dc07bbd763378bfdb598edb5b0836b6e4279f4221a52cf0fc6b40
output_sha256: 8507ddd6bb9a0579c3ca7705ee1ce418155d3227836d9acac24c784af65f28b4
hash_basis: LF-normalized bytes
---

# THM-2986 -- first-gap transfer-one checkerboard sign

**PROVED / VERIFIED-EXACT / INDEPENDENTLY HOSTILE-AUDITED.**

This theorem closes one structural sign debt exposed by THM-2982.  It does
not prove the remaining self-curvature sign, the wall invoice, no-return, an
all-width strict-ULC theorem, arbitrary radial coefficients, or GMC(2).

## Statement

For an integer M, put

    A_M(x)=(1-x^M, x-x^M, x^2-x^M)

and let R_M be the 3 by 3 matrix whose rows are A_M(2), A_M(3), A_M(4).
For k in {2,3,4}, let p_k be the corresponding column of R_M^{-1} and define
the cardinal response

    lambda_k(x)=A_M(x) p_k,
    lambda_k(j)=delta_jk  (j=2,3,4).

Let E=x d/dx and, for j unequal to k, set

    S_jk=(E^2+E)(lambda_k/lambda_j)|_(x=j).              (1)

For every integer M at least 3,

    sign S_23=+,  sign S_24=-,
    sign S_32=+,  sign S_34=+,
    sign S_42=-,  sign S_43=+.                           (2)

In the first-gap factorial resultant, let g_(jk,1) be the transfer-one
coefficient from the j-form toward the k-response coordinate.  Then

    g_(jk,1)=-(j-1)S_jk/2.                               (3)

Consequently, for every first-gap width M at least 6,

    sign g_(jk,1)=(-1)^(row(j)+row(k)),                  (4)

where row(2)=0, row(3)=1, row(4)=2.  In particular each of the three
opposite transfer-one products g_(jk,1)g_(kj,1) is strictly positive.

The lower boundary is sharp for this cardinal model: R_M is singular at
M=1,2.  It is nonsingular and all six signs in (2) are strict from M=3
onward.

## Proof

### 1. Cardinal reduction

Direct expansion gives

    det R_M=-2D_M,
    D_M=3*2^M-3*3^M+4^M-1.                              (5)

At M=1,2, D_M=0.  At M=3, D_M=6.  For M at least 4,
4^M>3*3^M, so D_M>3*2^M-1>0.

Let B_j and C_j be the first two subleading response rows and put

    b_jk=B_j p_k,  c_jk=C_j p_k.

The explicit first subleading factorial form gives

    g_(jk,1)=(j-1)(2b_jj b_jk-c_jk)/2.                  (6)

Differentiate A_M(x)p_k twice, use lambda_j(j)=1 and lambda_k(j)=0,
and apply the quotient rule.  The result is

    2b_jj b_jk-c_jk
      =-(E^2+E)(lambda_k/lambda_j)|_(x=j),               (7)

which proves (3).  This calculation uses only the first three factorial
slices and the cardinal normalization, not a full Macaulay determinant.

### 2. Six exact rational-exponential numerators

Write S_23=P_23/D_M^2, S_24=-P_24/D_M^2,
S_32=P_32/(2D_M^2), S_34=P_34/(2D_M^2),
S_42=-P_42/D_M^2, and S_43=P_43/D_M^2.  Inverting the
three by three response matrix and differentiating gives

    P_23 =
      6*16^M -12*12^M -3(M-1)(M+4)*8^M
      +9(M^2-M+4)*6^M +3(3M^2-9M-10)*4^M
      -24*3^M +3(M^2+15M-4)*2^M +24;                   (8a)

    P_24 =
      6*12^M -12*9^M -(M^2+M+8)*8^M
      +3(M^2-3M+14)*6^M +(3M^2-3M-22)*4^M
      -24*3^M +(M^2+13M+14)*2^M +4;                    (8b)

    P_32 =
      48*16^M +3(2M^2-40M+57)*12^M
      +18(M^2+2M-12)*9^M -288*8^M
      +9(2M^2+8M+33)*6^M +21*4^M
      -6(M^2-2M+6)*3^M -9*2^M +12;                     (8c)

    P_34 =
      (2M^2-16M+45)*12^M +6(M^2-10M+18)*9^M
      -48*8^M +3(2M^2+32M-123)*6^M +291*4^M
      -2(M^2+10M-54)*3^M -159*2^M +24;                 (8d)

    P_42 =
      (3M^2-35M+64)*16^M +(9M^2-15M+44)*12^M
      -108*9^M +(-9M^2+87M-320)*8^M +324*6^M
      +(3M^2-37M+148)*4^M -152*3^M -4*2^M +4;          (8e)

    P_43 =
      (3M^2-23M+36)*16^M +(9M^2-51M+92)*12^M
      +(-9M^2+123M-296)*8^M -108*6^M
      +(3M^2-49M+456)*4^M +16*3^M -244*2^M +48.        (8f)

Thus it remains to show P_jk>0.

### 3. Dominant-base tail certificate

For q(M)=aM^2+bM+c and M at least 1, use

    |q(M)| <= (|a|+|b|+|c|)M^2.                         (9)

For the first three rows below, the dominant coefficient is constant and
M^2(beta/B)^M decreases from the displayed threshold.  For the last three,
the dominant coefficient divided by M^2 increases from the displayed
threshold, while every (beta/B)^M decreases.

| pair | dominant term | tail N | dominant floor | exact positive margin |
|:---:|:---|---:|:---|:---|
| 23 | 6*16^M | 25 | 6 | 14081865922424343364772654133 / 39614081257132168796771975168 |
| 24 | 6*12^M | 25 | 6 | 1423058678816606772702161 / 13249474533898474022240256 |
| 32 | 48*16^M | 30 | 48 | 1334504361243610627652012697506433 / 5192296858534827628530496329220096 |
| 34 | (2M^2-16M+45)12^M | 18 | 5/4 | 86807696439009361 / 369768517790072832 |
| 42 | (3M^2-35M+64)16^M | 15 | 214/225 | 333379720232789933 / 32425917317067571200 |
| 43 | (3M^2-23M+36)16^M | 16 | 109/64 | 49905085875636783 / 288230376151711744 |

Here the margin is the dominant floor minus the sum of the coefficient-norm
bounds (9), after division by the dominant exponential and, in the last
three rows, by M^2.  All six margins are strictly positive.

For the constant rows, the worst base ratio is 3/4 and

    ((N+1)/N)^2*(3/4)<1

at each displayed threshold.  For the polynomial rows, if
p(M)=aM^2+bM+c is the dominant coefficient, then p(M)/M^2 increases once

    (-b)M(M+1)>c(2M+1),                                 (10)

and its forward difference continues increasing because
(-b)(M+1)>c.  The exact companion verifies these inequalities and the six
displayed margins.  It also checks P_jk>0 exactly on every finite prefix
3<=M<N.  Therefore all six P_jk are positive for every M>=3.  Equations
(5) and (8) prove (2), and (3) proves (4).

QED.

## Consequence and exact boundary

In the pure-power resultant second logarithmic jet, the transfer-one
contribution for an unordered pair {j,k}, with l the remaining degree, is

    -2l g_(jk,1)g_(kj,1).                               (11)

THM-2986 proves that (11) is strictly negative for every first-gap width.
Together with the direct transfer-two square formula and the structural
transfer-three zero, this settles every cross-edge sign in that local
Hessian.

It does not by itself settle the three self terms
`2beta_j-alpha_j^2`; THM-2988 proves their uniform negativity separately.
Even the resulting negative pure-resultant curvature still requires the
positive wall square-sum invoice and an interior no-return theorem before
yielding all-width strict ULC.

## Exact evidence

Reproduce with

    python 04-computation/gmc_first_gap_transfer_one_checkerboard_sign_thm2986.py --output .scratch/thm2986.normal.out
    python -O 04-computation/gmc_first_gap_transfer_one_checkerboard_sign_thm2986.py --output .scratch/thm2986.opt.out

The normal and optimized transcripts are LF-identical (938 bytes).  The
companion uses explicit optimization-stable guards, checks all six closed
formulas against independent rational response inversion for M=3..50
(288 equalities), verifies every finite prefix and tail certificate, and
prints record digest

    f91325c10019e93b4e55e17a4bafbb99fb2a646e8b087e6be4ef8fd5a67e1772.

Frozen LF hashes are

    script  f2e333bf4b6dc07bbd763378bfdb598edb5b0836b6e4279f4221a52cf0fc6b40
    output  8507ddd6bb9a0579c3ca7705ee1ce418155d3227836d9acac24c784af65f28b4

An independent hostile audit symbolically reinverted the generic response
matrix, reduced all six cleared numerators modulo `4^M-(2^M)^2`, and obtained
the displayed signed table.  A second factorial-slice extraction checked 66
transfer identities at eleven widths through `M=73`.  It also replayed the
finite prefixes, tail monotonicity and exact margins, scanned positivity
through `M=500`, and reproduced the normal, optimized and stored transcript
hashes.  No load-bearing defect remained.
