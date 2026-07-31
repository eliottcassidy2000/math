---
id: THM-2997
title: "First-gap wall-stripped all-width second-edge circuit positivity"
status: PROVED CANDIDATE + VERIFIED-EXACT / AWAITING INDEPENDENT HOSTILE AUDIT
source: codex-gmc-second-edge-circuit-2026-07-31
depends_on:
  - THM-2969-first-gap-wall-stripped-resultant-norm-core-atlas
  - THM-2973-first-gap-wall-stripped-norm-core-continuation-through-thirty-one
  - THM-2978-first-gap-wall-stripped-norm-core-at-thirty-two
  - THM-2982-first-gap-wall-stripped-norm-core-strict-ulc-through-thirty-four
  - THM-2989-first-gap-wall-stripped-all-width-leading-edge-positivity
related:
  - THM-2986-first-gap-transfer-one-checkerboard-sign
  - THM-2988-first-gap-self-curvature-negativity
script: 04-computation/gmc_first_gap_wall_stripped_all_width_second_edge_circuit_positivity_thm2997.py
output: 05-knowledge/results/gmc_first_gap_wall_stripped_all_width_second_edge_circuit_positivity_thm2997.out
script_sha256: 40959bd9e47fb9ea7bfd9b0ac98b6a303aa1b4b840047b9eef5b2ee214997d5e
output_sha256: 96b7c559f88dbf9e32df3dc796a881b76d2a06399c35cda1ac9e60b38eabb0ed
hash_basis: LF-normalized bytes
---

# THM-2997 -- first-gap wall-stripped all-width second-edge circuit positivity

**PROVED CANDIDATE + VERIFIED-EXACT / AWAITING INDEPENDENT HOSTILE AUDIT.**

This theorem proves one family-specific no-return step for the explicitly
encoded first-gap continuation wall invoice: the second normalized Newton
ratio is strictly larger than the first for every width `M>=34`.  It is
unconditional only at the proved actual width `M=34`; every actual width
`M>=35` still requires the continuation hypotheses below.  It proves neither
full ratio no-return nor full ULC, arbitrary-radial GMC, or GMC(2).

## 1. The circuit and its discrete-log meaning

Let

    N_M(n)=sum_(i=0)^d a_i n^i,  a_d>0,

and normalize its leading coefficients by

    h_k=a_(d-k)/(a_d binom(d,k)),  h_0=1,
    R_k=h_k^2/(h_(k-1)h_(k+1)).                       (1)

For every positive coefficient sequence,

    R_k>=R_(k-1)
      iff h_(k-2) h_k^3 >= h_(k-1)^3 h_(k+1)
      iff Delta^3(log h)_(k-2)<=0.                     (2)

Thus the first two ratio checks are the first two signs of the discrete third
derivative of the normalized log-coefficient sequence.  At `k=2`, (2) is

    R_2>R_1  iff  h_2^3>h_1^3 h_3.                    (3)

Put

    u=a_(d-1)/a_d,
    ell2=2a_(d-2)/a_d-u^2,
    ell3=3a_(d-3)/a_d-3u a_(d-2)/a_d+u^3.             (4)

Equivalently, `ell_j=j[n^-j]log(N_M/(a_d n^d))` for `j=2,3`, while
`u=[n^-1]log(N_M/(a_d n^d))`.  Then

    h_1=u/d,
    h_2=(u^2+ell2)/(d(d-1)),
    h_3=(u^3+3u ell2+2ell3)/(d(d-1)(d-2)).             (5)

Consequently (3) is exactly

    d(d-2)(u^2+ell2)^3
      >u^3(d-1)^2(u^3+3u ell2+2ell3).                  (6)

For

    x=-d ell2/u^2,  z=d^2 ell3/u^3,                   (7)

the difference in (6), after division by the positive factor `u^6/d^2`, is

    F(d,x,z)
      =(d-2)(d-x)^3-(d-1)^2(d^2-3dx+2z)
      =d^2(3x^2-2z-1)
       +d(-x^3-6x^2+3x+4z)+2x^3-2z.                  (8)

The last quadratic-in-`d` form is the key simplification.

## 2. Encoded continuation wall invoice

Use exactly the continuation package of
`THM-2989-first-gap-wall-stripped-all-width-leading-edge-positivity`.  Namely,
for every `M>=34`, let `R_M` be the degree `46M-26` first-gap pure resultant
supplied by the continued THM-2927 factorial Macaulay chart and the same flag
quotient as in the proved finite atlas.  Assume its complete positive wall is

    W_M=(n+1/2)^6 (n+1)^6 (n+2)^24
        product_(s=3)^floor(M/3) (n+s)^28
        product_(s=floor(M/3)+1)^floor(M/2) (n+s)^26
        product_(s=floor(M/2)+1)^floor(2M/3) (n+s)^24
        product_(s=floor(2M/3)+1)^(M-1) (n+s)^23
        (n+M)^20,                                      (9)

with one additional factor `n+(4M+1)/5` when `M=1 mod 10`, and that

    R_M=W_M N_M                                        (10)

up to a nonzero rational unit.  No new Smith, seam, matrix-only, or
pure-resonance wall may occur; the chart degree and flag identity must
continue.  Empty products in (9) are one.

This is an assumption, not an all-width construction of the actual core.

## 3. Statement

For every `M>=34` satisfying (9)--(10),

    R_2(N_M)>R_1(N_M).                                 (11)

At the actual proved width `M=34`, (11) is unconditional by THM-2982 and the
exact atlas.  If (9)--(10) hold for every `M>=35`, then the first-gap family
has

    1<R_1<R_2                                           (12)

at every width `M>=34`, where the first inequality is THM-2989.

The threshold is sharp in the frozen family: at `M=33`, `R_2<R_1`.

## 4. Exact third resultant jet

Set

    U=2^M, V=3^M, D=U^2+3U-3V-1.                     (13)

For the range in scope, `D>0`.  The exact degree-seven Macaulay construction
uses the monomial-priority 36-row chart for degrees `(2,3,4)`.  After the
three-by-three response adjugate, its pure leading matrix is the identity.
The ten degree-seven target monomials divisible by two of
`x_0^2,x_1^3,x_2^4` index the classical extraneous minor.  Hence

    log Res=log det(M_36)-log det(E_10).                (14)

For a normalized slice `I+tX_1+t^2X_2+t^3X_3`, both determinants use

    L_1=tr X_1,
    L_2=tr X_2-tr(X_1^2)/2,
    L_3=tr X_3-tr(X_1X_2)+tr(X_1^3)/3.                 (15)

The exact coefficient calculation first gives numerators over powers of the
response determinant `Delta=-2D`.  Their full removable factors are

    N_1=16D^2 P_1,
    N_2=16D^4 P_2,
    N_3=-32D^6 P_3,                                    (16)

so the reduced resultant log jets are

    L_1=P_1/D^2,
    L_2=P_2/(16D^4),
    L_3=-P_3/(128D^6).                                 (17)

The three exact polynomials `P_1,P_2,P_3` have respectively 27, 122, and 333
terms.  Their complete coefficient-table digest is

    cfb36557e1d54a0328a309375a948ace99c78e0688a54a014aef0906c1b90513.  (18)

This compact table is frozen directly in the companion.  The raw 36-minus-10
Macaulay rebuild independently reproduces every coefficient.  As a low-width
control, it gives exactly

    [n^-3]log Res at M=6
      =958351870363086969113/6204146484375000.          (19)

At `V=0,U->infinity`, (17) reduces exactly to

    L_1=23M^2-12M+66,
    L_2=-(184M^3-144M^2+47M+16758)/24,
    L_3=(92M^4-96M^3+47M^2+210669)/24.                (20)

## 5. Wall subtraction and the normalized box

Let `e,w_1,w_2,w_3` be the degree and first three power sums of (9), including
the half-root and possible mod-10 root.  Define

    d=46M-26-e,
    A=P_1-w_1D^2,
    B=P_2/8+w_2D^4,
    C=-3P_3/128-w_3D^6.                                (21)

Product logarithms give

    u=A/D^2, ell2=B/D^4, ell3=C/D^6,                  (22)

and therefore

    x=-dB/A^2, z=d^2C/A^3.                             (23)

The leading wall moments and core jets are

    e/M ->76/3, w_1/M^2->145/12,
    w_2/M^3->2551/324, w_3/M^4->1681/288,

    d/M->62/3, u/M^2->131/12,
    ell2/M^3->-2417/324, ell3/M^4->1631/288.           (24)

Thus

    x->599416/463347,
    z->12539128/6744273,
    3x^2-2z-1 ->21630685837/71563480803>0.             (25)

The all-width proof uses a rational box stronger than this asymptotic signal.
For every `M>=43`, the exact residue certificate proves

    A>0,
    129/100 <= x <=2,
    0<=z<=39/20.                                       (26)

On this box,

    3x^2-2z-1 >=923/10000,
    -x^3-6x^2+3x+4z >=-26,
    2x^3-2z >=-39/10.                                  (27)

The encoded degree is at least 701 for every `M>=34`, and

    923d^2/10000-26d-39/10
      >=271264123/10000>0 at d=701,                    (28)

with positive derivative thereafter.  Equations (8), (26)--(28) prove
`F>0` for every `M>=43`.

## 6. Complete residue certificate

Write `M=30T+r`, `0<=r<30`, and substitute

    U=2^r X, V=3^rY,
    X=2^(30T), Y=3^(30T).                              (29)

The five strict polynomial inequalities behind (26) are

    A>0,
    -100dB-129A^2>0,
    2A^2+dB>0,
    C>0,
    39A^3-20d^2C>0.                                   (30)

For each residue, the companion groups every polynomial by its exponential
base.  The exact uniform census is:

| inequality | expanded terms | bases | dominant / runner | largest tail start | finite prefixes |
|---|---:|---:|---|---:|---:|
| `A>0` | 27 | 9 | `(4,0)/(2,1)` | 2 | 0 |
| `x>=129/100` | 147 | 25 | `(8,0)/(6,1)` | 3 | 5 |
| `x<=2` | 147 | 25 | `(8,0)/(6,1)` | 2 | 0 |
| `z>=0` | 333 or 337 | 49 | `(12,0)/(10,1)` | 3 | 1 |
| `z<=39/20` | 435 | 49 | `(12,0)/(10,1)` | 3 | 11 |

In every row the runner/dominant base ratio is `3/4`.  If `q_0(T)` is the
dominant polynomial, `L` its leading coefficient, `A_-` its exact negative
lower-term tariff, and `C_abs` the absolute tariff of all lower bases, then
the transcript verifies at its recorded tail start `T_0` that

    P_r/(B^(30T)T^D)
      >= L-sum_(j<D) A_j/T^(D-j)
         -C_abs T^delta (3/4)^(30T)>0.                 (31)

The largest `delta` is two, and the common monotonicity bound is

    4(3/4)^30
      =205891132094649/288230376151711744<1.           (32)

All 150 grouped digests, exact tail margins, and 17 finite-prefix values are
frozen in the canonical transcript.  Their aggregate record digest is

    4e2ddd60fa34f1317733d76b1ec1bde842ff183b22861fa4ab93d2e8b62770b0.  (33)

The nine widths `M=34,...,42` are evaluated directly and exactly from (21).
They all have `F>0`; the ten-row boundary record including the negative
`M=33` control has digest

    191aee5b470ab0c438f058a6c4802fad6d435b54e28fbce09ea31e2372dd6f2e.  (34)

At the sharp first positive width,

    R_2-R_1 = 1.1691055537151569... * 10^-7 >0.        (35)

This completes (11).

## 7. Redundant cleared-circuit controls

The companion also clears (6) directly over `D^12`:

    H=d(d-2)(A^2+B)^3
      -A^3(d-1)^2(A^3+3AB+2C).                        (36)

For every residue, `H` has exactly 2,529 expanded terms and 169 bases, with
dominant `(24,0)`, runner `(22,1)`, and ratio `3/4`.  The dominant degree is
12 and every lower degree is at most 14.  The tail starts are `4` for
`r=0,1`, `3` for `r=2,3`, and `2` for `r=4,...,29`; all 32 finite prefixes are
positive.  The aggregate exact record digest is

    91873acd31b976bd6234a73d81e483701913c568300d19251c0a398baef4eb38.  (37)

As a deliberately uncancelled hostile control, residue four is rebuilt over
the original `Delta^24` denominator.  It has 9,369 terms and 625 bases,
dominant `(48,0)`, runner `(46,1)`, tail start three, and two positive finite
prefixes.  Its grouped digest is

    46e496fa09b968ef8b63254f8fe8dfe943d52bb29e4f1b4275b4070f3d8915ff.  (38)

The factored and uncancelled paths agree exactly.  This control also catches
the invalid discarded prototype: that prototype had only 34 distinct pivots
in its claimed 36-row selected chart and injected binary floats into 28 of 31
third-slice entries.  Neither defect is present in (14)--(18).

Therefore `H>0`, hence `F>0`, hence (6), and finally `R_2>R_1`.

QED.

## 8. Equality and failure boundary

Every certified inequality is strict.  Failure can enter only before this
certificate:

1. the continuation chart, flag quotient, or wall invoice changes;
2. a new Smith, seam, matrix-only, or pure-resonance wall appears;
3. one tries to retain the discarded non-pivot chart; or
4. one extrapolates from the first two ratios to an interior or global
   no-return statement.

The last warning is load-bearing.  Even the negative-real-rooted polynomial

    (x+1)^2(x^2+3x+1)=x^4+5x^3+8x^2+5x+1             (39)

is PF-infinity and strictly ULC, but its normalized ratios satisfy

    R_1=75/64, R_2=256/225, R_2<R_1.                  (40)

Thus real-rootedness, Hurwitz stability, PF-infinity, and ULC cannot replace
the family-specific third-log-jet calculation.

## 9. Exact evidence

Reproduce after installing the repository's exact Python dependencies with

    python 04-computation/gmc_first_gap_wall_stripped_all_width_second_edge_circuit_positivity_thm2997.py --output .scratch/thm2997.normal.out
    python -O 04-computation/gmc_first_gap_wall_stripped_all_width_second_edge_circuit_positivity_thm2997.py --output .scratch/thm2997.opt.out

Normal and optimized outputs are LF-identical: 80 lines and 103,555 bytes.
All checks use optimization-stable explicit guards.  Frozen LF hashes are

    script  40959bd9e47fb9ea7bfd9b0ac98b6a303aa1b4b840047b9eef5b2ee214997d5e
    output  96b7c559f88dbf9e32df3dc796a881b76d2a06399c35cda1ac9e60b38eabb0ed

Promotion still requires an immutable independent hostile audit of the raw
36-minus-10 Macaulay reconstruction, every residue digest and tail tariff,
the exact wall moments, the direct `M=33,34` boundary, and normal/optimized/
stored transcript equality.

**QED (candidate; independent audit pending).**

