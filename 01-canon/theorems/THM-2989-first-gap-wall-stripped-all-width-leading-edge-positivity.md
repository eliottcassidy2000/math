---
id: THM-2989
title: "First-gap wall-stripped all-width leading-edge positivity"
status: PROVED + VERIFIED-EXACT + INDEPENDENTLY HOSTILE-AUDITED
source: codex-gmc-all-width-sparse-edge-2026-07-30
audit: >
  An immutable independent audit rederived the leading-circuit normalization,
  the 3-by-3 determinant, all 90 symbolic and 720 direct wall-moment checks,
  and all thirty residue numerators.  It independently matched every grouped
  digest, tariff, tail margin, prefix, and the exact width-34 control, then
  replayed normal, optimized, and stored transcripts byte-for-byte.  No defect
  was found; the M>=35 conclusion remains explicitly conditional on the
  encoded continuation wall invoice.
depends_on:
  - THM-2969-first-gap-wall-stripped-resultant-norm-core-atlas
  - THM-2973-first-gap-wall-stripped-norm-core-continuation-through-thirty-one
  - THM-2978-first-gap-wall-stripped-norm-core-at-thirty-two
related:
  - THM-2982-first-gap-wall-stripped-norm-core-strict-ulc-through-thirty-four
  - THM-2986-first-gap-transfer-one-checkerboard-sign
  - THM-2988-first-gap-self-curvature-negativity
script: 04-computation/gmc_first_gap_wall_stripped_all_width_leading_edge_positivity_thm2989.py
output: 05-knowledge/results/gmc_first_gap_wall_stripped_all_width_leading_edge_positivity_thm2989.out
script_sha256: 5dba884e168f78f58a59911efe9620820b32acaf7fb5ebe4e7d383e3e4e20ce7
output_sha256: 6aa181ed1b3f59588b518b694c9f7c6ac3ff517b207faa987c6cb13babb23b77
hash_basis: LF-normalized bytes
---

# THM-2989 -- first-gap wall-stripped all-width leading-edge positivity

**PROVED + VERIFIED-EXACT + INDEPENDENTLY HOSTILE-AUDITED.**

This theorem proves the leading strict-ULC circuit for the explicitly encoded
continuation wall invoice.  For actual norm cores it is unconditional only on
the already proved range M=6..32.  Width 33 uses the separate THM-2982 finite
candidate, and every actual width M>=34 requires the continuation hypotheses
stated below.  There is no full-ULC, interior no-return, arbitrary-radial, or
GMC(2) conclusion.

## 1. The leading circuit

Let

    N_M(n)=sum_(i=0)^d a_i n^i,  a_d>0,

be a wall-stripped first-gap norm core.  Put

    u=a_(d-1)/a_d,
    ell2=2a_(d-2)/a_d-u^2.                              (1)

Strict ultra-log-concavity at the leading nontrivial index i=d-1 is

    a_(d-1)^2(d-1)
      >2d a_(d-2)a_d.                                   (2)

After division by a_d^2 and use of (1), (2) is exactly

    -u^2-d ell2>0.                                       (3)

Thus (3), not merely ell2<0, is the k=1 strict-ULC circuit.

## 2. Encoded continuation wall invoice

For every M>=34 let R_M be the first-gap pure resultant supplied by the
continued THM-2927 factorial Macaulay chart, with degree 46M-26 and the same
flag quotient used in the proved finite atlas.  Assume that its complete
positive wall divisor is

    W_M=(n+1/2)^6 (n+1)^6 (n+2)^24
        product_(s=3)^floor(M/3)       (n+s)^28
        product_(s=floor(M/3)+1)^floor(M/2) (n+s)^26
        product_(s=floor(M/2)+1)^floor(2M/3) (n+s)^24
        product_(s=floor(2M/3)+1)^(M-1) (n+s)^23
        (n+M)^20,                                       (4)

with one additional factor

    n+(4M+1)/5                                          (5)

when M is congruent to 1 modulo 10.  Empty products are one.

The hypothesis is that

    R_M=W_M N_M                                         (6)

up to a nonzero rational unit, with no new Smith, seam, matrix-only, or
pure-resonance wall beyond (4)--(5), and that the chart degree and flag
identity continue.  This package is called the **encoded continuation wall
invoice**.  THM-2989 proves a consequence of this package; it does not prove
the package itself for arbitrary M.

## 3. Statement

For every M>=34 satisfying the encoded continuation wall invoice,
the leading circuit (3) of N_M is strictly positive.

On actual proved cores, (3) is unconditionally positive for M=6..32 by
THM-2969, THM-2973, and THM-2978.  The separate THM-2982 candidate supplies
the exact actual-core checks at M=33,34.  Consequently, conditional on
promotion of that finite candidate and on (4)--(6) for every M>=35, the
leading circuit is positive for every first-gap width M>=6.

No assertion is made about the other d-2 ULC circuits.

## 4. Sparse proof

### 4.1 Resultant log jets

Set U=2^M and V=3^M.  The leading response determinant is

    Delta=-2(U^2+3U-3V-1).                              (7)

For M>=4, 4^M>3*3^M, hence the parenthesis in (7) is
positive.  In particular Delta is nonzero and Delta^8>0 for every width in
scope.

The exact top-three factorial-slice and adjugate calculation writes

    u(R_M)=A_M/Delta^4,
    ell2(R_M)=B_M/Delta^8,                              (8)

for explicit polynomials A_M,B_M in M,U,V.  Formula (8) is obtained using
only base coefficients and the two possible edge transfers in the local
resultant Hessian; no dense symbolic resultant is expanded.

Let

    e=deg W_M,
    w1=sum_s e_s s,
    w2=sum_s e_s s^2,                                  (9)

where s=1/2 is allowed and e_s is its wall multiplicity.  Product
logarithmic jets give

    d=46M-26-e,
    u(N_M)=A_M/Delta^4-w1,
    ell2(N_M)=B_M/Delta^8+w2.                           (10)

Therefore the numerator of (3), over the positive denominator Delta^8, is

    P_M=-(A_M-w1 Delta^4)^2
        -d(B_M+w2 Delta^8).                             (11)

### 4.2 Complete residue split

Write M=30T+r, 0<=r<30, and substitute

    U=2^r X,  V=3^r Y,
    X=2^(30T),  Y=3^(30T).                             (12)

The floors in (4) and the quartic indicator (5) are now polynomial data in
T.  For every residue, the exact companion expands (11) as

    P_r(X,Y,T)=sum_(a,b) q_(r,a,b)(T) X^a Y^b.          (13)

Every one of the 30 polynomials has 483 expanded terms and 81 grouped
exponential bases.  The complete grouped coefficient lists have separate
SHA-256 digests in the canonical transcript.  Their aggregate record digest
is

    92482d7309ff570ed71b8f1146a896b63259e76ff7f62eaebe42b8ef53675a93.

The complete 30-residue degree, first-wall-moment, and square-wall-moment
invoice has digest

    3dbcf0165edbe1ba8c2a0a1dd8660884cb9f5049008a3e018e5100b335eafeab.

These digests freeze all coefficients, not sampled evaluations.

### 4.3 Dominant-base certificate

For every residue, the dominant grouped base is

    B=2^16=65536,

and the runner-up is

    B'=2^14*3=49152,  B'/B=3/4.                         (14)

The dominant coefficient has degree 4 in T.  Every lower coefficient has
degree at most 5.  For each residue the transcript freezes:

- the full grouped-representation digest;
- the dominant polynomial, its leading coefficient, and its negative
  lower-coefficient tariff;
- the absolute coefficient tariff of all 80 lower bases;
- the exact dominant floor;
- the exact tail ratio and positive tail margin;
- the tail monotonicity bound; and
- every finite-prefix value and its digest.

If L_r is the dominant leading coefficient, A_r its negative tariff,
C_r the absolute lower tariff, and T_0 the recorded tail start, then for
T>=T_0

    q_dom(T)/T^4 >= L_r-A_r/T=:f_r,

and

    P_r/(B^(30T)T^4)
      >= f_r-C_r T(3/4)^(30T).                          (15)

The exact transcript verifies the right side at T_0 is positive.  It also
checks the common step bound

    2(3/4)^30
      =205891132094649/576460752303423488 <1,           (16)

so the lower majorant decreases thereafter.

The tail starts are:

- r=0: T_0=3, with the single finite prefix T=2;
- r=1,2,3: T_0=2, with no prefix;
- r=4: T_0=2, with the single finite prefix T=1; and
- r=5,...,29: T_0=1, with no prefix.

All two finite-prefix values are strictly positive.  Among all 30 exact tail
margins, the minimum occurs at r=5 and is

    35623413693023365408297094332961860518254695069409
    /59266408172243425004797336907050584018696273920000 >0.  (17)

Thus P_M>0 for every M>=34.  By (11), Delta^8>0, and the equivalence
(2)--(3), the leading strict-ULC circuit is positive.

QED.

## 5. Equality and failure boundary

Every certified inequality is strict; no residue or tail case attains
equality.  Failure can enter only before (11):

1. a new wall changes e,w1,w2;
2. the chart/resultant degree or flag quotient changes;
3. a new Smith or seam correction changes (6); or
4. one tries to infer an interior ULC circuit from the leading one.

THM-2986 and THM-2988 explain why the pure-resultant local curvature has the
correct sign, but they do not remove any of these wall or no-return
boundaries.

## 6. Exact evidence

Reproduce with

    python 04-computation/gmc_first_gap_wall_stripped_all_width_leading_edge_positivity_thm2989.py --output .scratch/thm2989.normal.out
    python -O 04-computation/gmc_first_gap_wall_stripped_all_width_leading_edge_positivity_thm2989.py --output .scratch/thm2989.opt.out

The normal and optimized transcripts are LF-identical: 43 lines and 58,202
bytes.  All checks use optimization-stable explicit guards.  The companion
also reconstructs the independent width-34 edge ratio

    27655303142730846948837614910195097534441250912778059699455304305433884593670175322013125
    /27642422226229352775679708685499047535273900449378640304403997574607551772198286144871066

and matches the frozen THM-2982 value exactly.

Frozen LF hashes are

    script  5dba884e168f78f58a59911efe9620820b32acaf7fb5ebe4e7d383e3e4e20ce7
    output  6aa181ed1b3f59588b518b694c9f7c6ac3ff517b207faa987c6cb13babb23b77

The independent hostile audit rederived the wall-invoice translation through
90 symbolic identities and 720 redundant direct samples, including the
half-wall and mod-10 quartic root.  It separately reconstructed all thirty
numerators (483 terms and 81 bases each), every grouped digest and tariff,
both finite prefixes, all tail margins and monotonicity bounds, the determinant
normalization, and the exact width-34 equality control.  Normal, optimized,
and stored transcripts match byte-for-byte.  No defect was found.

**QED.**


