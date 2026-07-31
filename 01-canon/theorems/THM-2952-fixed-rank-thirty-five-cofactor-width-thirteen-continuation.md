---
id: THM-2952
title: "Fixed rank-thirty-five cofactor width-thirteen/fourteen continuation"
status: >
  PROVED + VERIFIED-EXACT + INDEPENDENTLY HOSTILE-AUDITED.  The fixed
  rank-35 cofactor of PROVED THM-2949 is
  nonzero at every integer depth on all 66 normalized four-slot
  supports of width thirteen and all 78 of width fourteen.  Exact
  finite prefixes and shifted Gregory--Newton tails close first-window
  SFC(4) on these 144 new support types.  PROVED THM-2956 transfers
  the same nonvanishing atlas to its fixed lower-invoice
  21Q+9C+5F chart.  Width fourteen contains a
  sharp real-ray hostile: the (0,1,2,14) cofactor has 24 positive
  integer values, then 19 negative values, then a positive Newton tail
  from depth 43.  No width-fifteen, real-depth, arbitrary-width,
  SFC(5), NC2, or GMC(2) claim is made.
source: codex-gmc-fixed-cofactor-continuation-2026-07-29
audit: >
  An independent hostile read accepted the promoted THM-2949 import
  typing, fixed cofactor, shifted-Newton proof, width-fourteen real-ray
  boundary, THM-2956 corollary, and scope.  A separate exact serial
  check recomputed the two maximum-base families and three controls.
  Normal, optimized, and stored transcripts match after LF
  normalization with empty stderr.
depends_on:
  - THM-2949-fixed-rank-thirty-five-cofactor-newton-atlas
  - THM-2956-koszul-gale-fixed-fifth-compound-exchange
related:
  - THM-2925-general-width-terminal-pole-cancellation-and-macaulay-degree-law
  - THM-2947-conjugate-pair-corank-parity-and-one-minor-resultant-gate
script: 04-computation/gmc_fixed_rank_thirty_five_cofactor_width_thirteen_fourteen_thm2952.py
output: 05-knowledge/results/gmc_fixed_rank_thirty_five_cofactor_width_thirteen_fourteen_thm2952.out
script_sha256: a1b00ca16ca29eca50c389ebefbefa904bc70b3b5e6fc6791103fbfa8ca7ce23
output_sha256: 818c99cbf95823b4776982eb0c6d9331f44e64d31c4492345796ac0ee5805949
constructor_dependency_sha256: 9a1c7068e079e232dc97fd6eb925621aa74b3d636380a85995f8e0db8b30aa54
hash_basis: LF-normalized bytes
---

# THM-2952 -- fixed rank-thirty-five cofactor width-thirteen/fourteen continuation

**PROVED + VERIFIED-EXACT + INDEPENDENTLY HOSTILE-AUDITED.**

## 1. Statement

Let

```text
L:C[s] -> C,                              L(s^j)=j!.       (1)
```

Fix

```text
M in {13,14},                  0<a<b<M,                  n>=0. (2)
```

For every nonzero

```text
H=c_0 s^n+c_1 s^(n+a)+c_2 s^(n+b)+c_3 s^(n+M),          (3)
```

at least one of

```text
L(H),                    L(H^2),                    L(H^3), L(H^4) (4)
```

is nonzero.  Thus first-window SFC(4) holds on the

```text
C(12,2)+C(13,2)=66+78=144                              (5)
```

normalized support types of exact width thirteen or fourteen.
Together with PROVED THM-2949, this gives every normalized translated
four-slot support through width fourteen:

```text
C(14,3)=364 support types.                              (6)
```

## 2. The inherited fixed fifth-compound coordinate

Use the denominator-cleared ternary moment forms

```text
Q,                           C,                           F
deg 2,                       deg 3,                       deg 4       (7)
```

and degree-seven Macaulay map from THM-2949.  In the inherited square
chart select the global rows

```text
0,...,19; 21,...,29,35; 36,...,41,                     (8)
```

delete selected row-position `30` (global quartic row `36`, multiplier
`x2^3`) and target column `0` (monomial `x2^7`), and call the resulting
`35`-minor

```text
P_(M,a,b)(n).                                           (9)
```

It has `20` quadratic, `10` cubic, and `5` quartic rows.  The PROVED
THM-2925 degree law used by THM-2949 gives

```text
deg P_(M,a,b)
 <=20(M-1)+10(2M-1)+5(3M-1)
 =55M-35.                                               (10)
```

Every polynomial in the present exact bank attains equality.  Hence

```text
d_13=680,                             d_14=735.          (11)
```

This is exactly the same fixed cofactor as THM-2949, not an adaptively
chosen chart.

## 3. Exact shifted-Newton proof

For every pair `1<=a<b<M`, the companion:

1. evaluates `(9)` exactly at `0,...,d_M`;
2. interpolates the unique polynomial of degree at most `d_M`;
3. verifies direct determinants at the three independent depths

   ```text
   d_M+1,                         d_M+2,             2d_M+3; (12)
   ```

4. orients the polynomial to have positive leading coefficient;
5. finds a base `0<=B<=4M` such that

   ```text
   Delta^j P_(M,a,b)(B)>=0                    for 0<=j<=d_M,
   P_(M,a,b)(B)>0;                                      (13)
   ```

6. checks every earlier integer value exactly and finds no zero.

For every integer `n>=B`, the Newton expansion is

```text
P_(M,a,b)(n)
 =sum_(j=0)^d_M Delta^j P_(M,a,b)(B) binom(n-B,j)>0.    (14)
```

The finite prefix and `(14)` prove

```text
P_(M,a,b)(n)!=0                         for every n>=0. (15)
```

The exact census is

```text
width  families  degree  max B  base-zero  mixed prefixes
  13      66       680      35      22            12
  14      78       735      43      23            17.   (16)
```

The unique width-thirteen maximum is

```text
(a,b)=(1,12),                                  B=35,    (17)
```

and the width-fourteen maxima are

```text
(a,b)=(1,2),(1,13),                            B=43.    (18)
```

There are `432` direct outside-grid determinant controls.  The
per-family records contain the exact degree, base, complete low-depth
sign word, constant sign, and a digest of the entire shifted Newton
vector.  Their immutable exact digests are

```text
M=13:
c1e7eaf3a794773573a664f9e83c439817efa6795f4c690b2849b79b7acad923

M=14:
2158531e878df93c158a85d83de9b4f7b10a6d25b3c1143a71d9b2ff3721c3d8

joint:
02566dbe576277c3c1582ca778505353496c81d5d5530f9c430f43cc51d0f515. (19)
```

## 4. Projective consequence

At every physical depth the quadratic `Q` is positive definite.
Therefore the real projective common-zero locus is empty.  PROVED
THM-2949, through the conjugate-pair rank gap of THM-2947, gives

```text
one nonzero 35-minor
  ==> rank(Phi_7)>=35
  ==> Res(Q,C,F)!=0.                                  (20)
```

Equation `(15)` supplies that minor.  Hence `Q,C,F` have no common
nonzero projective zero.  Mean elimination gives `(4)`.

PROVED THM-2956 gives the universal identity

```text
q_200 P_(M,a,b)+c_300 P^opt_(M,a,b)=0,               (20a)
```

where `q_200(n)>0` on every physical factorial specialization and
`P^opt` is the fixed `21Q+9C+5F` cofactor.  Therefore `(15)` also
forces

```text
P^opt_(M,a,b)(n)!=0                     for every n>=0. (20b)
```

Its formal degree invoices are

```text
deg P^opt<=54M-35 =
 667  for M=13,
 721  for M=14.                                        (20c)
```

This exchange is a proved corollary, not an input to the present exact
bank, so the companion does not recompute the optimal determinants.

## 5. The width-fourteen real-ray hostile

For the support

```text
(0,1,2,14),                                            (21)
```

the normalized fixed cofactor has exact integer sign word

```text
P(n)>0,                           0<=n<=23;
P(n)<0,                          24<=n<=42;             (22)
```

and all its shifted Newton coefficients are nonnegative from base
`43`, with `P(43)>0`.  It therefore has a real zero in each of

```text
(23,24),                              (42,43),          (23)
```

but no nonnegative integral zero.  The certificate is a discrete
Pluecker thread, not positivity of one chart on the whole real ray.

## 6. Scope

Proved, after exact replay and promotion:

```text
slots:         four distinct translated exponents;
widths:        exactly thirteen and fourteen;
depths:        every integer n>=0;
window:        moments one through four;
certificate:   the fixed THM-2949 rank-35 cofactor.     (24)
```

Not proved:

```text
* width fifteen or any arbitrary-width continuation;
* fixed-cofactor nonvanishing for every real depth;
* a canonical or positive resultant factorization;
* a shifted moment window, SFC(5), NC2, or GMC(2).      (25)
```

## 7. Exact companion

The companion imports and LF-hash pins the promoted THM-2949 script,
then reconstructs all `144` families with exact FLINT arithmetic.
Run

```text
python 04-computation/gmc_fixed_rank_thirty_five_cofactor_width_thirteen_fourteen_thm2952.py
python -O 04-computation/gmc_fixed_rank_thirty_five_cofactor_width_thirteen_fourteen_thm2952.py
```

Both modes byte-match the stored transcript after LF normalization,
with empty stderr.

**QED.**
