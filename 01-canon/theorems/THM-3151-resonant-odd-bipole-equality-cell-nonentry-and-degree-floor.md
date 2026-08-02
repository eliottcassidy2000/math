---
id: THM-3151
title: "Resonant odd-bipole equality-cell nonentry and degree floor"
status: >
  PROVED + VERIFIED-EXACT + INDEPENDENTLY HOSTILE-AUDITED.  For every
  normalized nonsplit polynomial exact-square-prefix Faber bank with
  R=3k+2>=8 and D=4k+3, no balanced response in the equality cell N=2D,
  pole passport (D,D), enters the source chart.  Consequently every balanced
  entrant in that bank has N>=4D.  This is a chart theorem, not exclusion of
  unbalanced or nonpolynomial branches and not the planar Jacobian conjecture.
audit: >
  Two independent coefficient extractors verified the universal even/odd
  infinity fan, all retained lower rows, wall coprimality, zero-polynomial
  boundaries, resonance exponents, and the final order contradiction.
  Ordinary and optimized transcripts byte-match.
source: root/frontier-synthesis/2026-08-02
depends_on:
  - THM-2796-balanced-response-stieltjes-pade-normal-form-and-one-double-zero-classification
  - THM-2827-uniform-pole-order-faber-nonresonance-atlas-and-double-pole-exclusion
  - THM-3133-common-simple-zero-faber-exclusion-and-odd-bipole-boundary
related:
  - THM-3140-r8-odd-bipole-faber-chart-nonentry
  - MISTAKE-317
script: 04-computation/jc_all_resonance_equality_cell_nonentry_thm3151.py
output: 05-knowledge/results/jc_all_resonance_equality_cell_nonentry_thm3151.out
script_sha256: ade9667dc3cda461c89f597538ce48a4ccba1b29b30eac688945244be2557215
output_sha256: e5041706d1205d2ccf73cd25ca500868aa9d7dfbbbb92de0dfb184130342a698
independent_script: 04-computation/jc_all_resonance_equality_cell_nonentry_independent_audit_thm3151.py
independent_output: 05-knowledge/results/jc_all_resonance_equality_cell_nonentry_independent_audit_thm3151.out
independent_script_sha256: cc3e21a2326976716df0fdc6ee035fc24fe8ae47679d8a7184d98ac513b65411
independent_output_sha256: 67296119cbd2c594ecd2f68f6c4d87c75f830256ae93ec5b9bf18f8b493dc841
hash_basis: LF-normalized bytes
---

# THM-3151 -- resonant odd-bipole equality-cell nonentry and degree floor

**PROVED + VERIFIED-EXACT + INDEPENDENTLY HOSTILE-AUDITED.**

THM-3133 gives the sharp abstract balanced-response boundary `N>=2D` in a
resonant Faber bank, while THM-3140 excludes one explicit response at the
first equality cell.  The present theorem shows that the infinity mechanism
is response-independent and uniform in `R`: the whole equality cell is
empty, doubling the balanced chart floor.

## 1. Statement and response typing

Work over an algebraically closed field of characteristic zero.  Fix

```text
R=3k+2>=8,                 D=4k+3,                 k>=2.     (1)
```

Suppose a nonconstant balanced response with pole passport `(D,D)` enters a
normalized nonsplit polynomial exact-square-prefix Faber chart

```text
Q=E_(4R-2)+sum_(j=1)^(R-2) a_j E_(4j-2).                    (2)
```

THM-3133 removes response simple zeros.  In THM-2796's factor normal form,
the equality cell has

```text
S=1,                 e=D,                 h=2.               (3)
```

Let `T_p` be the squarefree monic quadratic with the two pole roots, and let
`E` be the squarefree double-zero polynomial.  The factors are pairwise
disjoint and

```text
deg(E)=D,
V=v T_p^(D+2),              M=E T_p,              v!=0.      (4)
```

In particular, `M` is squarefree.  Retain the source coordinates

```text
T=A_src^2/V,
d_F=C_0-B_src^2/(4V),
s_F=A_src B_src/(2V)-E_0,                                (5)
```

and the exact Faber sidecars

```text
Phi_Q=0,              Psi_Q in C,
K_Q=T H_Q,            A_src K_Q=lambda M,       lambda!=0,
H_Q in C[T,s_F,Td_F].                                  (6)
```

## 2. Divisor support and the resonant pole orders

Let `alpha` be a zero of `A_src` away from `V`, of order `a>=1`.  The three
arguments of `H_Q` in `(6)` are regular there and

```text
ord_alpha(T)=2a,
ord_alpha(K_Q)>=2a,
ord_alpha(A_src K_Q)>=3a>1.                              (7)
```

This contradicts the squarefreeness of `M`.  Thus every zero of `A_src` lies
at one of the two roots of `T_p`.

At either pole root,

```text
nu=ord(V)=D+2=4k+5,                 ord(M)=1.             (8)
```

The exact THM-2827 trichotomy excludes the polar and regular-`H` lanes.  Its
unique `R=3k+2` resonance has `delta=1` and forces

```text
a=ord(A_src)=1+(2k+1)=2k+2=(D+1)/2,
b=ord(B_src)>=2k+3=a+1.                                (9)
```

There are no other zeros of `A_src`, so for constants `c,c'!=0`,

```text
A_src=c T_p^a,
T=c'/T_p,
T_p^(a+1) divides B_src.                               (10)
```

The exponent identities

```text
2(a+1)-(D+2)=1,                 a+(a+1)-(D+2)=0          (11)
```

show directly from `(5)` that `d_F,s_F` are polynomials.

## 3. Universal infinity fan

Normalize `ord_infinity(x)=-1`.  Since `T=c'/T_p`, choose the local deck
coordinate `q` with `q^2=T`; then

```text
ord_infinity(T)=2,                ord_infinity(q)=1.        (12)
```

The Faber rows are extracted from

```text
[1+2d_F t^2+q t^3+(d_F^2-s_F)t^4]^(j-1/2)
 =[(1+d_Ft^2)^2+t^3(q-s_Ft)]^(j-1/2).                     (13)
```

Assume first that `d_F,s_F` are nonzero and put

```text
u=deg(d_F),                 w=deg(s_F),
Gamma=2w+2-u.                                               (14)
```

Direct expansion of the factorized expression in `(13)` gives every
`Phi_j` monomial the exponent form

```text
(q,d_F,s_F)=(2r,r-2ell,j-1-2r+ell),                         (15)
```

and every `Psi_j` monomial the form

```text
(q,d_F,s_F)=(2r,r-2ell,j-2r+ell),                           (16)
```

for `ell>=0` on the surviving support.  Relative to the pure-`s_F` face,
the infinity weight increment is exactly

```text
(r-2ell)Gamma+ell(4+3w).                                   (17)
```

This produces three exhaustive chambers.

### Chamber I: `Gamma>0`

If `w>=1`, the pure-`s_F` faces are unique and nonzero.  Their top-row
orders are

```text
ord(Phi_R)=-(R-1)w,             ord(Psi_R)=-Rw.             (18)
```

Every retained row `j<=R-2` is strictly higher in both channels.  The top
row cannot be cancelled.

### Chamber II: `Gamma=0`

Only `ell=0` survives on the wall.  With leading coefficients `d_0,s_0,q_0`
and

```text
z=d_0 q_0^2/s_0^2,                                            (19)
```

the two top faces are nonzero scalar multiples of

```text
s_0^(R-1) P_R(z),                   s_0^R Q_R(z),             (20)

P_j(z)=sum_r (-1)^r binom(j,2r+1)z^r,
Q_j(z)=sum_r (-1)^r binom(j,2r)z^r.                          (21)
```

The exact binomial identity

```text
Q_j(z)^2+z P_j(z)^2=(1+z)^j                                  (22)
```

proves `gcd(P_j,Q_j)=1`.  Indeed, a common root would be `z=-1`, but

```text
P_j(-1)=Q_j(-1)=2^(j-1)!=0.                                  (23)
```

At least one top flux therefore survives; in its channel, `(18)` again puts
it strictly below every retained lower row.  This coefficient sidecar, not
face uniqueness alone, is what avoids MISTAKE-317.

### Chamber III: `Gamma<0`

The maximal-`d_F` endpoint is unique.  When `R` is even it lies in `Psi_R`
and is a nonzero multiple of

```text
q^R d_F^(R/2).                                               (24)
```

When `R` is odd it lies in `Phi_R` and is a nonzero multiple of

```text
q^(R-1)d_F^((R-1)/2).                                        (25)
```

The chosen top endpoint is strictly below the same channel in every
`j<=R-2`, for both parities.  It is a genuine pole and cannot be cancelled
by `(2)`.

The three chambers exclude every nonconstant `s_F`.  If `d_F=0`, the
pure-`s_F` endpoint remains.  If `s_F=0`, deleting its monomials leaves the
parity endpoint `(24)` or `(25)`, which excludes every `deg(d_F)>=3`.
The same third chamber handles constant nonzero `s_F`.  Therefore

```text
s_F is constant (possibly zero),
d_F=0 or deg(d_F)<=2.                                      (26)
```

## 4. The order contradiction

Equations `(12),(26)` make each of `T,s_F,Td_F` regular at infinity.  The
ring statement in `(6)` therefore gives

```text
ord_infinity(H_Q)>=0,              ord_infinity(K_Q)>=2.     (27)
```

The alternative `H_Q=0` contradicts the nonzero right side of `(6)`.
But `(4),(6),(10)` require

```text
K_Q=lambda M/A_src=lambda' E/T_p^(a-1),

ord_infinity(K_Q)=2(a-1)-D=(D-1)-D=-1.                     (28)
```

Equations `(27),(28)` are incompatible.  No balanced `(D,D)` response enters
the normalized chart `(2)`.

## 5. Degree floor and uniqueness of the equality response

For any balanced entrant in the bank, THM-3133 gives

```text
s_resp=0,                   p_j=D delta_j.                   (29)
```

Since `N=2e` is even and `D` is odd,

```text
N/D=sum_j delta_j                                             (30)
```

is even.  Nonsplitting requires some `delta_j` odd.  At sum two, the only
partitions are the split single part `(2)` and the equality passport `(1,1)`
just excluded.  The next even sum is at least four, so every balanced entrant
satisfies

```text
N>=4D=4(4k+3).                                                (31)
```

There are no hidden abstract responses inside the equality cell.  After an
affine normalization of the two pole roots to `+1,-1`, THM-2796's first
integral becomes

```text
2(y^2-1)E'(y)-2D y E(y)=C!=0.                               (32)
```

Two monic degree-`D` solutions would differ by a polynomial of degree
`n<=D-1`; its leading term under the left side of `(32)` has nonzero
coefficient `2(n-D)y^(n+1)`.  Hence the monic solution is unique, recovering
THM-3133's explicit odd-bipole polynomial up to affine and scalar
normalization.  The chart proof above is stronger: it never uses those
coefficients or the affine normalization.

## 6. Exact controls and scope

The primary referee checks exact Faber supports through row `50`, both
parities, all lower rows, wall identity `(22)`, zero-polynomial boundaries,
and resonance arithmetic through `k=30`, for `18,320` gates.  The independent
audit reconstructs `62` flux rows by two separate expansions, checks `31`
wall gcds, and scans `2,288` hostile degree/parity/boundary chambers.  Run

```text
python3 04-computation/jc_all_resonance_equality_cell_nonentry_thm3151.py
python3 -O 04-computation/jc_all_resonance_equality_cell_nonentry_thm3151.py
python3 04-computation/jc_all_resonance_equality_cell_nonentry_independent_audit_thm3151.py
python3 -O 04-computation/jc_all_resonance_equality_cell_nonentry_independent_audit_thm3151.py
```

The theorem is confined to balanced responses entering the normalized
polynomial exact-square-prefix Faber chart.  Squarefreeness/disjointness,
characteristic zero, the missing `R-1` row, and
`H_Q in C[T,s_F,Td_F]` are load-bearing.  Nothing here excludes unbalanced,
nonpolynomial, or other chart branches, constructs a Keller map, or proves
`JC(2)`.  QED.
