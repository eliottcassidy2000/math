---
id: THM-3315
title: "Tournament cut-switching centered-coronal walk compiler"
status: >
  PROVED by exact algebra + VERIFIED-EXACT through order five.  For a
  tournament adjacency `A`, centered operator `C=2A-J`, and cut sign `d~-d`,
  put `P(z)=det(I-zC)` and `N_d(z)=d^T adj(I-zC)d`.  If `T^d` reverses every
  arc across the cut, then
  `det(I-2zA(T^d))=P-zN_d` and
  `G_(T^d)(2z)=N_d/(P-zN_d)`.  Thus the centered switching-class polynomial
  plus one signed-observer numerator compiles the complete total directed-walk
  sequence and a degree-at-most-`n` constant-coefficient recurrence.  Uniform
  labelled cut-cube averages are `E N_d=nP-zP'` and
  `E det(I-2zA(T^d))=(1-nz)P+z^2P'`.  A sharp transitive-triangle hostile
  proves that centered spectrum plus cut cardinality is insufficient.  The
  interface does not restore Hamiltonian/path-cover, SCC, join, or
  substitution-run data.
audit: >
  The proof is the matrix determinant lemma and Sherman-Morrison after exact
  conjugation `D(I-2zA(T^d))D=I-zC-zdd^T`.  The companion independently uses
  integer permutation determinants and matrix powers to exhaust all 1,099
  labelled tournaments of orders one through five and 16,933 cut signs modulo
  global negation.  It checks centered conjugation, observer moments,
  numerator truncation, direct ordinary adjacency determinants, direct walks,
  coefficient convolution, recurrences and cut-cube averages in 252,811
  moment cells.  Normal and optimized outputs byte-match the frozen transcript;
  the source has zero assertion nodes and zero floating literals.
source: root/creative-synthesis-recover/2026-08-03
depends_on:
  - THM-1415-switching-is-the-canonical-star-quotient
  - THM-3181-tournament-half-grid-reciprocity-and-repeated-join-recurrence
  - THM-3248-q4-paired-owner-stirling-compiler
related:
  - THM-1862-order-join-reduction-principle
  - THM-1936-signed-redei-join-multiplicative
  - THM-2501-switching-fourth-moment-signed-c4-and-gram-energy
  - THM-3202-c3-repeated-join-moving-jet-formula-and-cfinite-obstruction
script: 04-computation/tournament_switching_centered_coronal_walk_compiler_20260803.py
output: 05-knowledge/results/tournament_switching_centered_coronal_walk_compiler_20260803.out
script_sha256: 1b68983c9ddbf41c06936287a1184554a889bb0876485c53c2951c03243e2e54
output_sha256: 5ecfb818d9f2edd1401333b792e8aab38c005a1e6801909180aaf8325fa079de
semantic_sha256: 282a0c09458e9b539e349135f3800cd8bfc850971b4e8782fc14de75ea36c1be
hash_basis: LF-normalized bytes
---

# THM-3315 -- tournament cut-switching centered-coronal walk compiler

**PROVED by exact algebra + VERIFIED-EXACT through order five.**

This theorem gives a closed-form sequence compiler for an intrinsic
tournament operation.  Cut switching preserves the centered spectrum but
moves the ordinary all-ones observer.  One signed-observer numerator is the
exact sidecar needed to recover the target total-walk sequence.

## 1. Cut switching and centered adjacency

Let `T` be a labelled tournament on `n` vertices, with row-to-column adjacency

```text
A_ij=1 if i -> j,       A_ii=0.                             (1)
```

Put

```text
C=2A-J.                                                     (2)
```

Thus `C_ii=-1` and the off-diagonal signs record the intrinsic arc
orientations.  For `d in {+1,-1}^n`, let `D=diag(d)` and identify `d` with
`-d`.  Reversing every arc between the two sign parts gives a tournament
`T^d` satisfying

```text
C(T^d)=DCD,              2A(T^d)=DCD+J.                    (3)
```

Define

```text
P(z)=det(I-zC),
U_d(z)=d^T(I-zC)^(-1)d=N_d(z)/P(z),                         (4)
N_d(z)=d^T adj(I-zC)d.                                     (5)
```

Here `P` has degree at most `n`, `N_d` has degree at most `n-1`, and both are
unchanged by the gauge `d -> -d`.

For the switched tournament define the total-walk generating function

```text
G_(T^d)(x)=1^T(I-xA(T^d))^(-1)1
          =sum_(k>=0) w_k(T^d)x^k.                          (6)
```

## 2. Exact determinant and coronal identities

Conjugating by `D` and using `D1=d` gives

```text
D(I-2zA(T^d))D=I-zC-zdd^T.                                 (7)
```

The matrix determinant lemma yields

```text
det(I-2zA(T^d))
 =P(z)(1-zU_d(z))
 =P(z)-zN_d(z).                                             (8)
```

Sherman-Morrison gives

```text
G_(T^d)(2z)
 =d^T(I-zC-zdd^T)^(-1)d
 =U_d(z)/(1-zU_d(z))
 =N_d(z)/(P(z)-zN_d(z)).                                    (9)
```

The rational identities initially hold in the rational-function field;
`(5)`, `(8)` and `(9)` clear denominators, so they are polynomial and formal-
series identities without a numerical invertibility assumption.

Consequently `(P,N_d)` is a complete constant-coefficient compiler for the
total-walk sequence of the selected switch.  `P` records the centered
switching class; `N_d` records the placement of the selected cut observer.

## 3. Coefficient and finite-moment compiler

Write

```text
U_d(z)=sum_(k>=0)u_kz^k,        u_k=d^TC^kd,
s_k=2^k w_k(T^d).                                            (10)
```

Equation `(9)` is equivalent to

```text
s_0=u_0=n,
s_k=u_k+sum_(i=0)^(k-1)u_i s_(k-1-i),       k>=1.           (11)
```

If `P(z)=sum_i p_i z^i`, then the first `n` signed moments determine `N_d`:

```text
[z^m]N_d(z)=sum_(i=0)^m p_i u_(m-i),        0<=m<n.         (12)
```

Cayley-Hamilton makes all later coefficients of `P U_d` vanish.  Put

```text
Q_d(z)=P(z)-zN_d(z)=sum_(j=0)^n q_jz^j.                     (13)
```

Then

```text
sum_(j=0)^n q_j s_(k-j)=0,                  k>=n,           (14)
```

with zero padding if the actual degree is smaller.  Dividing by `2^k`
immediately gives a constant-coefficient recurrence for `w_k`.

## 4. Uniform cut-cube averages

Average uniformly over the `2^(n-1)` sign vectors modulo global negation.
This is a labelled cut-cube average, not an isomorphism-class average.  Since

```text
E_d[dd^T]=I,                                                 (15)
```

one has

```text
E_d U_d(z)=tr(I-zC)^(-1).                                   (16)
```

Jacobi's identity gives

```text
tr(I-zC)^(-1)=n-zP'(z)/P(z).                                (17)
```

Therefore

```text
E_d N_d(z)=nP(z)-zP'(z),                                    (18)
E_d det(I-2zA(T^d))=(1-nz)P(z)+z^2P'(z).                   (19)
```

Thus the mean numerator and mean ordinary characteristic denominator collapse
to centered spectral data, even though individual target sequences do not.

## 5. Sharp hostile

Let `T` be the transitive tournament `1 -> 2 -> 3`, including `1 -> 3`.
Its centered denominator is

```text
P(z)=1+3z+6z^2+4z^3.                                      (20)
```

Switching at the source and switching at the middle vertex are singleton cuts
of equal cardinality.  The source switch remains transitive, whereas the
middle switch is the directed triangle:

| switched vertex | `N_d(z)` | `det(I-2zA(T^d))` | total walks |
|---|---|---|---|
| source | `3+6z+4z^2` | `1` | `3,3,1,0,0,...` |
| middle | `3+6z+12z^2` | `1-8z^3` | `3,3,3,3,3,...` |

Hence

```text
centered spectrum or switching class + cut cardinality
does not determine the ordinary total-walk sequence.         (21)
```

The missing coordinate is exactly `N_d`.  This hostile is minimal: below
order three there is no directed cycle.

## 6. Exact audit and loss ledger

The finite companion exhausts

```text
orders 1 through 5,
1099 labelled tournaments,
16933 switch cases,
252811 observer-moment cells.                                (22)
```

For every case it independently verifies `(3)`, the `d~-d` gauge, direct
matrix-power moments, numerator truncation, direct determinants, direct walks,
the convolution `(11)`, the recurrence `(14)`, and the averages `(16)--(19)`.

The exact interface is

```text
source:      centered switching class plus signed cut observer (C,d);
map:         diagonal conjugation followed by rank-one J recovery;
target:      total directed-walk OGF and recurrence of T^d;
preserved:   centered characteristic polynomial P;
sidecar:     signed-observer numerator N_d, or first n moments;
lost:        fixed all-ones observer when N_d is forgotten;
not restored: Hamiltonian/path-cover, SCC, join, substitution-run data.    (23)
```

In particular, the theorem gives no converse to order-join or signed-Rédei
multiplicativity.  It supplies an exact sequence compiler for cut switching,
not a universal tournament invariant.
