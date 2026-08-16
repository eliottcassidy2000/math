---
id: THM-3484
title: "Ternary weighted-determinant minimal recurrence and cubic Fourier degree drop"
status: >
  PROVED + VERIFIED-EXACT + INDEPENDENTLY AUDITED.  The explicit
  period-three degree-seven word defined below
  has minimal characteristic polynomial
  (x-1)^8(x^2+x+1)^7 and hence minimal recurrence order 22, two below the
  naive order 24.  The drop occurs because all three residue polynomials have
  the same leading coefficient, lowering both nontrivial cubic Fourier
  colours to degree six.  Every residue lane has natural and harmonic
  coefficient 1/3.  By THM-3482 this is the canonical private-gradient
  determinant word; no relation-current or LRC claim follows.
source: codex-2026-08-15-ternary-spectral-recurrence
audit: >
  self-contained integer and Q(zeta_3) arithmetic; exact degree decomposition;
  order-22 recurrence through k=5000; generating numerator; two shortened-
  factor hostiles; rational Berlekamp-Massey reconstruction; exact nonzero
  22x22 Hankel determinant and prime-power factorization; security, semantic,
  and normal/optimized/stored replay gates; independent derivation audit of
  every polynomial lane, Fourier component, minimality argument, recurrence,
  generating numerator, hostile, harmonic coefficient, and Hankel determinant
depends_on:
  - THM-3482-private-count-gradient-weighted-spectral-closure-without-absolute-h1-flux
related:
  - THM-3473-three-times-p-eight-owner-private-sheet-partition-and-irredundancy
  - THM-3455-berggren-q-spine-cap-seven-atom-sieve-and-fibonacci-rank-spectrum
script: 04-computation/lrc_ternary_weighted_determinant_minimal_recurrence_thm3484.py
output: 05-knowledge/results/lrc_ternary_weighted_determinant_minimal_recurrence_thm3484.out
script_sha256: e0bf5f055b65a508c572530fceae70888e12cafed249d840ea0508cf8c4204fd
output_sha256: d6f68554cea8243670e84a4d7ac32d610553dbb582091b0e77ee4a7b2af83ea3
semantic_sha256: 0c9cac55b66f8a8aa4241ad728e3336dd000cff25544f3df7eef9af120012518
hash_basis: LF-normalized bytes
---

# THM-3484 -- the ternary determinant word has order 22, not 24

**PROVED + VERIFIED-EXACT + INDEPENDENTLY AUDITED.**

The recurrence proof is self-contained for the explicit sequence below.
Proved THM-3482 identifies that sequence with its private-gradient weighted
determinant word.

## 1. The three polynomial lanes

For `k>=0`, define `D_k=P_r(k)` when `k==r (mod 3)`, where coefficient
tuples below are written constant term first:

```text
P_0(k)=-2048k^4+4096k^5+8192k^6-16384k^7,

P_1(k)= 2048k^4+4096k^5-24576k^6-16384k^7,

P_2(k)=-256k+3072k^2-15360k^3+40960k^4
       -61440k^5+49152k^6-16384k^7.                (1)
```

These are the three factored determinant lanes proved in THM-3482, expanded.
Here `(1)` may also be read as a self-contained sequence definition.

The naive quasi-polynomial bound says that any period-three degree-seven word
is killed by

```text
(x^3-1)^8,                                          (2)
```

an order-24 characteristic polynomial.  The present word has two fewer
degrees of freedom.

## 2. Cubic Fourier decomposition

Let `omega^2+omega+1=0`.  Define

```text
Q_j(k)=1/3 sum_(r=0)^2 omega^(-jr) P_r(k),
                                                     j=0,1,2.       (3)
```

The cubic indicator identity gives

```text
D_k=sum_(j=0)^2 omega^(jk) Q_j(k).                  (4)
```

All three polynomials in `(1)` have leading coefficient `-16384`.  Therefore
the degree-seven coefficients cancel in both nontrivial Fourier colours.  An
exact calculation in the basis `(1,omega)` gives

```text
deg Q_0=7,    lead(Q_0)=-16384,

deg Q_1=6,    lead(Q_1)=32768/3+24576 omega,

deg Q_2=6,    lead(Q_2)=-40960/3-24576 omega.        (5)
```

Both nontrivial leading coefficients are nonzero.  Equivalently, the triple
of degree-six coefficients

```text
(8192,-24576,49152)                                 (6)
```

is not constant, so its two conjugate nontrivial Fourier transforms cannot
vanish.

## 3. Minimal characteristic polynomial

A sequence `lambda^k R(k)` with `deg R=d` is annihilated by
`(x-lambda)^(d+1)` and by no smaller power when the leading coefficient of
`R` is nonzero.  Distinct exponential colours are linearly independent.
Equations `(4)--(5)` therefore give the minimal characteristic polynomial

```text
chi(x)=(x-1)^8(x-omega)^7(x-omega^2)^7
      =(x-1)^8(x^2+x+1)^7.                          (7)
```

Its degree is

```text
8+7+7=22.                                           (8)
```

This proves both the order-22 recurrence and its minimality.  The two-degree
drop from `(2)` is exactly the shared leading term, one cancellation in each
nontrivial cubic colour.

There is an independent rational certificate.  The first `22x22` Hankel
matrix

```text
H=(D_(i+j))_(0<=i,j<22)
```

has determinant

```text
det H=-2^382 3^191 5^22 7^8 61^7 !=0.              (9)
```

Thus no recurrence of order at most 21 can generate the word, while `(7)`
supplies one of order 22.

## 4. Generating function and exact recurrence

With `G(z)=sum_(k>=0)D_k z^k`, the minimal denominator is

```text
Q(z)=(1-z)^8(1+z+z^2)^7
    =(1-z)(1-z^3)^7.                                (10)
```

Its coefficient tuple is

```text
(1,-1,0,-7,7,0,21,-21,0,-35,35,0,
 35,-35,0,-21,21,0,7,-7,0,-1,1).                  (11)
```

Consequently, for every `n>=22`,

```text
sum_(j=0)^22 Q_j D_(n-j)=0,                         (12)
```

where `Q_j` denotes the `j`th entry of `(11)`.  The exact numerator
`N(z)=Q(z)G(z)` has coefficient tuple

```text
(0,-34816,-338432,-28657152,-335106048,-313495296,
 -3294224640,-9788725248,-4813545216,-26592198912,
 -36296777728,-11030363648,-39322132992,-27167668224,
 -5002493952,-11956875264,-3768176640,-338608896,
 -502246656,-40216576,-239360,-186624).             (13)
```

Deleting one copy of the trivial factor first fails at index `21`, with
residual `-180592312320`.  Deleting one copy of the quadratic cubic factor
first fails at index `20`, with residual `55897620480`.  These are cheap
hostiles for both multiplicities in `(7)`.

## 5. Tournaments, ternary colour, and harmonic subsets

THM-3482 obtains `(1)` from two four-vertex coactivity packets.
Each `K4` has six undirected edges and contributes a degree-three tree sum;
the forced bridge contributes degree one.  Thus the geometric factorization

```text
six-edge K4  +  six-edge K4  +  one bridge,
degree three + degree three + degree one             (14)
```

becomes the degree-seven word whose ternary Fourier colours are analysed
above.  This is an exact connection between tournament-sized-four packets,
their six pairwise observables, and a period-three recurrence.  It is not an
identification with a three-child tree.

For each `r`, the subset

```text
A_r={k>=1:k==r (mod 3)}                              (15)
```

has natural density and harmonic coefficient `1/3`:

```text
sum_(k<=N,k in A_r) 1/k=(1/3)log N+O(1).            (16)
```

Hence `(1)` is a three-coloured subset decomposition of the harmonic series,
and `(7)` measures how much recurrence data survives after the colours are
interlaced.  Arbitrary subsets of the natural numbers still define harmonic
subseries, but they need not possess a natural density or a harmonic
coefficient; periodicity is the decisive sidecar in `(16)`.

## 6. Exact companion and scope

Run from the repository root:

```bash
python 04-computation/lrc_ternary_weighted_determinant_minimal_recurrence_thm3484.py
python -O 04-computation/lrc_ternary_weighted_determinant_minimal_recurrence_thm3484.py
```

The standard-library companion derives `(10)--(13)`, audits `(12)` through
`k=5000`, computes the exact `Q(omega)` degrees and leading coefficients,
reconstructs the recurrence by rational Berlekamp--Massey, and independently
computes and checks the Hankel determinant `(9)`.  It includes both shortened-
factor hostiles, security gates, and a frozen semantic digest.

The independent immutable-package audit re-expanded all three THM-3482 lanes,
rederived the cubic Fourier decomposition and its degrees, proved minimality
from the distinct exponential colours, checked all recurrence and numerator
coefficients, reproduced both shortened-factor failures and the exact
`22x22` Hankel factorization, and verified the three harmonic coefficients.
It found and repaired two stale status/scope sentences left from THM-3482's
former provisional stage; no mathematical or computational defect remained.

This theorem proves no property of an LRC relation current, no physical edge
assignment, no bispectrum nonvanishing, no Jacobian statement, and no
LRC(14) conclusion.  Equation `(7)` is the minimal recurrence of THM-3482's
canonical private-gradient determinant word.
