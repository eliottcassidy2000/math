---
id: THM-3490
title: "Full Fourier-Dirichlet residue tomography of periodic-polynomial words"
status: >
  PROVED + VERIFIED-EXACT + INDEPENDENTLY AUDITED.  The complete table of
  character-twisted Dirichlet residues
  is exactly the Fourier coefficient array of a periodic-polynomial tail and
  recovers every lane coefficient.  Its minimal recurrence is the lossy
  projection retaining only each character row's highest occupied rung.
source: codex-2026-08-16-full-fourier-dirichlet-residue-tomography
audit: >
  independently audited by death-star-2026-08-16 after integration: the
  Hurwitz-residue normalization, Fourier inversion, exact information order,
  declared-period refinement, finite-prefix loss, demodulate-deflate
  reconstruction, rational descent, and THM-3484/THM-3487 scopes all pass;
  exact companion checks 82 split-field packets, direct Berlekamp-Massey
  recovery, alternating/composite-period/prefix hostiles, normal/optimized
  replay, and AST security.  The audit also repaired a stale script hash.
depends_on:
  - THM-3485-periodic-polynomial-fourier-jordan-recurrence-classification
  - THM-3486-critical-harmonic-transform-of-periodic-polynomial-words
related:
  - THM-3364-cyclotomic-boolean-clocks-berggren-t4-xor-and-crt-phase
  - THM-3487-two-twenty-four-state-fibonacci-bundles-cycle-type-obstruction
script: 04-computation/periodic_polynomial_fourier_dirichlet_residue_tomography_thm3490.py
output: 05-knowledge/results/periodic_polynomial_fourier_dirichlet_residue_tomography_thm3490.out
script_sha256: 52a58263a2ea00af72a1038d1d64c955e411b0f308ba7fa94621439f00d79d43
output_sha256: 72f7b97f3acabeef953c1184a35904b1fa28d0494cb9f4d0e2100c85aa3ebc3e
semantic_sha256: 39daa9a75e84e0fcb28141fa9079c3f616c9aa5d0689934cc239b4d911fc713e
hash_basis: LF-normalized bytes
---

# THM-3490 -- every Fourier/Jordan rung is a Dirichlet residue

**PROVED + VERIFIED-EXACT + INDEPENDENTLY AUDITED.**

THM-3486 identifies one scalar: the top trivial Fourier coefficient is the
critical harmonic mass.  Keeping every character and every possible pole
height gives a lossless transform of the whole periodic-polynomial tail.
The minimal recurrence is then a strict, lossy projection of that table.

## 1. Setup

Fix `p>=1`, a primitive `p`th root `zeta`, and lane polynomials

```text
P_r(t)=sum_(m=0)^d c_(r,m)t^m,          0<=r<p,          (1)
```

where zero top coefficients are allowed and not all lanes are zero.  Define

```text
a_n=P_r(n) when n=r mod p,                               (2)
Q_j(t)=1/p sum_(r=0)^(p-1) zeta^(-jr)P_r(t).             (3)
```

Thus THM-3485's Fourier inversion is

```text
a_n=sum_(j=0)^(p-1) zeta^(jn)Q_j(n).                    (4)
```

For each character `j`, form the twisted Dirichlet series

```text
D_j(s)=sum_(n>=1) zeta^(-jn)a_n/n^s,     Re(s)>d+1.      (5)
```

The period and the primitive-root/character gauge are part of the input.

## 2. Meromorphic residue formula

For `r in {0,...,p-1}`, let `r+` be its representative in `{1,...,p}`;
in particular `0+=p`.  Splitting (5) into arithmetic progressions gives

```text
D_j(s)=sum_(r=0)^(p-1) sum_(m=0)^d
 zeta^(-jr)c_(r,m)p^(m-s) HurwitzZeta(s-m,r+/p).         (6)
```

This finite expression meromorphically continues `D_j` to the plane.  Its
only possible poles are simple and lie at

```text
s=1,2,...,d+1.                                          (7)
```

Hurwitz zeta has residue one at argument one.  At `s=m+1`, the scale factor
in (6) is `p^(-1)`.  Therefore

```text
boxed:
R_(j,m):=Res_(s=m+1) D_j(s)
        =1/p sum_r zeta^(-jr)c_(r,m)
        =[t^m]Q_j(t).                                   (8)
```

This proves the residue tomography formula.

## 3. Exact information content

Fourier inversion in the lane coordinate gives

```text
c_(r,m)=sum_(j=0)^(p-1) zeta^(jr)R_(j,m).               (9)
```

Hence the labelled table

```text
R=(R_(j,m))_(0<=j<p,0<=m<=d)                            (10)
```

recovers every lane polynomial and therefore the entire
periodic-polynomial word.  The transform is lossless on the declared tail.

The sharp lost coordinate is a finite prefix.  Changing finitely many
`a_n` adds a Dirichlet polynomial to each `D_j`; that contribution is entire
and has no residues.  Thus (10) classifies periodic-polynomial **tails**, not
transient initial data.

Replacing `p` by a nonminimal multiple inserts redundant zero character
rows and relabels surviving roots.  The table depends on the declared period;
the recovered sequence and its minimal tail recurrence do not.

## 4. Recurrence is a lossy projection

For each character row set

```text
d_j=max{m:R_(j,m)!=0},                                  (11)
```

with the row omitted when all its residues vanish.  Substitution of (8) into
proved THM-3485 gives

```text
boxed:
chi_a(x)=product_(j:row j nonzero)(x-zeta^j)^(d_j+1).   (12)
```

Thus the pole table determines the minimal shift polynomial.  The converse
is false: (12) retains only the support and highest occupied height of each
row, forgetting all residue values and lower rungs.  Already `a_n=n` and
`b_n=2n` have the same minimal polynomial `(x-1)^2` and different residue
tables.  More generally, infinitely many words have any fixed nonzero
support-depth profile.

The correct information order is

```text
periodic-polynomial tail
  <-> labelled full residue table
   -> Fourier/Jordan support-depth profile
   -> minimal recurrence polynomial.                   (13)
```

Only the first arrow is reversible after the period and character gauge are
fixed.

## 5. Real-variable reconstruction without continuation

At the top degree, proved THM-3486 applied after character demodulation gives
for every `j`

```text
sum_(n<=N) zeta^(-jn)a_n/n^(d+1)
  =R_(j,d) H_N+C_j+O(1/N),                              (14)
```

and hence

```text
R_(j,d)=lim_(N->infinity) 1/log(N)
                  sum_(n<=N) zeta^(-jn)a_n/n^(d+1).     (15)
```

After recovering the degree-`d` row, subtract

```text
n^d sum_j R_(j,d)zeta^(jn)                              (16)
```

from `a_n` and repeat at exponent `d`.  Iterating downward gives a
demodulate--measure--deflate algorithm that reconstructs every row of (10)
using critical logarithmic limits alone.  Formula (6) packages all of those
steps into one meromorphic object.

THM-3486 is exactly entry `(j,m)=(0,d)` of this table.  Ordinary harmonic
mass is one coordinate, not the whole spectral object.

## 6. Sharp hostiles and specializations

### Alternating top block

For

```text
a_n=(-1)^n n^d                                           (17)
```

at period two,

```text
R_(0,d)=0,             R_(1,d)=1.                       (18)
```

The ordinary critical harmonic coefficient sees nothing, while the
nontrivial row recovers the full `(x+1)^(d+1)` block.  This is the cheapest
proof that the trivial residue is insufficient.

### Composite-period refinement

The period-four presentation of the parity word has residue-degree profile

```text
(-infinity,-infinity,0,-infinity);                       (19)
```

only the order-two character survives.  Its recurrence equals that of the
minimal period-two presentation.  Primitive order-four rows are exactly zero,
so refining the period creates no spurious factor.

### THM-3484's ternary word

Its degree-seven lane vector is constant `(-16384,-16384,-16384)`.  Hence

```text
(R_(0,7),R_(1,7),R_(2,7))=(-16384,0,0),                 (20)
```

while both nontrivial rows have nonzero degree-six residues.  Equations
(11)--(12) recover the proved profile `(7,6,6)` and order `22` without a
Hankel fit.

### Rational descent

For rational lane coefficients the rows lie in `Q(zeta)` and Galois acts by
permuting character labels.  Grouping rows by root order turns (12) into the
cyclotomic factorization of THM-3485.  The labelled table is gauge-dependent;
the rational recurrence polynomial is not.

## 7. Subsets of the harmonic series and the 24-state bundles

This degree-zero specialization recovers THM-3364's proved finite cyclotomic
transform of an eventually periodic Boolean support.  For a periodic Boolean
word, `d=0`, and (8) says that the residues at `s=1`
are exactly the full discrete Fourier transform of its address indicator.
The trivial residue is only its density.  All character residues together
recover the subset within one period.

Apply this carefully to THM-3487.  After declaring a uniform 24-state address
enumeration, the trivial residue of a shift-invariant union gives only the
half-grade or quarter-grade in THM-3487 (29).  The complete 24-character bank
recovers which cycle union was selected, even when two unions have the same
cardinality.  A cyclic reindexing multiplies nontrivial rows by phases, so an
address origin remains a load-bearing gauge.

This is not a time average along one dynamical orbit.  It also does not
classify arbitrary subsets of `N`: nonperiodic subsets may have no density,
and their reciprocal subseries may converge or diverge.

## 8. Cross-frontier boundary

The useful operation-level lesson is to demodulate before taking a critical
residue.  It parallels the live LRC endpoint warning that a trivial or coarse
character can identify a forced bridge that a refined character separates.
But the objects differ: (5) is an index-word Dirichlet transform, not an LRC
endpoint response or a Jacobian divisor norm.

Likewise, the residue table is a spectral `H^0`-type address observable.  It
does not by itself construct THM-3487's seam class in `H^1(C_6;V4)`, much less
the desired LRC word-current to Jacobian-flux `H^1` map.  Character labels,
coefficient local systems, ancestry, and physical-current typing must still
be supplied.

This theorem proves no tournament classification, full Berggren-tree
transplant, LRC bispectrum nonvanishing, Jacobian counterexample
classification, or case of LRC(14).

## 9. Exact companion

Run

```bash
python -B 04-computation/periodic_polynomial_fourier_dirichlet_residue_tomography_thm3490.py
python -B -O 04-computation/periodic_polynomial_fourier_dirichlet_residue_tomography_thm3490.py
```

The companion works independently in the split field `F_2521`, whose
multiplicative group supports every period `1..9`.  On `82` structured
packets it checks exact Fourier inversion, extracts (11)--(12), and compares
the prediction with direct Berlekamp--Massey recovery.  It separately checks
the alternating, composite-period, seven-term-prefix, and THM-3484 hostiles.
Normal and optimized replays agree.  These finite controls support but do not
replace the proof of (6)--(12).

The companion verifies the finite controls; the all-period statement follows
from the Hurwitz-residue and Fourier proof in Sections 2--4.
