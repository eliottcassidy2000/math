---
id: THM-3103
title: "Low-depth diameter-forty fixed Macaulay atlas"
status: >
  PROVED + VERIFIED-EXACT + INDEPENDENTLY HOSTILE-AUDITED.  One fixed
  36-by-36 degree-seven Macaulay minor is nonzero modulo 1000003 for every
  four-slot factorial support {n,n+a,n+b,n+M} with 0<=n<=3 and
  0<a<b<M<=40.  This exactly certifies first-window SFC(4) on all 39520
  cells.  Of these, 38640 lie beyond the diameter-twelve all-depth atlas.
  The result is finite and modular: it proves rational nonvanishing, not a
  determinant sign, a continuous-depth statement, or an all-depth theorem.
source: root-low-depth-wide-macaulay-atlas-2026-08-02
audit: >
  An independent hostile audit accepted the rational-to-modular implication,
  denominator-unit bound, 20Q+10C+6F fixed-row typing, projective evaluation
  argument, exact universe and transverse-new-cell counts, and the positive-
  depth THM-3096 radial corollary.  It independently compared 64 additional
  deterministic cells with the pinned rational inclusion-exclusion
  constructor.  Normal, optimized, and stored runs agree; both LF hashes and
  repeated-exponent zero controls pass.
depends_on:
  - THM-3096-physical-support-complete-intersection-and-arbitrary-radial-pair-radical
related:
  - THM-2921-diameter-four-nonconsecutive-macaulay-newton-closure
  - THM-2947-conjugate-pair-corank-parity-and-one-minor-resultant-gate
  - THM-2949-fixed-rank-thirty-five-cofactor-newton-atlas
  - THM-2952-fixed-rank-thirty-five-cofactor-width-thirteen-continuation
  - THM-3097-translated-support-monge-compactification-and-cofinite-bad-set-induction
script: 04-computation/gmc_low_depth_diameter_forty_macaulay_atlas_thm3103.py
output: 05-knowledge/results/gmc_low_depth_diameter_forty_macaulay_atlas_thm3103.out
script_sha256: 06e0b04c4332c2fea189979bd1d7e9fbebf1c9750c9255731975161ff7203960
output_sha256: e96b7b58921bcdf185b33b33836362a0ad438afbd94357c6c3caecbe77979ce3
hash_basis: LF-normalized bytes
---

# THM-3103 -- low-depth diameter-forty fixed Macaulay atlas

**PROVED + VERIFIED-EXACT + INDEPENDENTLY HOSTILE-AUDITED.**

The existing finite SFC(4) atlases are deepest in the translation direction:
THM-2949 proves every depth for diameter at most twelve.  The complementary
finite question is whether shallow supports remain good after their internal
gaps become much larger.  A single unchanged Macaulay chart answers that
question through diameter forty and depths zero through three.

## 1. Exact finite statement

Let

```text
L(s^j)=j!,                       f_j=s^j/j!,             (1)
```

and fix integers

```text
0<=n<=3,                  0<a<b<M<=40.                  (2)
```

For every nonzero polynomial

```text
H=c_0 f_n+c_1 f_(n+a)+c_2 f_(n+b)+c_3 f_(n+M),         (3)
```

at least one of

```text
L(H),                  L(H^2),                  L(H^3), L(H^4) (4)
```

is nonzero.  Equivalently, the physical first-window resultant of every
support in `(2)` is nonzero.

There are exactly

```text
4 sum_(M=3)^40 binom(M-1,2)=4 binom(40,3)=39520         (5)
```

support-depth cells.  THM-2949 already covers the `880` cells with `M<=12`;
the other

```text
4[binom(40,3)-binom(12,3)]=38640                        (6)
```

are a new, transverse enlargement of the exact atlas.

## 2. The fixed chart

Because `L(f_j)=1`, eliminate the first moment by

```text
c_3=-(x+y+z).                                           (7)
```

After division by the positive base moment `L(f_n^r)`, let

```text
Q_(n,a,b,M), C_(n,a,b,M), F_(n,a,b,M)                  (8)
```

be the resulting ternary forms of degrees `2,3,4`.  Form their degree-seven
Macaulay rows in the inherited order

```text
21 quadratic rows, 15 cubic rows, 10 quartic rows.       (9)
```

Use the same selected global rows as THM-2921:

```text
0,...,19; 21,...,29,35; 36,...,41.                     (10)
```

Thus `(10)` gives one square `36`-row chart against all ternary monomials of
degree seven.  Call its determinant `Delta_(n,a,b,M)`.

The exact companion proves

```text
Delta_(n,a,b,M) != 0 mod 1000003                       (11)
```

for every cell `(2)`.  The prime exceeds every factorial argument in the
construction:

```text
4(n+M)<=172<1000003.                                   (12)
```

Consequently every coefficient denominator is a unit modulo the prime.  If
the rational determinant were zero, its reduction would be zero; hence
`(11)` proves

```text
Delta_(n,a,b,M) != 0 in Q.                              (13)
```

This is a one-way exact modular certificate, not a probable-prime heuristic
and not a claim that every nonzero rational minor stays nonzero modulo this
prime.

## 3. Why one full-rank minor proves SFC(4)

Suppose `(4)` vanished.  Equation `(7)` would give a nonzero projective point

```text
[x:y:z] in P^2
```

common to `Q,C,F`; the eliminated vector cannot be zero unless all four
original coefficients vanish.  Evaluation at this point annihilates every
degree-seven multiple of the three forms.  The degree-seven Macaulay row
space would therefore have rank less than the `36`-dimensional target space,
contradicting `(13)`.  This proves `(4)`.

No conjugate-pair parity inference is needed here: the displayed maximal
minor itself has full rank.  Conversely, this chart may vanish outside the
declared universe even when another chart or the resultant survives.  The
theorem asserts only the implication certified by `(11)`.

## 4. Arbitrary-radial two-charge corollary

For the positive-depth subatlas

```text
1<=n<=3                                                   (14)
```

all four support exponents are positive.  THM-3096 therefore applies
literally.  For every such support and arbitrary complex radial coefficients
on it, the two-charge family

```text
P=alpha conjugate(Z)+Z C(Z conjugate(Z))                (15)
```

has moment-null variety exactly

```text
{alpha=0} union {C=0},                                  (16)
```

and moments through degree eight give the sharp pair-radical certificate
with exponent seven.  This covers `3 binom(40,3)=29640` arbitrary-radial
support-depth cells.  The `n=0` row proves factorial SFC(4), but is not put
into `(15)` because THM-3096's polynomial radial support is strictly
positive.

## 5. Exact evidence and boundaries

The companion is dependency-light and auditable in four ways.

1. Its fast constructor expands the original four coefficient variables and
   only then makes the substitution `(7)`.
2. A hash-pinned THM-2921 constructor instead uses normalized rising
   factorial tensors and inclusion-exclusion.  The two constructors agree
   coefficientwise in `15` representative forms and on `5` selected minors.
3. Three repeated-exponent controls make the same `36`-minor exactly zero,
   showing that the computation is not a tautological rank assertion.
4. Every one of the `39520` residues is folded into depthwise and global
   immutable SHA256 digests.

Run

```text
python 04-computation/gmc_low_depth_diameter_forty_macaulay_atlas_thm3103.py
python -O 04-computation/gmc_low_depth_diameter_forty_macaulay_atlas_thm3103.py
```

Both modes must equal the stored transcript after LF normalization.

The theorem proves no determinant sign, real-depth positivity, all-depth
width-forty statement, SFC(5), arbitrary support theorem, or closure of the
Gaussian Moment Conjecture.  Its role is a substantially wider exact finite
core for THM-3097/THM-3099 compactification and zero-escape mechanisms.

**QED.**
