---
id: THM-3606
title: "Exponent-two three-by-four scalar and singleton gate atlas"
status: >
  PROVED + FINITE-EXACT + VERIFIED-EXACT / AWAITING INDEPENDENT HOSTILE
  AUDIT.  Among the 149 connected oriented three-by-four additive-support
  words of THM-3603, exactly 38 admit an integer realization, a scalar fibre,
  a scalar-arm address, and both-zero-or-same-strict-sign weights on every
  singleton fibre.  There are exactly 178 such anchor/orientation schemes.
  Thirty schemes have a rectangle-exposed scalar double and reduce to the
  forbidden two-by-two sector, leaving 148 schemes on 31 words.  Only four
  size-nine words remain, each with one unexposed scalar scheme.  These are
  necessary weight gates only: no multi-address cancellation, Darboux pair,
  planar Jacobian counterexample, or proof of JC(2) is claimed.
source: kps-s188 / THM-3603 coefficient-gate continuation, 2026-08-21
audit: >
  Author exact audit.  The standard-library companion hash-pins THM-3603,
  solves every nonhomogeneous word-cone branch by rational RREF and a complete
  one-parameter congruence search, and independently scans every positive
  five-gap vector in [1,16]^5.  The global, generic-scaled, and bounded word
  sets agree.  Ordinary and optimized runs are byte-identical to the stored
  transcript after LF normalization.  Independent hostile audit remains
  pending.
depends_on:
  - THM-3603-three-by-four-additive-support-collision-cones-and-fibre-cut-atlas
  - THM-3592-universal-exponent-two-three-by-three-weight-darboux-nonentry
  - THM-3569-danielewski-two-by-three-weight-darboux-nonentry
related:
  - THM-3583-universal-exponent-two-two-by-four-weight-darboux-nonentry
  - THM-3600-danielewski-arm-plane-atlas-singular-shear-and-no-filling
script: 04-computation/jc2_three_by_four_scalar_singleton_gate_atlas_thm3606.py
output: 05-knowledge/results/jc2_three_by_four_scalar_singleton_gate_atlas_thm3606.out
script_sha256: 9d558f5637c5cc3573214fe00817bd6acb4c749b85a9046644e2390bd1d0ad91
output_sha256: 58b1160bd8831139fc5a6d4eb5102c876df4b47a823293080413c8f6c3f995b4
semantic_sha256: c21d12829dc415dded4fcd2347c23137991904da598cbe4cd8ad4e8844f5c1f0
bounded_semantic_sha256: 68437764f3eb632545e70224b265f1e74b01a30826d2cac5c6658e150358ff07
hash_basis: LF-normalized bytes
---

# THM-3606 -- exponent-two three-by-four scalar and singleton gate atlas

**PROVED + FINITE-EXACT + VERIFIED-EXACT / AWAITING INDEPENDENT HOSTILE
AUDIT.**  This is the first coefficient-independent weight gate after the
additive classification of THM-3603.  It is global over positive integral
weight gaps.  It deliberately stops before the coupled polynomial coefficient
equations.

All rings are over `C`.  Let `Sigma` be squarefree with `deg Sigma>=2` and use
the exponent-two Danielewski grading

```text
A_Sigma=C[b,c,e]/(c^2 e-Sigma(b)),       wt(b,c,e)=(0,1,-2).       (1)
```

## 1. The finite input and the two gates

After translating the two supports independently, write

```text
A={0,X,X+Y},
B={0,U,U+V,U+V+W},                  X,Y,U,V,W in Z_(>0).           (2)
```

For initial weights `p,q`, the actual weights of the three pieces of `P` and
four pieces of `Q` are

```text
p_i=p+a_i,                   q_j=q+b_j.                            (3)
```

The fibre at total support value `s` is

```text
F_s={(i,j):a_i+b_j=s}.                                             (4)
```

THM-3603 proves that every putative `3 x 4` Darboux pair has one of 149
connected oriented fibre words, and that every corresponding word cone has
dimension at most two.  We now impose the following necessary conditions from
THM-3592.

1. A fibre which supplies the scalar bracket must contain an address with

   ```text
   (p_i,q_j)=(-2,1) or (1,-2).                                    (5)
   ```

2. If a non-scalar fibre consists of a single address, its bracket must vanish.
   For retained nonzero coefficient polynomials this requires

   ```text
   p_i,q_j both strictly positive, both strictly negative, or both zero.     (6)
   ```

   Strict opposite signs are impossible by logarithmic derivatives, and one
   zero weight would make the weight-zero coefficient removable.

A **scheme** is a choice of scalar fibre, one anchored address in that fibre,
and one of the two orientations in `(5)`, together with at least one positive
integral five-gap realization satisfying `(6)` at every singleton fibre.

## 2. Exact affine reduction

The fixed numbers `-2` and `1` make this an affine, not projective, problem.
For example, in orientation `(-2,1)` anchored at `(i,j)`, a singleton
`(k,l)` has weights

```text
(-2+a_k-a_i, 1+b_l-b_j).                                          (7)
```

Because both supports are increasing, condition `(6)` becomes a finite list
of integral linear alternatives.  Away from the cross quadrant these are
ordinary inequalities.  Two representative cases are

```text
k<=i, l<j:       b_j-b_l>=2,
k> i, l>=j:      a_k-a_i>=3.                                     (8)
```

The quadrant `k>i,l<j` has precisely the two possibilities

```text
a_k-a_i=1 and b_j-b_l>=2,       or
a_k-a_i=2 and b_j-b_l=1.                                        (9)
```

The other orientation is the reflected table.  In particular, every
nonhomogeneous alternative in `(9)` fixes an interval distance to `1` or `2`.
Since a word cone has dimension at most two, adjoining one such equality
leaves an affine line or point.

The companion solves these residual systems exactly over `Q`.  On a line it
writes

```text
gaps=alpha t+beta,                                                  (10)
```

intersects all inequalities to obtain an integer interval for `t`, and checks
one complete residue period, the least common multiple of the denominators in
`alpha,beta`.  This is exhaustive: integrality repeats with that period and
the inequalities define an interval.  If no soft alternative occurs, three
times the positive integral chamber sample satisfies every threshold in
`(8)`.  Thus every scheme is decided globally; no bounded stabilization is
used for the theorem.

The exact result is

| `|A+B|` | surviving words | schemes |
|---:|---:|---:|
| 6 | 1 | 12 |
| 7 | 7 | 66 |
| 8 | 26 | 90 |
| 9 | 4 | 10 |
| **total** | **38** | **178** |

The two scalar orientations occur 90 and 88 times.  The scalar fibre is a
double in 161 schemes and a triple in 17.  Every scheme has an exact witness
whose largest gap is at most 12.  The transcript lists every word, every
scheme digest, and a least stored witness.

## 3. The projective-quotient trap

Testing only the primitive chamber representative retains 18 words, not 38.
Testing three times every representative retains exactly the global set of
38.  This discrepancy is not a numerical defect: common dilation preserves
the fibre word but does not preserve `(7)`, because the scalar-arm anchor fixes
absolute weights.

Therefore a projective additive atlas is not by itself a classifier for an
affine scalar gate.  The absolute scale is a required sidecar.  This is the
same logical type of quotient loss that forced THM-3592 to retain orientation,
although the lost coordinate here is dilation rather than reversal.

## 4. Rectangle-exposed scalar fibres are impossible

Let the scalar fibre be the double

```text
F={(i,j),(k,l)},                i!=k, j!=l.                         (11)
```

Call it rectangle-exposed when the two cross-corners `(i,l)` and `(k,j)` are
both singleton fibres.  Their bracket rows vanish by `(6)`.  Retain only the
two `P` pieces indexed by `i,k` and the two `Q` pieces indexed by `j,l`.
The two diagonal addresses in `(11)` still give the complete scalar row, while
the cross-corners give zero.  The truncation would therefore be a `2 x 2`
Darboux pair, contradicting THM-3569.  Hence every rectangle-exposed scalar
scheme is impossible.

Exactly 30 of the 178 schemes are rectangle-exposed:

| `|A+B|` | exposed schemes | unexposed schemes |
|---:|---:|---:|
| 6 | 0 | 12 |
| 7 | 0 | 66 |
| 8 | 24 | 66 |
| 9 | 6 | 4 |
| **total** | **30** | **148** |

Seven of the 38 words lose every scheme.  Thus a putative exponent-two
`3 x 4` Darboux pair must lie among exactly 148 schemes on the following 31
oriented words:

```text
W001 W002 W003 W004 W005 W006 W007 W008
W009 W010 W011 W012 W015 W016 W022 W030
W031 W032 W033 W034 W036 W037 W040 W041
W042 W043 W045 W072 W073 W141 W149.                  (12)
```

These labels are the canonical ordering in the hash-pinned THM-3603 atlas.

## 5. The four sharp size-nine remnants

Every surviving size-nine word has exactly one unexposed scheme, and its
scalar fibre is a one-fibre projection cut:

| word | scalar anchor/orientation | gap witness | fibre word |
|---|---|---|---|
| `W072` | `F_3:11,+-` | `(9,3,3,6,6)` | `00|01|02=10|11=20|03=21|12|22|13|23` |
| `W073` | `F_3:11,+-` | `(9,6,6,3,12)` | `00|01|02=10|11=20|12|03=21|22|13|23` |
| `W141` | `F_5:12,-+` | `(6,9,12,3,6)` | `00|10|01|02=20|11|03=12|13=21|22|23` |
| `W149` | `F_5:12,-+` | `(3,9,6,6,3)` | `00|10|01|11|02=20|03=12|13=21|22|23` |

Simultaneous reversal pairs `W072` with `W149` and `W073` with `W141`, but the
four orientations are retained because regularity is reversal-sensitive.
Deleting the scalar collision disconnects a projection graph in each case;
this is not yet a legal deletion of polynomial pieces.  It identifies the
next exact target: combine the unique bridge with singleton common-base
factorizations and the all-arm Hermite--Pade divisibility of THM-3600.

## 6. Independent hostile control and verification

Independently of the affine solver, the companion scans every positive
five-gap vector in `[1,16]^5`, without primitive quotienting.  Among 6,724
vectors having connected projection graphs, 1,288 vectors and 3,790 scalar
placements survive the two gates.  Their 38 fibre words agree exactly with
the global set.  Counts of surviving vectors by sumset size are

```text
6:16,             7:54,             8:1064,             9:154.    (13)
```

This is a hostile representation check, not the proof of globality.  Reproduce
the exact atlas with

```bash
python3 04-computation/jc2_three_by_four_scalar_singleton_gate_atlas_thm3606.py
python3 -O 04-computation/jc2_three_by_four_scalar_singleton_gate_atlas_thm3606.py
```

Both runs are byte-identical to the stored transcript after LF normalization.
The companion hash-pins the THM-3603 implementation and pins both the global
scheme ledger and the bounded control by independent semantic digests.

## 7. Boundary

This theorem imposes only the scalar-arm and singleton-row necessary gates,
plus the lawful rectangle-exposed `2 x 2` reduction.  It does not decide
multi-address differential cancellation, possible variation of the active
scalar address among roots of `Sigma`, nonalternation across arms,
Hermite--Pade regularity, existence of a Darboux pair, or the planar Jacobian
conjecture.

**QED.**
