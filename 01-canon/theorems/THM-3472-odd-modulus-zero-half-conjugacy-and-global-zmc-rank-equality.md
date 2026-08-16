---
id: THM-3472
title: "Odd-modulus zero/half conjugacy and global ZMC-rank equality"
status: >
  RESERVED / PROVISIONAL PROOF CANDIDATE + VERIFIED-EXACT / INDEPENDENT
  AUDIT PENDING.  For every odd q>=3, the full zero-mode-cochain rank equals
  the literal half-twist cover rank, including infinity.  Consequently the
  full odd cap-seven atlas has exact ranks 4,6,7 and >7 with no rank five;
  its annotated word has minimal period 729664650 and exact natural/harmonic
  coefficients.  No even-modulus, endpoint-current, or LRC(14) conclusion
  follows.
source: codex-2026-08-15-odd-layer-conjugacy
audit: >
  self-contained sheet-bijection, sign normalization, primitivity, divisor
  descent, and dilation proof; exact 10787000-cell conjugacy, 2613750-cell
  dilation, weighted-CRT density, minimal-period, dependency, semantic,
  security, and normal/optimized replay gates; independent audit pending
depends_on:
  - THM-3405-common-centre-gcd-gauge-and-boolean-half-twist
  - THM-3453-global-literal-half-twist-cap-seven-support-classification
related:
  - THM-3464-u-spine-q123-rank-eight-break-and-divisor-layer-certificate
  - THM-3469-three-times-p-half-twist-eight-owner-cover-boundary
script: 04-computation/lrc_odd_zero_half_conjugacy_global_rank_thm3472.py
output: 05-knowledge/results/lrc_odd_zero_half_conjugacy_global_rank_thm3472.out
script_sha256: 039bf8871f04be15ead6ac8725033c81ff9782371adeffe9a01f3926ba126ab0
output_sha256: b52a1d1d6767db9110991f40c0d4beadd3d49d870c449ceb5a48cae75d5f4269
semantic_sha256: d8bc9ad4a49f954ec1c76db01a7506a5965f0dc5f58881bb7315de402d151221
hash_basis: LF-normalized bytes
---

# THM-3472 -- odd-modulus zero/half conjugacy and global ZMC-rank equality

**RESERVED / PROVISIONAL PROOF CANDIDATE + VERIFIED-EXACT / INDEPENDENT
AUDIT PENDING.**

The proof and deterministic companion pass their internal gates.  This file
remains outside the proved dependency graph until an independent audit checks
the universal layer conversion and immutable package.

## 1. Two ranks and the inheritance problem

For `Q>=2`, define the fixed-zero and literal half-twist masks

```text
Z_(Q,s)={j in Z/QZ:14 dist_Q(sj,0)<Q},                (1)

B_(Q,r)={ell in Z/QZ:
  14 dist_(2Q)(r(2ell+1),0)<2Q}.                      (2)
```

Owners are transverse, sign-normalized, and distinct.  Let `rho_H(q)` be
THM-3453's minimum number of masks `(2)` covering every sheet.  Let
`rho_ZMC(q)` be THM-3405's minimum number of selected modes at one common
centre with complete cochain zero.  Both ranks take values in the positive
integers together with infinity when no transverse cover exists.

Every half-twist cover is a zero-cochain cover at the fixed common centre
`1/(2q)`, so

```text
rho_ZMC(q)<=rho_H(q).                                  (3)
```

THM-3405 first reduces an arbitrary zero-cochain family to a primitive
fixed-zero or half-twist cover on some `Q|q`.  Oddness supplies the missing
conversion back to the literal half layer.

## 2. Odd fixed-zero masks conjugate into the half layer

Assume `Q` is odd.  The affine sheet map

```text
phi_Q(ell)=2ell+1 mod Q                                (4)
```

is a permutation.  Canonicalize a fixed-zero owner to
`1<=s<=(Q-1)/2`; replacing `s` by `Q-s` does not change `(1)`.  Then `2s<Q`
is a transverse half owner and

```text
B_(Q,2s)(ell)=Z_(Q,s)(phi_Q(ell)).                    (5)
```

Indeed, for every integer `a`,

```text
dist_(2Q)(2a,0)=2 dist_Q(a,0),                        (6)
```

so the two strict inequalities are identical after `(4)`.  Distinct
canonical `s` give distinct doubled owners.  Moreover, for every owner family
`S`,

```text
gcd(Q,{2s:s in S})=gcd(Q,S),                          (7)
```

because `Q` is odd.  Thus the conversion preserves primitivity as well as
cardinality and coverage.  This direction is all that is needed; `(5)` does
not assert that every odd half owner arises from the same fixed-zero chart.

## 3. Equality of the global ranks

Let `q>=3` be odd and suppose a zero-cochain cover uses `k` owners.  By
THM-3405's active-gcd reduction, some primitive divisor `Q|q` carries a
`k`-owner cover in one of the two Boolean layers.

- If it is already the half layer, retain it.
- If it is the fixed-zero layer, apply `(4)--(7)` to obtain a primitive
  half-twist cover on `Q` with at most `k` owners.

Write `q=lambda Q`.  THM-3405 dilation is visible directly in `(2)`:

```text
B_(lambda Q,lambda r)(ell)
  =B_(Q,r)(ell mod Q),                                 (8)
```

because cyclic numerator distance scales by `lambda`.  The divisor cover
therefore pulls back to a half-twist cover of all `q` sheets with at most `k`
owners.  Hence `rho_H(q)<=rho_ZMC(q)`.  Together with `(3)`, this proves

```text
rho_ZMC(q)=rho_H(q) for every odd q>=3.                (9)
```

The equality includes infinity.  It is an equality of minimum grades, not an
identification of ancestry: as THM-3464 shows at `q=123`, distinct primitive
quotient layers may attain the same minimum.

## 4. The full odd cap-seven atlas

THM-3453's all-modulus half-twist atoms have ranks

```text
rank 4:  8,9;
rank 5:  10,12;
rank 6:  11,15,23,25;
rank 7:  13,14,29,38,51,68,148.                       (10)
```

Only odd atoms can divide odd `q`.  Equation `(9)` therefore promotes the
literal classification to the full zero-cochain rank:

```text
rho_ZMC(q)=4  iff 9|q;

rho_ZMC(q)=5  never;

rho_ZMC(q)=6  iff 9 does not divide q and
  one of 11,15,23,25 divides q;

rho_ZMC(q)=7  iff no lower-rank atom divides q and
  one of 13,29,51 divides q;

rho_ZMC(q)>7  otherwise,                              (11)
```

for every odd `q>=3`.  Equation `(11)` does not classify the exact rank
inside the final `>7` stratum.

## 5. Exact subsets of the harmonic series

Define an annotated word on all positive integers by writing `0` at even
indices and `4,6,7,>7` at odd indices according to `(11)`, with `>7` placed
at the out-of-domain index `q=1` by convention.  Its minimal period is

```text
P=2*3^2*5^2*11*13*17*23*29=729,664,650.              (12)
```

The powers `3^2` and `5^2` retain the atoms `9` and `25`; the factor `17`
retains `51`.  A weighted CRT product has only `576` local valuation states.
In one period the census is

| stratum | count | natural density / harmonic coefficient |
|---:|---:|---:|
| even marker `0` | `364,832,325` | `1/2` |
| rank 4 | `40,536,925` | `1/18` |
| rank 5 | `0` | `0` |
| rank 6 | `64,859,080` | `4/45` |
| rank 7 | `31,171,360` | `283376/6633315` |
| rank `>7` | `228,264,960` | `691712/2211105` |

For example, after excluding rank four, the odd density is `4/9`; excluding
the four rank-six atoms leaves `16/45`, so rank six has density `4/45`.
The final-stratum density factors as

```text
1/2 * 3088/3825 * 10/11 * 22/23 * 12/13 * 28/29
  =691712/2211105.                                    (13)
```

Rank seven is the difference between `16/45` and `(13)`.  The companion gives
a changing shift witness for every prime factor of `(12)`, proving
minimality:

```text
(2,1,>7,364832326,0),
(3,3,>7,243221553,4),
(5,25,6,145932955,>7),
(11,7,>7,66333157,6),
(13,13,7,56128063,>7),
(17,51,7,42921501,>7),
(23,17,>7,31724567,6),
(29,29,7,25160879,>7).                                (14)
```

Consequently each displayed density `delta` is also the coefficient in

```text
sum_(q<=N,q in stratum) 1/q=delta log N+O(1).          (15)
```

Conditioned on odd `q`, the coefficients are twice those in the table.  This
is a complete cap-seven subset-of-the-harmonic-series theorem; it does not
turn the unclassified `>7` stratum into one exact rank.

## 6. Sharp boundary and information contract

Oddness is load-bearing for the proof.  At `Q=8`, the map `(4)` reaches only
the four odd sheets, while `Z_(8,4)` occupies four even sheets and transports
to the empty half mask.  This refutes use of `(4)` as an even-modulus sheet
conjugacy.  It does not assert that rank equality itself fails at every even
modulus; no even conclusion is made.

```text
source:      a primitive fixed-zero or half-twist cover on Q|q
target:      a literal half-twist cover on the ambient odd q
maps:        ell -> 2ell+1, owner doubling/sign, divisor dilation
preserved:   strict mask incidence, cover size, active gcd, zero cochain
destroyed:   original layer label and primitive-divisor ancestry after grading
sidecars:    odd modulus, canonical owner sign, active divisor Q
hostile:     Q=8 nonbijective sheet map
```

The theorem gives no endpoint current, relation-residue coefficient,
bispectrum, physical LRC row, decrement, or LRC(14) conclusion.

## 7. Exact companion

Run from the repository root:

```bash
python 04-computation/lrc_odd_zero_half_conjugacy_global_rank_thm3472.py
python -O 04-computation/lrc_odd_zero_half_conjugacy_global_rank_thm3472.py
```

The standard-library companion checks every `40,200` owner row and
`10,787,000` sheet cells for odd `3<=Q<=401`, including active-gcd
preservation; `7,650` divisor-dilation rows and `2,613,750` cells at scales
`3,5,7`; the exact `576`-state CRT census, all density reductions, the eight
minimal-period witnesses, dependency pins, AST/security, and a frozen
semantic digest.  It uses explicit exceptions under `-O` and performs no file
write, dynamic evaluation, subprocess, or network action.
