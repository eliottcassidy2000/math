---
id: THM-2973
title: "First-gap wall-stripped norm-core continuation through width thirty-one"
status: >
  PROVED + VERIFIED-EXACT + INDEPENDENTLY HOSTILE-AUDITED.
  For supports {n,n+1,n+2,n+M}, 27<=M<=31, the two THM-2969 Macaulay
  charts have the same explicit common flag and the same primitive
  wall-stripped pure-resultant core. Every core is dense, coefficient-positive,
  and PF2. At M=31 the correction support must include the simple quartic
  wall (31,25). M=32 and arbitrary radial GMC(2) are outside scope.
source: codex-gmc-first-gap-norm-core-continuation-2026-07-30
depends_on:
  - THM-2969-first-gap-wall-stripped-resultant-norm-core-atlas
  - THM-2964-general-pure-factorial-moment-resonance-ladder
related:
  - THM-2942-macaulay-extraneous-flag-factor-and-pluecker-mutation
  - THM-2943-width-seven-eight-two-chart-macaulay-resultant-closure
script: 04-computation/gmc_first_gap_wall_stripped_norm_core_continuation_thm2973.py
output: 05-knowledge/results/gmc_first_gap_wall_stripped_norm_core_continuation_thm2973.out
script_sha256: 3539d1fcd7dd317659aab8296a86ab928999a7efff237e3f00571481a80cdad3
output_sha256: 7146c4bb062d453de50af3285a9d8f7764184b4716fb3d690ee2eed78d5e7bd8
hash_basis: LF-normalized bytes
---

# THM-2973 -- first-gap wall-stripped norm-core continuation through width thirty-one

**PROVED + VERIFIED-EXACT + INDEPENDENTLY HOSTILE-AUDITED.**

## Statement

Let `M` be an integer with `27<=M<=31`, let `n>=0`, and let `H` be a
nonzero polynomial supported on

```text
{n,n+1,n+2,n+M}.                                          (1)
```

Use the original and stable-mutated full Macaulay charts of THM-2969,
translated to normalized support `(0,1,2,M)`. Their specialized determinants
have the common flag

```text
H_M(n)=(n+M) product_(r=3)^floor(M/2) (n+r).               (2)
```

Let `q_M,c_M,B_M,E_M` be the pure quadratic, cubic, local Smith, and seam
factors in the notation of THM-2969. For `M<=30` retain its correction
support. At `M=31` enlarge the local Smith correction by the simple quartic
wall

```text
(M,r)=(31,25).                                             (3)
```

After exact division of the common chart gcd by

```text
q_M^5 c_M B_M E_M,                                        (4)
```

the common gcd leaves, up to a nonzero rational unit, the primitive
wall-stripped pure-resultant core `N_M(n)`. Equivalently, both individually
divided charts contain `N_M`, with their calibrated cofactors `U_M,V_M`.
For every `27<=M<=31`:

1. `N_M` is dense and every coefficient is strictly positive;
2. its coefficient sequence is PF2;
3. the calibrated original/mutated chart pair is coprime in `Q[n]`;
4. the identity `R_M=W_M N_M`, together with coefficient positivity, implies
   that `L(H),L(H^2),L(H^3),L(H^4)` cannot all vanish.

The exact core atlas is

| `M` | degree `N_M` | SHA-256 of primitive coefficient tuple |
|---:|---:|:---|
| 27 | 556 | `c00583d9eded605bd99fed73892a358e91224a7d4086a19c0a55f0b394285b7f` |
| 28 | 577 | `9a6b8b02f59f2d1f31a7bc262863c4461aa8d1d0c89c1e4b6901237a68558591` |
| 29 | 599 | `98f7a7599769664a781edb7fa1f61bd5b4a71092c3fef345de74c66094a42e18` |
| 30 | 617 | `ba6fc2329242f5f222295ec261b36d3a57e353206716667d67ee26e38b17c7a1` |
| 31 | 639 | `c9ffa50924c9d0298fd4995b1343fb200d3da0c38e1fcbe07347d97b8908e81c` |

Thus the finite first-gap nonvanishing range of THM-2969 continues from
`M=26` through `M=31`. No statement for `M=32`, arbitrary width, a fixed-prime
gate, or arbitrary radial coefficients is asserted.

## Proof

### 1. Full charts and complementary wall division

The companion imports the immutable THM-2969 engine by its LF hash and, for
each width, reconstructs both full determinant interpolants from exact integer
evaluations. It uses three outside-grid paired depths per width, hence ten
interpolants, fifteen paired controls, and thirty individual determinant
evaluations. The full determinant degrees remain the THM-2942 values and the
specialized chart resultants agree exactly.

Equation `(2)` is recovered as the primitive chart gcd. The THM-2969 identity

```text
H_M divides B_M E_M,             W_M=(B_M E_M)/H_M          (5)
```

continues at all five widths after applying `(3)`. Since `R_M=W_M N_M`, the
raw common factor is

```text
q_M^5 c_M H_M W_M N_M = q_M^5 c_M B_M E_M N_M.         (5a)
```

Thus `q_M^5 c_M H_M W_M` is precisely the removed divisor `(4)`, not the
whole raw common factor. Exact polynomial division leaves `N_M`, the same
primitive core as the corresponding wall-stripped pure resultant.

### 2. The new wall and its exact invoices

The pure integer-wall census on `1<=r<=M` is empty at `M=27,28,30`. At
`M=29` there is the simple cubic wall `r=22`. At `M=31` the quadratic wall
`r=21` and quartic wall `r=25` occur simultaneously. Each has valuation one
in its pure coefficient. No assertion about unrelated complex roots of the
high-degree pure coefficient polynomials is made.

In the order

```text
(q,c,B,E,H,W,raw gcd,removed divisor,core),                (6)
```

the exact local invoices are

```text
(29,22): (0,1,19,4,0,23,24,24,0),
(31,21): (1,0,19,4,0,23,28,28,0),
(31,25): (0,0,20,4,0,24,24,24,0).                         (7)
```

The last row is the load-bearing repair. Without `(3)`, the old correction
support leaves one unremoved quartic factor and the `M=31` audit correctly
fails. With `(3)`, raw and removed orders agree and the core order is zero.
Hence the THM-2964 wall is an invoice mutation, not a zero or sign failure of
the primitive core.

### 3. Degree law, positivity, and coprimality

The corrected degree is

```text
deg N_M
 =23M-2 floor(M/3)-2 floor(M/2)-floor(2M/3)-3-e_M,
e_M=1_(M in {11,12,21,31}).                               (8)
```

Adding `e_M` back, the six-width baseline difference remains exactly `124`
for sources `M=21,...,25` and targets `M=27,...,31`. Direct primitive
normalization gives the five degrees and hashes in the table. Every
coefficient is positive, there are no internal zeros, and every adjacent PF2
minor is nonnegative. Exact gcd computation of the two calibrated chart
polynomials has degree zero at every width. This last coprimality check is a
separate chart-calibration statement; it is not the moment implication.

For real `n>=0`, coefficient positivity gives `N_M(n)>0`. The complementary
factor `W_M` is a product of positive linear factors on this ray, so the
identity `R_M=W_M N_M` gives `R_M(n)!=0`. By the defining property of the
ternary resultant, the denominator-cleared quadratic, cubic, and quartic
factorial-moment forms have no common projective zero. Restoring the
mean eliminated in the construction proves that the first four factorial
moments cannot all vanish, exactly as in THM-2969 Section 3.

## Exact evidence

The canonical reproduction is

```text
python 04-computation/gmc_first_gap_wall_stripped_norm_core_continuation_thm2973.py --output .scratch/thm2973.normal.out
python -O 04-computation/gmc_first_gap_wall_stripped_norm_core_continuation_thm2973.py --output .scratch/thm2973.opt.out
```

The stored transcript is 2,064 UTF-8/LF bytes. Its core/flag record digests
are

```text
core=cf5513fec8a47aea1b090dca3fa1d931d2a07704cdd576cd20626b197f409fb3,
flag=92a1d2eb0a13263f16884fdd796d4e910e048a74522ecb9f985121b72bceead3.
```

The frozen LF hashes are

```text
script=3539d1fcd7dd317659aab8296a86ab928999a7efff237e3f00571481a80cdad3,
output=7146c4bb062d453de50af3285a9d8f7764184b4716fb3d690ee2eed78d5e7bd8.
```

The independent audit rederived the full common-gcd identity
`G=q^5 c R H=q^5 c B E N`, all three local wall invoices, the corrected
degree law and six-width baseline, positivity/PF2 nonvanishing, and the
ternary-resultant implication for the first four moments. It separately
confirmed that calibrated-chart coprimality is a chart statement rather than
the moment argument, replayed normal and optimized execution against the
stored transcript, and recovered the displayed hashes.

The computation is finite evidence for `(1)` only. In particular, it neither
computes `M=32` nor proves that no later pure wall changes `(8)` again.

**QED.**
