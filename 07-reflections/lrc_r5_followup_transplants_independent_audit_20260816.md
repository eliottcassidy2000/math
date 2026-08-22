# Independent audit of the fixed-root and root-difference follow-up transplants

**Status: ACCEPT, SCOPED FINITE-EXACT.**  The two candidate packages are
mathematically sound in their stated linear scope: a frequency-independent
channel map followed by one common `C_13` circulant, equivalently one scalar
projective multiplier at each drift frequency.  The audit constructs no
`U_full` common-ancestry relation, Boolean coupling, physical current,
absolute `H^1` class, bispectrum, scalar-row exclusion, or LRC(14) theorem.

The independent companion is
`04-computation/lrc_r5_followup_transplants_independent_audit_20260816.py`,
with stored output
`05-knowledge/results/lrc_r5_followup_transplants_independent_audit_20260816.out`.
Neither candidate script is imported or edited.

## 1. Inheritance and hostile baseline

The closest proved source mechanism is THM-2594, whose exact
`N(u,q,ell,theta)` table is formed on one Boolean ancestry base before every
marginal and DFT.  Its canonical hostile is the fixed-absolute-root table:
it also fires, so affine theta slaving is not the unique source of
nonvanishing (MISTAKE-295).  The closest target mechanism is THM-3514's
rank-four `K4 x F_13` endpoint Walsh bank.  THM-3518 is the least-used
load-bearing sidecar here: it pins the target's `13^-3` normalizer and the
positive left-sheet phase `zeta^(q_t a)`, while proving that the endpoint
edges remain gradients with zero cycle pairing.

The audit therefore asks only whether specified linear images of the proved
THM-2594 source curve are projectively equivalent to the proved endpoint
Walsh curve.  It does not identify their supports.

## 2. Clean-room reconstruction

The verifier hash-pins and reruns

```text
THM-2594 source:
  lrc14_stage2_theta_contraction_opus_20260728.py
  09c43af0a0a56c7a0833bbfd13ed6a96bc5a7a3718aa1bc6b77a144bde101a06

THM-3514 atom source:
  lrc_ufull_owner_boundary_k4xf13_endpoint_factorization_independent_audit_20260816.py
  f89be10c65bb77270199f9399b155d5a2c82c0da121b3e8589fe3c1f7e9824fc
```

It obtains the full four-index THM-2594 table from the proved constructor,
not from either follow-up.  Independently, it materializes all `169`
unguarded endpoint atom tables and reconstructs `q_H` and `q_q5` twice:

1. direct summation over all `tau` slices; and
2. contraction against a separately formed translated guard kernel.

The routes agree entrywise.  Reversing the required left-sheet phase changes
all `52` active target buckets.  The reconstructed target has rank four and
all `52` Walsh-frequency entries nonzero.

## 3. Package A: fixed absolute root

For the slaved and absolute-root tables, the audit independently obtains

```text
physical support:             12/91, 18/91
C7 x C13 spectral support:    91/91, 91/91
rational row ranks:           3, 3
marked-union rational rank:   4
marked intersection:          2.
```

The same ranks occur in the certified split field, and all `91` spectral
coordinates survive independently at split primes `547`, `1093`, and
`2003`.  The candidate residue-table digests agree exactly; the audit also
retains the unreduced integer numerators for the rational-rank calculation.

The projective ledgers are

| source | source rank | wedge rank | nullity | source-annihilator dimension | excess |
|---|---:|---:|---:|---:|---:|
| slaved seven rows | 3 | 12 | 16 | 16 | 0 |
| absolute seven rows | 3 | 12 | 16 | 16 | 0 |
| marked union, fourteen rows | 4 | 16 | 40 | 40 | 0 |

All `56` marked fourth-sidecar selections have rank four and zero nullity.
Among the `127` nonempty binary absolute sidecars, `126` have rank four and
zero nullity.  Mask `127` is the unique zero sidecar, leaving rank three,
nullity four, and exactly four annihilators.  This independently explains
the exceptional mask by the vanishing `ell=0` row.

All twelve drift dilations and all thirteen common translations retain zero
excess.  A stronger hostile allows the slaved and absolute banks a relative
phase before applying each common dilation.  Across all `12*13=156` gauges,
the union rank varies through `3,4,5,6`, with signatures

```text
(3,12,44,44,0), (4,16,40,40,0),
(5,20,36,36,0), (6,24,32,32,0).
```

Thus `rank=4` and intersection dimension two are statements in the
candidate's marked common gauge, not invariants under independently moving
the two source origins.  The no-go is stronger: every relative gauge is
still annihilator-only.

## 4. Package B: typed root difference

The common-base marginal

```text
B(ell,s)=sum_(u,theta) N(u,u-s,ell,theta)
```

independently has exactly `72/91` nonzero physical cells, zero `s=0`
column, all `91/91` Fourier coordinates nonzero, and rational/split-field
rank six.  Retaining the three theta windows gives rank eight.

The corresponding projective ledgers are

| source | source rank | wedge rank | nullity | source-annihilator dimension | excess |
|---|---:|---:|---:|---:|---:|
| root difference, seven rows | 6 | 24 | 4 | 4 | 0 |
| theta-resolved, twenty-one rows | 8 | 32 | 52 | 52 | 0 |

For owner frequency `k=0` and either source convention `zeta^(-u)` or the
literal target-aligned convention `zeta^(+u)`, the ranks are

```text
rank B_0=6, rank B_1=6, rank(B_1-B_0)=6,
rank(B_0 union B_1)=9.
```

The first three systems have `(rank,equation rank,nullity,annihilator,excess)
=(6,24,4,4,0)`; the union has `(9,36,20,20,0)`.  The audit then exhausts all
twelve nonzero owner harmonics and all `12*12=144` coupled choices of owner
harmonic and drift dilation.  Every owner union has rank nine and signature
`(9,36,20,20,0)`.  Hence neither a Fourier-sign choice nor a generator
change repairs the connection.

All eight folded selections—mode zero plus one member from each of
`{1,6},{2,5},{3,4}`—have rank four and projective rank sixteen, hence zero
nullity.  Independently filtering all `7^4=2401` labelled allocations gives

```text
common exact kernels:             0
common amplitude/translation:    0
kernel classes in every case:    4.
```

This confirms both the folded and unrestricted allocation completeness.

## 5. Projective formulation audit

For source columns `X_k` and target columns `Y_k`, the candidate condition is

```text
M X_k = lambda_k Y_k.
```

The verifier constructs this relation in three independent forms:

1. all six pairwise wedge equations at every frequency;
2. three equations using the nonzero first target coordinate as anchor; and
3. the full `52`-equation augmented system in the entries of `(M,lambda_k)`.

The first two ranks agree in every tested bank.  For each load-bearing
augmented system, the projection of the nullspace to all thirteen
`lambda_k` coordinates has rank zero, and every nullspace basis vector
annihilates the complete source span.  Thus “excess zero” is not an artifact
of eliminating projective scalars.

This condition is necessary and sufficient for the declared common
circulant ansatz after excluding the annihilator solutions.  It does not
test an arbitrary `13 x 13` right operator, a frequency-dependent channel
map, or a nonlinear relation before Fourier marginalization.

## 6. Sign, torsor, and denominator ledger

- The target's `+a` phase is forced by direct `tau` reconstruction; the
  opposite sign disagrees in all active buckets.
- Reversing the drift DFT is multiplier `12` in the exhausted dilation
  family.  Reversing the septimal DFT only permutes source rows and leaves
  the arbitrary fixed channel-map test unchanged.
- Common drift translations multiply every source column by one nonzero
  frequency phase and are absorbed by `lambda_k`.  Package A's possible
  *relative* source-origin phase is separately exhausted in all 156 gauges.
- Under a fixed regular-torsor action, equivariant gauges are translations.
  Dilations are the broader generator-change hostile, not an additional
  fixed-action equivariance theorem.
- The THM-2594 entries have one common denominator `DENC`.  It is a unit at
  the certified prime and the three independent split primes.  The verifier
  explicitly compares normalized and unnormalized projective systems; both
  package signatures are unchanged.

No omitted sign, affine gauge, denominator, or projective-equation
formulation changes either scoped no-go.

## 7. Scope recommendation

Accept both packages as **FINITE-EXACT representation obstructions** with
the following precise scope:

```text
excluded:
  one frequency-independent channel map followed by one common C13
  circulant/projective multiplier, for every enumerated marked bank,
  sidecar, sign, and affine/generator gauge;

not excluded:
  arbitrary right operators, nonconvolutional address maps,
  frequency-dependent geometrically justified connections,
  nonlinear Boolean support relations, and common-stalk couplings formed
  before the source or endpoint marginal.
```

The next lawful object remains the pre-marginal support relation on
`(u,s,ell,theta; a,d,C,D)`, not another fit between aggregate spectra.

## 8. Reproduction

Run from the repository root:

```text
python -B 04-computation/lrc_r5_followup_transplants_independent_audit_20260816.py
python -B -O 04-computation/lrc_r5_followup_transplants_independent_audit_20260816.py
```

The pinned semantic digest is

```text
69e4dd6d71290663a560a5809c9fc509a0649891fd28b403a69d92b822046fc5.
```

The LF-normalized file hashes are

```text
script  13743ff9964210d26c2c3f6ae3ec73767116f46c7ed10b272d3715a0f6863856
output  11c54fbbb29d5a5c4610400de7ac0ce09e816933bef63964e0604cb87715c6d4
```

Normal and optimized executions of the independent verifier reproduce the
stored output byte-for-byte.  Normal and optimized executions of each
untouched candidate also reproduce its own stored output byte-for-byte.
