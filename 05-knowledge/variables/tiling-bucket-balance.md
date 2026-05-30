# Variable: tiling bucket balance

**Symbol:** `bucket_balance_q,M(b)`
**Type:** integer identity / quotient-matrix row constraint
**Defined in:** THM-346; Lean half-line core THM-348; unordered Lean layer THM-350; Boolean-mask specialization THM-351; transport-row checksum THM-352; concrete staircase specialization THM-353; quotient gap rows/columns THM-355

## Definition

Let `q : Q_m -> B` be a quotient map from tilings to buckets, and let `M`
be a set of nonzero Hamming masks. For bucket `b`, define:

```text
bucket_balance_q,M(b) =
  |q^{-1}(b)| * |M| - 2*self_b(M) - incident_cross_b(M).
```

Here `self_b(M)` counts unordered mask-lines with both endpoints in bucket
`b`, and `incident_cross_b(M)` counts unordered mask-lines with exactly one
endpoint in bucket `b`.

THM-346 proves:

```text
bucket_balance_q,M(b) = 0
```

for every quotient `q`, mask family `M`, and bucket `b`.

Lean now proves the oriented finite-set version (THM-348):

```text
|selfHalf_b(M)| + |crossHalf_b(M)| = |q^{-1}(b)| * |M|.
```

The full unordered tiling formula is this half-line identity plus the
mask-involution pairing `(x,u) <-> (x xor u,u)` on internal lines. THM-350
proves the abstract unordered algebra and the finite orbit-parity theorem:
fixed-point-free involutive moves force `selfHalf_b(M)` to have even
cardinality and therefore satisfy the unordered balance. THM-351 closes the
Boolean-mask specialization by proving that xor by any nonzero mask is a
fixed-point-free involution.

## Values at Small n

The S5 script checked this for the merged tournament quotient `G_n/Z_2` and
the even-graph quotient `E_n`.

| n | quotients | mask families | violations |
|---|---|---|---:|
| 3..6 | `G_n/Z_2`, `E_n` | `d=1`, middle, `d=m`, all waggly | 0 |

For complement-tiling (`|M|=1`) the parity corollary was also checked:

```text
incident_cross_b == |q^{-1}(b)| mod 2.
```

It had 0 violations for `n=3..6` in both quotients.

## Equations Involving This Variable

- `2*self_b(M) + incident_cross_b(M) = |q^{-1}(b)|*|M|` (THM-346).
- `|selfHalf_b(M)| + |crossHalf_b(M)| = |q^{-1}(b)|*|M|`
  (THM-348, Lean).
- If `Even |selfHalf_b(M)|`, then
  `2*internalLineCount_b(M) + |crossHalf_b(M)| = |q^{-1}(b)|*|M|`
  (THM-350, Lean).
- If each selected move is a fixed-point-free involution, then
  `Even |selfHalf_b(M)|` and the same unordered balance follows without a
  separate evenness assumption (THM-350, Lean).
- If the space is a finite Boolean cube and each selected mask is nonzero,
  the selected xor moves are fixed-point-free involutions and the unordered
  balance follows directly (THM-351, Lean).
- If `transportHalf_b,c(M)` counts oriented half-lines from bucket `b` to
  target bucket `c`, then

```text
sum_{c != b} |transportHalf_b,c(M)| = incident_cross_b(M).
```

  Combining this with THM-351 gives the Boolean-cube row checksum

```text
2*internalLineCount_b(M) + sum_{c != b} |transportHalf_b,c(M)|
  = |q^{-1}(b)|*|M|.
```

  This is THM-352 in Lean and is the preferred matrix-audit form for
  quotient transport computations.
- In the concrete staircase tiling cube `StTiling n`, the same checksum is
  now formalized for all nonzero masks, single-tile masks, the complement
  mask for `n >= 3`, and the finite good-cut quotient
  `goodCutBucket : StTiling n -> Fin (n+1)` (THM-353, Lean).
- If a quotient bucket has empty fiber, every transport row out of that bucket
  and every transport column into that bucket is empty (THM-355, Lean).  For
  `goodCutBucket`, the finite image is `{0} ∪ {2,...,n-1}` when `n >= 3`, so
  buckets `1` and `n` are certified zero rows/columns.
- If `|M|=1`, then `incident_cross_b(M) == |q^{-1}(b)| mod 2`.
- Normalized escape plus neutrality:

```text
escape_b(M)   = incident_cross_b(M) / (|q^{-1}(b)|*|M|)
neutral_b(M)  = 2*self_b(M) / (|q^{-1}(b)|*|M|)
escape_b + neutral_b = 1.
```

## Relationships

- Related to [projection defect](projection-defect.md): projection-defect
  profiles classify where cross-lines go after two quotient maps are compared.
- Related to `G_n/Z_2`: its bucket cross-lines decompose into SC-SC spine,
  SC-NS ribs, and NS-NS sea.
- Related to `E_n`: the same balance holds for the even-graph quotient, giving
  a dual audit whenever merged-tournament buckets are computed.

## Sources

- `01-canon/theorems/THM-346-tiling-quotient-bucket-balance.md`
- `01-canon/theorems/THM-348-finite-bucket-halfline-balance.md`
- `01-canon/theorems/THM-350-finite-unordered-bucket-balance-layer.md`
- `01-canon/theorems/THM-351-boolean-cube-mask-bucket-balance.md`
- `01-canon/theorems/THM-352-quotient-transport-row-checksum.md`
- `01-canon/theorems/THM-353-staircase-tiling-transport-checksum.md`
- `01-canon/theorems/THM-355-quotient-gap-transport-vanishing.md`
- `01-canon/theorems/THM-345-merged-bucket-parity.md`
- `04-computation/lean/TournamentH7/TournamentH7/BucketBalance.lean`
- `04-computation/lean/TournamentH7/TournamentH7/StaircaseBucketTransport.lean`
- `04-computation/tiling_quotient_bucket_balance_s5.py`
- `05-knowledge/results/tiling_quotient_bucket_balance_s5.out`
- `05-knowledge/results/lean_bucket_balance_kind_pasteur_2026-05-30-S1.out`
- `05-knowledge/results/lean_verify_unordered_kind_pasteur_2026-05-30-S2.out`
- `05-knowledge/results/lean_boolcube_bucket_balance_opus_2026-05-30-S1.out`
- `05-knowledge/results/lean_verify_boolcube_bucket_balance_opus_2026-05-30-S1.out`
- `05-knowledge/results/lean_staircase_bucket_transport_codex_2026-05-30.out`
- `05-knowledge/results/lean_bucket_gap_transport_codex_2026-05-30.out`
- `07-reflections/merged-tiling-bucket-constraints.md`
- `07-reflections/unordered-bucket-balance-orbits.md`
- `07-reflections/boolean-cube-balance-as-checksum.md`
- `07-reflections/staircase-transport-is-boolean-transport.md`
- `07-reflections/fiber-gaps-and-residue-boundaries.md`

## Tags

#tiling #bucket #merged-metagraph #even-graphs #projection-defect #engineering-checksum
