# Variable: tiling bucket balance

**Symbol:** `bucket_balance_q,M(b)`
**Type:** integer identity / quotient-matrix row constraint
**Defined in:** THM-346; Lean half-line core THM-348; unordered Lean layer THM-350

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
proves the abstract unordered algebra in Lean once `selfHalf_b(M)` has even
cardinality, and proves the partner-map closure/non-self-pairing facts needed
to derive that evenness from a fixed-point-free involution.

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
- `01-canon/theorems/THM-345-merged-bucket-parity.md`
- `04-computation/lean/TournamentH7/TournamentH7/BucketBalance.lean`
- `04-computation/tiling_quotient_bucket_balance_s5.py`
- `05-knowledge/results/tiling_quotient_bucket_balance_s5.out`
- `05-knowledge/results/lean_bucket_balance_kind_pasteur_2026-05-30-S1.out`
- `05-knowledge/results/lean_verify_unordered_kind_pasteur_2026-05-30-S2.out`
- `07-reflections/merged-tiling-bucket-constraints.md`
- `07-reflections/unordered-bucket-balance-orbits.md`

## Tags

#tiling #bucket #merged-metagraph #even-graphs #projection-defect #engineering-checksum
