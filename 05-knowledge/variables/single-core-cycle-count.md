# Variable: single-core odd-cycle count

**Symbol:** `r_core(s)`
**Type:** integer statistic of a binary signature
**Defined in:** opus-2026-05-29-S11, HYP-1758

## Definition

Let a tournament T have a core vertex v such that T-v is transitive. Write the
source-to-sink order of T-v as `1,...,m`, and define a binary signature `s` by:

- `s_i = 1` if `v -> i`
- `s_i = 0` if `i -> v`

Then the directed odd cycles of T are exactly the cycles
`v -> i -> ... -> j -> v` where `i < j`, `s_i=1`, `s_j=0`, and an even number
of interior vertices from `{i+1,...,j-1}` is chosen.

Thus

```text
r_core(s) = sum_{i<j, s_i=1, s_j=0} 2^{max(j-i-2, 0)}.
```

The exponent convention gives weight 1 when `i,j` are consecutive.

## Values and Gaps

The S11 target search over all signatures of length `m <= 16` found:

| target r | status up to m=16 | interpretation |
|---:|---|---|
| 3 | absent | complete-core H=7 obstruction |
| 10 | absent | complete-core H=21 obstruction |
| 31 | first appears at m=7 | n=8 H=63 unlock |

At `m=7`, the signatures with `r_core=31` include `1001100` and `1100110`,
exactly the two THM-344 H=63 classes.

## Equations Involving This Variable

- If T has a single cycle-family core v and T-v is transitive, then
  `|V(Omega(T))| = r_core(s)`.
- Since all odd cycles share v, `Omega(T) = K_{r_core(s)}` and
  `H(T) = I(Omega(T),2) = 1 + 2*r_core(s)`.

## Relationships

- Related to [Omega(T)](omega-graph.md): gives the complete-conflict case
  `Omega(T)=K_r`.
- Related to [H(T)](hamiltonian-paths.md): `H=1+2r_core(s)` in the
  single-core family.
- Related to HYP-1758 and INV-191.

## Sources

- `04-computation/omega_extreme_fingerprints_s11.py`
- `05-knowledge/results/omega_extreme_fingerprints_s11.out`
- `05-knowledge/results/single_core_signature_targets_s11.out`

## Tags

#single-core #complete-omega #forbidden-H #signature
