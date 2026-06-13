# HYP-1785: Deletion-residue rank filters first Omega obstructions

**Status:** OPEN engineering/structural hypothesis; S355 feature layer added.
**Session:** codex-2026-05-30 S355

## Statement

For odd-cycle conflict graphs, the first real-rootedness and ghost-cycle
obstructions are better filtered by deletion-residue rank than by raw odd-cycle
count alone.

Given the old-projection deletion residue

```text
R_v(T) = {directed odd cycles C of T : v notin C},
```

record the max-loss vertex `v*` and the rank

```text
rank_res(v*) = max {k : alpha_k(R_v*) > 0}.
```

Exact kills have `rank_res=0`; dangerous near-kills have small nonzero
residue size but `rank_res>=2`.

## Evidence

The S355 probe compares the same examples that motivated HYP-1779/HYP-1780:

| Example | Max-loss residue |
|---|---|
| H=63 single-core signatures `1001100`, `1100110` | `alpha=(0,0)`, `rank_res=0`, exact kill |
| THM-025 n=9 real-root failure | `alpha=(2,1)`, `rank_res=2`, `I(R,2)=9` |
| Paley T7 | broad residue `alpha=(20,1)`, not near-kill |
| Interval T7 | broad residue `alpha=(16,2)`, not near-kill |

This separates three cases that scalar `H`, total odd-cycle count, and even
max deletion loss blur together:

1. empty exact kill (`rank_res=0`);
2. small dangerous near-kill (`keep=2`, `rank_res=2`);
3. broad fiber/disjointness residue (`keep>=16`).

## Predictions

1. Non-real-rooted `I(Omega(T),x)` examples at `n>=9` should be enriched in
   the stratum `0 < keep_v* <= 2` and `rank_res(v*) >= 2`.
2. H=63-style complete-Omega unlocks should sit in the exact-kill stratum
   `rank_res(v*)=0`, even when total odd-cycle count is large.
3. HYP-408-style ghost-cycle failures should prefer low-size nonzero residues
   with `rank_res>=2` over exact kills or broad residues.
4. Paley/Interval contrasts should be modeled as fiber/disjointness residues,
   not deletion-rank near-kills.

## Implemented Features

`04-computation/tournament_tda.py` now emits:

- `scc_count`, `scc_defect`, `scc_largest`, `scc_is_strong`;
- `omega_min_deletion_residue_cycles`;
- `omega_min_positive_deletion_residue_cycles`;
- `omega_near_kill_vertices`;
- `omega_near_kill_rank2_vertices`;
- `omega_max_loss_residue_alpha_1`;
- `omega_max_loss_residue_alpha_2`;
- `omega_max_loss_residue_rank`;
- `omega_max_loss_residue_I2`.

## Related

- HYP-1779: support-residue calculus.
- HYP-1780: residue-rank stratification.
- THM-025, THM-344, THM-354.
- `05-knowledge/variables/residue-rank.md`.
- `04-computation/residue_rank_probe_s355.py`.
- `05-knowledge/results/residue_rank_probe_s355.out`.
