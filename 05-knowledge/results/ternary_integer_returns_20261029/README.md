# Integer return certificates, elliptic positivity, and high-rank families

Research packet, 2026-10-29. New results are scoped in their individual files;
no major open conjecture is claimed solved.

- [INTEGER_RETURNS](INTEGER_RETURNS.md): universal entry into the cheaper
  ternary queue, both integer clocks, and the exact five expanding return types.
- [RENEWAL](RENEWAL.md): same-source failure-word composition and three exact
  weighted recurrences; the queue-reset payment is still open.
- [check_returns.py](check_returns.py), [checks.json](checks.json): exact finite
  reduction, two native executions, hostile sources and recurrence checks.
- [LEAST_SOURCES](LEAST_SOURCES.md): a pointwise source-lift recurrence;
  exact least all-growing sources through forty blocks, and the first failure
  of the Mersenne-only candidate at fifteen blocks.
- [FRUIT_RETURNS](FRUIT_RETURNS.md): a sharp thirteen-shift real-circle cover,
  first-return clocks `4,9,13`, and a rational-lifting/positivity steering rule
  with sharp coefficient budget `76=4+6*12` on the named cyclic subgroup.
- [RANKS_AND_CONNECTIONS](RANKS_AND_CONNECTIONS.md): rank21's heuristic scope,
  optimal fifteen/sixteen/sixteen-prime witness banks, and a genuine Fano/F4
  projective map explaining the controller probabilities `1/7` and `1/21`.
- [AUDIT](AUDIT.md): inheritance, independent algorithm paths, scope and open
  obligations. No independent-agent audit of this new packet is claimed.

The concurrent Session97 [seven-prefix proof](../session97_overnight_arithmetic_carriers_20261029/controllers/FIVE_EXCEPTIONAL_CYLINDERS.md)
is preferred over this packet's larger finite proof of the same five classes.
Its exact rank/isogeny certificates are also retained as dependencies. The
new extensions are linear division clocks, weighted failure renewals, least
source recurrences, effective elliptic returns, and the observer/projective
bridges; overlapping results are not counted twice.

## Inheritance and live concepts

Anchor: replace the preceding auxiliary Haar almost-sure statement by an
integer-specific certificate, then retain the actual carry across returns.
Niche: the September fruit cubic's positivity obstruction and August local
rank certificates. Wildcard: interpret high-rank families and tournament gaps
through their actual counting mechanisms rather than identifying their labels.

| Live concept | Closest mechanism / hostile / missing coordinate |
|---|---|
| Cheaper ternary grammar | [Prior moments](../affine_tree_recurrences_20261028/MOMENTS.md); `3 -> 5`; queue-reset debt |
| Source refinement | [Prior tail bound](../affine_tree_recurrences_20261028/REFINEMENT.md); nested null source; exact `(J,D,C)` |
| Positive fruit solutions | [Fruit subgroup tree](../prime_shells_20260921_elliptic.md); all children of `9G` fail positivity; return clock |
| High-rank elliptic curves | [August rank-30 certificate](../elliptic_rank30_optimal_local_observer_certificate_codex_20260821.out); point incidence alone does not prove independence; local character matrix |
| Family construction | [August Mestre sidecar](../elliptic_mestre_six_tuple_sidecar_rule30_scout_codex_20260821.out); zero-sum sextuples need another equation; independent section data |
| Tournament fibres | THM-4812, [all-order exclusion](../../../01-canon/theorems/THM-4812-presentation-fibres-7-and-21-are-absent-at-every-order-orbit-product-lemma-vertex-transitive-lemma-closed-window.md); bounds alone leave a finite window; full orbit/port compatibility |

The corrected near misses are the distinction between stopping and descending,
the old extra decimal zero in two fruit coordinates, and interpreting a rank
lower bound or a heuristic threshold as an unconditional exact rank. The
least-used inherited sidecars here are the Mestre coefficient relation and
the torsion/free-coefficient separation on the fruit curve.

Reproduce the core certificate from the repository root:

```sh
python3 -B 05-knowledge/results/ternary_integer_returns_20261029/check_returns.py
python3 -B -O 05-knowledge/results/ternary_integer_returns_20261029/check_returns.py
```

The output is deterministic. Its embedded source digest pins the checker;
the proof explains exactly which finite check closes the universal reduction.

The additional standard-library scripts are `check_fruit_returns.py`,
`check_least_sources.py`, `check_rank_banks.py`, and `check_projective_bridge.py`.
Run each normally and with `-O`; their deterministic JSON receipts are stored
alongside them. The rank checker reads the three named, tracked curve snapshots
in the incoming Session97 `ranks` directory. No live network is needed to replay
the proofs. The raw 12.9 MB leaderboard download is not committed.
