# Integer return certificates, elliptic positivity, and high-rank families

Research packet, 2026-10-29. New results are scoped in their individual files;
no major open conjecture is claimed solved.

- [INTEGER_RETURNS](INTEGER_RETURNS.md): universal entry into the cheaper
  ternary queue, both integer clocks, and the exact five expanding return types.
- [RENEWAL](RENEWAL.md): same-source failure-word composition and three exact
  weighted recurrences; the queue-reset payment is still open.
- [check_returns.py](check_returns.py), [checks.json](checks.json): exact finite
  reduction, two native executions, hostile sources and recurrence checks.

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
