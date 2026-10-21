# Two-shore exclusions and prime-center refinement

## Scoped reconciliation at the 5872 incoming boundary

[Independent reconciliation and attribution](session22-incoming-center-review/RECONCILIATION.md), [portable reproduction](session22-incoming-center-review/REPRODUCTION.md), and the separate [Aut-orbit-color corollary](session22-incoming-center-review/ORBIT_COLOR_COROLLARY.md) preserve the proof and exact consumer controls.

The incoming [combined cutoff](../deficit_six_centers_20261021b/COMBINED_CUTOFF.md), published at f4353d7da3e141ddad7d90f13c40ac1ded239552, combines the prior general closing-cycle obstruction with the independently reviewed odd-shore argument and its order10 certificate. It proves degree R6 <=35 and leaves ten possible residual orders. The two-shore insertion theorem below independently excludes order14, leaving exactly the following nine potentially unknown coefficient positions:

    {12,15,18,20,21,24,25,28,35}.

This is a necessary support list, not a classification or evaluation of those coefficients. The safe order18 recurrence onset at this boundary is n>=48, not an earliest-onset claim. Incoming work receives credit for the35 cutoff; this packet contributes the exclusion of14. Frozen45/20 and24-position materials below remain honest intermediate derivations, not the latest global frontier.

Standalone intermediate consequence: degree R6 <=45, with20 potentially unknown coefficient positions and safe order18 recurrence onset n>=58. The degree47 and24-position results are valid intermediate bounds, retained with their own proofs and review lineage rather than silently rewritten. No complete center coefficient census or earliest recurrence onset is claimed.

- [Analytic two-shore insertion exclusion](session22-q2structural/PROOF.md): q=2 with odd shore size m>=5 gives deficit at least2m-2, excluding deficit-six centers at orders10 and14.
- [Exact marked order10 universe](session22-order10/README.md) and [independent review](session22-q2-review/ORDER10_REVIEW.md): extra finite evidence, not a census of all order10 tournaments.
- [Intermediate24-position support reduction](session22-q2structural/SUPPORT_REDUCTION.md)
- [Final prime-circulant refinement](session22-primecenters/PRIME_CUTOFF.md) and [separate independent review](session22-prime-review/REVIEW.md): excludes prime orders37,41,43,47; composite39 and45 remain possible.
- [Inherited degree47 packet](../odd_shore_deficit_centers_20261021/README.md)

The prime connectivity argument uses the cited Cauchy–Davenport theorem only for primes; a literal order9 composite hostile is retained. The complete low-order coefficients through9 are explicitly inherited, not freshly re-enumerated here. Parameter support means possible indices, not existence or values of their coefficients.

## Later reviewed incoming refinement

The separately committed incoming work at 359ddd20b63b00e24bc64204e1e7e78ed5abd5e4 excludes orders12 and25 and independently corroborates14: [order12](../deficit_twelve_20261021c/README.md), [order25](../center_transitivity_20261021c/README.md), [order14](../deficit_fourteen_20261021c/README.md). The final reviewed cutoff therefore leaves seven possible unknown positions {15,18,20,21,24,28,35}, with degree35 and safe onset48. The [final synthesis](../session22_tournament_tail_synthesis_20261021/README.md) routes that later frontier and its independent reviews. The nine-position intersection above is this packet's scoped intermediate contribution, not a claim of latest sharpness or priority over the separate incoming14 census. Frozen proofs are not rewritten.

## Reproduce

From this packet directory:

```sh
(cd session22-q2structural && python3 -B check.py && python3 -B -O check.py && python3 -B support_check.py)
(cd session22-order10 && python3 -B check.py && python3 -B -O check.py)
python3 session22-q2-review/independent.py
python3 -O session22-q2-review/independent.py
(cd session22-primecenters && python3 -B check.py && python3 -B -O check.py)
python3 session22-prime-review/independent.py
python3 -O session22-prime-review/independent.py
python3 session22-incoming-center-review/check_intersection.py --repo-root ../../..
python3 -O session22-incoming-center-review/check_intersection.py --repo-root ../../..
python3 session22-incoming-center-review/orbit_color_control.py
python3 -O session22-incoming-center-review/orbit_color_control.py
```

Producer checkers overwrite receipts with actual clock fields. Replay in a temporary copy to preserve published evidence; compare scientific fields excluding start_utc/end_utc/wall_seconds. The prime independent driver has a transparent path-only adaptation to use final producer receipts instead of copied review snapshots; portable-paths.patch and portable-source-pins.json record exact changes and both hashes. Original review source pins remain provenance, and adapted source is not claimed byte-identical. Other selected files are unchanged; copied audit execution directories and exploratory files are excluded.

Publication baseline: 336b76f137a6bc4843a9b1cb702132ad41574d58. The mathematical incoming cutoff is explicitly359ddd20; later incoming work is not claimed reviewed here. The manifest pins all included files except itself. Separate proof and census scope are maintained throughout.
