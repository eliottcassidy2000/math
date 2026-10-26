# Native alignment-moment checks

Final mathematical checker: `check_alignment_moments.py`.

- Normal: PASS, 24,741 explicit checks.
- Optimized `-O`: PASS, 24,741 explicit checks.
- Outputs `check.normal.json` and `check.optimized.json` are byte-identical.
- Complete native universe: k=1,...,7, every one of 2^k alignments and every one of 2^k simultaneous-subset guards at each k; 21,844 actual guard evaluations.
- Seed rigidity: all 7!=5,040 permutations for X and separately for tau(X).
- Additional exact controls: source/target rooted identifications, strongness, all-black colour, degree-separated wrapper rigidity, Boolean triangular rank, independent GF(2) elimination, every truncation rank, every proper full marginal, exact native moment evaluation versus superset transform, and Boolean inversion.
- No inherited producer module is imported. All check conditions use explicit exceptions and survive Python optimization.
- The all-order theorem is proved in REPORT.md; finite replays are controls, not an all-order enumeration.

Commands, from this directory:

```
python3 -B check_alignment_moments.py > check.normal.json
python3 -B -O check_alignment_moments.py > check.optimized.json
cmp check.normal.json check.optimized.json
```

Independent peer audit is in the separate `session70-bridge-audit` lane. That audit reconstructs seeds from free masks, checks maps/rigidity, native guards through k=6, the integral inverse and ranks modulo 2,3,5. Final independent prose review accepted REPORT.md hash da137bfa8ad50378a29f963b87e7b7950b8878c4cfea79590c06783c0d4d4792 and code hash e3836182474bf9a336c148f73ab1fb24640bef79cd31b077545592da3e734d88, with no mathematical scope edits. Its independent normal/-O receipts are byte-identical PASS with 78,653 checks. Parent should include those frozen report and receipts for separately audited publication.

Source snapshots are listed by git blob and SHA-256 in SOURCES.json. Native construction remains pinned at fb872272125510c317b8114dd3f697b433c595f9; incoming fresh-re-marking records are pinned separately at 7d5c27c1f53ef90f2da1efad4793b90f53a75cd5.
