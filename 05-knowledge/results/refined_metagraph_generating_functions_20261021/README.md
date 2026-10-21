# Refined metagraph generating functions

Status: PROVED from the explicitly inherited finite-atom and SCC laws; independent review PASS. FINITE-EXACT controls are not a tournament census.

See the frozen [proof](session23-refined/PROOF.md), [source guide](session23-refined/README.md), and [independent review](session23-refined-review/REVIEW.md). The extensions are multivariate marked rationality, exact support-dependent denominator laws, the all-r power-of-three denominator family, and a raw automorphism-order non-D-finiteness obstruction. Exact recurrence orders 9 and 18 are inherited controls, not new results. Unknown strong-atom and center coefficients, recurrence onset, and general odd-degree cancellation are not settled by these results. No global literature-priority or connectivity claim is made.

All source and review bytes are unchanged, in their original sibling layout. Source pins retain original repository provenance. The top-level manifest records exact lengths and hashes of every payload file, including the original source manifest.

## Replay

From this packet directory, with Python 3 and SymPy:

```sh
python3 session23-refined/check.py > /tmp/refined-producer.json
python3 -O session23-refined/check.py > /tmp/refined-producer-O.json
python3 session23-refined-review/check.py > /tmp/refined-review.json
python3 -O session23-refined-review/check.py > /tmp/refined-review-O.json
cmp /tmp/refined-producer.json session23-refined/normal.json
cmp /tmp/refined-producer-O.json session23-refined/optimized.json
cmp /tmp/refined-review.json session23-refined-review/normal.json
cmp /tmp/refined-review-O.json session23-refined-review/optimized.json
```

Independent controls: 1,280 checks and 11,435 marked formal words. Synthetic atom signatures remain explicitly synthetic, except for the separately identified certified small tournament signatures. No canonical or navigation files are changed in this checkpoint.
