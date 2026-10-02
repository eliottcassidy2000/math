# Exhaustive search for edge-multiset resolving sets of Q_7

Write-up: [`05-knowledge/results/edge_multiset_dimension_q7_20261002.md`](../../05-knowledge/results/edge_multiset_dimension_q7_20261002.md)
(method, proofs of the normal-form and weight lemmas, validation, results, reproduction).

Result (FINITE-EXACT; a blind independent audit on 2026-10-02 found no false claim and no gap): Q_7 has no
edge-multiset resolving set of size k <= 12, so 13 <= edim_m(Q_7) <= 19.

Layout:
- `src/`: C sources. `edimsearch.c` (the search), `orbreps.c` (Aut(Q_n)-orbit representatives by orderly
  generation), `canon.c` (canonical form and stabilizer order), `repcheck.c` (brute-force check of
  representative files), `anneal.c` (simulated annealing).
- `py/`: the resumable driver, the validator (tiling, C leaves = DP leaves, Python re-verification), the
  independent DP leaf count, the pure-Python checker, Burnside counts, the C-versus-Python cross-check, the
  normal-form test, table and annealing summaries.
- `runs/*.sh`: the run scripts used. Their `cd` line was changed to this directory; `checks.sh` also reads its
  negative control from `data/`, and one comment line of `anneal.sh` changed.
- `runs/d7_e0/`: records of the Q_7 run, k = 1..12 (`chunks.jsonl`, `summary_k*.json`, `validation_k*.json`).
- `runs/d6_e0/`, `runs/d6_e1/`: the Q_6 reproduction (k = 1..15 with engine 0; k = 14, 15 with engine 1).
- `logs/`: driver log for k <= 12, validation, verification-mode, check, DP, `orbreps`, `is_canon` self-test and
  annealing logs, and `provenance_sha256.txt` (sha256 of the production binaries, sources and representative files).
- `data/`: Burnside counts and the negative control of `runs/checks.sh`. The representative files
  (`data/q5`, `data/q6`, about 44 MB for a <= 10) are regenerated in under a minute and not stored.

Quick start (from this directory):

```bash
gcc -O3 -march=native -Wall -o src/edimsearch src/edimsearch.c
gcc -O3 -march=native -Wall -o src/orbreps src/orbreps.c
gcc -O3 -march=native -Wall -o src/canon src/canon.c
mkdir -p data/q5 data/q6
./src/orbreps 5 32 data/q5 && ./src/orbreps 6 10 data/q6
for k in $(seq 1 10); do python3 py/driver.py 7 $k --workers 2 --target 60 --engine 0 --run runs/fast_d7; done
python3 py/validate.py 7 runs/fast_d7 $(seq 1 10)
```

The growth runner `04-computation/experiments/edge_multiset_dimension_growth_20261002_run.py --q7` runs this
fast subset (plus Q_6, k <= 13) in a temporary directory.
