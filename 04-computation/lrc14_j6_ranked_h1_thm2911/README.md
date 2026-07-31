# THM-2911 finite ranked-H1 computation package

This directory contains the exact finite rank-four/H1 computation and its
branchwise composition with the proved THM-2905 `G5` and THM-2904
hostile-centre pivot certificates.

The canonical artifacts are under

```text
05-knowledge/results/lrc14_j6_ranked_h1_thm2911/
```

They comprise 32 optimized shard outputs, the all-survivor Hunter output,
an ordinary/optimized shard-replay attestation, the promoted full THM-2904
hostile-centre ledger, and the locked three-route verifier output.
Duplicate ordinary shard files are intentionally omitted.

## Reproduction

Set:

```bash
PKG=04-computation/lrc14_j6_ranked_h1_thm2911
OUT=05-knowledge/results/lrc14_j6_ranked_h1_thm2911
PIVOT=04-computation/lrc14_j6_all_hard_ranked_h1_hunter_pivot_census_codex_20260729.py
```

Reproduce one finite-H1 shard, using four workers only for shard zero and
one worker for every other shard:

```bash
python3 -O "$PKG/scout.py" \
  --workers 4 \
  --max-cutoff 15000 \
  --max-combinations 50000000000 \
  --scope all \
  --shard-index 0 \
  --shard-count 32

python3 -O "$PKG/scout.py" \
  --workers 1 \
  --max-cutoff 15000 \
  --max-combinations 50000000000 \
  --scope all \
  --shard-index 1 \
  --shard-count 32
```

Change `--shard-index` through `31` and compare each complete output with
`$OUT/all_sNN_of32.out`. The verifier hash-pins all 32 files, their
aggregate semantic digest, and the independent ordinary-mode replay
attestation.

Reproduce the 24 nonstar Hunter repairs:

```bash
python3 "$PKG/hunter_all_24_star_survivors.py"
python3 -O "$PKG/hunter_all_24_star_survivors.py"
```

Both modes must match
`$OUT/hunter_all_24_star_survivors.out` byte for byte. The locked verifier
also joins the exact 24 `(body, rank, pivot)` identities in the shard
exceptions, the pinned source target battery, and the Hunter output.

Reproduce the THM-2904 ledger:

```bash
python3 "$PIVOT" --workers 8 \
  --ledger /tmp/thm2904.ordinary.ledger \
  > /tmp/thm2904.ordinary.out
python3 -O "$PIVOT" --workers 8 \
  --ledger /tmp/thm2904.optimized.ledger \
  > /tmp/thm2904.optimized.out
cmp /tmp/thm2904.ordinary.out /tmp/thm2904.optimized.out
cmp /tmp/thm2904.ordinary.ledger /tmp/thm2904.optimized.ledger
cmp /tmp/thm2904.ordinary.ledger \
  "$OUT/thm2904_hostile_centre.ledger.out"
```

The two summary outputs must also match the canonical THM-2904 output at
`05-knowledge/results/lrc14_j6_all_hard_ranked_h1_hunter_pivot_census_codex_20260729.out`.

Run the locked aggregate verifier:

```bash
python3 "$PKG/verify.py"
python3 -O "$PKG/verify.py"
```

Both modes must match `$OUT/locked.out` byte for byte.

## Scope and architecture

`scout.py` reconstructs literal apex carriers and actual excluded prefixes,
uses the strict discrepancy cutoff to compute complete finite H1 cores, and
checks the ordered earliest-label pivot. It closes 5,975 branches directly
and sends 24 branches to the Hunter audit.

`hunter_two_star_exceptions.py` supplies exact overlap and maximum-tree
machinery. `hunter_all_24_star_survivors.py` audits all 24 exceptions and
repairs all 54 hostile five-sets.

`verify.py` does not rerun those interval scans. It verifies every raw
artifact, reconstructs the 5,999 finite-H1 and 2,964 `G5` branch keys from
the hardened THM-2905 ledger, and recovers the 279 hostile-centre pivot
keys from the promoted THM-2904 ledger. Every key retains

```text
(body, gate size, rank, apex, excluded prefix).
```

The exact branch atoms are

```text
pair-cap exception only       52
G5 only                       55
G5 and finite H1            2909
finite H1 only              2875
finite H1 and pivot          215
pivot only                    64
three-route union           6118
```

The 52 pair-cap exceptions are disjoint from the three-route union and are
a hostile control, not silently treated as closures.

Whole-root recomposition gives:

```text
proved baseline through THM-2904     88
finite-H1 roots                     132
finite-H1/G5 roots                  134
three-route roots                   138
three-route intersection             45
three-route additions                 93
proved union                         181
residual                            3251
```

The pivot route increases the current union by exactly one root beyond the
finite-H1/G5 composition:

```text
(1,2,3,5,6,8,13).
```

Its rank-one branch is finite-H1 closed (cutoff `6930`) while its rank-two
branch is pivot-only. This is a genuine H1/pivot composition, not a
single-route closure.

The 181 rank-four/H1 branches with cutoff above 15,000 remain outside this
finite route. Failure to enter any route is not evidence for a cover, and
this package does not prove LRC(14).
