---
id: HYP-1995
status: OPEN
source: codex-2026-06-01-S522
related:
  - HYP-1990
  - HYP-1986
  - THM-387
  - THM-384
---

# HYP-1995: Exact gap races are a compactified wall ledger

## Statement

HYP-1990 should be sharpened from "some LS race reaches LL" to:

```text
Every primitive speed set has at least one wrap-around reset whose exact
post-reset wall trace reaches the compactified LL wall before a short-left
loss or the next reset.
```

The first hit of a successful race is naturally a closed wall hit.  Open LL
intervals may follow, but the proof object is the threshold wall where the
right gap first reaches `1/n` while the left gap is still long.

## Evidence

`lrc_gap_race_swap_s522.py` traces the exact wall sequence after every
wrap-around reset, using the THM-387 directed fiber flow rather than the
S518 no-swap approximation.

Bounded primitive windows match the full closed-LL audit with zero
race/full mismatches:

```text
n=3, max_speed=30: 277/277 exact reset races reach closed LL
n=4, max_speed=20: 997/997
n=5, max_speed=15: 1325/1325
n=6, max_speed=12: 786/786
```

In all four windows the first successful hit is `wall_LL`; there are no
`open_LL` first hits, because an open LL cell is entered through its closed
threshold wall.  Initial segments through `n=20` are exactly the wall-ledger
extremals: open LL measure is zero, but compactified LL wall hits persist.

The S522 Tournament Analysis over selected rows is transitive, with
fingerprints:

```text
score_hist=((0,1),(1,1),(2,1),(3,1),(4,1),(5,1),(6,1),(7,1))
c3=0
SCCs=(1,1,1,1,1,1,1,1)
Hamiltonian_paths=1
```

The top rows are gate/open-measure rows (`n14_gate_7`, `n18_gate_8`), while
initial segments remain wall-only lower in the order.

## Predictions

1. A proof of HYP-1990 should count or force compactified reset-wall hits,
   not open LL intervals.
2. The exact wall ledger should factor through arithmetic reset classes:
   resetter speed, nearest-left runway, and the next threshold wall.
3. Initial segments should be the universal equality model for the ledger.
4. THM-369 sieve completeness and CRT/gate decompositions should be used to
   guarantee at least one reset wall with enough left runway.

## Sources

- `04-computation/lrc_gap_race_swap_s522.py`
- `05-knowledge/results/lrc_gap_race_swap_s522.out`
- HYP-1990
- THM-387
- THM-384
