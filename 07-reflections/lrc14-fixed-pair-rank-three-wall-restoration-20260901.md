# Fixed-pair rank-three wall restoration

Date: 2026-09-01

Status: **FINITE-EXACT, TWO IMPLEMENTATIONS / EXPLORATORY ONLY.**  This note
does not promote reserved THM-4333, prove arbitrary third-outsider entry, or
prove LRC(14).

## Scoped object

Use THM-4326's fixed thirty-label pool `P`, choose a nine-label body, fix the
two outsiders `(50,70)`, and let a third outsider `s` range over the 737 typed
values in `1<=s<=769`.  The exact wall decomposition assigns to each
outsider-safe cell the subset of pool labels that fail there.

The rank-three truncation keeps only failure masks of size at most three.
For a chosen body its retained mass is therefore a lower bound for the full
safe mass, but it deliberately discards every rank-at-least-four cell.

## Exact result

Two independent implementations agree:

- `729` third outsiders have minimum rank-three mass strictly above `4/63`;
- none is equal;
- exactly `s=45,55,65,90,110,130,220,260` fall below;
- those eight rows have respectively
  `1,12,422,1729,32756,5549,36,6` subthreshold bodies, `40,511` total;
- restoring all wall cells makes every one of those bodies strictly exceed
  `4/63`;
- the smallest full mass is at `s=110`, body mask `0x00384583`, and equals
  `1021183/13113800>4/63`.

The primary path uses a custom exact weighted maximum-coverage
branch-and-bound.  The clean-room path reconstructs the rational walls by a
different floor classifier, proves every minimum `OPTIMAL` with OR-Tools
CP-SAT, exhausts every subthreshold solution single-threadedly, and evaluates
its full mass directly.  Their key fields agree exactly.

## Connection and loss ledger

```text
source:       exact outsider-safe wall cells for P union {50,70,s}
target:       weighted failure hypergraph, first truncated at rank three
map:          cell -> failed pool-label mask and exact rational width
preserved:    exact body avoidance for every retained mask
destroyed:    every rank>=4 cell, cyclic address, endpoints, and owner phase
sidecar:      the complete all-rank mask table for each s
decisive test: restore full mass for every rank-three-subthreshold body
```

The useful signal is negative and positive at once.  Rank three alone is not
a uniform certificate: eight outsiders and 40,511 bodies cross below the
target.  But every crossing is repaired by the exact information that the
truncation threw away.  Any cofinal or arbitrary-pair theorem must control
that restoration uniformly; the bounded scan itself supplies no such
transfer.

The now-proved THM-4333 overlaps but has a different quantifier structure: it
proves a rank-three pair-safe surplus above `2/27` throughout the exact
THM-4231 residual-pair universe and derives a cofinal third-tail cutoff.  This
probe instead fixes `(50,70)`, inserts the third outsider before measuring,
and restores every discarded rank on the bounded range `s<=769`.  It is a
finite diagnostic beside THM-4333, not an independent audit of that theorem
or an arbitrary-third-tail extension of it.

## Reproduction

Primary:

```text
python3 04-computation/lrc14_rank3_wall_hypergraph_probe_20260901.py \
  --scan-through 769
```

Frozen output: `05-knowledge/results/lrc14_rank3_wall_hypergraph_probe_20260901.out`.

SHA-256: `b12c41e5d4ca2a3130dd0127dec85a0fc6dbd0365aaa245624fc897bd72e7d70`
and `f22d022ba93c906b5c519ff981db7db4998873995839182ad85de40bc06673d4`.

Independent:

```text
python3 04-computation/lrc14_rank3_wall_hypergraph_probe_20260901_independent.py
```

Frozen output: `05-knowledge/results/lrc14_rank3_wall_hypergraph_probe_20260901_independent.out`.

SHA-256: `cc7eceb3d6d5f40f02251a93735275edd70c1dde17da3653fe91dd31a1b5dc03`
and `16fa3be0a10add8dab8625380ae3c6a2445e3aef5165df28eb9b25e42a5ea476`.
