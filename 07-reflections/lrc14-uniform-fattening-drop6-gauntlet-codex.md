# LRC14 uniform fattening: the drop-6 core is the current worst witness

**Source:** codex-2026-06-17-S1.  Prompt: work creatively on an LRC(14) proof.  This note targets
OPEN-Q-108, the uniform fattening lemma needed by the current singular-series proof skeleton.

## Result

The exact rational gauntlet `04-computation/lrc14_uniform_fattening_gauntlet_codex.py` tests
12-speed cores `C` by computing

`G_C={t in [0,1): ||v t|| > 1/14 for all v in C}`.

The smallest core in all tested families is the AP drop-6 core

`C6={1,2,3,4,5,7,8,9,10,11,12,13}`,

with

`meas(G_C6)=7/858 = 0.008158508...`.

This matches the decoupling floor that kind-pasteur's THM-523 proof skeleton already singled out,
but here it is tested directly as a 12-core extremal candidate.  The run found:

- AP drop-one cores: drop `6` is the unique minimum; next is drop `12` with `426/35035`.
- Two AP deletions plus one replacement, `w<=180`: `13026` exact cores, zero with `meas(G_C)=0`, and best positive `3859/420420` at delete `(6,10)`, add `20`.
- Random primitive 12-cores with values `<=90`: `3000` exact cores, zero with `meas(G_C)=0`; best was much larger, about `0.0926`.
- Single-swap greedy stress moves the sporadic core `{1..11,13}` directly to the drop-6 core by swapping `6 -> 12`, then stops.

## New subtarget

OPEN-Q-108 can be sharpened into a concrete extremal statement:

> Prove either that every primitive 12-core has `meas(G_C) >= 7/858`, with equality at the AP
> drop-6 core, or prove at least that every non-AP coordinated-growth core has
> `meas(G_C) >= 7/858`.

This is stronger than the uniform fattening lemma, but it is much more rigid: a single exact core
and a single rational value replace an existential lower bound.

## Tournament Analysis

The script uses a speed-load tournament.  For a pair of speeds `(a,b)`, remove both and compare how
much of `G_{C\{a,b}}` is covered by `D_a` versus `D_b`; orient toward the more load-bearing speed,
with ties resolved by the speed-order Hamiltonian path.

For both the AP drop-6 core and the best random core, the speed-load tournament is transitive:
score histogram `{0:1,...,11:1}`, `0` directed 3-cycles, SCC sizes all `1`, and exactly one
Hamiltonian path.  Scale-by-2 has `0` edge flips, as it should.  Switching the threshold from
`1/14` to `1/13` flips `7` edges for the drop-6 core and `4` for the random core.

This is useful but also limiting: runner-level load order is not where the obstruction hides.
The next quotient should use safe components, endpoint events, q-grid obligations, deleted AP
positions, residues mod 14, pair-sum wall crossings, or proof-obligation packets.  The speed
tournament preserves pairwise load contribution to `meas(G_C)`, but it destroys endpoint alignment
and component adjacency, exactly the data OPEN-Q-108 probably needs.

## Honest status

This is evidence, not a proof.  It does not settle coordinated growth with three or more large
arithmetically related speeds.  Its value is that the surviving proof target is now smaller and more
testable: explain why the visible worst core is `AP_13 \ {6}` and why coordinated growth cannot make
the `1/14` safe set thinner than `7/858`.

Cross-links: OPEN-Q-108, THM-523, THM-522, HYP-2561, HYP-2566, HYP-2568.
