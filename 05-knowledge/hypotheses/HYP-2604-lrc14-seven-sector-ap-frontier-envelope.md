# HYP-2604: LRC(14) Seven-Sector AP-Frontier Envelope

**Status:** OPEN proof program with a useful refutation of the too-strong
AP-extremal conjecture.

Script: `04-computation/lrc14_sector_ap_extremal_codex.py`  
Output: `05-knowledge/results/lrc14_sector_ap_extremal_codex.out`

## Claim

For the seven-sector cover set

`S7(E) = {x : {e*x mod 1 : e in E} hits all seven fixed sectors}`,

the strong conjecture

`meas(S7(E)) <= meas(S7({0,1,...,k-1}))`

is the right shape only for the dangerous sector-cap rows `8 <= k <= 11`.
For `k=12,13` it is false, but the failures are near-AP tail rows with large
slack against `cap_k`.

The corrected proof target is therefore an **AP-frontier envelope**:

1. prove AP extremality for `8 <= k <= 11`;
2. prove a coarse AP-rich envelope for `k=12`;
3. ignore `k=13` at this layer, since `cap_13=1` makes the sector cap
   tautological;
4. carry the continuous sector result into the colored finite-denominator
   placement layer of HYP-2593/HYP-2595.

## Exact Scout

The script uses an exact integer common-refinement engine for sector endpoints
`m/(7e)`, avoiding slow rational sorting while keeping exact measures.

Widened primitive boxes found no AP-beater in the dangerous range:

| `k` | exact primitive box | max `meas(S7)` | maximizer | AP beaters |
|---:|---:|---:|---|---:|
| 8 | `maxE<=18`, `31788` rows | `481/1470` | `{0,...,7}` | `0` |
| 9 | `maxE<=17`, `24309` rows | `2447/5880` | `{0,...,8}` | `0` |
| 10 | `maxE<=16`, `11440` rows | `8899/17640` | `{0,...,9}` | `0` |
| 11 | `maxE<=16`, `8008` rows | `3419/5880` | `{0,...,10}` | `0` |

The strong conjecture fails in the easy-slack range:

| `k` | box max | maximizer | cap margin |
|---:|---:|---|---:|
| 12 | `11381/17640` | `(0,1,2,3,4,5,6,7,8,9,10,12)` | `3739/17640` |
| 13 | `12029/17640` | `(0,1,2,3,4,5,6,7,8,9,10,12,15)` | `5611/17640` |

So the false strong conjecture fails in the harmless direction: it produces a
better envelope target rather than a sector-cap counterexample.

## Proof Route

The live continuous S3 proof route now has three nested formulations.

- HYP-2603: prove `meas(S7(E)) <= cap_k`.
- THM-532: write `meas(S7(E)) = M7(k) + corr(E)` and bound the relation-height
  correction.
- HYP-2604: prove the AP-frontier envelope, using the exact AP cap checks for
  `k=8..11` and coarse slack for `k=12,13`.

This route is attractive because it turns the low-height residual from an
open-ended relation-tail enumeration into a discrete rearrangement problem:
among one-parameter coimages `x -> (e*x)_e`, interval-like offset sets appear to
maximize full sector occupancy exactly where the cap margin is smallest.

## Caveat

This does not prove LRC(14).  It only improves the continuous seven-sector
carrier.  The finite CRT placement layer remains separate: after finding a
continuous reservoir, one still has to prove the colored grids
`x=(b+14t)/(14V)` hit it, as in HYP-2593 through HYP-2595.

## Tournament Analysis

Vertices are proof routes, not runners:

`AP_frontier_envelope`, `relation_height_tail`,
`colored_resonance_placement`, `two_runner_removal_or_global_witness`,
`single_deletion_arcwidth`, `local_gap_compression`.

The pairwise observable is current closure power plus exact stress evidence
minus known false-target risk.  The resulting tournament is transitive
(`score_hist=0..5`, no directed `3`-cycles, singleton SCCs, one Hamiltonian
path), led by `AP_frontier_envelope`.  The quotient preserves the continuous
sector-cap proof obligation and destroys CRT color, finite denominator, and
original speed labels.

See also HYP-2603, THM-532, HYP-2599, HYP-2595, HYP-2593, and OPEN-Q-108.
