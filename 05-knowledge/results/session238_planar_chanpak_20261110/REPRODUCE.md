# Reproduction and exact scopes

Run from the repository root with Python3; only the standard library is needed.
The combined runner is `python3 05-knowledge/results/session238_planar_chanpak_20261110/replay.py`.
It runs all six programs ordinarily and under `-O`, compares stdout exactly,
and requires each `--inject-failure` run to fail under both modes.

1. `cyclic_cost/check.py`: all skew sets for primes3,5,7,11; subset DP;
   independent permutation enumeration only through7; ordinal and mixed-cycle
   controls. The11 constant is reproduction, not a new record or needed input.
2. `cyclic_cost_audit/check.py`: independent literal DFS of30 circulants on
   3,5,7,9 vertices and208 invariant tournaments across six mixed profiles;
   equality characterization, lcm normalization and P7 element orders.
3. `04-computation/experiments/session238_graph_groups/verify.py`:17 regular-action
   graph cases, including C3xC3; physical orientations, completion orbits,
   selected supports and extra triangles; actual Hamilton DP in the stated
   smallest cases only.
4. `graph_groups_audit/check.py`:73 independent cases, including multiple
   vertex orbits; all16 smallest completions checked by separate cycle-packing
   and Hamilton-path algorithms.
5. `planar_actions/verify.py`:11,415 exact controls,1,596 positive words through
   weight14,66 physical tournament instances; commutators, inverse signs,
   padded cone-hostile outputs and unique binary decoding.
6. `planar_actions_audit/verify.py`:8,343 independent controls,1,596 words and48
   physical instances, using actual odd-cycle packings/root-path enumeration
   rather than the producer's endpoint DP. Includes padded hostile states.

Paths in items1,2,4,5,6 are relative to this packet. No finite result is described
as an unrestricted tournament census. All-order conclusions follow from the
proofs and their explicitly named inherited dependencies.

`REPLAY.json` records program exit codes, durations and hashes of matching
success outputs and hostile stderr. `SHA256SUMS` covers this packet, its external
experiment script and the three canon summaries, excluding itself. Check it
with `sha256sum -c 05-knowledge/results/session238_planar_chanpak_20261110/SHA256SUMS`.
The separate audit reports retain their own reviewed-target hashes.

Documentation check: `python3 agents/check_docs.py`.
Whitespace check applies only to explicit owned paths; unrelated old CRLF
checkout changes are preserved and never staged.
