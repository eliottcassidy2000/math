# THM-4333 rank-three packet

This packet freezes the rank-three retained-mass census on THM-4231's exact
`181,194`-pair remainder.

Contents:

- `pair_rank3_2over27_screen.cpp`: primary wide event sweep, rank-three
  aggregation, degree screen, and exact submodular branch-and-bound;
- `pair_rank3_2over27_screen_all.csv`: one screen row for every residual pair;
- `pair_rank3_2over27_exact_v2.csv`: exact minima for precisely the `67,198`
  inconclusive screen rows;
- `verify_pair_rank3_packet.py`: closed universe, hash, identity, positivity,
  extremum, and cofinal-cutoff verifier;
- `independent_pair5070_flat_audit.cpp` and `.out`: literal-midpoint-wall,
  unpruned audit of all `C(30,9)` bodies on the canonical `(50,70)` hostile;
- the two `.out` files recording the primary and verifier summaries;
- `SHA256SUMS`: raw-byte packet manifest.

Quick verification from the repository root:

```text
python 05-knowledge/results/lrc14_rank3_surplus_third_tail_thm4333/verify_pair_rank3_packet.py
python -O 05-knowledge/results/lrc14_rank3_surplus_third_tail_thm4333/verify_pair_rank3_packet.py
```

The full primary replay consumes
`05-knowledge/results/lrc14_rank2_wall_graph_complete_pair_closure_thm4326/thm4231_remainder181194.csv`
and writes the two CSV ledgers. Compile with C++20 and `-I.` from the repository
root. Full O2 and O3 replays produced byte-identical ledgers. The independent
audit imports neither the event sweep nor the branch-and-bound.

The packet proves a finite retained-mass lower bound. It does not identify the
minimum full safe mass, cover outsider pairs outside the frozen remainder,
provide arbitrary-row entry, or prove LRC(14).
