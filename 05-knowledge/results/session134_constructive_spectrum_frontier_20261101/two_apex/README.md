# Two-apex actual-word recovery

**PROVED:** twelve Hamiltonian continuation counts, with suffixes of length
at most three, recover the entire two-apex word modulo leading source 0s.
This gives an exact, polynomial-time decode-and-replay membership test for
the inherited 17-coordinate actual orbit. **Prescribed-H completeness OPEN.**

Read [THEORY.md](THEORY.md). Run either:

```sh
python3 check_word_recovery.py
python3 -O check_word_recovery.py
```

Both use only the standard library and work from any working directory.
All validity tests are explicit exceptions, not optimized-away assertions.
`word_recovery.py` exposes `probes`, `extract_codes`, `recover_word`,
`full_profile`, and `recover_full_profile`. Optional `find_count(h)` makes a
complete all-order decision within the two-apex language, using the cited
Busch/SCC spine bound; it may be expensive and does not settle completeness. Inputs are exact integer values;
nonintegral/negative extracted codes and nonreplaying candidate profiles are
rejected.

The original 36-state transfer, 17-dimensional linear representation, integer
lattice and coverage through 99,999 are inherited. The contribution here is
the explicit observable binary codes, rigidity modulo source padding, and
exact actual-orbit decoder. The exploration files outside `PUBLIC_FILES.txt`
are not publication payload or proof dependencies.
