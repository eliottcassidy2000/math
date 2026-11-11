# Reproduce the structural controls

Python3 standard library only. From any directory:

```
python3 /path/to/repo/05-knowledge/results/session273_tournament_foundations_20261111/reproduce.py
```

The wrapper runs four programs normally and with `-O`, requires byte-identical
stdout, parses it as JSON and compares the saved result. The optimized replay
checks packaging consistency; assertions in normal runs supply checks, and
`REVIEW.md` records separate mathematical inspection. All-order statements are
proved in the notes, not extrapolated from these controls.

- Orbits: all labelled tournaments through5; signed-orbit conditional formula
  versus direct completions, exact transition normalization and hostile.
- Microscopic: subset Hamiltonian DP and exact triangle-chain sample counts.
- Interfaces: physical subset-DP cells and endpoint controls; larger word
  counts use SCC signatures and proved product laws, not full26-vertex DP.
- Clock:900 exact triple-gap pairs,24 finite-language controls, rational Markov
  identities,2400 margin checks and an equal-class/different-future hostile.

Repository `python3 agents/check_docs.py` also passed. Only this session's
packet and new navigation page are published; unrelated checkout normalization
artifacts are not staged. SHA256SUMS covers final packet files except itself.
