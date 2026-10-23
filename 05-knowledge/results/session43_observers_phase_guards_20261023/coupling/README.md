# Session43 restricted graph/tournament coupling

Read THEOREM.md for the rigorous statements and limits. Run `python3 check.py` and `python3 -O check.py` from this directory (Python 3 with SymPy). No network or repository mutation is performed. The output is deterministic JSON and should match receipt.json and receipt.optimized.json. explore.py supplies explicit ordinary and affine permutation actions, orbit labels, and the incidence-based chain; its direct entry point is optional exploratory output.

The receipts retain exact incidence matrices and characteristic polynomials, not numerical eigenvalue approximations. For soft guards, word_count means total integer weight (accepted words have weight 2, rejected weight 1). All finite decisions use rational arithmetic and explicit exceptions, not assertions.

Scope: quotient coupling, not native metagraph adjacency. Conditioning, coarse observations, marked actions, and restrictions of the old linear operator are explicitly distinguished. No Git writes were performed by this research lane.
