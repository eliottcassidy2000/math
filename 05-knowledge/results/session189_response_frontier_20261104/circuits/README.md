# Session 189: constant-region port budgets

The main result is an exact strengthening of Session 188's consumed-fanout
budget. For a certified constant region C in a live single-output binary
DAG, let a count incoming gate edges and let kappa count components of the
outside gate graph. Then

    Q = E+a-kappa+1 = beta(G)-beta(G-C) >= E,
    sum Q_j <= Phi(initial)-Phi(final).

The proof keeps the physical port incidences in a multigraph quotient and
supplies graph-minor branch sets for every allowed simplification. It also
splits the existing deletion certificate exactly:

    |C|+|B| = Q+lambda+kappa-1.

Q is consumed topology, **not additional gate deletion**. A four-gate Boolean
hostile deletes exactly two gates while Q=1. Conversely, universal shielding
preserves any eligible core through arbitrarily many exact two-gate steps
with Q=0; an explicit XOR-grid family retains Phi=(m-1)^2 and an m-by-m grid
minor. Minimum-circuit claims, arbitrary resynthesis, and uncharged affine
coordinate changes are outside scope.

Files:

- `PROOF.md`: all-size proofs, model, physical-port warning, sharp hostile,
  grid shielding, and source attribution.
- `verify.py`: standard-library exact verifier, with no removable assertions.
- `RECEIPT.json`: deterministic complete finite-check counts and retained
  named hostile/selected shielding traces.
- `run-normal.txt` and `run-optimized.txt`: identical successful outputs.

Reproduce from this directory:

    python3 verify.py
    python3 -O verify.py

Finite topology controls cover every live two-external-input distinct-pin
binary topology through six gates (64,938 closed regions), and every live
repeated-pin topology through five gates (29,202 closed regions). They check
cycle-space image/kernel ranks, the equality pattern, exact identities,
sequential inequalities, and branch-set certificates. The truth-table control
covers 6,240 live circuits with at most two gates on three inputs, all 16
binary tables and repeated unordered pins, and 162,240 proper coordinate
restriction cases. There are 11,712 nonconstant-output cases with a nonempty
constant region. Grid shielding covers m=2,...,9 and k in {1,2,3,8,16,32}.

The finite tests corroborate the all-size proofs. They do not constitute a
formal proof-assistant certificate or exhaustive large-circuit validation.
