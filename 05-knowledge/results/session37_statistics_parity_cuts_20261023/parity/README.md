# Session37 parity lane

1. `EXTERIOR_TOURNAMENT_FOCK.md`: full exterior representation of labelled tournaments, distinct exterior-admissible graph basis, genuine symmetric/exterior connected-component occupation, homological sign test, exact connected balance and dyadic recurrence. Main deliverable.
2. `SUPER_PBW_AND_SIGNED_CURL.md`: separate vertex-parity super-PBW interpretation of SCC words, required odd-square correction, signed curl-necklace criterion, alternating-group orbit refinement, and Hamiltonian-supertrace obstruction.
3. `check_exterior_statistics.py` / `exterior_statistics.json`: standard-library exact permutation and graph-orbit checks; optimized receipt matches.
4. `check_super_curl.py` / `super_curl.json`: standard-library weighted-word checks, exact rational free Lie-super span generation, and direct labelled curl orbit checks; optimized receipt matches.

Run each Python checker with ordinary Python and with `-O`. No third-party packages are required. All finite checks use explicit exceptions rather than removable assertions. Source checkout was read only, with intentional deliverables kept outside it.

The two parities are not interchangeable. The first construction's components are graph orientation lines graded by edge number. The second starts with tournament SCCs and grades by vertex number. Neither identifies a tournament with a physical boson/fermion system, and neither supplies a natural operation-preserving tournament/graph set bijection.
