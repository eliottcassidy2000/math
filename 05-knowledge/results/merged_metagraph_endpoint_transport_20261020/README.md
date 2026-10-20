# Endpoint-score transport and complement leakage

A research packet for the coloured, fixed-Hamiltonian-path complement metagraph of converse-merged tournaments.

Read `endpoint-transport-and-complement-leakage.md` for the proofs, conventions, equality cases, counterexample and limitations. `AUDIT.md` records the independent proof and computation audit. `DEPENDENCIES.json` fixes the mathematical-source and independent-classifier fingerprints and tested runtime versions.

## Main statements

- The coloured triangle-current and folded-score row at a node is prescribed by its Hamiltonian-path endpoint-score inventory.
- Triangle-current rigidity is equivalent to regularity of the source and sink strongly connected components. This is different from path-rigidity (`H=|Aut|`) and pure-blue status.
- For every order at least four, all triangle-current-rigid nodes form an independent set. Each colour has at most half its half-edge mass at those nodes.
- Regular nodes can meet only the almost-regular score shell, with an exact endpoint-path formula for return mass.
- Complement leakage has an exact covariance, rank and equality description; a colour decomposition separates within-colour and between-colour losses.
- An order-five hostile preserves node colour degrees, loop counts and all global joint triangle/score distributions, yet fails the endpoint transport law. Its support changes, but no previously unsupported pair is introduced.

These are necessary arrangement rules, not a sufficient classification of realizable incidence. This is not the full single-arc-flip metagraph of THM-588, and no prize-problem closure or inherited spectral gap is claimed.

## Reproduction

Run from this directory with Python 3. The first three programs use only the standard library. The fourth also requires SymPy; the tested version is recorded in `DEPENDENCIES.json`.

    python3 -u probe.py
    python3 -u verify.py > verification.out
    python3 -u verify_hostile_and_colour.py > hostile_and_colour_audit.out
    python3 -u independent_transport.py > independent_transport.out

`probe.py` performs an independent score-refined canonical classification through order seven and writes `n3.json` through `n7.json`.

`verify.py` enumerates Hamiltonian paths and automorphisms, checks current and score-profile transport, the source/sink characterization, the independent-set and colour-capacity laws, and exact leakage ranks through order seven. It writes `audit_n3.json` through `audit_n7.json`.

`verify_hostile_and_colour.py` checks the hostile and the full rational colour-covariance matrix identity through order six. It also compares all tilings through order seven with the included independently generated census projections, explicitly translating tile order and bit gauge.

`independent_transport.py`, using only the separate exhaustive-permutation classifier in `independent_classifier.py`, reconstructs and checks all nodes through order six without importing the primary implementation. Its comparison with `n*.json` is additional, not required to construct its own objects.

All tests use exact integer or rational arithmetic. The leakage rank certificate uses reduction modulo a prime dividing no fibre mass, together with the proved matching upper bound. No floating-point eigenvalue test is used.

## Files

- Main proof, `AUDIT.md`, and `DEPENDENCIES.json`.
- Five primary census files, five detailed endpoint audits, and five compact independent census projections.
- Primary and independent Python programs with exact output receipts.
- `hostile_n5.json` gives the small explicit reassignment.

The order-seven triangle-leakage fraction is retained exactly in `audit_n7.json`; its decimal in the paper is only a readable approximation.
