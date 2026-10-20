# Session 17 finite-field surface addendum

This is a separate extension of the marked F5/icosahedron bridge. It does not
replace or modify the frozen `session17_prime_shapes_20261020` packet.

- `FINITE_FIELD_SURFACES.md`: uniform proof, source/target and loss contracts,
  literature boundary, and hostile examples.
- `verify_surfaces.py`: exact small-field checker, using only Python's standard
  library.
- `receipt.json`: 4,544 passing checks at orders 3, 5, 7, and 9, with a
  characteristic-two negative control.
- `LITERATURE.json`: checked primary-source locator and imported scope.
- `SHA256SUMS`: hashes of the other five payload files, including this README.

Main result: the determinant-±1 triangle complex over Fp is a closed oriented
surface of genus 1+(p²-1)(p-6)/24. Over Fq, q=p^r, every vertex link splits into
q/p p-cycles, and normalization gives q(q²-1)/(p(p²-1)) copies of that prime-field
surface. In particular F9 normalizes into 30 tetrahedron boundaries and F25 into 130
icosahedron surfaces. The graph carrier is classical, so no novelty claim is made.

Status: proofs complete; independent review requested before publication.
Target: `05-knowledge/results/session17_prime_surfaces_20261020`.

Reproduce with `python3 verify_surfaces.py`; compare its output with
`receipt.json`. All checks also run under `python3 -O`.
