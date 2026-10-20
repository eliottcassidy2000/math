# Odd blue spectra: a class-level explanation

The odd-order split-spectrum obligation is closed by a proved
hollow-frame toggle-cube theorem, not by extrapolation.

Read `ODD_BLUE_SPECTRUM_AND_HOLLOW_CUBES.md` for the proof, mechanism,
choice convention, exact triangle response and failure boundaries.
`VALIDATION.md` summarizes the two independent audits.

Core outputs:

- `verification.json`: complete n3/5/7/9 counts and assertions.
- `pairing_n*.json`: every SC class, selected partner and auxiliary cube.
- `n5_auxiliary_squares.json`: exact diagram-ready two-square data.
- `n7_nonadjacency_witness.json`: the purity/mass/adjacency hostile.
- `audit.md`: first independent proof audit.
- `audit_independent/`: second independent proof and computational audit,
  including a direct n9 B=39 counterexample to the historical 3-power cap.

Reproduce the standalone main census with `python3 verify_odd_pairing.py`.
For cross-checks against the separately published full census, supply
`--inherited ../merged_metagraph_black_quotient_bridges_20261020`.
Then `python3 build_n5_squares.py` regenerates the exact diagram data.
The optional failed-switch probe takes the same path via `--census`.

The preserved invariant is blue degree, not the original blue graph.
