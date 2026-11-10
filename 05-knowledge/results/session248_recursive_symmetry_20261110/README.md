# Session248: recursive tournament symmetry

**COMPLETE, 10 November 2026.** Five connected proved-in-note mechanisms
are ready. No independent canon audit is claimed; inherited results are
credited inputs, not re-verified. Chan–Pak and LRC(14) remain open.

[Full synthesis](SYNTHESIS.md).

- [Growth oscillator](OSCILLATOR.md): an exact creation/annihilation
  representation, with one vacuum for every newly born invariant mode.
- [Births and compatible replacements](birth/PROOF.md): exact symmetry
  exchange and an explicit one-birth C3-to-C5 replacement, plus the
  automorphism-weighted birth/death laws.
- [Hierarchical modules](modules/PROOF.md): scale-resolved symmetry loss,
  actual 3-adic actions, two rigid binary observer types at every depth,
  and symmetry-derived arithmetic charge that survives its group.
- [Metagraph spectrum](bridges/PROOF.md): exact vertex-retention spectrum,
  invariant births, and the graded birth/death commutation relation.
- [Observer code banks](observers/PROOF.md): sharp physical thresholds for
  distinguishing ensemble means, including a recursive 2k-sector bank
  that needs k/4 exceptional vertices.

These are different carriers of symmetry. Individual automorphisms,
marked hierarchy actions, invariant function spaces, and ensemble mean
symmetries are deliberately kept separate.

Small new-formula probes accompany each lane. Their finite scope does not
replace the written all-order proofs or assert a general spectrum result
for the repository's older tiling-line metagraph.

## Reproduction and scope

From this directory, run:

    python3 birth/probe.py
    python3 modules/probe.py
    python3 bridges/check_resampling.py
    python3 observers/probe.py

Python3 is sufficient except that the exact characteristic-polynomial check
in the bridge script uses SymPy. Stored outputs record the finite universes
and named controls. The new oscillator proof is algebraic; it adds no
empirical claim or inherited-census dependency.
