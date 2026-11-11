# Wildcard bridge: arithmetic clocks as a finite-state tournament process

[Proof and source map](PROOF.md) · [exact controls](controls.json) · [reproduction](clock_lift.py)

The old LRC half-turn and planar chirotope lenses produce a useful three-way split:

* **Local arithmetic:** every speed triple has an explicit three-edge probability law, a one-character perturbation of uniform randomness. Cyclic occupancy lies sharply between `1/6` and `1/2`.
* **Global capacity:** real rank-two circle tournaments admit at most `2^(n-1)` classes, a vanishing fraction of A000568. This thin class family is not thereby ruled out as an H-spectrum generator.
* **Faithful dynamics:** for each fixed speed row, an explicit finite hidden Markov model predicts all scale-indexed tournament-word probabilities. Its word entropy is exactly one bit per scale. Rooted dyadic signs recover actual phases and LRC margins to error at most `2^(-K)`.

The present four-vertex class can have either next H=1 or next H=5: a concrete obstruction to collapsing the physical state to H or to a tournament isomorphism class. Finite-state size depends on the speed differences; no uniform proof cutoff follows. Exact equality witnesses need a boundary sidecar.

This is a composed structural bridge, not a proof of LRC(14), Chan–Pak, or planar JC. Standard square-wave Fourier and Markov-partition facts are credited and reproved in this setting rather than claimed as new mathematics.
