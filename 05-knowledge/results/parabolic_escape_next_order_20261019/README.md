# Next-order quadratic cusp phase

2026-10-19. A separate extension of the initial quadratic cusp phase packet.

## Result

For the unique positive threshold parameter epsilon_n with x_n(epsilon_n)=R, where R>1/2 is fixed, put delta_n=sqrt(epsilon_n). The theorem proves

n = pi/delta_n + C_R + (pi/3)delta_n + O_R(delta_n^2 log(1/delta_n)).

The proof gives a controlled endpoint/middle matching argument, a parameter expansion through order n^-4, and the first correction to adjacent escape-band widths. It resolves the explicitly open hypothesis (H) in the initial checkpoint. It does not calculate a numerical universal constant in the O_R term or establish the next coefficient.

## Files and reproduction

- `parabolic_escape_next_order_20261019.md`: full analytic proof and its scope.
- `verify_parabolic_escape_next_order_20261019.py`: modest numerical controls.
- `parabolic_escape_next_order_20261019.out`: successful captured output.
- `MANIFEST.sha256`: file integrity hashes.

Run `python verify_parabolic_escape_next_order_20261019.py`. Dependencies: Python 3 and mpmath 1.3.0; tested with Python 3.12.14. Runtime is a few seconds. No input files, network, or external state changes are needed. The verifier uses a phase interval established in the initial packet, and all calculations printed here are explicitly labeled NUMERICAL. The proof does not depend on them.

The [initial companion note](../parabolic_escape_precision_20261019/parabolic_escape_precision_20261019.md) defines and proves the convergent constant C_R, establishes the leading bounded error and threshold phase, and gives the phase interval certificate. The next-order proof imports those established results rather than restating their full proofs.

## Review and attribution

The analytic argument has received mathematical review and the verifier passed. A [separate independent audit](review/independent_phase_review.md) has passed. The [review packet](review/README.md) includes a self-contained degree-eight phase enclosure proof and checker. No worldwide novelty claim is made: the leading pi-law and bounded-error refinement are established prior literature; this is an independently derived, new-to-session elementary next-order proof.

The local multiplier period independently explains the pi/3 coefficient, but the proof does not infer a global transit theorem merely from that local complex expansion. The quantitative matching estimates are essential.
