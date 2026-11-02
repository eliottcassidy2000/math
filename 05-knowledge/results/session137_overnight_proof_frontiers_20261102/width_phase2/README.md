# Width-three local universality — Phase2

**PROVED; accepted by two independent all-order reviewers.**

The fixed grammar with seed K1, ordinal joins, and
`B(T,U)=C3[T,U,K1]` realizes every parity-compatible residue modulo every
integer M, with arbitrarily large strong witnesses. Consequently local spectral width
is exactly3 for every positive odd target other than1. This is scalar local
universality, not exact integer completeness: inherited `w(13)=w(19)=5` and
`w(31)=6` remain unchanged.

The construction uses actual positive contexts, two tracked polynomial
evaluations, dyadic/3-primary resets, a shear, and a safe ten-state mod5
controller. It explicitly avoids a singularity at13 and never infers CRT
synchronization from separate prime-power coverage.

See `PROOF.md`. Run `python3 -B local_universality.py` and the corresponding
`-O` command. Receipts match: 1,950 actual polynomial gate checks, 226 complete
positive-word replays, 24,688 simultaneous modular landings, and direct
vertex-level checks of the bilinear foundation and anchor seeds. Both
independent reviews and reproducible checkers are included under `independent/`.

The separate `operators.py` / `probe.json` exploratory orbit search is not a
dependency of the theorem; reported budget exhaustion is not an obstruction.
Phase1 files remain frozen in the neighboring `width` directory.
