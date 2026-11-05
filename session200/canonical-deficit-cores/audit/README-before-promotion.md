# Session 200: canonical physical deficit cores

Proof candidate and reproducible implementation, awaiting independent audit.

- `PROOF.md`: sharp endpoint/height dichotomy; canonical sink-core projection;
  exact physical dyadic-ray joint fibers; constructive C1 target compiler;
  quantitative obstruction to uniformly small deficit.
- `verify.py`: complete finite inventories, minimum J costs, formula and
  hostile controls, literal graph tests.
- Frozen dependencies: `CARMICHAEL-DEPENDENCY.md`, `JOINT-DEPENDENCY.md`.

Reproduce:

    python3 verify.py --cap 5001 > VALIDATION_NORMAL.json
    python3 -O verify.py --cap 5001 > VALIDATION_OPTIMIZED.json
    cmp VALIDATION_NORMAL.json VALIDATION_OPTIMIZED.json

The full {7,21} coverage theorem remains open. No external conjecture or
Diophantine-definability theorem is assumed.
