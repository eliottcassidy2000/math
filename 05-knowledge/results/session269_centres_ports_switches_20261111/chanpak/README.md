# Chan–Pak: optimal-exponent dyadic control inside rigid strong chains

New: every odd residue modulo 2^k, k>=3, has a target-directed physical realizer on exactly k²−k+7 vertices, with all realizers strong, rigid, score-identical and having the minimum n−2 cyclic triangles among strong tournaments. Equal-cost Moon-chain switches implement triangular bit control. Partition entropy proves quadratic order is necessary within all Moon interval chains. The infinite digit limit is an exact 2-adic isometry and uniform measure parametrization.

This is modular control, not an exact-integer spectrum theorem. Every output avoids divisibility by 7, despite approximating 7 to arbitrary 2-adic precision. The exact Chan–Pak conjecture is all positive odd counts except the known universal holes 7 and 21. The stronger family-specific hostile is 63: general tournaments attain it, but no chain in this architecture does. Chan–Pak remains open.

- `PROOF.md`: full new proofs, inherited mechanisms, boundaries and architecture-relative lower bound.
- `reproduce.py`: exact residue decoder, literal adjacency constructor and independent subset-DP controls.
- `results.json`: passing tests and exact target-7/target-63 examples.
- `results-optimized.json`: byte-identical output under `python3 -O`; checks use explicit exceptions.

Reproduce: `python3 reproduce.py > results.json`.
One target: `python3 reproduce.py --target 7 --bits 8 --graph`.

Source pinned at `1ca9c5f83b55bcac7dacc0a6f5def38e02f245e7`; source theorem THM-4744 is cited with its full filename inside the proof. No external novelty claim and no re-audit of inherited canon.

Attainable-target hostile: `python3 reproduce.py --target 63 --bits 8 --graph` yields H=12644053665875775, congruent to 63 modulo 256. This is not an H=63 witness.
