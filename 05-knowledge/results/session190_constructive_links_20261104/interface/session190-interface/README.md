# Six-port physical parity collision

**Proved:** a ring with six forced subdivided boundary spokes and the same ring with a hub have different full profile parities but identical Hamiltonian parity in every permitted planar disk exterior. Their integer Hamiltonian counts differ by exactly twice the exterior's all-ones partial-forest total. This gives a physical measurement of the one bit lost at six ports. The parity equality also extends to arbitrary graph exteriors; the exact extraction identity remains planar.

`PROOF.md` contains the construction, exact gluing hypotheses, response-lattice theorem including the integral-kernel term, and the signed-Walsh comparison. `verify.py` independently checks 377 physical exterior cases, the complete profiles and embeddings, and the precision claims. `RESULTS.json` and its optimized replay are identical. `FROZEN_MANIFEST.json` gives hashes for independent audit.

This is not a planar-to-tournament bijection, a new general Smith-normal-form theorem, or a claim that the probabilistic Walsh squaring loss is a linear modular kernel. It credits the Session 190 planar and bijection packets as inputs.
