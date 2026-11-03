# Independent all-order proof audit for the forbidden values 7 and 21

**PROVED previously established exclusion / FINITE-EXACT independent replay.**
Every finite tournament has a Hamiltonian-path count different from7 and21.
This was already proved in [THM-1370](../../../01-canon/theorems/THM-1370-h-spectrum-omits-7-21-all-n.md);
the full assertion that every other positive odd value occurs remains OPEN.

This session supplies a compact [self-contained proof](PROOF.md), including
the strong-deletion lemma, and a fresh standalone census. No imported
atlas, prior transfer matrix, automorphism-normalized count or source
selector is used to generate the finite exclusion certificate.

- The complete Hamilton-frame cover has2,8,64,1024,32768 presentations at
  orders3..7, with independent path-count implementations agreeing on every
  one. A separate Python permutation audit checks all532 classes at orders1..7.
- Strong minima3,5,9,15,25 and absence of21 at order6 establish the finite
  head. Elementary strong deletion and SCC multiplication establish all orders.
- The reconstructed56 minimal induced obstructions agree exactly with the
  older [Session133 classification](../session133_fractran_spectrum_refinement_20261101/spectrum/README.md).
  Their orders are6,7,8,9 with multiplicities29,13,13,1. This is recertification.
- A [marked core-counting refinement](PORTS.md) retains one-port and joint
  two-port orbits and gives exact rational series throughout H<=21.

Reproduce from the repository root:

```powershell
python -X utf8 -B 05-knowledge/results/exclusion_7_21_20261103/reproduce.py
```

The C++17 compiler is required; the Python verifier uses only the standard
library. [Audit, recovery and boundaries](AUDIT.md) identify the genuinely
open continuation. The proof audit is by this agent with independent
implementations, not a separate-agent or Lean verification.
