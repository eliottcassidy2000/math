# Session43 cycle and geometry packet

## Useful results first

- Recovered source operation: prime-to-exponent map from the small-prime part of |Monster|/|Baby Monster|, not Collatz.
- Positive arithmetic repair: h_N(b)=max{k:b^k|N} creates the literal cycle4→2→5→3→7→4. States13 and4 have the same image2, so their fold is dynamically legitimate, although not modular arithmetic.
- Exact whole-divisor census:143 bases enter the cycle,5616 terminate at1; Φ_N(z)=5616z+120z²+16z³+4z⁴+2z⁵+z⁷.
- Complete prescribed-height realizability theorem: combine lower prime-exponent bounds, then test one upper witness per base. Unique least N under divisibility when feasible.2240 independently checked systems.
- Separate reversible completion: the(5,2) monodromy passport becomes a connected degree7 double-star Belyi map with S7 monodromy. It realizes a cycle analogy, not a Monster representation.
- Sphere domain and hyperbolic(10,2,7) branch orbifold coexist. This gives a precise distinction rather than declaring one geometry wrong.

## Files and replay

`PROOF.md` contains definitions, proofs, boundaries, and source attribution.
`PROVENANCE.json` pins repository sources and external references.
`verify_cycle_geometry.py` checks the whole divisor census by independent repeated divisibility, finite monodromy groups, passports, and rational critical values.
`verify_realization.py` independently compares the exact finite constraint criterion against an exhaustive exponent box.
`verification.json` and `realization_verification.json` are normal execution receipts; optimized executions byte-match them.

Run both scripts with ordinary Python and with `python -O`. Both use only the standard library, explicit exception checks, exact integers, and rational arithmetic. No assertions or float literals occur. No Git writes were made by this worker. Native tournament-port realization remains unproved.
