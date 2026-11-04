# Session 188: weighted prime-factor closure

The main result is an exact local compiler from required prime valuations to
the least compatible exponent length. Its closure solves simultaneous
lower/upper valuation intervals without exponent enumeration and exhibits a
concrete nondistributive lattice of arithmetic profiles.

* `PROOF.md`: proofs, exact arithmetic examples, conditional OPN overflow
  certificates, and explicit local/global limits.
* `verify.py`: standard-library verifier with direct independent geometric-sum
  and bounded-search controls; all checks remain active under optimization.
* `normal.json`, `optimized.json`: identical PASS receipts.

Reproduce with:

    python verify.py > normal.json
    python -O verify.py > optimized.json
    cmp normal.json optimized.json

Representative result: for ordinary source 3, valuation profiles (2,0) and
(0,1) on targets (11,23) are individually realizable, while their least joint
realization is (3,1), not (2,1). Separately, prime targets (4561,368089,71)
give an explicit nondistributive M3 sublattice through orders (15,21,35).

41,040 valuation boxes were compared with exhaustive physical lengths through
180 across 190 configurations. These are implementation checks; the proofs
establish the unrestricted local results. The theorem does not ensure no
outside prime factors, global column/exponent compatibility, or control of all
possible odd-perfect supports. No OPN nonexistence or new general bound is
claimed.
