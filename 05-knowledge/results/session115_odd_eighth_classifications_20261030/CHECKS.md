# Verification and reproducibility

Run `python run_checks.py` from this packet. It executes sixteen scripts normally and with `python -O`, requires successful exits and compares output bytes. Python 3 and SymPy are needed; most checkers use only the standard library. Explicit exception gates remain active under optimization. The combined receipt records actual runs and output hashes.

## Proof review and finite universes

- **Norm/Boolean phases:** 11,796 producer gates, including sixteen normalized gauges, 128 zero-constant gauges, all 256 Boolean functions in dimension3, linear-coordinate actions and uniform-defect controls through dimension7. The [independent phase audit](audit_phases/AUDIT.md) adds 66,103 gates, exhausts all Boolean functions in dimensions3 and4, and verifies the dimension boundary. Written proofs establish the arbitrary-dimensional characterization.
- **Word-address transports:** 89,001 producer gates, all 21,845 multiplier residue pairs through depth8, full multiplication tables through depth5 and 504 first-defect cases. The independent literal shortcut-word decoder adds 44,313 gates through depth7, including the signed odd-coset relation and first-failure formula.
- **Finite groups:** 222 producer gates and 152 [independent group gates](audit_groups/AUDIT.md), covering all 45 odd subgroups of GL₃(2), their affine groups, actual element orders, centralizer weights and Sylow-D₈ controls. The universal group proof explicitly checks the hypotheses of the cited Burnside transfer criterion; finite examples are not a replacement for that theorem. No arbitrary switching-group realization is claimed.
- **Harmonic LRC:** 7,578 producer gates over 274 preserving rows, 138 global maximum rows and affine-period/character controls. The [independent audit](audit_lrc/AUDIT.md) adds 8,531 gates: 390 literal-grid rows, 94 full maxima, 21 period rows and 99 character identities. All maximizing points and their actual endpoint owners are checked. The lower-clearance chamber remains outside the classification.
- **Collatz reset2 and constants:** respectively 155,094 and 23,442 producer gates. [Independent review](audit_collatz/AUDIT.md) adds 50,894 literal affine/entry/digit gates. These verify the stated all-stall bounds and equality classifications; parameter-Haar means are not ordinary integer averages.
- **All reset letters:** 329,135 producer gates, including 320 exact maximal formulas, 40,960 entry-bound and 40,960 equality tests, 100,464 literal intersection tests and 140,700 mixed-radix comparisons. The [separate audit](audit_collatz/GENERAL_RESET_AUDIT.md) adds 82,904 gates, with reset letters through41 and high-valuation controls101,501,2501, literal depths through4 and all entries of clock at most4. The proof, not that finite range, covers every positive entry and all a>1,L>=0.
- **Cross-field synthesis:** 1,136 producer gates with exact cyclotomic arithmetic, actual affine orders, named character readings and all sixteen odd-residue multiplier pairs. The [independent synthesis audit](audit_synthesis/AUDIT.md) adds 5,923 gates. The raw OR versus normalized-gauge observer shift is explicitly checked.

## Concurrent incoming audit

The final ancestry preserves f128982d824452a919602687554ffea9ef33155c. Its conditional-source energy and phase decoder receive their original credit. The [narrow independent audit](audit_incoming/AUDIT.md) reviews E1–E3 and the multiplicity-sensitive decoder. It adds **numerical**, not exact-arithmetic, supplements: 93 finite conditional Gram universes, 7,914 matrix entries and 48 remainder/decoder controls, with maximum floating error 1.12e-16. This sixteenth script is explicitly distinguished from the preceding exact checkers. Normal/optimized byte agreement does not turn floating calculations into a proof certificate. The written conditional-Haar proof was reviewed separately; the source's full checker was not replayed here.

## Important repaired and retained boundaries

A prepublication synthesis draft identified the correct relative gauge class but omitted its observer shift. The published statement is p_OR(0)=p_s(111)=7/8, whereas the normalized representative has p_s(0)=3/8. No unshifted identification is used.

The general-reset extension keeps the full intersection modulus. The Lucas rule remains exact in the 2,5,5,... radix, but the other radix chains require their own no-borrow rule. The a=2 bound below1/2 and constant-controller bound5/8 are not generalized to arbitrary reset letters by analogy.

The positive LRC weight can become nonreal under a unitary. Positivity-preserving transport and invariance of the specific weight have different congruence gates. Similarly, group centralizer odd-order weights need not be reciprocal Sylow orders.

No Lean build, CI pass, literature-priority conclusion or major conjecture closure is claimed. Every complete classification names its family and retained data.
