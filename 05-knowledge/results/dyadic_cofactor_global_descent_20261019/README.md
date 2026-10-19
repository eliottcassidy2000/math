# Global descent for the dyadic cofactor map

The [proof](COFACTOR_MAP_GLOBAL_CLASSIFICATION.md) classifies every orbit of the auxiliary map F(n)=n-p+2c-1, where n=cp, c is its exact power-of-two factor, and p is odd. Each orbit reaches a fixed point c(2c-1), the four-cycle 8→22→14→10→8, or the three-cycle 16→46→26→16. Every nonfixed start outside {2,4,8,12,16} has a smaller iterate within six steps; 48 shows that six is sharp.

F agrees with the aliquot sum only when p is an odd prime. Both nontrivial auxiliary cycles violate that guard. Consequently no nontrivial aliquot cycle of any length lies entirely in this class, and a genuine aliquot trajectory that stays in the class forever reaches an even perfect number. This strengthens the [earlier A3/D3 three-cycle theorem](../dyadic_aliquot_three_cycle_exclusion_20261019/README.md). It gives no conclusion about an aliquot trajectory after it leaves the class, or about general Collatz convergence.

## Reproduce

Python 3 and SymPy are required. Run from this directory:

- `python -B verify_cofactor_descent.py`: compare stdout with `verification.json`.
- `python -B -O verify_cofactor_descent.py`: compare stdout with `verification-optimized.json`.
- `python -B review/audit_cofactor_global.py`: compare stdout with `review/cofactor_global_audit_results.json`. This checker also writes that receipt next to itself. Keep assertions enabled; do not use `-O` for the independent audit.

The [independent audit](review/COFACTOR_GLOBAL_INDEPENDENT_AUDIT.md) checks the exhaustive proof cases and well-ordering argument. Its finite controls include two million even starts, 150,558 structured and large-integer checks, a separate orbit scan, and native-prime guard controls. These computations supplement the proof; they do not supply its universal quantifier. Historical priority has not been established.

`SOURCE_PROVENANCE.json` records the reviewed hashes and the status, presentation, and portability updates. Mathematical content and test logic are unchanged.
