# Verification

**PASS:** `python3 -B verify.py` and `python3 -B -O verify.py` produced byte-identical `normal.json` and `optimized.json`.

The verifier is standard-library-only and imports no inherited research implementation. Exact finite-state universes are stated in `PROOF.md`; count claims do not come from a bounded tournament-order census.

- 519 clock moduli: all1..512 plus625,1024,1331,2048,3125,6561,14641.
- 70,588 direct clock-period comparisons.
- All 1099 labelled tournaments through order5 checked for prime-function signatures/parity at3,5,7,11.
- 35,168 Psi-function/signature comparisons: all 1099 labelled tournaments through order 5 at every modulus 1..32.
- 225 literal graph/subset-partition DP checks of P/C/D against the polynomial operators.
- Exhaustive actual-operation quotient orbits:6,48,686,322100 parity-tagged states at3,5,7,11.
- Literal polynomial replays of the positive zero-factor recipes.
- Explicit actual H5/H11 polynomials, disjoint-support productH55, and complementary idempotent projectors modulo5.
- Hostiles: wrong dyadic period; erasing p-adic transients; ordinary-function equality at modulus16.

Independent review completed: [ALL_MODULI_AUDIT.md](../independent/ALL_MODULI_AUDIT.md) and [ALL_MODULI_PADDING_AUDIT.md](../independent/ALL_MODULI_PADDING_AUDIT.md). Both accepted the frozen main proof and verifier; independent full-function orbits, permutation/cut profiles, and sharp universal padding checks passed.
