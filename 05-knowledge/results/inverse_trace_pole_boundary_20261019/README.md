# Inverse trace-pole boundary gate

**Current global scope, 2026-10-20 (CITED).** [Domrina's 2000 theorem](https://www.mathnet.ru/eng/im273) excludes complex polynomial Keller maps of generic/function-field degree four. The repository's [primary-PDF source record](../../reference/CORE-PAPERS-PLANAR-JC.md) confirms that degree convention. Quartic A4 is therefore only a residual of this local method, not a globally open Keller case. The local pole models and conditional trace gates remain valid; the live global program must address higher-degree families, without claiming a particular first open degree.

19 October 2026. **PROVED scoped synthesis; independent review PASS.**

For a genuine complex polynomial Keller map of generic degree d, one actual finite inverse sheet survives at every affine target divisor. Combining that fact with Newton identities gives:

- In a primitive source coordinate x with P monic in the transverse coordinate, the pole divisors of Tr(x^k), 1≤k≤d−1, are exactly the nonproper-value set.
- The pole order of an inverse coordinate is exactly recoverable from those trace moments; maximal delayed detection forces a root-of-unity pole packet.
- In even degree d, polynomiality of the first d−2 moments forces alternating monodromy. In degree four, S4 therefore requires a pole in Tr(x) or Tr(x²). A4 remains only as this method's residual; Domrina globally excludes degree four.
- A detected trace pole propagates to all higher pole orders along a polynomial dual-derivation chain.

The [proof and scope](TRACE_POLE_BOUNDARY.md) include polynomial non-Keller and rational symplectic hostiles, an exact denominator compiler, the even discriminant-order proof, and a local triple-packet discrepancy model. No unconditional trace regularity, Keller counterexample, or Jacobian-conjecture resolution is claimed. Historical priority is not claimed.

Run:

    python3 verify_trace_pole_boundary.py
    python3 -O verify_trace_pole_boundary.py

Both runs pass 808 exact checks with identical stdout. Python 3 and SymPy are required. Finite controls supplement the all-degree proofs.


## Independent reproduction

The [independent review](review/REVIEW.md) passes the stated local hypotheses. Its dated description of quartic A4 exclusion as open is superseded by the global scope correction above; unconditional trace regularity and the Jacobian conjecture remain open. Run `python3 -B review/independent_checks.py` and repeat with `python3 -B -O review/independent_checks.py`; compare stdout with the normal/optimized receipts and the generated JSON with `review/independent-results.json`. The copied author verifier and its receipts are included in the review directory as well. `SOURCE_PROVENANCE.json` records the original reviewed/published hashes and this dated scope correction; Mathematical test logic is unchanged; status, presentation, and package-relative source-hashing updates are documented.
