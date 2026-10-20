# Finite dilation and power-equivariance repair

**Status: PROVED in the scopes below; independently reviewed.** This is the corrected successor to [THM-4725](../../../01-canon/theorems/THM-4725-no-uniform-dilation-rule-for-prime-circulants-signed-set-burnside-and-the-p-mod-4-dichotomy-with-the-paley-rule.md), as it stood at commit `127698803818637f07ae89c32ac5d04c59161292`.

## Results

1. Every fixed finite signed dilation stencil fails to give an even-circulant type bijection at infinitely many primes p = 9 mod 16. The repaired proof uses positive rational radicals, correct Kummer correspondence, a total-positivity/norm obstruction, and Chebotarev. Negative dilation constants complement input bits.
2. Signed-set Burnside and the Paley Gamma-equivariant construction survive. The set-theoretic equivariance dichotomy includes the exceptional prime p=5: existence at p=3 mod4 and p=5, failure at p=1 mod4 with p>=13.
3. Parameter-free multiplicative formulas are power-equivariant, not automatically multiplier-equivariant. The obstruction applies if their outputs descend to a well-defined bijection on multiplier types, via the quotient power action.
4. The particular Paley construction uses parity and is not uniformly first-order definable in the ordinary multiplicative language. Canonical equivariance is distinct from first-order definability.

**Scope:** This does not rule out arbitrary multiplier-equivariant algorithms or every quantified multiplicative rule with named constants. It does not prove that only additive constructions remain. It is not the earlier affine F2-ideal obstruction.

Read [the repaired proofs and full audit](REPAIR_AND_AUDIT.md), [review record](audit/REVIEW.md), and [primary sources](PRIMARY_SOURCES.md). The exact [independent output](independent_controls.out) reproduces every Gamma-orbit row through p=37, with signed-stencil and small-prime hostiles.

Reproduce from this directory with `python3 independent_controls.py` (Python 3 and SymPy). The script exits nonzero on an assertion failure. Check file integrity with `sha256sum -c SHA256SUMS`.
