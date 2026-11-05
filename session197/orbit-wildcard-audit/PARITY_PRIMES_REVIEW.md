# Independent review: parity-sensitive primes

**PASS**, no corrections requested, 2026-11-05.

Reviewed frozen files in `session197-orbit-wildcard/`:

- PARITY_PRIMES.md: cb894848471c7ba99b2dd8a80906c54bf56e0e24d3f3afa7e9eabf5f9984d91c
- verify_parity_primes.py: cd20da7d808c37c8218a76ff63cb04de464fcfe2a87ff30602a9e091bb781ccb
- parity-primes-normal.json: 8b07b3fe250533f50e6b65eb05ad11cf1d6e747d91fa6ad20d600561e7704b6d

The unconditional gcd follows from the unitriangular Stirling matrix and the already-reviewed positive context resolution. The a-to-beta change of basis is integer unitriangular in the reverse index order, so it preserves the exact gcd. The first D odd or positive-even clone values have matrix V_e with entries b^e(b^2)^t and determinant product b^e times product_(i<j)(j-i)(j+i). No prime greater than 2D-1 divides this determinant. Invertibility after localizing at such a prime therefore preserves the minimum coordinate valuation, proving the stated equality for both parities and the unrestricted gcd.

The adjugate argument also proves the fully integral bound g | G_e | g*Delta_e. Since g | G_all | G_e, the relative parity divisor G_e/G_all divides Delta_e. The edge case D=1 is correct: Delta_e=1 and there are no parity-specific prime factors.

The inherited clone-matching identity converts 2D-1 to 2nu_r+1. Its geometric size bound uses nu_r<=floor((q-1)/2), giving q-1+2D<=2q. These are physical certificate sizes for universal divisibility, not inverse-spectrum witness bounds.

Normal and optimized verifier outputs agree byte for byte. The controls exhaust 486 scaled formal coefficient vectors and check the three named physical templates and their matching degrees. Formal vectors are correctly labelled algebraic controls rather than asserted graph realizations. The all-order prime valuation theorem is proved by the determinant argument, independently of the finite controls.
