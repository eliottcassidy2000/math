# Independent audit of THM-4711 and the HYP-9330 (B),(C) routes

**PASSED for the explicit all-prime claims below.** Independent audit agent
`detector_audit`, 2026-10-18. The audit derives the identities independently,
checks the submitted proof, and runs its own exact programs. It does not
promote the rest of HYP-9330 or infer an all-prime statement from a finite bank.

Audited candidate:
`01-canon/theorems/THM-4711-legendre-wronskian-and-binomial-cube-lucas-carry-divisibility.md`.
The parent owns promotion, navigation updates, commit and push. This agent
has not changed those files or used the parent's verification script.

## Approved mathematical scope

For every prime p >= 7:

1. For u_n=binom(2n,n)^2, F_p is squarefree and gcd(F_p,V_p)=1.
2. For u_n=binom(2n,n)^3, monic gcd(F_p,V_p)=rad(F_p). This includes
   the singular factor 1-64x at p=3 mod4.
3. The squarefreeness and coprimality conclusion holds for all six Zagier
   sporadic sequences with the six parameter triples in the theorem.

Here V_p is the actual integer first-block quotient
sum_{k<p}((u_{p+k}-u_1u_k)/p)x^k, not an arbitrarily chosen companion.
The proof uses known first-block Lucas congruences for the sporadic
sequences, cited to Straub Theorem1.2/Malik-Straub; it does not assume the
unproved detector slope or general HYP-9330.

## Adversarial proof checks

- The carry-free binomial shift is valid modulo p^2. In the omitted
  k>(p-1)/2 range both central binomials have p-valuation one at p>=7.
  Raising to powers2 or3 makes the divided difference zero modp.
- The Legendre recurrence has no unaccounted terminal coefficient:
  m+1/2=0 modp kills both the squared factor and its first derivative.
- The Wronskian primitive has degree p, not degree belowp. Its derivative
  vanishes, but the primitive is 1-z^p. The final right side is
  (1-z)^(p-1), which is nonzero at every nonsingular root of F.
- Clausen is differentiated in the correct parameter direction: the
  lower parameters have slopes2 and1, whose sum is3. The tangent is
  exactly 6(H_{2k}-H_k). An arbitrary deformation or a higher derivative
  would not automatically be the same Lucas carry.
- Euler's exponent is identically1/2 throughout that deformation. It
  contributes no logarithmic tangent; its square contributes the factor
  1-z to both F and its companion in the p=3 mod4 branch.
- The quarter-factor tail is integral: after its first p-factor, its
  period coefficients have valuation2 and its harmonic companion at
  least1. Reduction through degree p-1 therefore leaves polynomials
  of degree at most(p-1)/2, so there is no hidden z^p remainder.
- The strengthened quarter Wronskian has degree2d<p and equals
  (1-z)^(2d). Vandermonde gives B(1)=binom(2d,d), a unit. Thus B is
  squarefree, gcd(B,C)=1, and the cube gcd is exactly its radical.
- The six-family recurrence argument preserves both ends. At k=1 the
  extraneous u_{p-1} term is multiplied by cp^2. At k=p the u_{2p}
  and u_p terms are multiplied by4p^2 andp^2. After division byp all
  vanish modulo p. At degree p+1 both the homogeneous and forcing
  tails vanish. No unspecified term beyond the first block remains.
- For h=1-ax+cx^2, h(F^2+x(FG'-F'G)) has zero derivative and degree
  at most2p. Hence it is Q^p, degQ<=2. Both simple roots of h are
  roots of Q, and the constant1 fixes Q=h. At a simple h-root a zero
  of F of order r<p would have nonzero indicial coefficient
  alpha*h'(alpha)*r^2. This excludes the singular roots exactly.
- The six discriminants 81,-27,64,125,16,1 and b,c are units forp>=7.
  In particular the all-prime quantifier has not silently dropped a
  prime in the declared range.

The submitted proof needed a header scope correction after the six-family
extension: the earlier Franel exclusion was removed. The remaining Apéry
exclusion means zeta(3), since ZagierD/Apéry zeta(2) is included. No false
mathematical statement was promoted, and no new MISTAKE entry is required
for the audited mathematics.

## Independent computations

Reproduction from the repository root:

```bash
python3 04-computation/cm_detector_sporadic/audit_bc_20261018l/audit_identities.py 397
python3 04-computation/cm_detector_sporadic/audit_bc_20261018l/audit_zagier.py 397
```

The saved outputs are `identities_397.json` and `zagier_397.json`.

- 77 odd primes3..397: direct arbitrary-precision binomial blocks verify
  the Legendre Wronskian, quarter and three-quarter Clausen factorization,
  its actual log companion, quarter Wronskian, exact cube gcd, and the
  alternative quadratic pullback companion identity. The unchanged
  existing `ordinary_detector_and_mechanism.py` function `mechanism`
  independently checks the stated B/C predicates on the integer sequence.
- 450 Zagier prime-family cases: every prime7..397 in each of six
  families verifies all coefficients of LF=0 and LG=-MF, including
  the terminal coefficients, and the complete Frobenius Wronskian.
  The existing `mechanism` checks squarefreeness and coprimality.
- The second checker generates integer sequence values by exact division;
  independent binomial-sum formulas reproduce every n=0..100 for
  Franel, ZagierC, Apéry zeta(2), and ZagierE.
- Hostile controls: atp7 the Legendre Wronskian primitive is the
  nonconstant derivative-zero polynomial1-z^7. Removing the singular
  Euler factor from the cube companion makes its value atz1 equal6,
  so that root is lost. Atp5 ZagierD has h=(1+2x)^2 and
  F=(1+2x)^4; squarefreeness and coprimality both fail. The discriminant
  hypothesis is therefore substantive.

The computation is finite evidence for the formulas; the all-prime proof
is the algebra and degree/valuation argument audited above.

## Inheritance and literature boundary

The existing result `05-knowledge/results/overnight9_20260906_jets_deuring.md`
§3 already uses a divided-reflection companion and a degree-bounded
Wronskian. Its companion is different from the Lucas carry V; it is an
inherited mechanism, not the exact new target.

The classical ingredients were checked in the primary reference:

- Clausen: https://dlmf.nist.gov/16.12.E2
- Euler transformation: https://dlmf.nist.gov/15.8.E1
- Straub: https://arxiv.org/html/2301.12248v1
- Caruso–Fürnsinn–Vargas-Montoya–Zudilin:
  https://arxiv.org/html/2510.23298v1

The CFVZ square factorization alone does not identify the carry tangent.
For Apéry zeta(3), the specific matching deformation and its complete
truncation law remain unproved here. A simple root of a separate singular
prefactor is not protected merely by an expression2BC+DF. This is a
precise obstruction to that shortcut, not a refutation of HYP-9330(C).

Nothing in this audit proves a uniform Collatz bound or six-colour forcing.
