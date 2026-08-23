# Independent audit: cubic Russell seed, quotient cocycle, and affine density wall

**Status:** FINITE-EXACT INDEPENDENT AUDIT, 2026-08-23.  The audit passes the
named scout without a mathematical repair.  It verifies two finite cells and
one named boundary chart; it neither constructs a Keller pair nor changes the
status of `JC(2)`, which remains open.

The independently implemented companion is
[`jc2_russell_cubic_seed_constant_density_independent_audit_20260823.py`](../04-computation/jc2_russell_cubic_seed_constant_density_independent_audit_20260823.py),
with frozen output
[`jc2_russell_cubic_seed_constant_density_independent_audit_20260823.out`](../05-knowledge/results/jc2_russell_cubic_seed_constant_density_independent_audit_20260823.out).

## Audit verdict

The following claims in the primary scout were rederived by a separate route
and agree exactly:

1. Taylor columns at normal grades `2,3,7` all give

   ```text
   M=[[-1,       2u],
      [-(3/2)u, 3u^2-1]],                    det M=1.
   ```

   Its inverse has polynomial entries.  The quadratic-to-cubic carry in the
   second coordinate is `-(3/2)eta_2`, and every fresh cubic pair is in the
   image of `M`.

2. For a general cubic source addition

   ```text
   X=x+xi(u)x^3,                 U=u+eta(u)x^3,
   ```

   direct symbolic differentiation verifies the full chain-rule cocycle

   ```text
   J(A_0 o Phi,C_0 o Phi)=J(A_0,C_0)(Phi) J(Phi),

   J(Phi)=1+3xi x^2+eta' x^3+3(xi eta'-xi' eta)x^5.
   ```

   Thus the coefficient quotient forgets the source determinant.  For the
   stronger binary statement, the audit uses the two representatives in the
   same formal `I_1` orbit:

   ```text
   Phi=id:                              J=1+(3/2)x,
   Phi=(chi(x),u),
   chi=(2/3)(sqrt(1+3x)-1):             J=1.
   ```

   This explicitly proves that “Jacobian is constant” does not descend to the
   bare **formal** quotient.  The algebraic `chi` is not a polynomial or a
   rational function, so this sentence must not be broadened to polynomial
   source equivalence.

3. The raw affine equations were solved through coefficient matrices, not
   `solve`.  Both the `E_1` and subsequent `E_2` matrices have full rank and
   force

   ```text
   p=3/4,       q=(9/8)u,       r=-9/8,       s=-(27/16)u.
   ```

   The resulting fixed-seed precomposition by

   ```text
   chi_3=x-(3/4)x^2+(9/8)x^3
   ```

   has first surviving debt `E_3=135/16`.  Hence the raw cell with all four
   coefficient functions affine is empty.

4. Triangular reversion of

   ```text
   (1+(3/2)chi)chi'=1
   ```

   independently gives

   ```text
   c_2=-3/4,              c_3=9/8,              c_4=-135/64.
   ```

   The derivative contribution of `c_4x^4` is `-135x^3/16`, exactly the
   missing cancellation.  The exact solution has the simple linear radicand
   `1+3x`, so its square root is not in `k(x)`.  This is the normalized
   THM-3846 Catalan/Hensel gate, used only as a consistency check; the audit
   also verifies the nonsquare directly by its odd prime valuation.

5. Direct substitution verifies both monic recovery equations

   ```text
   u^3+(2-3A)u+2C=0,
   x^3-(2/3)x^2+(8/9)x-(8/9)chi_3=0.
   ```

   Thus the unique two-row near-candidate is on the integral side of the
   hidden-coordinate divide.

6. The parameterization

   ```text
   p=3/4+2u eta,
   q=(9/8)u+(3u^2-1)eta
   ```

   is exhaustive because `gcd(2u,3u^2-1)=1`.  With affine relative cubic
   residuals `alpha,beta`, `E_2=0` bounds polynomial `eta` by degree four.
   Recomputing `E_3` from the universal Jacobian gives the successive leading
   gates

   ```text
   12h_4^2, 12h_3^2, 12h_2^2, 12h_1^2, 12h_0^2.
   ```

   They force `eta=0`; the remaining full-rank `E_2` system returns the raw
   continuation, and `E_3=135/16` again.  The precise scope is therefore
   “the complete polynomial `E_1` fibre with **affine relative cubic
   residuals**,” not every polynomial or rational cubic correction.

7. On the named exceptional quadratic chart, numerator/denominator
   valuations give

   ```text
   ord_t(u x^3)=-4,             ord_t(x^3)=-3,
   lc_t(u x^3)=lc_t(x^3)=1/(8e^3).
   ```

   Consequently a nonzero affine slope cannot be cancelled by its constant
   term; after the slope vanishes, the constant term still has a pole.  This
   holds separately in both target coordinates, so regularity on this chart
   forces all four affine cubic coefficients to vanish.  Cubic additions also
   leave the already nonconstant `E_1` row unchanged.  This effectivity claim
   is only for the displayed pole chart and raw coordinate corrections; it
   does not classify decaying rational coefficients, other boundary charts,
   or target recombinations.

## Wording audit

No correction is required.  One local logical distinction should remain
visible in any promoted synthesis:

```text
nonconstant source determinant  => the density function is not quotient data;
constant versus nonconstant      => use the identity/Catalan pair above.
```

The primary reflection contains both ingredients, so its conclusion is
sound.  The exact safe scope is:

```text
FINITE-EXACT:
  raw p,q,r,s affine cell is empty;
  full polynomial E1 fibre with affine relative residuals is empty;
  affine cubic corrections are ineffective on the named nonintegral chart.

OPEN / NOT TESTED:
  rational corrections with decay;
  unrelated pole charts;
  nonmonogenic multi-control orders;
  arbitrary polynomial cubic coefficient functions;
  JC(2).
```

## Reproduction and hashes

```bash
python3 -B 04-computation/jc2_russell_cubic_seed_quotient_constant_density_scout_20260823.py
python3 -B -O 04-computation/jc2_russell_cubic_seed_quotient_constant_density_scout_20260823.py
python3 -B 04-computation/jc2_russell_cubic_seed_constant_density_independent_audit_20260823.py
python3 -B -O 04-computation/jc2_russell_cubic_seed_constant_density_independent_audit_20260823.py
```

All four runs byte-match their frozen outputs.  The independent audit performs
`66` optimization-safe exact checks.

```text
primary script SHA-256      3d67cc6674c0f432c4018d9cd03a0f646a9f9599b4cb4c00335c160049854a49
primary output SHA-256      8883668f8c8666bc30adcc1a89fc7e95bd80c295041b05f0e8c13096d0641d8f
primary semantic SHA-256    ad0debdf06adf279ee0f7e84c9e89b8257d2b14575f394c00b519a2d2f9121c9

audit script SHA-256        c833a718f4ec614a02d144af35a463b2468595a0da555381d99b98cbf8b52411
audit output SHA-256        43579b7be2daf109ae96d9007f85b58ce1f93d2a0c4422d7cce7edae2b7b5a27
audit semantic SHA-256      62b5bbf8a6906ce830a14d3dd17dca368d2b888aeb30a86310cde367374863f3
```
