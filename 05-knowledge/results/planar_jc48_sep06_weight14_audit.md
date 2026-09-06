# Independent audit of the compensated weight-22 transport

**Status: PASS / analytic and exact independent reconstruction.** Root
review, September 6, 2026. The producer source and its frozen output were
read completely. THM-4308 and THM-4438 were read in full, including their
finite-depth and field boundaries. The new result is relative to the
inherited boundary existence theorem and does not reclassify that boundary.

## Analytic audit

1. The displayed low rows agree with THM-4308 equation (4) at Phi=0,
   Delta=896/15 and K=-32/5. The definition of the depth modules and their
   monomial spanning law are exactly THM-4308 equations (18)–(19).
2. The correction starts at row12. Hence unknown background rows4 and
   above cannot affect a bracket coefficient below15 or a defect
   coefficient below16. In the cubic defect, the other background factor
   has nonnegative valuation, so it cannot lower this bound.
3. Quadratic bracket variation starts at row23. Quadratic and cubic
   defect variation start at rows24 and36. This proves finite additive
   transport for every scalar, independently of the formal parameter's
   size. There is no infinitesimal-to-finite inference left unpaid.
4. All eight correction rows obey the declared degree caps. More strongly,
   all 27 supplied monomials obey the actual depth bounds and reconstruct
   every projected coefficient through row15, with no earlier rows.
5. THM-4438's bracket polynomial P and conic depth condition Q are distinct.
   Its section j=-145P/(24C) exists throughout its characteristic-zero G_m.
   The new nonzero constant beta=235202/27945 permits w=-j/beta everywhere.
   The source after translation has weight at most22, and the unchanged
   late source coefficients and conic parameter remain on the same chart.
6. Translation and its inverse preserve every defining equation in the
   stated finite solution space, so the inherited A10 terminal fibers are
   isomorphic. A zero of j over an extension field gives the identity
   translation and creates no singular locus. Rational no-response
   irreducibility is inherited only over Q.
7. THM-4438's minimal weight24 remains correct in its frozen-prefix packet.
   The new packet moves that prefix; no claim that22 is globally minimal
   has been made. Infinite lifts, full depth membership, polynomial
   termination, chart entry, JC(2) and DC(2) are all outside the conclusion.

## Independent exact path

The separate
[Fraction verifier](../../04-computation/planar_jc48_sep06_weight14_audit.py)
uses standard-library sparse polynomials with rational coefficients. It
imports neither SymPy nor the producer. Its only producer input is the
SHA-pinned immutable coefficient certificate in the frozen output; the
low background rows are independently transcribed from THM-4308.

The verifier reconstructs the complete bracket and cubic-defect variations,
all 27 depth generators and every projected coefficient, the omitted
compensator response, finite-action higher-order terms, and the exact
replacement sign. It also verifies that changing background row4 is
actually visible at the first untested bracket/defect rows. This provides
an explicit limit to the universal-background assertion.

All **93** referee checks pass normally and with optimization, producing
identical [output](planar_jc48_sep06_weight14_audit.out). The **169** producer
checks were independently rerun normally and with optimization; both
streams equal its frozen output. No research dependency program is imported
or rerun by the independent verifier.

```
python3 -B 04-computation/planar_jc48_sep06_weight14_audit.py
python3 -B -O 04-computation/planar_jc48_sep06_weight14_audit.py
source 03784e94c50a07772b0a1dce3a7961a3c0c2f51fd666073801c02f87e23808d0
output 1a1d5ca435e313273f5f406d395dc5ce7e7e10b17fa9466a3a974ee2e37c1929
```

No mathematical repair was needed. The accepted result is exact finite
transport with an inherited positive family, not a new full Keller pair.
