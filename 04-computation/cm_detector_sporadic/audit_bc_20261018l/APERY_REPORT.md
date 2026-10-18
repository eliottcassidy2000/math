# Supplemental independent audit: Apéry carry descent

**PASSED — PROVED for every prime p>=7.** Independent audit agent
`detector_audit`, 2026-10-18. This supplement closes the Apéry zeta(3)
obstruction left open by the first audit in `REPORT.md`. It does not
promote other cases of HYP-9330(C).

Audited candidate: `/tmp/apery_descent_candidate_20261018l.md`, intended
for insertion into THM-4711. The parent verified that it is now reproduced
unchanged except heading numbering in
[THM-4711 §8](../../../01-canon/theorems/THM-4711-legendre-wronskian-and-binomial-cube-lucas-carry-divisibility.md#8-apéry-zeta3-the-reflected-companion-descends).
The asserted consequence is
rad(F_A) divides V_A for the actual integer Lucas carry quotient,
including roots over the algebraic closure and the simple singular
prefactor roots. No claim is restricted to rational CM parameters.

## Exact object checked

Let H be the Franel period, and use normalized actual companions

```
G=V_Fr/2-(V_Fr(0)/2)H,
G_A=V_A/5-(V_A(0)/5)F_A.
```

Set t=x(1-8x)/(1+x), sigma=(1-8x)/(8(1+x)), K=(1+x)^(p-1),
eps=1 or-1 according as p=1 or5 mod6, and

```
B=G+eps*K*G(sigma)-eps*G(1/8)*H.
```

The independently audited rational identities are

```
F_A(t)=H^2/K,             G_A(t)=H*B/K.
```

Subtracting the displayed constant multiples of the periods leaves the
root-vanishing conclusion unchanged, so no supercongruence for the
constant coefficient is needed in the proof.

## Proof audit

1. The period identity follows directly from the characteristic-zero
   identity f_A(t)=(1+x)h_Fr(x)^2 and the known Lucas equations,
   F_A=f_A^(1-p), H=h_Fr^(1-p). CFVZ Lemmas5 and7 give the quadratic
   fixed field and H=eps*K*H(sigma), respectively.
2. The Franel ordinary differential equation is invariant under
   y(x) -> (1+x)^(-1)y(sigma(x)). The independent symbolic program
   computes both normalized derivative-coefficient residuals and obtains
   zero. In characteristicp, K differs from this gauge by the
   differential constant1+x^p.
3. The symmetric-square equation for products of two Franel solutions
   pulls back to the full Apéry third-order operator after gauge1+x.
   All four derivative-coefficient residuals vanish independently.
   The nonzero multiplier is
   x^3(x+1)^2(8x-1)^3/(8x^2+16x-1).
   This verifies the entire operator, not just one matched solution.
4. Formal logarithms can be interpreted as differential symbols with
   prescribed derivative. Define logt=logx+logsigma; since
   x*sigma=t/8, its derivative is t'/t. No convergence of a characteristic-p
   logarithm series and no choice of a numerical log8 is assumed.
   The normalized analytic companion is fixed afterwards at t=0.
5. Adding the reflected logarithmic Franel solutions gives Hlogt+B.
   The symmetric-square product and the gauge1/K yield the Apéry
   logarithmic solution (H^2/K)logt+HB/K. The characteristic-p gauge
   differs from1+x by the differential constant1/(1+x^p).
6. K*K(sigma)=1 and the twisted reflection fixes both H and B.
   Therefore HB/K is sigma-invariant and belongs to F_p(t).
   Clearing the linear denominator of sigma shows degB<=p-1 exactly.
   HB/K has poles only at-1 and infinity, each of order at mostp-1.
   They are the two simple poles of t. The descended rational function
   is thus a polynomial R(t) of degree<=p-1. This is the decisive
   global degree argument; a congruence modulo x^p would not suffice.
7. B(0)=0 fixes R(0)=0. The actual Apéry block difference satisfies
   the same inhomogeneous logarithmic recurrence in degrees1..p-1.
   At k=1 the term outside the block has coefficientp^3 before
   division byp, and vanishes modulo p afterwards. Each new unknown
   has coefficient k^3, a unit. This determines R=G_A through their
   entire degree range. No Frobenius-power ambiguity survives belowp.
8. Every finite root of F_A has a lift to a point different from-1
   and infinity, where K is nonzero. The first identity forces H=0;
   the second forces G_A=0. Ramification introduces no exception.
   Restoring the subtracted multiple of F_A proves rad(F_A)|V_A.

No substantive proof repair was required. The direct parameter-deformation
route remains unavailable here; the complete pole-divisor descent supplies
an alternative proof of the same target.

## Independent exact reproduction

```bash
python3 04-computation/cm_detector_sporadic/audit_bc_20261018l/audit_apery.py 397
python3 04-computation/cm_detector_sporadic/audit_bc_20261018l/audit_apery_operators.py
```

Saved outputs: `apery_397.json` and `apery_operators.json`.

Every prime7..397 (75 primes) passed the complete cleared-denominator
period and actual carry identities, both reflection-equivariance identities,
and the unchanged existing `ordinary_detector_and_mechanism.py` function
`mechanism`. These are coefficientwise polynomial checks, not tests only
at finite-field evaluation points. Independent integer recurrence values
for Apéry agree with the binomial sum at every n=0..100.

The hostile reciprocal-gauge control fails atp7, as it should: replacing
H^2/K by K*H^2 does not satisfy the cleared polynomial identity. The
correct sign of the gauge exponent is derived above. The inspected HTML
rendering of the proof of CFVZ Theorem2 displays a reciprocal discrepancy;
this audit relies on their stated period identity and Lucas equations,
which unambiguously give H^2/K, rather than that displayed intermediate
line. No assertion about the source PDF's typography is needed.

The finite bank supports the implementation; the all-prime conclusion
rests on the complete differential, field-descent, pole and recurrence
proof. Domb, Almkvist–Zudilin, other unresolved families, and the split-prime
slope law remain outside this audit's approved scope.
