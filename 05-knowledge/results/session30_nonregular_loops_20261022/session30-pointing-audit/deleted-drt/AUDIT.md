# Independent audit of the deleted-DRT formal constants

Status: PASS, conditional on the previously independently certified five/26-core
kernel definitions. This audit concerns those finite resummed predictors, not
the exact Hamiltonian count, and makes no novelty claim.

Source audited: session30-deleted-drt/DELETED_DRT_FORMAL_CONSTANTS.md and its
symbolic scripts. No finite tournament sample or numerical extrapolation is
used in this audit.

## Analytic checks

1. E=(J+ddᵀ)/n is an orthogonal projection: 1 and d are orthogonal and each has
   squared norm n. S² acts as−I on this subspace and−(n+1)I on its complement.
   Therefore the stated scalar G,H,O resolvents follow. For same-sign unequal
   vertices, SE has zero entry and E has entry2/n. For opposite-sign vertices,
   G and H have zero off-diagonal entry. The count of ordered unequal same-sign
   pairs is n(n−2)/2. These facts justify the exact scalar X1 reduction.

2. Endpoint signs are correct: R−=R+ᵀ, so summing a leaf at the first endpoint
   switches the leading +zd/n row sum to−zd/n. Q leaf sums are O(n⁻²) because
   S²1=−1. In contrast, R leaf sums are O(n⁻¹). Summing these leaves before
   applying entrywise estimates is essential and is correctly done.

3. After terminal elimination every weight-two core has E=v+2 factors, with v
   remaining summed indices. The audit script checks this and that no edge
   connects two terminal vertices for all26 supplied cores. Every error-factor
   replacement costs another n⁻¹. A nonloop diagonal contribution imposes one
   equality of previously free indices and likewise costs n⁻¹.

4. Removing even repetitions of an S edge produces only diagonal-collision
   errors of order O(n⁻³) after normalization. For a remaining simple edge,
   conditioning on all other indices produces bounded functions f,g and the
   rigorous estimate |fᵀSg|≤n√(n+1). The v−2 other summations then bound the
   normalized contribution by O(n⁻⁵ᐟ²). This argument does not assume random
   vertices or probabilistic independence. It accommodates repeated vertex
   labels and the bounded d weights.

5. When no signed edge remains, odd d powers sum to zero and even powers sum
   to n. Thus the survival criterion used by the finite contraction is valid.
   All constants in these finite calculations can be uniform on bounded
   z-compact sets, since the exact resolvent denominators stay away from zero
   for sufficiently large n.

## Independent finite arithmetic

The script audit_constants.py imports the certified tensor data but no producer
code. It separately contracts the leading local factors and verifies the four
survivors:12z⁶,14z⁶,3z⁴,4z⁴. Their sum is7z⁴+26z⁶.

It also independently applies the Euler differential polynomials directly to
exp(X0) and exp(X0)A1. This avoids reusing the producer's Bell-polynomial or
mixed-shift implementation and checks those normalizations independently.
It obtains, with u=1/n:

- Y1=−3u/2+13u²/3;
- Y2=35u²/24;
- κ1/e=1−5u/2+20u²/3;
- κ2/e=1−5u/2+65u²/8.

Independently expanding the stated exact rho formula gives
log rho=1−u+u²/3+u³/3+O(u⁴), hence rho/e=1−u+5u²/6−u³/6+O(u⁴), consistent
with both predictor conversions.

The resulting first- and second-order y coefficients3/2,−13/3 and−139/24 are
correct for these predictors. The second-order analytic remainder is the proved
O(n⁻⁵ᐟ²) in κ2, or O(n⁻³ᐟ²) in y2, not a stronger remainder inferred from a
truncated symbolic polynomial.

## Limit of endorsement

The finite-kernel construction itself was audited separately and is an explicit
input here. More importantly, none of this bounds the omitted infinite forest
gas uniformly at z=1. The exact Hamiltonian count's uniform asymptotic expansion
remains open. The source document preserves this distinction correctly.
