# Independent review: two-exponent source-flow classification

2026-10-19. **ACCEPTED ANALYTICALLY** after reading the full producer proof. The algebraicity locus and rational-map locus are correct for the specified original-source scalar-time flow, over characteristic-zero fields, with m>=2, k>=0 and nonzero time. The stated distinction between rationality of a pair and separate coordinate rationality is maintained.

## 1. Source flow and universal obstruction

The source derivation raises t-order by at least m-1, so the scalar-time exponential and its inverse exist. The invariant equation has positive degree in either coordinate when the other nonzero coordinate is adjoined, including k=0. Thus one algebraic coordinate forces both algebraic, justifying individual transcendence outside the algebraicity locus.

The dilation degree is a+1, where a=k-2m. For a!=-1, the restricted double-series and injective time-substitution argument correctly promotes algebraicity at one nonzero scalar time to algebraicity of the full formal-time pair. For a!=0, the generic-initial-constant trace argument then forces exactness of w dv on w^a=v^m(1+v). The valuation at v=-1 gives full Kummer degree in either sign of a, and the trace is taken after enlarging the constant field, so no unsupported exactness descent is needed.

Character projection reduces exactness to

    a v(1+v)R'+[m+(m+1)v]R=a v(1+v).

Poles away from 0,-1 are impossible; poles of positive order r at those points require ar=m or ar=1. This excludes all finite poles when a>0 does not divide m, or when a<0. The polynomial-degree argument and the displayed linear coefficient comparison then yield precisely the necessary divisor conditions. In the negative case the possible degree is (m+1)/b>=2, with b=-a; b=1 is excluded by its nonzero logarithmic residue. These arguments apply to all integers in the stated range, not merely checked examples.

## 2. The exceptional scaling weight

For a=-1 the source series U depends only on u and solves the autonomous equation with vector field u^m(1+u). The algebraicity descent from K(u,E) to K(u) is valid using E transcendental over K((u)) after the source substitution. The partial fraction expansion with c=(-1)^m is correct: the finite inverse-power sum has exponents at most -2, and its primitive is rational.

Consequently the displayed h and r satisfy dh=dr/r. If h is nonconstant, the finite-field norm to the rational h-field gives an impossible nonzero constant rational logarithmic derivative. If h is constant, r is constant; its series constant term is one. Then r=1 immediately forces U=u, contradicting the nonzero lambda u^m first displacement. This excludes every nonzero scalar time at a=-1 without trying to apply the nonresonant dilation argument.

## 3. Constructive positive-divisor family

For a>0 dividing m, write r=m/a and z=xu^r. Direct substitution confirms

    E=z^a(1+u),
    x=z^(m+1)/(E-z^a)^r,
    t=(E-z^a)^(2r+1)/z^(2m+a+2).

These give genuine birational field coordinates. The derivative is delta(z)=z^(m+a+2)/(E-z^a)^r. Setting s=1/z gives delta(s)=-1/[s^a(Es^a-1)^r]. Hence the polynomial primitive -integral s^a(Es^a-1)^r ds has degree m+a+1 and derivative one along the flow.

Applying the existing source exponential to this rational source primitive is legitimate in K(x)((t)): the derivation raises every integer Laurent order, also for rational x-coefficients. Its polynomial time-translation identity proves algebraicity at each scalar time directly. It does not select an arbitrary algebraic root or rely on a converse to the trace obstruction.

## 4. Constructive negative-divisor family

For a=-b with b>=2 properly dividing m+1, put r=(m+1)/b>=2 and z=x/u^r. The inverse formulas

    u=1/(Ez^b-1),
    x=z/(Ez^b-1)^r,
    t=(Ez^b-1)^(2r-1)/z²

are correct, as is delta(z)=-1/[z^(b-2)(Ez^b-1)^(r-2)]. The polynomial primitive has degree m-b. This gives scalar-time algebraicity by the same source-Laurent reasoning. Its degree is positive because r>=2.

The additional specialization safeguard is valid: after removing the maximal common tau-lambda coefficient factor from a cleared polynomial relation, domain cancellation preserves the identity, and coefficientwise evaluation in K[x,tau][[t]] yields a nonzero specialized relation. The explicit primitive identities already prevent a degenerating scalar-time claim, but this supplementary argument is sound.

## 5. Rationality of the pair

The polynomial-translation lemma is correct over any characteristic-zero coefficient field. A rational R with P(R)=P+lambda can have no finite poles, hence is polynomial; degree comparison makes it affine. The translation case forces zero shift by the next-to-leading coefficient, while an affine map with nontrivial root-of-unity multiplier has finite order, contradicting iteration by a nonzero additive lambda. Thus no such rational R exists when degree(P)>=2.

Applying this over K(E) rules out rational coordinate pairs in every positive-divisor case and all negative-divisor cases except degree m-b=1. The latter equation combined with b dividing m+1 and b>=2 gives exactly b=2,m=3,k=4. The argument does not overclaim separate nonrationality or exact degree for every algebraic coordinate.

## 6. The extra rational map and its scope

For (m,k)=(3,4), z=1/(x³t²) and delta(z)=-1. Writing D=1-lambda x³t² and B=(1+u)D²-u, substitution of z-lambda into the birational inverse gives exactly

    X=xD/B², T=tB³/D².

Both denominators are units in the source completion. The inverse rational map is the opposite-time expression. Since D is coprime to t and u, and B is congruent to -u modulo D, the denominator D² cannot cancel from T at nonzero time. This proves nonpolynomiality of the pair. The a=0 rational family is independently established by the common-u formula and its first-coordinate coprimality argument.

All rational identities are understood on their rational domains; source-completed Jacobian one is not being promoted to global polynomial regularity. The theorem consistently retains its two-exponent scope and leaves individual degrees in most algebraic cases open. No mathematical correction was required in the reviewed draft.

## 7. Final verification and independent replay

The accepted final proof SHA-256 is `3f0500d2c3321feac096ae0f427f5c3906e6e4d05789d3595f327c8669499b9c`. The companion checker source was inspected and independently run normally and with `python -O`. Both runs completed successfully on 372 controls and matched each other and the staged canonical JSON byte for byte. The reproduced JSON SHA-256 is `58d6092aa8bca9e8f3b3f22ef22b450200fa51cd8d8540969f6b15dc4fd67c9b`. The controls correctly test representative parameter identities and the extra rational example without claiming that primitive degree equals a field-extension degree. The universal divisor loci and rationality classification rest on the analytical proof. Review is complete with no unresolved mathematical blocker.
