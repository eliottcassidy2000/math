# An exact formal singular-part decomposition for the quadratic unfolding

19 October 2026. **PROVED, formal:** the decomposition and its explicit convergent scalar coefficients below. **NOT CLAIMED:** convergence of the remaining formal series, or a real transit-time error estimate from formal algebra alone. This is an extension of [parabolic_cotangent_invariant_barrier.md](parabolic_cotangent_invariant_barrier.md), Section 4.

## The correct formal space

Work in R=C[[q,delta]], with both q and delta of degree one, and use the formal map

    F(q,delta)=(q+q²+delta²,delta).

It is tangent to identity in the two variables. Its pullback logarithm is a filtration-raising derivation which fixes delta, hence is

    log F*=v(q,delta) partial_q.

For a nonzero numerical delta this construction is not a point-adic logarithm at q=0. All formal assertions are in the joint (q,delta)-adic ring.

Put h=q²+delta². Since substitution difference A(q+h,delta)−A(q,delta) is divisible by h for every polynomial A, and the ideal (h) is closed in R, the logarithm series applied to q is divisible by h. Its leading term is h. Therefore

    v=hW,   W∈R, W(0,0)=1.                                (1)

The symmetry delta↦−delta implies W∈C[[q,delta²]].

## The theorem

There are unique A,B∈C[[delta]] and U∈C[[q,delta]] such that

    1/v = [A(delta)+B(delta)q]/(q²+delta²)+U(q,delta).       (2)

The two scalar series actually converge near zero and have the explicit formulas

    A(delta)=i delta[1/log(1+2i delta)−1/log(1−2i delta)],
    B(delta)=1/log(1+2i delta)+1/log(1−2i delta).             (3)

Both are even; the apparent singularities at delta=0 are removable. Their first terms are

    A=1+delta²/3−19delta⁴/45+863delta⁶/945+O(delta⁸),
    B=1−delta²/3+3delta⁴/5−275delta⁶/189+O(delta⁸).          (4)

The interpretation of (2) is an identity in the localization R[1/h]. The regular remainder U is necessarily divergent: even its restriction U(q,0) does not converge near q=0, as proved below.

**Proof.** Weierstrass division of the formal unit W^(-1) by the distinguished polynomial q²+delta² gives uniquely

    W^(-1)=A+Bq+(q²+delta²)U.

This proves the existence and uniqueness of (2). Substitution q=±i delta is valid as a continuous map R→C[[delta]], because both variables have degree one. The sections q=±i delta are fixed by the formal flow: v vanishes on them by (1). Its derivative with respect to q along either section is therefore the exponential of the linearized vector field. At time one this yields

    exp(v_q(±i delta,delta))=1±2i delta.

The linearized vector-field coefficient has zero constant term, so formal exponential/logarithm inversion gives

    v_q(±i delta,delta)=log(1±2i delta).                    (5)

Since h_q(±i delta)=±2i delta, equation (1) gives

    W(±i delta,delta)=log(1±2i delta)/(±2i delta),

where the quotients are formal units. Evaluating the division identity therefore gives

    A+i delta B=2i delta/log(1+2i delta),
    A−i delta B=−2i delta/log(1−2i delta).

Solving these two equations gives (3). The standard branches of the logarithms at 1 show convergence and removable singularities; the formulas also show the parity. Expanding them proves (4). ∎

**Divergent remainder, with cited analytic input.** At delta=0, equation (2) reads 1/v(q,0)=q^(-2)+q^(-1)+U(q,0). If U(q,0) converged, then v(q,0)=q²/[1+q+q²U(q,0)] would converge. But v(q,0) is the iterative logarithm of q+q², which diverges by the explicit classical result recalled in the Aschenbrenner–Bergweiler primary text cited in the main note, page 3. Thus U(q,0), and therefore U as a joint series, is divergent. The scalar singular coefficients A and B converge even though the remaining series does not.

## A fully determined singular clock, with a remainder still to control

For a nonzero numerical delta near zero, a multivalued analytic primitive in q of the rational singular part in (2), whose coefficients A and B converge, is

    [A(delta)/delta] atan(q/delta)
       +[B(delta)/2]log(q²+delta²).                        (6)

For positive real delta the divergent full-passage contribution of (6) is

    pi A(delta)/delta
       =Re[2pi i/log(1+2i delta)]
       =pi/delta+(pi/3)delta−(19pi/45)delta³+… .            (7)

The equality with the real part follows directly from conjugating the two logarithms in (3). The logarithmic coefficient B tracks the total residue of the two local time forms, while A tracks their opposite imaginary parts. The corrected clock in the escape proof sets A=B=1; its next required singular correction is exactly pi delta/3.

There is a particularly useful low-order check. Weighted expansion of the logarithm gives

    v=h[1−q+(3/2)q²+delta²/6]+O((q,delta)^5).

Thus, in the localized weighted sense,

    1/v=(1+q)/h−1/2+delta²/(3h)+terms of higher weight.

The improved elementary clock adds −q/2+(delta/3)atan(q/delta) to the earlier corrected clock. This independently produces the same pi delta/3 contribution as the multiplier period.

A separate [matched-tail proof](../parabolic_escape_next_order_20261019/parabolic_escape_next_order_20261019.md) establishes the corresponding next-order transit term. That proof, and not (2) alone, controls the difference between formal time, real endpoint coordinates, and integer crossing. Formula (7) does not assert that the coefficient −19pi/45 is the entire next coefficient of the threshold expansion: endpoint-dependent even-order terms and remainder estimates must also be controlled. Nor does (2) assert a common exact analytic logarithm for numerical delta across the whole bottleneck.

## Verification and limits

`verify_formal_unfolding.py` independently computes the joint formal pullback logarithm through total degree ten, divides its inverse unit by q²+delta², and compares the resulting A and B with (3). It checks the fixed-section derivative identity (5) and the displayed low-order modified equation. The general proof is the formal division and linearization argument above; a finite computation does not prove convergence or real asymptotics.

This is a concrete bridge from the quadratic recurrence to its formal Hamiltonian cotangent family: multiplication by p gives the joint formal Hamiltonian pv, and the residue data A,B converge while the regular remainder is necessarily divergent. The distinction between a convergent singular coefficient and a convergent complete primitive is essential.
