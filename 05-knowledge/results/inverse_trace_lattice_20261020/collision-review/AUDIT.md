# Independent audit: collision gap and singular graph locus

20 October 2026. **PASS, within the stated actual polynomial Keller-map hypotheses.**

Reviewed addendum: `session11-jacobian/COLLISION_GAP_AND_SINGULAR_LOCUS.md`.

SHA-256: `6d6c29f5b1f879376a1c584971562e0a478e4f3bce69172ef01058e1b60bed59`.

This bounded audit is separate from the finite-selector and source-line audits. Their proofs and receipts were preserved.

## Proof review

### 1. Exact gap-unit certificate

The divided differences U,V are polynomials, including on the diagonal. Their diagonal values are P_y,Q_y. Because J=P_xQ_y−P_yQ_x is constant and nonzero, the numerator defining Ω vanishes identically on s=t and is divisible by Δ=s−t. The resulting identity gives Δ(−Ω/J)=1 in the full collision quotient, with no reducedness assumption. Hence the saturation identity I:Δ^∞=I and the absence of diagonal points are correct.

As an independent structural check, this identity also gives (Δ)+(U,V)=(1). The full graph fiber product has ideal (ΔU,ΔV), which is consequently the intersection (Δ)∩(U,V). Thus its diagonal and collision pieces split as disjoint closed pieces. This agrees with the addendum's direct finite-map proof. A unit certificate supplies no positive analytic lower bound on the gap.

### 2. Finiteness and exact singular image

Monicity of P in y makes U monic in t of degree deg_y(P)−1. Therefore T is finite over the source plane (a,s), and composition with the finite graph map G makes T finite over W. When that degree is zero, U=1 and T is empty. Passing to T_red preserves finiteness.

The graph map is finite birational, with normal source A², and is therefore the normalization of W. A smooth point of W has exactly one normalization preimage. Conversely, if its fiber is a singleton, finiteness confines all nearby preimages to one source inverse-chart neighborhood, and the Keller inverse function theorem makes the complete local graph image smooth. Thus singular points are exactly the graph points with more than one source preimage. Their source x-coordinates agree, producing precisely the ordered collisions parametrized by T. The reduced collision scheme maps finitely and surjectively onto the reduced singular locus. The manuscript correctly means finite over W, not finite over C or zero-dimensional.

### 3. The R_Z criterion and specialized singularities

Differentiation of R(P,Q,x)=0 gives the two displayed chain-rule identities. The coefficient matrix in R_u,R_v has determinant J. At a source point, R_Z=0 therefore forces R_u=R_v=0; the reverse implication is immediate. Surjectivity of G transfers this to all graph points, proving the set equality Sing(W)=V(R,R_Z). It is not asserted as an equality of unreduced ideals.

On Z=a, the derivatives of R_a are R_u,R_v. If these vanish at a section point, any source preimage and the first chain-rule identity force R_Z=0. Conversely, all graph partials vanish on the singular locus. Because the source-line packet proves R_a reduced, its singular set is exactly Sing(W)∩{Z=a}. The reducedness requirement is retained rather than inferred merely from irreducible support.

### 4. Every parameter occurs; no properness claim

The previously verified one-line injectivity theorem implies that a noninvertible Keller map has a collision on every affine source line x=a. Hence the parameter map from T, and equivalently from its finite image Sing(W), is surjective on all complex points. This is stronger than dominance. Each parameter fiber is finite because a finite birational line normalization has only finitely many singular-image identifications and finite fibers over each. These facts do not make the projection proper or finite; the note explicitly avoids that inference.

### 5. Boundary exclusion and the second unit

At H=0 the graph equations reduce to R=AZ+B and R_Z=A. Their simultaneous vanishing would force H=A=B=0, which the frozen gate excludes. Thus Sing(W) lies in D(H), and no actual ordered collision lies over the nonproperness curve. The polynomial H(P(a,s),Q(a,s)) has no zero on the collision scheme's complex support. The Nullstellensatz therefore makes its class a unit in the full quotient, including any nilpotents. The same argument makes H a unit on the reduced singular locus. No explicit inverse formula for this second unit is claimed.

## Scope and external input

The graph and gap arguments use actual polynomial source coordinates, finite G, primitivity, monicity, and the Keller determinant. Boundary exclusion additionally uses the global low-trace gate. The all-parameter nonemptiness uses Gwoździewicz's one-line injectivity theorem, recorded as Theorem 2.7 in [Jelonek, arXiv:2011.03472v3](https://arxiv.org/html/2011.03472v3). That exact version and theorem were independently checked in the immediately preceding source-line audit.

The result gives a necessary collision/singularity configuration, not its general exclusion. It does not turn a dominant projection into a surjective one, a unit into an analytic bound, or a surjective quasi-finite projection into a finite map. No Jacobian-conjecture resolution or new degree exclusion is asserted.

## Bounded exact controls

The independent verifier passes 42 exact controls under normal and optimized Python, with byte-identical output and JSON. It checks the gap identity and chain-rule equations on polynomial Keller automorphisms; an actual finite immersed nodal graph family with collisions in every parameter but non-Keller Jacobian; exact singular-image equations; gap comaximality; a dominant projection missing zero; and a unit in a nonreduced algebra. The nodal family also exhibits extra smooth points of V(R,R_Z) when the Keller condition is dropped.

The producer's verifier was independently replayed from a staging copy: 192 checks pass in normal and optimized Python, with byte-identical receipts matching the producer's JSON. Its abstract surjective nonproper curve is correctly parametrized by t≠0, a=t²−t, s=t+t⁻¹. Every a has a nonzero root t, while t→0 gives bounded a and unbounded s; this verifies the precise missing properness implication. No hostile is claimed to be a Keller realization.

Only this review staging directory was written. No producer proof, earlier audit, repository checkout, Git state, or remote was changed.
