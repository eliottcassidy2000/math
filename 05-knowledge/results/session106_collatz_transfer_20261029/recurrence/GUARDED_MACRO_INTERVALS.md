# Why the joint-source carrier closes in one dimension

**PROVED elementary generalization.** This extends the interval assembly in `SOURCE_INTERVAL_TRANSFER.md`, not the already inherited exact affine lift. It applies beyond3x+1 and supplies no convergence estimate.

Suppose a current cell has original source n=r+Mq and endpoint y=t+2Aq, with A,t odd positive, M positive even, and q in an integer interval [L,U] with 0<=L<=U (allowing U=infinity). A finite menu of authenticated guarded affine macros has entries

    y = g_i mod 2^(d_i+1),
    y_new=(alpha_i*y+c_i)/2^d_i,

where d_i>=1, alpha_i positive odd, g_i odd, and the guard guarantees a positive odd endpoint. Authentication must also establish any claimed native itinerary; rational divisibility alone is not an itinerary proof. No shared3x+1 multiplier is needed for the interval theorem.

For P=2^d, source splitting and the new state are

    s=((g-t)/2) A^(-1) mod P,
    q=s+Pq',
    (r',M',A',t')=(r+Ms, MP, alpha*A,
                   [alpha*(t+2As)+c]/P).

The guard proves exact division and oddness. The parameter interval is transformed by the same ceiling/floor formulas as before. This is the incoming arbitrary-block least-source lift, written in joint source/endpoint coordinates.

## Conjunctive affine observations close on two extrema

Allow any finite conjunction of rational affine inequalities involving the **immutable source and current endpoint**, for example

    u*n+v*y >= w.

After clearing positive denominators and substitution, each is

    (u*M+2*v*A)q >= w-u*r-v*t.

Each condition is a lower half-line, upper half-line, empty set, or all integers. Equalities are two inequalities; strict integer inequalities can be adjusted by one. Thus the conjunction with an interval is exactly an interval, computable by the largest lower bound and smallest upper bound. Under a macro substitution q=s+Pq', every earlier affine condition remains affine, and its already stored interval transforms exactly. Induction proves all finite observation histories close on the same two bounds.

This is the positive structural reason that the prior guard history can be forgotten here: **all historical inequalities act on one common scalar source parameter**. Their intersection is order-convex, and only its two extrema matter. The retained tuple still keeps the joint source/endpoint alignment. Forgetting that alignment would change coefficients and right-hand sides of later guards. This contrasts with tournament attachment ports, whose joint compatibility need not be ordered or convex. It is an analogy about sufficient continuation state, not a map between tournament objects and Collatz orbits.

## Exact boundaries

* A disjunctive guard need not preserve an interval. On q>=0, the condition `q<=1 OR q>=3` keeps0,1,3,4,... but excludes2. No single interval represents it. Taking its convex hull adds an unauthenticated source. A finite union of intervals is a safe repair, possibly with growing component count.
* Additional congruence guards must be handled by cylinder splitting/reparameterization, not interval intersection alone: even q are not one interval in q. After writing q=2q', interval closure resumes on that refined cylinder.
* For a menu with pairwise disjoint guards, child source masses add without double counting. If the guards exhaust the parent domain, the sum equals the parent mass; otherwise omitted branches retain their missing mass and cannot be declared terminated. Overlapping macros require a disjoint refinement or explicit proved weights disintegrating each source unit before treating their weighted sum as source mass. Tagging preserves branch-incidence identity but does not turn duplicated source atoms into a source partition or justify naive additive source probabilities. A finite menu of growing returns is not the complete native transition domain.
* The theorem is not a globally minimal-state theorem. It proves a sufficient exact carrier for the declared affine-conjunctive continuation language; a smaller state could suffice for a weaker observation. It does not bound integer sizes, cell counts, or the running time to reach an unbounded horizon.

The independent audit checker tests the macro equations for all five inherited growing-return words, on multiple source parameters, and a two-sided mixed source/endpoint conjunction. The all-macro scope follows from the substitutions above, not from those finite checks.
