# A six-plus-four split extends equal-profile arbitrary-ambient separation to order ten

**PROVED, subject to the explicit exact certificates retained here:** three actual named child tournaments with a common complete ordered-cover profile and common order10 are separated from opposite cyclic quotient by the canonical root query in every fixed ambient tournament. Target child interiors may change while preserving the common profile. Combined with the previously certified orders1 through9, the conclusion holds for every common child order at most10. No result beyond10 is claimed.

## 1. A small-input Hamilton bound

The known small extremal Hamilton-path counts are

    max_(|A|=4) H(A)=5,        max_(|A|=6) H(A)=45.       (1)

These are not new extremal values. PROVENANCE.md records their published numerical provenance. For a self-contained input, `maxima.cpp` independently enumerates all64 labelled four-tournaments and all32,768 labelled six-tournaments, computes every exact Hamilton count by subset endpoint recursion, and verifies(1). It retains the complete histograms and explicit attaining codes. The encoding is the lexicographic list of pairs i<j, with bit1 meaning i->j. No isomorphism reduction is used.

For every ten-name tournament A, split each directed Hamilton word after its first six names. If S is its first-six-name set, its prefix and suffix are Hamilton paths in A[S] and A[V-S]. This map is injective into pairs of such paths. The connecting arc is a further restriction that may be ignored for an upper bound. Hence

    H(A)<=sum_(|S|=6) H(A[S])H(A[V-S])
         <=binom(10,6)·45·5
         =47,250.                                     (2)

This is an elementary positional splitting argument using the independently checked old small inputs. It is not a claim that47,250 is the exact ten-name maximum.

## 2. The continuous order-ten coefficient certificate

The classical odd-cycle-cover expansion writes every actual ten-name profile as

    q=t^(10)+w2 t^(2)+w4 t^(4)+w6 t^(6)+w8 t^(8),
    w_s>=0,        q1=1+w2+w4+w6+w8.

Together with(2), the lower weights lie in a simplex of total weight at most47,249. No parity or integrality restrictions on the weights are required by the following certificate; it works on the entire continuous simplex.

The complete port-augmented ordered-cover contrast for the source cyclic quotient A1->A2->A3->A1 and the reversed quotient, with a=A1,c=A2, is the exact inherited functional

    D_h(q)=(h/3)L_q[K E0^(h-1)/(1-XYZ)^(h+1)].           (3)

Normalize w_s=47,249 x_s and introduce the nonnegative slack coordinate with total sum1. Homogenize the cubic to degree three by the same exact procedure as the earlier bounded-cone certificate. Every one of its1,085 coefficients is nonpositive, across all31 coordinate polynomials h=2,...,32. `verify.py` retains and recomputes every coefficient, including zeros. Thus D_h(q)<=0 for every higher coordinate. The remaining coordinates vanish by degree/support.

The scalar positive-square identity gives D1(q)=-P(q)<0 for every nonzero finite nonnegative q. Therefore every actual order-ten common profile yields componentwise domination of the augmented-core ordered covers, strict at the Hamilton coordinate.

The weaker adjacent-pair bound used in `../equal_nine_ambient/PROOF.md` does not provide this order-ten coefficient certificate. That packet's reported positive coefficient concerns a larger continuous simplex, with Hamilton allowance113,400. The smaller valid domain47,250 in(2) has the certificate sign. There is no conflict, and no counterexample to actual domination was claimed at that earlier boundary.

## 3. The actual physical query in an arbitrary ambient

Let Q be any fixed tournament with distinguished root r, and let T=C3[A1,A2,A3] replace r. Each child has the same complete profile q and ten fixed physical names. Let U reverse the cyclic quotient on those sets, allowing any internal target replacements that keep profile q. All exterior names and arcs of Q remain unchanged.

Attach x->y. Within the old core use a=1_(A1),c=1_(A2), where a_v=[v->x], c_v=[v->y]. For every exterior z, give both ports the root incidence: z->r means z->x,z->y, while r->z means x->z,y->z. The augmented core is then a common module with collapsed quotient Q.

Exact module-run substitution supplies coefficients c_h>=0 and c1=H(Q)>0. For this single actual named query,

    F_(Q[T])-F_(Q[U])=sum_h c_h D_h(q)
        <=-H(Q)P(q)<0.                                 (4)

This is a strict integer inequality in every fixed strong or non-strong ambient quotient of any size. It does not require a nonzero modular residue. The details of the same physical lift are preserved independently in `../equal_small_ambient/PROOF.md`.

## 4. Scope and reproduction

This is an equal-complete-profile cyclic-quotient result at the certified child orders, not reconstruction of arbitrary30-name tournaments. The source and target must have the same named child sets, common child profile and fixed ambient context. Changes involving exterior names or different quotient replacements are outside the claim.

Run `python verify.py > normal.json` and `python -O verify.py > optimized.json`, then compare. The driver uses only the standard library and an already available C++17 compiler, recompiles the complete small extremal check in a temporary directory, and recomputes the integer-polynomial certificate without importing producer modules. Explicit exceptions and process failures remain active under optimized execution. The receipt pins all proof and source inputs.

The known extrema, positional split, classical profile cone, and nonpositive-coefficient method carry no novelty claim. This packet certifies their stated response-theoretic application.
