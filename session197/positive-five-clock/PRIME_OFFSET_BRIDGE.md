# Prime-offset descent in the arbitrary-module quotient

This is a separate structural consequence of the audited arbitrary-module identity. It is not a proof of unrestricted prime-spectrum coverage.

## Theorem

Fix an odd prime q>=5. In the seven-block quotient of `PROOF.md`, with arbitrary nonempty P,U,R,V and T=Tr_m, an eventual pure ray

    H_m=alpha*3^m+q,  alpha>0,

occurs if and only if:

1. P and R are both Tr2;
2. one of U,V has Hamilton count q;
3. the other is transitive, of some positive order v.

Writing W for the count-q module, the amplitude is exactly

    alpha=18 B_W(2) 2^v.

Thus such a prime-offset ray exists precisely when q itself is already realizable as a tournament Hamilton count. Its construction cannot produce a prime offset without an input module realizing that prime.

## Proof

Set delta_X=B_X(2)-3H(X). The audited cut bound gives B_X(2)>=(|X|+1)H(X). Together with the explicit one- and two-vertex cases this shows that delta_X is never zero, that its only negative value is-1 at K1, and that delta_X=1 occurs exactly at Tr2.

The all-order coefficient identity is H_m=A3^m+B2^m+C, with A>0 and

    C=h_U h_V delta_P delta_R.

An eventual pure ray forces B=0 and C=q. Since h_U,h_V are positive integers, delta_P and delta_R have the same sign.

If both deltas are negative, they are both-1, so P=R=K1. But then

    B=2(b_V h_U+b_U h_V)>0,

a contradiction.

If both are positive, the primality of q forces exactly one of h_U,h_V,delta_P,delta_R to equal q and the other three to equal1.

If h_U=q or h_V=q, then both deltas are1. Hence P=R=Tr2, the other Hamilton count is1, and that module is transitive. Both binary terms vanish, and the amplitude is exactly the stated value.

If instead delta_P=q, then R=Tr2 and h_U=h_V=1. The coefficient B reduces to

    3 b_V (4h_P-b_P).

Vanishing forces b_P=4h_P and hence delta_P=h_P=q. The cut bound now gives |P|<=3. But every tournament on at most three vertices has H at most3, contradicting q>=5. The case delta_R=q is symmetric.

This proves necessity. Conversely, the three stated conditions give B=0,C=q and the stated positive amplitude directly. The eventual identity is then exact for every positive m.

The assertion H=1 iff transitive follows from ordered SCC factorization and the fact that every nontrivial strong component has at least three Hamilton paths, or from directed-triangle insertion.

## Exact q=5 specialization

Every tournament W with H(W)=5 has one strong four-vertex H5 component, with any number s of singleton SCCs before and after it. Indeed, primality leaves only one nontrivial SCC; Busch's bound excludes a strong core of order at least5, the cyclic triangle has count3, and the strong four-tournament is unique up to isomorphism.

That four-vertex core has bin polynomial x^4+4x^2. Bin polynomials multiply under ordered joins, since each labelled bin uniquely concatenates its path pieces from successive SCCs. Consequently

    B_W(2)=32*2^s.

Combining with the theorem, the complete set of amplitudes of positive-constant5 pure rays in this quotient is

    1152*2^k,  k>=0.

Every k occurs by placing k singleton SCCs in W and taking the other module K1. The main packet constructs the least amplitude in this **specified quotient and fixed constant5 class**, namely1152. It does not assert a global least-amplitude theorem for arbitrary tournament mechanisms.

The strong-order bound and four-vertex classification are standard inputs already checked in the Session197 low-count packets. The separate verifier exhausts the small physical module profiles, all labelled H5 four-vertex cores and their isomorphisms, and the coefficient factorization. Those finite checks support the argument; the prime-offset equivalence itself is all-order.
