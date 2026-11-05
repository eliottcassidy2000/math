# Short physical source certificates for low-cycle scalar membership

Session198. **PROVED, using the named complete source reductions:** for each fixed c in {1,2,3}, deciding whether a positive binary integer N is H(T) for some strong tournament with exactly c Hamilton cycles belongs to NP. This does not prove polynomial-time membership, cofinite coverage, or universal odd coverage.

Let F(D,E,U,d)=D−E+U, as in `PROOF.md`. For every physical source D<=U and E>=0, so F<=2U. Formula(1) of that proof gives, under every permitted source move,

    F_new−F_old>=U_old>=F_old/2.

Starting from the singleton with F=1, any source recipe of length t therefore satisfies

    (3/2)^t <= F <= 2U.                            (1)

In particular every source with U<=B has a physical recipe of length at most floor(log_(3/2)(2B)). This applies to every forward recipe, not only a chosen shortest recipe.

For a gadget step M_L, D_new=gD−(L mod2)E+LU>=LU since D>=E and g>=1. Also U increases at every step and D_new<=U_new. Consequently L<=B in any recipe whose final U<=B. Each index L has O(log B) binary digits. A source recipe can therefore be encoded in O(log^2 B) bits using O(log B) operation tags and binary gadget indices.

The C1 inverse uses a source with F=N, hence U<N. The complete C2 reduction uses two sources with U<=floor((N−1)/6), apart from its explicit H15 exception. The complete C3 reduction uses two sources with U<=floor(N/3), or the sharper floor(2N/15) bound proved in the separately reviewed source-coverage packet. Thus B<=N is sufficient for all three languages. C1 needs one source recipe; C2/C3 need two and one of finitely many template tags. Singleton boundaries and the C2 H15 exception have constant descriptions.

A verifier checks the length and index bounds, starts from(1,1,1,1), computes each physical forward source step, checks the final source domains, and evaluates the appropriate exact scalar form. It accepts exactly when the result is N. The physical source grammar and gluing theorems ensure that every accepted recipe describes an actual graph of the claimed cycle count. Completeness follows from the proved C1,C2,C3 normal forms. No unrestricted formal profile is accepted without a source recipe.

For bit complexity, the gadget constants g(L),Q(L),B(L) are computed by their binary recurrences in O(log L) stages. Their bit lengths are O(log L): for example g(L)<=2L, B(L)<2L (with L=1 handled separately), and Q(L)<2B(L), all directly inductive. A verifier may reject any intermediate profile with U>N; genuine witnesses never require it. Products before this check have O(log N) bits up to a fixed constant factor. There are O(log N) source moves, followed by fixed-size bilinear arithmetic, so verification is polynomial in the binary input length. This proves NP membership.

The actual graph can also be reconstructed with O(log^2 N) vertices: every O adds one vertex, every M_L adds1+bit_length(L−1) vertices, and the final template contributes only a constant number. The NP verifier does not count Hamilton paths of an arbitrary graph of that size; it verifies the arithmetic construction certificate whose exact counting identities have been proved. Confusing those two verification tasks would be a scope error.
