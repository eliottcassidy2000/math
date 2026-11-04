# Arbitrary original-module families at an equal-profile cyclic root

**PROVED, using the companion decomposition theorem:** Let T=C3[A1,A2,A3] on fixed physical names, where the three nonempty children have identical complete ordered-path-cover profiles. Every original-module XOR family containing the whole vertex set changes the full named two-port response. The family may contain arbitrarily many nested or crossing proper modules. Two specified actual queries suffice to distinguish the original T from any such result: (0,0) and (V1,V2).

This extends the equal-profile two-level reversal cube to arbitrary depth and overlap inside its children. It does not prove arbitrary full-response reconstruction, injectivity between every two states of an arbitrary original-module orbit, or the arbitrary strongly connected ambient lift.

## 1. Proper modules cannot cross the cyclic children

Every proper module M of C3[A1,A2,A3] is contained in one child. If M met Ai and Aj but omitted a vertex of the third child Ak, that vertex would see the two occupied children in opposite directions, contradicting modularity. Thus M would have to contain all of Ak. Any omitted vertex of Ai then sees Aj and Ak in opposite directions, and similarly for Aj. Hence meeting two children forces M to be the whole tournament.

Consequently, in a family containing the whole root exactly once after parity reduction, only that root operation toggles the quotient triangle. All proper selected sets act inside children. The resulting U therefore has the form

    U=C3^op[B1,B2,B3],

on the same physical child sets. No proper operation changes the inter-child quotient again. These are simultaneous reversals of modules of the ORIGINAL T; no legality of an arbitrary crossing sequential ordering is assumed.

## 2. Two-query separation

The companion `module_decomposition/PROOF.md`, Theorem 1 and Section 5, proves for any original-module XOR that H(U)>=H(T), with equality implying equality of every ordered-cover profile on every physical induced subset.

If H(U)>H(T), the actual query (a,c)=(0,0) separates because F(0,0)=H. Otherwise apply hereditary equality to each fixed child set Vi:

    q_h(Bi)=q_h(Ai)=q_h       for every i,h.

The physical block-uniform query a=V1,c=V2 depends only on these three profiles and the quotient orientation, by the maximal-child-run bijection. The unconditional equal-profile theorem in `covered_response/PROOF.md` therefore gives

    F_U(V1,V2)-F_T(V1,V2)=P(q)>0.

The strict sign is the positive-square identity in `covered_positivity/PROOF.md`. It does not require identical interior orientations, isomorphic children, transitivity, log-concavity, or a bound on descendant depth. These two cases exhaust all families, proving the theorem.

For a laminar family the first case cannot occur, because current modularity is preserved through every selected reversal and all scalar counts are invariant. Thus the single root-uniform query always separates such a root-containing family, independently of every descendant operation.

## 3. Scope and actual controls

The conclusion concerns the entire cyclic core T, or an ordered earlier/later SCC padding to which exact localization applies. For a selected cyclic module in an arbitrary strongly connected outer context, the scalar-increasing case is already handled by the decomposition theorem, but the scalar-equality case still needs an ambient lift of the root-uniform separator. The separate `covered_ambient` reduction records that open sign problem; it is not assumed here.

The checker enumerates every root-containing original-module family of C3[TT3,TT3,TT3], including crossing adjacent intervals inside a child. It constructs each XOR target on the same nine names and checks the specified physical queries by literal predecessor subset DP. It also preserves selected twelve-name deep laminar controls. These finite checks are not the proof of the unbounded statement, which is Sections 1–2 with the explicitly cited decomposition and positive-square theorems.

## 4. Complementary corollary: all transitive child sizes at least two

The same argument covers T=C3[TT_a,TT_b,TT_c] for every a,b,c>=2, whether or not their profiles agree. Equal sizes use the theorem above. For unequal sizes, choose a child i whose other two cyclicly ordered child sizes differ. The companion `covered_quotient/PROOF.md`, Sections 2–4, proves the strictly signed aggregate start-count contrast

    S_i(T)-S_i(U)
      =sum_(k>=1) q_i(k)[q_j(k)q_l(k-1)-q_j(k-1)q_l(k)] !=0

whenever U has the converse quotient and the same individual child profiles as T. That identity depends only on child profiles, not on their physical interior orders. Hereditary profile equality in the scalar-equal case therefore gives exactly the required hypothesis, even after arbitrarily many nested or crossing descendant reversals.

Its actual endpoint battery uses (e_u,0) and (0,e_u) for every physical u in child i. The identity

    2S_u=F(e_u,0)-2F(0,e_u)+F(0,0)

shows that at least one of these 2|Ai| queries separates in the scalar-equal case. Add (0,0) for the scalar-increasing case. Thus a fixed battery of at most 2|Ai|+1 actual queries separates the original T from every root-containing original-module XOR family. No minimal-query claim is made.

The checker additionally enumerates all 128 root-containing families on C3[TT2,TT3,TT3], retaining every selected query value and the complete endpoint battery in scalar-equal cases. This corollary inherits the same whole-core and ordered-SCC ambient restriction; it does not settle arbitrary unequal nontransitive child profiles.
