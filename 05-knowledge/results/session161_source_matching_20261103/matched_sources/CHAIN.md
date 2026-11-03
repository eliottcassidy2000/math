# An exact two-step chain with matched physical sources

**Scope.** The constructions below give specific strong tournaments and individual edge reversals, not a free addition rule or a global spectrum theorem. A shared intermediate is the same labelled graph, not merely another tournament with the same count.

Let Z₀ be the transitive tournament on 0<1<2<3<4. Obtain Z₁ by reversing 0–2; then obtain Z₂ by reversing 0–3. Thus the two steps each reverse precisely one physical edge, and all other labels and edges remain fixed. For j=0,1,2,

    F_Zj(r)=r⁵+2j r³;      (H(Zj),F_Zj(2))=(1+2j,32+16j).

Proof: vertex 4 is a sink. In Z₀, Z₁, Z₂ the first four vertices contain, respectively, 0, 1, 2 triangles, and no other odd cycle. The odd-cycle packing path-cover formula (S159's [cycle profiles](https://github.com/eliottcassidy2000/math/blob/7589463e54e20585aff6a2e68019ee4a4af582cd/05-knowledge/results/session159_constructions_20261103/cycle_profiles/NOTE.md)) yields r⁴+2j r² for these four vertices, and the final sink multiplies this by r. In particular both edge changes have the *identical full response* 2r³.

Let G(X) be the triangle substitution X→a→b→X and T_k the k-vertex transitive tournament. Both identities

    F_{T_k⊕X}(r)=r^k F_X(r),      H(G(X))=H(X)+F_X(2)

are given in S159's [positive profile operator](https://github.com/eliottcassidy2000/math/blob/7589463e54e20585aff6a2e68019ee4a4af582cd/05-knowledge/results/session159_constructions_20261103/profile_operator/NOTE.md). The tournament G(T_k⊕Zj) is strongly connected, including at j=0, and has k+7 vertices. The triangle wrapper and transitive padding add exactly the same named vertices and edges at each j, so the two source transitions still flip exactly 0–2 and then 0–3 within the old block. The middle graph agrees literally in both steps. The common formula is

    H(G(T_k⊕Zj)) = 1+32·2^k+j(2+16·2^k),
    common jump = 2+16·2^k,                         j=0,1,2.

In particular:

- k=0: **33→51→69**, consecutive +18 jumps on seven vertices.
- k=3: **257→387→517**, consecutive +130 jumps on ten vertices.

The response 2r³ also equals that of the earlier balanced P=1,R=3 pair in [PROFILE](https://github.com/eliottcassidy2000/math/blob/main/05-knowledge/results/apex_wedge_constructor_20261103/PROFILE.md); its raw source has H=131,F₂=1612, so the strong wrapped starting count is 131+1612·2^k with the same jump 2+16·2^k. For example at k=0 it yields 1743→1761. This is equality of full edge-response profiles under arbitrary shared contexts assembled from fixed ordinal sums and G; the source counts and graphs remain different. It gives several actual matched starting counts, not all prescribed starts.

Run `python3 chain_check.py` and `python3 -O chain_check.py`. The controls enumerate odd cycles, compute full ordered slot profiles independently by subset/path convolution, and count Hamilton paths of all three literal wrapped graphs at k=0 and k=3. The checker explicitly verifies the physical edge changes and strong connectivity. As elsewhere, finite controls do not establish unlimited chain length or global optimality.
