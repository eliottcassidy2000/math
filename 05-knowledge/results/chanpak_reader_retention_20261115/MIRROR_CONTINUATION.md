# The mirror unit switch carries a complete slot difference

The [five-site mirror switch](../chanpak_tangent_survey_20261115/MIRROR_NORMS.md) was known to change the ordinary count by 2 for a transitive exterior, or generally by twice the square of the exterior's count. A background follow-up suggested the full polynomial difference. We prove it here and check literal covers independently. It puts a sharp bound on what an unmarked reader can do to this exact unit move.

Let $L$ be any painted tournament (the paint marks outgoing neighbors of an external apex), and use the ordinal spine $L\Rightarrow C\Rightarrow R$ where $R$ is the marked-complement converse mirror of $L$. In the transitive five-site center use $(u,1-v,0,v,u)$ with $u,v\in\{0,1\}$. Fix $v$, and call the one-apex tournaments $U_0,U_1$ according to $u$. Write $\Pi_z(X)$ for the slot polynomial of $X$, evaluated as the ordered $z$-slot path count at positive integers. Then

$$\Pi_z(U_1)-\Pi_z(U_0)=2z^4\Pi_z(L)^2.$$

**Proof.** In the odd-cycle expansion any changed cycle contains the apex; no other cycle in that packing can intersect it. The switch only changes arcs from the apex to sites 1 and 5, so a changed cycle must use exactly one of them as a first or last spine endpoint. A cycle using both of them as endpoints never works: their marks agree. Reflect the entire spine and reverse the traversal to pair a gained cycle with a lost cycle whenever its other endpoint is not site 3. All arcs on the spine change to their converses under reflection, and marks away from site 3 and the two switched endpoints complement under reflection. Thus the reflected object is present at the opposite switch setting. Disjoint cycles and unused vertices correspond too, so their full packing contributions cancel, even if an internal vertex of the changed cycle is site 3.

Site 3 is unmarked in both settings. A changed cycle with this other endpoint must run from the apex to site 1, then to site 3 and back to the apex. The only possible intermediate spine vertex is site 2, and including it makes an even cycle. The remaining gain is exactly one triangle, with weight $2z$; the unused sites 2, 4, 5 contribute $z^3$, and no cycle outside the apex can leave its ordinal component. The left and right profiles multiply and are equal by converse. This proves the identity.

## Exact consequence under every unmarked pointed reader

Set $V=L\Rightarrow L^*\Rightarrow TT_4$, with $L^*$ the converse. Its slot polynomial is $z^4\Pi_z(L)^2$. The visit expansion for any tournament $S$ with one marked substitution slot is linear in the nonempty cover coordinates, which are the binomial-basis coordinates of the slot polynomial. Therefore

$$ H(S[U_1])-H(S[U_0])=2 H(S[V]). $$

This is an equality of integer **differences**, not an assertion that a formal signed object can be repeated as a tournament. The step stays $+2$ exactly when both $L$ and $S$ are transitive: then $S[V]$ is transitive, while if either contains a triangle, so does $S[V]$, and the odd-cycle formula gives strictly more than one path.

For a strong nonsingleton $S$ the marked vertex is on a triangle. One way to see this is that absence of such a triangle gives the ordinal separation of the in-neighbors, marked vertex and out-neighbors in the matching note. Thus the reader has at least one two-visit template and at least three one-visit templates, the latter because $H(S)\ge3$. Its step is at least $2(3H(V)+P_2(V))=2(H(V)+\Pi_2(V))$, i.e.

$$H(S[U_1])-H(S[U_0])\ \ge\ 2H(L)^2+32\Pi_2(L)^2.$$

This is sharp for $S=C_3$, whose one- and two-visit template counts are 3 and 1. At a transitive binary envelope of length $\ell$ the sharp new step is $2+32\cdot4^\ell$. So the old unit move cannot simply be fed through an unmarked strong reader. A selective reader is outside this obstruction.

Run python verify_mirror_continuation.py. The verifier uses literal adjacency and subset convolution, including nontransitive envelopes. The full-profile linearity is the inherited unmarked substitution argument; the explicit polynomial difference and sharp unit-step obstruction are the new parts.
