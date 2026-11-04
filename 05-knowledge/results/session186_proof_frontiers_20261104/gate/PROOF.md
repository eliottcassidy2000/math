# Three-state clock, named twin separation, and a scalar-context hostile

**PROVED** the all-order statements below, with **FINITE-EXACT** checks separately described. These statements connect Session185's three-state gate to the named two-port observer. They do not prove universal reconstruction, remove the prime twelve-vertex residual, or exhibit any full-response collision.

## Conventions and inherited inputs

Outgoing-neighbor masks use fixed names. The named two-port response F_T(a,c) adjoins x->y with bit a_u=[u->x], c_u=[u->y]. C_uv(T) counts ordered spanning covers by two nonempty directed paths, first ending at u, second starting at v, including absent-edge gradient directions. H counts directed Hamilton paths. q_k counts ordered spanning covers by k nonempty directed paths.

The inputs are the explicit law in [Session185's gate](../../session185_proof_frontiers_20261104/alternative/PROOF.md); C2 in [the eight-core algebra](../../two_port_eight_core_20261103/algebra/PROOF.md); and T4 in [the twelve-principles proof](../../two_port_twelve_principles_20261103/PROOF.md). C2 gives a named ordered-cover sign recoverable linearly from the full two-port table. T4 proves F_T-F_(T^op) is cubic with vanishing pure quadratic parts and the exact 1+2n+n^3 battery. This cubic assertion concerns **whole converse**, not arbitrary internal module reversal.

Use Session185's seeds

 Q0=(16,5,33,39,46,3), Q1=(16,5,49,39,42,3),
 Q2=(18,4,49,39,42,3).

Replace vertex5 by the transitive chain 5->6->...->(r+4), obtaining T_j(r). All are strong and successive states differ by a single old edge. The inherited exact law is

 H(T_j(r))=4·3^r+9·2^r-3+2j,
 H(Qj[M])=(27+2j)H(M)+21q_2(M)+4q_3(M).

## G1. The gate clock supplies an exact named-response separator

**Theorem.** For r>=3, put a=5,b=6,w=7 and

 L=C_bw-C_aw+C_wa-C_wb.

Then

 L(T_j(r))=8·3^(r-2)+18·2^(r-2)-6+4j >0.             (G1)

Both the whole converse T_j(r)^op and the internal module converse Qj[TT_r^op] have the negative of this value. Thus their contrast differences from T_j(r) are twice the positive quantity in (G1). Neither pair has equal or positively proportional named two-port tables. In particular the whole-converse cubic difference is nonzero; the cubic battery must detect it.

**Proof.** The consecutive pair {a,b} is homogeneous, with a->b. Contract it to z, obtaining Qj[TT_(r-1)]. C2's swap-and-contract bijection gives

 L=C_zw+C_wz.

In this contraction {z,w} is again a homogeneous consecutive pair. Joining the two paths at their named boundary and contracting that block is a bijection from each of the two cover sets to Hamilton paths of its contraction. The actual mutual arc is irrelevant because the boundary is removed in the cover. Hence both C_zw and C_wz equal H(Qj[TT_(r-2)]). Apply the inherited gate law to obtain (G1).

Reversing the entire tournament transposes C and negates L. Reversing only the module changes a->b to b->a, so C2 gives the negative sum of the contracted cover counts. The remaining consecutive pair again contracts, now to Qj[TT_(r-2)^op]. Every q_k is converse-invariant by reversing each path while keeping the ordered slots, so the uniform-run expansion gives H(Qj[TT_(r-2)^op])=H(Qj[TT_(r-2)]). This proves the exact negative value in the internal case as well. C is linearly recoverable from the full named two-port table, and positive scaling preserves the sign of L, proving separation. The last cubic assertion follows only for the whole-converse case from T4. QED.

For r=3 the three L values are54,58,62; for r=4 they are138,142,146. The exponentially growing H clock is therefore an explicit growing response separator. It is not evidence of a cancellation of the three-path response tensor.

## G2. Uniform port attachments are intrinsically blind to module converse

**Theorem.** For every nonempty M, every quotient Q, and every number of new ports whose incidences are uniform across M, replacing M by M^op leaves every resulting Hamilton count unchanged.

**Proof.** Reversal of each path gives q_k(M)=q_k(M^op), retaining the k ordered slots. Absorb all uniformly attached ports into the exterior quotient. Its occupied-run expansion is sum_k c_k q_k(M), with the same integer skeleton coefficients on both sides. QED.

Thus the full response's nonuniform attachments are essential. For TT_r, G1 certifies the information lost when those attachments are removed. The scalar formula and its3^r/2^r clocks alone never certify equality of the full response.

## G3. Every scalar substitution context can miss the existence of a pair

Consider the actual five-name modules

 A=(4,1,2,7,15), B=(10,0,3,6,15).

A is the ordered join of singleton4, singleton3, and the directed triangle on0,1,2. B is the ordered join of singleton4, the directed triangle on0,3,2, and singleton1. A has exactly the homogeneous pair {3,4}; B has no homogeneous pair. They have the same complete profile

 (q_1,q_2,q_3,q_4,q_5)=(3,42,162,240,120).           (G3)

Consequently H(Q[A])=H(Q[B]) for **every** quotient and root, and for every collection of exterior ports with uniform incidences to the module. Nevertheless their actual named two-port values at(a,c)=(0,1) are15 and9. This is an explicit scalar-context collision, **not a full-response collision**. Since the seeds Qj are strong, Qj[A] and Qj[B] are also strong scalar-context controls of order10, with equal H=1611+6j. To lift the separating query, give each port the same incidence to every quotient exterior vertex as the replaced root, retaining(a,c)=(0,1) on the module. The actual extended module has seven vertices. Applying the gate law, or direct subset DP, gives ambient query values21015+30j and17793+18j. Thus these three strong ambient controls are explicitly separated too.

**All-context proof of profile equality.** For an integer t>=0, let P_M(t) count assignments of vertices to t named colors, with a Hamilton path on each nonempty color class, and the empty class contributing1. Then

 P_M(t)=sum_k binom(t,k)q_k(M).

For an ordered join X->Y, every induced color-class Hamilton path is the concatenation of its X and Y paths, so P_(X->Y)(t)=P_X(t)P_Y(t). This product is insensitive to the order of the components. A and B consist of one directed triangle and two singleton components, hence have identical P, and therefore identical q_k by independence of the binomial polynomials. More explicitly,

 P_singleton(t)=t, P_C3(t)=t^3+2t,
 P_A(t)=P_B(t)=t^5+2t^3.

Expanding t^m=sum_k k!S(m,k)binom(t,k) gives (G3). The module law proves the all-context claim. Checking the ten pairs directly gives the stated pair sets; literal Hamilton DP on the two seven-vertex port extensions gives15 and9. QED.

The pair-free property here is only absence of a two-vertex module; B still contains its triangle as a proper module. No prime residual counterexample is claimed.

## Reproduction

`python verify.py` and `python -O verify.py` compute every ordered C entry from independently computed start/end subset DP tables. They check all three states, r=3,...,7, the explicit law, H equality, whole converse, and internal module converse:15 positive controls and30 reversed controls, without using the C2 identity in computing C.

`python hostile.py` and its optimized replay independently enumerate color assignments for every q_k of A and B, inspect every pair, and count the literal two-port extensions, including all three strong ambient substitutions. `python query.py` and its optimized replay give explicit whole-converse two-port queries for each state at r=1,...,4. In every case(a,c)=(1,0) already separates; these are finite controls, not an extrapolated formula. All gates use explicit exceptions rather than removable asserts. Normal and optimized JSON outputs agree.

The exploratory file explore.py is not used as proof or validation. No install, external literature novelty claim, or complete larger-class enumeration is needed.
