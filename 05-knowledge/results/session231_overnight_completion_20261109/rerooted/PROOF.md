# Rerooted grafts, right-suffix obstruction, and a smaller interval

Status: PROVED identities and structural arguments; FINITE-EXACT named-port checks. Independent audit PASS in ../rerooted_audit/AUDIT.md. No unrestricted completeness or literature-priority claim.

Inheritance: Session231 fibers/SHORT_INTERVALS.md supplies the audited deletion family L_j with fixed (h,q,A,G,S,E), h=S=E=45, and T_j=T_0+2j. Its star closure has vertices L_j,v,r and count D_j=2A+h-T_j=D_k-2j. Session229 dynamic/PROOF.md supplies the seven-coordinate ordinal composition law and definitions; inverse/PROOF.md supplies the all-gadget scalar graft formula. The sibling ../placement/PROOF.md proves same-r collapse. Hostile: a changed root does not preserve the pure-D four-state fiber automatically; physical sidecars, not scalar gcd, control amplification.

## 1. Apex identity and the other closing root

For any nonempty cut-labelled tournament K, add one apex r beating exactly the cut M. Let (h,q,A,G,T,S,E) denote K's seven-state. Then

                 H(K+r)=A-T+h-E.                    (1)

A counts the ordered covers whose second path starts in M. Subtracting T removes those whose nonempty first path ends in M. Joining remaining covers through r counts exactly paths in which r is not last, including an empty first path. The remaining h-E paths end at r. This establishes (1) without a transitivity assumption.

Let F_j=A-T_j=D_j-A-h. Root the original star closure at v. Its physical state is

                 (D_j,F_j,D_j,F_j).                (2)

Deleting v leaves the apex extension of L_j, giving F_j by (1) since E=h. A full path ending at v is any full path on L_j+r, because every other vertex except r beats v; r cannot be the endpoint of such a path, as the last selected singleton of L_j is a global sink in L_j+r. Thus E_v=d_v=F_j. For U_v, the root-starting bin either is just v, contributing F_j, or begins v,r. Removing v from the latter bin gives an ordered cover on L_j with one nonempty selected-start path, counted by A. Cases with r in the other bin and v alone were already counted; hence U_v=F_j+A+h=D_j. More directly, the v,r bin case includes the option of r alone after v, contributing h, and selected-start nonempty continuation, contributing A.

For a gadget (H,Q), the inherited scalar row gives

                 N=H D_j+(Q-H)(A+h).               (3)

This is literally the original star closure with L_j⇒B⇒v, where r beats all B. So an arbitrary v-graft is an arbitrary all-selected suffix block. Repeated same-v grafts replace H,Q by their products, and unit spacing requires every gadget transitive. The unit outputs are D_j+(A+h)(2^m-1), already selected singleton suffix padding.

## 2. Arbitrary root changes in the exterior suffix

Fix the entire protected L_j, including its trailing selected singleton. Start with exterior suffix R={v}, cut empty on v, and distinguished apex r. Allow arbitrary nonempty gadget grafts at r OR at any vertex of R, including newly added vertices, with arbitrary changes among these root choices. No graft may target a vertex of L_j.

At every stage the physical tournament minus r is exactly L_j⇒R, with an arbitrary induced tournament R and arbitrary cut on R. This follows because every protected vertex beats every new gadget vertex, regardless of which allowed root is selected. Rooting at r appends a block to R; rooting within R changes R internally. The cut on original v remains unselected.

The seven-state composition law shows that all coordinates of L_j⇒R except T are fixed in j, while its T changes by 2H(R) per move. By (1), the exact output progression is therefore

                 N_j=N_0-2H(R)j.                  (4)

This is an actual progression with all k+1 witnesses, not a gcd statement. For k≥1, consecutive-odd spacing is equivalent to H(R)=1, hence R transitive. The standard characterization is that a nontransitive tournament has a nontrivial strong component; a spanning directed cycle in it supplies at least three Hamilton paths, and ordinal SCC counts multiply.

If R is transitive, remove terminal cut-selected vertices: each is a global sink in the full tournament and does not affect the path count. At least one unselected vertex remains, since original v is unselected. The final unselected vertex serves as the closing v, and preceding R vertices form an ordinary binary singleton suffix. Thus even arbitrary root changes and arbitrary gadgets in this entire exterior-suffix region add no unit-window placements beyond singleton padding of the original L_j. The argument does NOT cover grafts at the protected leading/trailing ports or core vertices.

## 3. A genuinely smaller physical interval from deletion

Although (2) does not give a new unit-preserving graft, its deletion coordinate gives a smaller exact interval. Delete v from the original graph, then delete L_j's trailing selected singleton. The latter is a global sink after v is removed. The resulting graph has exactly 12k+11 vertices and count

                  F_k-2j,  0≤j≤k,                 (5)
                  F_k=D_k-A_k-45,
                  A_k=7516·4096^k+343260 J_k+45,
                  J_k=(4096^k-1)/4095.

Thus k+1 consecutive odds occur, with F_0=2303 and F_1=9274153, compared with original D_0=9909 and D_1=40403039. This baseline cannot be reached by adding vertices to the corresponding original source, since insertion injects its Hamilton paths into those of every supertournament. This is a smaller-family placement, not an all-target selector.

The reduced graphs are strong. The apex beats the initial selected singleton, which reaches every later deletion vertex. If j<k, the final Y block ends at an unselected singleton, reached by every earlier vertex and beating r. If j=k, the final deletion block is the strongly connected eight-vertex C, containing unselected vertices, so every vertex of C reaches r. Therefore every vertex reaches r, and r reaches every vertex, in both cases.

Before removing the trailing global sink, the rooted state at r is especially simple:

                 (F_k-2j,0,A_k+45,45).             (6)

Deletion has h=45. No full path ends at r because the trailing selected vertex is a global sink. Root-start covers consist of r alone with any deletion Hamilton path, plus the selected-start covers counted by A. Hence U=A+h.

A single arbitrary gadget at this root gives

                 N_j=H(F_k-2j)+(Q-H)A_k+45H.       (7)

For a transitive m-vertex gadget, m≥1, this is F_k-2j+(2^m-1)A_k+45. At m=1 it exactly recovers D_k-2j, and larger m match the original r-graft padding family. Thus the smaller m=0 source is the new placement; its first graft returns to the old family. The next section gives a further certified deletion.

## 4. Bare-core intervals on 12k+9 vertices

A further reduction, suggested by the parent research thread and checked independently here, removes the initial s1,s0 as well. Define B_j=X^j⇒C⇒Y^(k−j), and add apex r beating its cut. These are strong tournaments on 12k+9 vertices with counts

                 Z_k−2j,  0≤j≤k,
                 Z_k=F_k−2(q_B−G_B),
                 q_B=1440·4096^k,
                 G_B=577+2383200 J_k,
                 A_B=878·4096^k+171630 J_k.         (8)

Equivalently, Z_k=(13617152·4096^k+1602510k+717259)/24843. Their first upper endpoints are 577, 2245227, 9196067477, and 37667091739327.

Proof: the cut state of prefix s1⇒s0 is (h,q,A,G,T,S,E)=(1,4,2,1,0,1,0). Thus prefixing B changes A to 2q_B+A_B and T to 2G_B+T_B, leaving h and E unchanged. By (1), its apex count differs from the bare count by 2(q_B−G_B). The prefix-plus-B apex count is F_k−2j, since it is the reduced graph of Section3. To see q_B and G_B are independent of j, strip the prefix and trailing selected singleton from the known full L interface: q_L=8q_B and G_L=4q_B+4G_B+45. Similarly A_L=4q_B+2A_B+45. Substitution gives the formulas in (8), proving all-order unit spacing without a new unverified carry identity.

For strongness, if j=0 then the initial core C is strong and has a nonempty cut, so r reaches C and every later vertex. If j>0 the first X bit is selected, so r reaches that initial vertex and every later vertex. If j<k the last Y bit is unselected, so every vertex reaches it and then r. If j=k the final C is strong and contains unselected vertices, so every vertex reaches r through C. Strongness of C is also verified by literal directed reachability.

The full physical rooted state at r is

                 (Z_k−2j,45−e_j,A_B+45,45),
                 e_j=0 if j<k, and e_k=14.         (9)

Indeed e_j is the cut-end count of B_j, determined by its last block; all preceding Hamilton multipliers are1. Formula (9) follows by the apex endpoint and root-start-cover arguments of Section3. Hence the first k witnesses form a pure-D rooted fiber, but the last witness changes E by−14. An arbitrary common graft across the full bare interval cannot be assumed to preserve its width. This is an exact new lower baseline family, not an all-target coverage theorem.

## 5. Hostile finite named-port controls

The exact ordered-block transfer in explore.py tracks actual Hamilton-path counts and physical four-states (D,E,U,d) for all named ports. It is independent of the seven-coordinate composition formula. For a Hamilton path, remove apex r and assign all ordinary vertices to the two ordered gaps around r. For U at a≠r, hold both a and r out: r lies either after a in the root-starting path, or in the other path. In either case there are three ordered gaps. Each ordinal block distributes its vertices among gaps; literal subset endpoint DP inside the block gives weights. An empty/last-endpoint-good/last-endpoint-bad flag for each gap enforces exactly the adjoining special-vertex arcs. This computes at most three-bin block covers, never 2^(12k+13) subsets of the full graph.

At k=0 all 13 full rooted states agree with independent full-graph subset DP. At k=1 the two physical witnesses give these differences:

- r: Δ(D,E,U,d)=(-2,0,0,0).
- v: (-2,-2,-2,-2).
- trailing selected port: (-2,0,10,12), so ΔJ=8, ΔV=-2 and ΔN=12H-2Q, a multiple of four. It cannot equal ±2.
- leading selected, next unselected, and all eight named core vertices: both ΔJ and ΔV are negative, with ΔJ<-2. Since H≥1 and Q≥2H, every nonempty gadget gives ΔN<-2.

Thus no other named port preserves even the k=1 consecutive pair under any common nonempty gadget. This is FINITE-EXACT over the two sources, but universal over gadgets by the scalar inequality. It is deliberately NOT generalized to all k: at k=2 some leading-port differences change sign, and the trailing port has (-2,0,-4,-2) on the first move and (-2,0,10,12) on the last. The stored k1.json and k2.json retain every full sidecar.

## 6. Reproduction and remaining escape routes

Run python3 -I -B verify.py and python3 -I -B -O verify.py, and compare with VERIFY.json and VERIFY.optimized.json. These verify (5),(6),(8),(9), graph strongness and sink removal through k=3 using exact transfer, then all k=0 original rooted states and the deleted-v state by independent literal subset DP. explore.py k emits all named-port sidecars. check.py verifies the original k=0 states alone. Explicit failures are used, not optimizable assertions.

Open: grafts into the protected block, sequences involving such grafts, interleaving dilated progressions, and an all-target baseline selector for either the original D or smaller F family. The suffix obstruction and finite named-port failures do not rule those out.


## Independent audit qualification

The standalone audit confirms the theorem and every named-port difference using full-subset DP, including both 25-vertex source graphs. A terminal bare-apex graft breaks the indexed progression, but not always its unordered consecutive set: at k=6 with H=1, all seven values are still consecutive. For k>=7 the final value duplicates an interior one. This is detailed in ../rerooted_audit/AUDIT.md and is not an all-target repair.
