# Late incoming work: all-size edge roles and aligned ownership

Incoming commit d45865c3323af6d1b3ce44525d6d4789b22778d7 arrived during final publication checks and is preserved as an ancestor. Its [frontier](../metagraph_owned_boundary_frontier_20261023.md) incorporates the earlier Session40 checkpoints in turn. The results below are credited to that incoming packet rather than presented as newly derived here.

## 1. Coarse square-sum single-edge status is now closed at all sizes

The incoming mistakes-ledger repair notices an elementary operation missing from its draft: chord insertion into a Hamiltonian cycle. In any finite simple undirected graph with a Hamiltonian cycle and at least three vertices, every edge is included in some Hamiltonian path and omitted in another.

For a cycle edge, cutting it omits it, while cutting a different edge includes it. For a chord uv, orient the cycle, remove u--succ(u) and v--succ(v), and add uv. The two remaining paths join into one spanning path containing the chord. Any path obtained by cutting the original cycle omits the chord.

This proof is independent of enumeration. Our added checker constructs both witnesses for every cycle edge and possible chord at orders3..64, and includes a fixed-endpoint hostile. It does not re-prove the external square-sum cycle theorem.

**CITED input:** Robert Gerbicz's own 21January2018 comment in [OEIS A090461](https://oeis.org/A090461) states that square-sum graphs have Hamiltonian cycles for every n≥32. The comment was directly re-read during this integration; the full gluing construction is not re-audited here.

**Consequence:** every fixed positive distinct a,b with square a+b is optional in Q_n at EVERY n≥max(32,a,b). Combined with the exact lower-order census, the coarse single-edge role is classified for all n. The finite-only scope in our original square-sum packet remains valid, but this later incoming argument supplies the stronger extension. It does not settle exact usage frequencies, prescribed endpoints, or simultaneous requirements on several edges. Those are the next real questions.

## 2. Exact aligned source-target refinement is stronger than source-only data

The incoming [aligned-port proof](../aligned_port_refinement_20261023/PROOF.md) retains actual alignment of adjacency matrices A and B and the actual path P. Its identity is

A+B−(J−I)=F_P−F_P^T.

Joint equitable refinement into the SAME owned cells is discrete within floor(n/2) rounds: divergence first isolates the endpoints, then singleton incidence peels the path inward from both ends. The only nonnegative doubly stochastic matrix commuting with both A and B is the identity. These are all-size proofs, but neither finds B or the alignment from compressed source data.

The separate agreement matrix D=A∘B still needs endpoint in/out-degree guards and a global no-cycle test. With those guards, nilpotence, det(I−D)=1 and an(n−1)-step start-to-end walk are equivalent native-path tests. A path plus detached directed cycles is the hostile if that global condition is omitted.

This complements, rather than supersedes, the Session40 endpoint-marked class decoder and the five-vertex source-only hostile. Complete aligned source AND target information is a strictly stronger input.

## 3. The actual order31 source and a genuine Sierpinski interface

The incoming [owned-endpoint packet](../endpoint_overlap_ports_20261023/README.md) treats the particular D(D(P7)) source, not all order31 tournaments. Among its930 ordered pairs used as marks,672 cases close intrinsically,252 need one extra actual owner, and six need two. Not every tested ordered pair is asserted to be Hamiltonian-feasible. The separate actual-port control has932 distinct Hamiltonian paths on510 endpoint pairs.

Shared-row equality resolves the difficult seven-cells; affine negative-loop fixed points then resolve some triple systems. Genuine residual C3 or order21 automorphisms prove why actual extra owners, not merely another unmarked statistic, are needed in the remaining cases.

The incoming [boundary-recursion packet](../hopf_sierpinski_port_recovery_20261023/REPORT.md) gives a genuine three-boundary affine-F3 section carrier with64states: empty,27singletons, and36three-point affine graphs. Its triangular gluing operation is closed and all64^3 compositions are checked.

There is a decisive applicability boundary: every one of the420 reached eight-triple systems has an empty global parallel-section relation. Some nevertheless have residual C3 symmetry, while others are intrinsically resolved. The centralizer of the holonomy action distinguishes these possibilities; the empty section relation does not. Therefore this compact compatibility interface cannot be substituted for full attachment or wiring data on that source. Existence of a common fixed owner and the group of commuting symmetries are different predicates.

The same incoming work lifts the actual path's Pascal clock to modulus2^s using our already-published binomial valuation: for n≥2 and s≥1, I+F has exact order2^(s−1+ceil(log2 n)). This clock, the64state boundary interface, the native electrical mod8 pairing, and the Hopf Chern question remain separate typed objects.

## 4. What changes in the final frontier

Unrestricted square-sum single-edge roles are no longer an open target. Work should instead retain prescribed endpoints and simultaneous edge incidence. On tournaments, the next target is economical querying of actual target incidence into owned source cells, without assuming the unknown target is already supplied. The new regular-block coefficient formulas and displacement-necklace repair provide additional exact carriers, but do not erase that input distinction.
