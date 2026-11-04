# Vertex-transitive extreme children survive original-module families

## 1. Every module of a vertex-transitive tournament is vertex-transitive

All tournaments and groups here are finite. A module is a nonempty set seen uniformly by every outside name. Two modules overlap when they intersect but neither contains the other. A strong module overlaps no module.

**Lemma 1.** Every module of a regular tournament is regular, and every module is strong.

For a module M, every vertex of M has the same number of out-neighbors outside M. Regularity of the ambient tournament therefore makes its internal outdegree constant, proving the first assertion.

Suppose modules M,N overlap. Put A=M\N, B=M intersection N, C=N\M; all three are nonempty. For a in A, b in B, c in C, the module conditions imply either A->B->C with A->C, or the complete reverse order. Indeed a sees all of N identically, and c sees all of M identically. These directions are uniform across each of the three sets. In particular M is an ordered union A->B or B->A. Such a union cannot be regular: the average internal outdegree in its first side is (|first|-1)/2+|second|, whereas the second side's average is (|second|-1)/2. Their difference is |M|/2>0. This contradicts the first assertion. Thus no two modules overlap.

**Lemma 2.** Every module M of a vertex-transitive tournament D is vertex-transitive on its induced physical names.

Vertex transitivity implies regularity, so M is strong by Lemma 1. For any u,v in M, choose an automorphism g of D with g(u)=v. The image g(M) is a module of the same cardinality as M and intersects M at v. Strongness forces containment, and equal finite cardinalities force g(M)=M. Restricting g gives an automorphism of D[M] taking u to v. This includes singleton and whole modules.

No claim that every regular tournament is vertex-transitive is used. Finiteness and regularity have distinct roles in the argument.

## 2. Marked-profile preservation under balanced module converse

For a tournament D and physical name u, let E_u(h),S_u(h) count ordered spanning covers by h nonempty directed paths whose first slot respectively ends or starts at u. Let q_h(D) count all such ordered covers. A tournament is endpoint-balanced if E_u(h)=S_u(h) for every u,h. Vertex transitivity gives E_u(h)=S_u(h)=q_h(D)/|D|.

**Lemma 3.** If M is an endpoint-balanced module of an arbitrary tournament D, conversing only M preserves q_h(D), E_u(h) and S_u(h) for every physical u and every h.

Decompose each ambient cover path into maximal M-runs and exterior pieces. With the exterior skeleton and all run positions fixed, ordinary run fillings are counted by the complete ordered-cover profile q(M), which is invariant under converse. If u is outside M, its marked endpoint condition belongs to the exterior skeleton, so ordinary q(M) still suffices. If u is in M, its run is the terminal run of the designated path for E, or its initial run for S. Among all M-run slots this is one prescribed position. Permuting internal cover slots transports it to the first slot bijectively, so its filling count is E_u^M(t) or S_u^M(t), without a multiplicity t. Converse exchanges these two equal numbers. Summing over all skeletons proves every asserted equality. Empty exterior pieces and singleton runs are included.

This lemma concerns one marked physical endpoint, not a factorization of joint two-marked data. It does not assert invariance of a complete two-port response.

**Theorem 4.** Let D be vertex-transitive. Simultaneously XOR-converse any family of modules of the original D, including D itself if desired. The resulting D' has exactly the same q_h, E_u(h) and S_u(h) at every physical u as D.

By Lemmas 1–2 the original modules form a laminar family and each is vertex-transitive, hence endpoint-balanced. Perform the selected converses in increasing inclusion order. A selected module is still a module: earlier selected sets are contained in it or disjoint from it. Earlier contained operations preserve all its marked profiles by Lemma 3, inductively, so it is still endpoint-balanced when its own turn arrives. That converse again preserves every marked profile of the whole D. The final edge orientations are exactly the requested XOR, since the edge toggles commute. This proves the theorem even when the intermediate tournament is not known to be vertex-transitive.

The same induction shows that every original module's own marked profiles are preserved by selected descendants; selected ancestors merely whole-converse it and therefore still preserve its balanced profiles. No sequential-legitimacy assumption is imposed on a family that has not first been proved laminar.

## 3. Alternating-hub root-family addendum

Let m=2r+1>=3 and Q=W_m, with transitive chain 0->...->2r-1 and a hub beating even chain indices and losing to odd ones. Inflate the hub by A=C3. Inflate the two extreme chain names by arbitrary nonempty vertex-transitive tournaments D_0,D_(2r-1). Every other chain child is an arbitrary nonempty tournament. All physical names and original child sets are fixed.

Select any family of modules of this original substitution T containing the whole root with odd multiplicity, and XOR-converse the family to obtain V.

**Theorem 5.** The full named two-port responses of T and V differ. A fixed battery depending only on T, of at most m+20 actual queries, detects every such root-containing family.

As W_m is prime, every proper original module is contained in one child. The quotient is conversed exactly once. The original-module monotonicity theorem gives H(V)>=H(T). The scalar query (0,0), already present in the quartic battery below, detects strictness. Assume equality. Its hereditary equality conclusion gives identical complete q profiles on every corresponding child.

The hub C3 has no proper nontrivial modules. If that child is selected, its orientation is restored after the root converse. By Theorem 4, all proper selected operations inside each vertex-transitive extreme preserve every physical E/S profile there; the root converse also preserves them because they are balanced. Choose any physical i,j in the two extremes, and a directed edge p->q in the original hub. The two marked extreme profiles therefore equal their original balanced profiles, while all unmarked middle profiles agree by scalar equality. The exact marked-profile restoration in cyclic_child_disjoint/PROOF.md turns

    [a_p a_i c_q c_j](F_T-F_V)

into precisely the positive fixed-interior W coefficient from that proof. Its lower bound is

    2 E_i(1) E_j(1) product_(interior chain t) H(D_t)>0.

Here E_i(1)=H(D_0)/|D_0| and E_j(1)=H(D_last)/|D_last| are positive integers. The sixteen actual finite-difference queries therefore contain a separator. The argument uses equality of the marked profiles, not an unproved assertion that changed child interiors are identical.

If the hub child is not selected, its orientation is conversed. Queries uniform on every other whole child restore its target q profile to the converse of the original profile. The target count is then exactly the whole-converse count of T. Split the three hub names into singleton blocks and apply the exposed-anchor battery from returns_next/UNIFORM.md at p. Its actual residue-two separator applies.

There are m+2 compressed blocks, so the inherited safe anchor bound is m+4. Combined with sixteen quartic queries this is at most m+20; overlap may reduce it. No extra scalar query is needed because (0,0) is already in the quartic list.

## 4. Scope and dependencies

Theorem 4 is a marked-profile invariance result, not two-port invariance. Its conclusion remains useful when an intermediate child is not known to retain vertex transitivity. Endpoint balance implies regularity by the marked next-to-last-coordinate score formula, so Lemma 1 also makes its modules laminar. However, endpoint balance alone has not been shown here to make every submodule endpoint-balanced. That weaker hypothesis is therefore not substituted for vertex transitivity.

Theorem 5 strengthens the all-C3 root-family result to vertex-transitive extreme children and arbitrary middle children. It is detection against the original source for families containing the root. It does not claim full orbit injectivity for arbitrary vertex-transitive children: a family supported entirely within one such child may lie outside the current reconstruction results. It does not claim arbitrary prime quotient separation, arbitrary strong outer-ambient lifting, or invisibility of cubic/quartic coordinates. The whole displayed core, plus separately justified ordered SCC padding, is the scope.

Dependencies: cyclic_child_disjoint/PROOF.md for the exact disjoint quartic and all-order W positivity; module_decomposition/PROOF.md for hereditary scalar equality; returns_next/UNIFORM.md for the singleton-anchor battery; alternating_prime_quotient/PROOF.md for W primeness. The module and marked-profile lemmas above are proved here rather than inferred from finite tests.
