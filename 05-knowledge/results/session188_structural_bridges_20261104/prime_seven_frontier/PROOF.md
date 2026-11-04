# Original-module collisions require highest prime quotients of order at least seven

## 1. Exact action and reviewed inputs

Fix a named tournament T. A family F consists of modules of this ORIGINAL T. Its simultaneous action toggles each edge according to the parity of selected modules containing both endpoints. Repeated selections cancel modulo two. The family may overlap, nest or cross; its members need not remain modules in an arbitrary sequential ordering. Let U be the resulting tournament on exactly the same names.

Write F_T(a,c) for the complete physical two-port Hamilton response with x->y and a_v=[v->x], c_v=[v->y]. Equality means equality for every actual attachment mask, not only uniform queries, a low-degree truncation, or an isomorphism class.

This addendum uses the following separately certified inputs:

- `../prime_five_frontier/PROOF.md`: any nontrivial original-module full-response collision can have only highest changed prime quotients of order at least six. It incorporates the complete arbitrary-child/fixed-ambient prime-three and prime-five closures and the linear-node exclusion.
- `../prime_six_reconnaissance/PROOF.md`: every six-name quotient has nonzero Hamilton endpoint parity, by a complete exact table of all32,768 labeled tournaments. Consequently its quotient converse with arbitrary nonempty named children and matching corresponding complete profiles is separated by a single uniform physical endpoint query in every fixed ambient.
- `../module_decomposition/PROOF.md` and `../highest_quotient_reduction/PROOF.md`: the original strong-module tree gives a common recursive modular partition of source and target. Scalar Hamilton equality implies equality of every ordered-cover profile on every fixed induced subset. Highest changed quotient nodes are disjoint, all changed edges lie inside them, and a changed prime quotient is its converse.
- `../monotone_dual_queries/PROOF.md`: core converse-complement query duality and monotonicity under coordinatewise child-profile increases, with exterior port incidences fixed.

Here strong means a module overlapping no other module; it is not a synonym for strongly connected. The original strong-module tree need not be the target's canonical strong tree. The six-name finite certificate and this integration each require their own independent review before the frontier is changed.

## 2. Root-family consequence and its limits

Let a fixed source core be Q[A_1,...,A_m], where Q is prime of order m in{3,5,6}, and every child is nonempty on a fixed physical set. Fix any rooted ambient R. Allow any family of original modules contained inside this core, selecting the whole core oddly. No selected operation involves an exterior name.

Each proper original module lies inside a root child. The target root quotient is therefore Q^op; target child profiles are coordinatewise no smaller than their source profiles by hereditary original-module monotonicity. Internal converse caused by a selected ancestor does not change a child's complete profile. No equality between distinct source children is required.

For m=3,5 use the reviewed arbitrary-child fixed-ambient baseline separator; for m=6 use the nonzero endpoint-parity coordinate in the complete six-name table. Let z be this child-uniform baseline query. At every old exterior vertex, BOTH ports imitate the original core root, on both sides. Pair z with Jz=(1-c,1-a) on CORE bits only, keeping the same exterior bits in both queries. Their exact baseline ambient contrasts are negatives. Fix the member with negative source-minus-baseline contrast. Nonnegative run substitution makes target profile increases unable to cancel it.

Thus the fixed pair{z,Jz} detects EVERY such root-containing original-module family relative to this source; one baseline-sign-selected member works for all its targets in the specified fixed ambient. The sign selection may depend on that ambient. This statement neither covers root-even families by itself nor asserts pairwise injectivity between arbitrary family outputs. It permits crossing descendants, but it does not permit operations outside the core.

The global collision argument below does NOT assume all operations lie inside one core or apply this root-family assertion directly to exterior changes. Instead it uses scalar equality and exact restoration of the other highest nodes.

## 3. Excluding a highest changed six-child quotient

Suppose U differs from T but F_T=F_U. The physical query a=c=0 has x,y before every old name, so F_T(0,0)=H(T) and H(T)=H(U). The hereditary theorem gives

    q_h(T[S])=q_h(U[S]) for every fixed old-name subset S and every h.

The reviewed prime-five frontier now says every highest changed original strong quotient is prime of order at least six. Suppose one such node M has exactly six original children A_0,...,A_5. Their physical sets remain modules in source and target; their complete corresponding profiles agree. The target quotient is the converse of the original six-name Q.

Form V by changing T[M] to U[M], leaving every exterior name and arc as in T. The exterior incidences of M are unchanged. The six-block theorem supplies a child-uniform physical endpoint query on T and V, with both ports imitating M at the exterior, such that

    F_T(a,c)-F_V(a,c)=2 modulo4.                         (1)

Descendant changes inside A_i do not invalidate (1): it explicitly allows arbitrary target interiors with the same corresponding complete profiles. One may equivalently restore such child interiors using profile substitution, but no literal equality of interiors is assumed.

For every other highest changed node N, the complete source and target profiles agree. It is disjoint from M and from the other highest nodes. The old cross arcs between M and N are uniform; both ports have M's incidences to N and hence see N uniformly. Therefore N remains a module in the actual augmented graph for this SAME query.

Replace each such N in V by its U interior. Exact occupied-run substitution preserves the queried Hamilton count because its complete q profile agrees. Disjointness and unchanged boundaries keep every replacement valid. All changed edges lie inside the highest nodes, so these replacements end at U. Hence F_V(a,c)=F_U(a,c), contradicting (1).

This also covers multiple highest nodes, nodes with changed descendants, and another highest node with an unresolved larger prime quotient. There is no assumed generic scalar transfer or cancellation principle for arbitrary integer gaps; the argument uses a concrete physical six-block query and exact equal-profile restoration.

## 4. The sharpened necessary frontier

Combining Section3 with the reviewed order-at-least-six frontier proves:

**Corollary.** If T differs from an original-module XOR target U but their complete named two-port responses coincide, EVERY highest changed original strong quotient is prime of order at least seven.

In particular, if every prime quotient in T's original strong-module tree has order at most six, EVERY original-module XOR target different from T is separated by the full named response. Linear nodes are already excluded; the only possible prime orders below seven are3,5,6, with the no-prime-four fact and the linear two-child convention supplied by the earlier frontier.

This is a necessary condition on a still-possible collision, not an existence claim at order seven or an exclusion of larger prime quotients. It is source-relative and concerns simultaneous action of modules of the original T. It is not unrestricted tournament reconstruction and is not extended to arbitrary sequences whose later modules were not original modules.

## 5. Controls and frozen scope

The checker pins the reviewed proof/data inputs and builds two source-relative multiple-highest-node examples in a strong cyclic ambient. One has a six-child prime node with a changed transitive descendant and a separate prime-five node. The other has a six-child prime node and a separate prime-eight node whose endpoint parity is zero. It computes the entire original module/strong-node structure, checks the actual original-module XOR target, identifies the exact highest-node antichain, verifies every required child and other-node complete profile equality, and counts all actual source/intermediate/target physical queries.

Those focused controls demonstrate the restoration hypotheses, not a census from which the corollary is extrapolated. The complete named cubic/quartic banks remain retained in the pinned component packets. Both normal and optimized modes run the same explicit checks.
