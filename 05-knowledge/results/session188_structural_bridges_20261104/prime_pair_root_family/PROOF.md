# Root-containing families with a marked pair child

## 1. Hypotheses and fixed battery

Let T=Q[D_i], where Q is prime of order m>=3, the D_i are arbitrary nonempty named tournaments, and D_z={p,q} is exactly a pair p->q. Suppose the exterior R=T-D_z has an ordered module X->Y, where X,Y are nonempty disjoint unions of complete other child blocks and |X union Y|>=3.

Let U result from simultaneous XOR reversal of any family of modules of the original T containing the whole root with odd multiplicity. Repetitions are reduced modulo two. There is no bound on the number, nesting depth or overlap of selected proper modules inside other children.

**Theorem.** F_T and F_U differ. A fixed battery depending only on T and the displayed ordered module works for every such family:

- the sixteen queries extracting [a_q a_X c_p c_Y];
- the compressed whole-converse anchor battery at p, with every D_i for i!=z kept uniform and p,q kept as separate singleton blocks, including its endpoint query.

The scalar query (0,0) is already among the sixteen quartic queries. There are k=m+1 compressed blocks. The anchor battery has at most k+2=m+3 queries, so the union has at most m+19 queries. For odd m, k is even and the explicit anchor construction has k main queries plus one endpoint query, improving this safe bound to m+18. Overlap between the two lists may reduce it further. These are detection bounds against the original source, not pairwise reconstruction of every target.

The two actual ports have fixed x->y. Their bits a_v,c_v mean old v->x and old v->y. No physical vertex is identified; 'uniform' restricts its port attachment only.

## 2. Proper modules and descendant parities

Every proper module of T is contained in one child. To see this, the set of occupied child indices of a module meeting at least two children is a module of Q: each unoccupied index sees all occupied pieces uniformly. Primeness forces occupation of every child. A prime tournament of order at least three is strong, so every quotient vertex has both an in-neighbor and an out-neighbor. Any omitted name in a partially occupied child would see occupied names in these two neighbors differently. Thus no name can be omitted and the original module is the whole root.

Consequently the root operation reverses Q exactly once, and all other selected operations are internal to the original D_i. Write U=Q^op[D_i'].

The original pair D_z has no proper nontrivial submodule. Hence precisely two cases occur. If D_z itself is selected, its edge is reversed twice, so D_z'=D_z. If it is not selected, D_z'=D_z^op. This explicitly retains the descendant contribution from the selected root.

## 3. Scalar strictness and the two equality branches

The original-module decomposition theorem gives H(U)>=H(T). If strict, the actual scalar query (0,0) separates. If equality holds, that theorem gives equality of all induced-subset Hamilton counts and therefore all ordered path-cover profiles. In particular,

    q_h(D_i')=q_h(D_i) for every i,h.                    (1)

On any query uniform on a child, maximal child-run compression shows its Hamilton count depends only on the full q profile of that child. Thus each other D_i' can be restored to D_i, or to D_i^op, without changing the actual query count. This statement is about entire uniform child blocks, not arbitrary endpoint marginals.

**Pair selected.** The pair interior is unchanged. On the quartic battery, every other child is uniform, so (1) makes the target response identical to Q^op[D_i] with original interiors fixed. The marked-pair and ordered-module identities give

    [a_q a_X c_p c_Y](F_T-F_U)
      = C_YX(R)-C_XY(R)
      >= H(R/(X union Y))
          [q_2(X)H(Y)+H(X)q_2(Y)+q_2(X)q_2(Y)] > 0.

The exact sixteen-query contrast therefore separates. All degree-four terms are retained.

**Pair not selected.** The pair interior is reversed. On any query uniform on all other children, (1) makes F_U exactly equal to F_(T^op), since every other child may be restored to its converse profile. Split the pair into singleton blocks {p},{q}. The uniform compression theorem in ../returns_next/UNIFORM.md applies at the uncovered anchor p. Its explicit anchor battery contains an actual whole-converse contrast 2 modulo 4. Therefore the same query separates T and U.

These two branches exhaust the scalar-equality case. No sequential legitimacy is assumed for overlapping original modules; the decomposition theorem is used before the profile restoration.

## 4. Critical prime quotient consequence

If Q is any critically indecomposable tournament with at least five vertices, every D_i is nontrivial, and at least one child is a pair, the critical-quotient synthesis supplies the required ordered exterior module for any chosen pair root. Thus every original-module XOR family containing the whole root is visible to the full named response, without a bound on proper descendant operations. Since the critical quotient order m is odd, a fixed battery of at most m+18 queries suffices.

The sharper cyclic, staircase and alternating-hub local size conditions can replace the all-nontrivial-child hypothesis whenever their ordered exterior module has at least three physical names. The C3 case is included by the same argument when its two exterior children have total order at least three.

For laminar selected families, scalar equality is automatic; this does not permit omitting the scalar check for arbitrary crossing families. It is already present at no extra cost in the quartic battery.

## 5. Scope and dependencies

This theorem concerns the whole displayed core, or ordered earlier/later SCC padding under the existing exact localization theorem. Arbitrary strong outer ambient lifting is NOT proved: the positive quartic branch is an exact integer result without a guaranteed nonzero residue modulo eight. The whole-converse branch alone does have modular transfer, but that does not establish transfer of the combined theorem.

The result distinguishes the original T from each root-containing target. It does not assert injectivity among all pairs of targets, arbitrary prime quotient reconstruction, all original-module families omitting the root, or legal sequences using modules that were not original modules.

Dependencies: ../module_decomposition/PROOF.md (scalar strictness and hereditary equality); ../returns_next/UNIFORM.md (compressed exposed-anchor separator); ../marked_pair_quotient/PROOF.md and ../ordered_module_skew/PROOF.md (the exact quartic and positive bound); ../critical_quotient_synthesis/PROOF.md (classification-based critical scope).

## 6. Verification

The checker enumerates root-containing original-module XOR families for all three five-name critical quotients with pair children, a staircase case with a TT3 child and overlapping selected intervals, a cyclic third-child control, and a larger cyclic quotient under its sharper local size hypothesis. It computes every retained query as a full integer Hamilton count on the augmented physical tournament. The scalar-strict, pair-selected and pair-not-selected branches are checked separately. Sources, dependencies and complete retained family records are hash-bound. Both normal and optimized modes must agree.
