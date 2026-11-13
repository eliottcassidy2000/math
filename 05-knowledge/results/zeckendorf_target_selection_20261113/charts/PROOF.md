# Active ancestors, constrained marks, and a physical symmetry selector

**PROVED:** weighted greedy selection on a fixed tree, the decorated maximum-spanning-tree identity, exact constrained counting/ranking, and a root-child lift to actual tournament automorphisms. **FINITE-EXACT:** the small controls. These results select attainable symmetry targets under declared port constraints; they do not select arbitrary Hamiltonian-path totals or prove Chan–Pak.

## 1. What is inherited and what has to be added

[Session312's cluster chart](../../session312_cluster_pythagorean_symmetry_20261113/CLUSTER_CHART.md) puts positive minors on ternary leaf addresses. Their 3-adic valuations recover common-prefix lengths; their signed leading units retain cyclic child orientation. Its farthest-leaf greedy theorem is for a homogeneous tower. [Session311's marking budgets](../../session311_recursive_symmetry_20261113/MARKING_BUDGETS.md) gives the corresponding unconstrained extrema and no-gaps exponent interval.

**Incoming priority, origin/main47580b0cca:** [Session313 PREFIX_SELECTION](../../session313_zeckendorf_target_selection_20261113/PREFIX_SELECTION.md) independently proves the same nonnegative weighted-tree greedy theorem, including every tie, the active-rotation interpretation, and exact attainable-weight recursion. These parts below are concurrent recoveries, not advances beyond that incoming result. The added interfaces here are the weighted common-ancestor maximum-spanning-tree formula, original-position adjacency and role quotas, owner-consistent dictionary choices, full rank lifts, and the declared root-child birth compiler. Session313 also supplies the same colored interior-gap mechanism, which we retain as a shared hostile control.

Here the leaves can already have different physical colors, some ports may be unavailable, and a chosen subset may have color quotas or forbidden adjacent positions. These are genuinely different selectors. The live coordinates are ancestry / oriented child classes / active rotations / physical port eligibility / boundary carries / exact residual target.

The orthogonal source is [boundary palette compression,§§1–2](../../boundary_palette_compression_20261011.md): keep the named boundary relation rather than just a count of colors, and retain lifts of every surviving state. The [multi-source completion theorem,§2](../../session45_injective_multisource_moves_20261023/injective-completion/THEOREM.md) similarly separates an actual resource choice from a coupled guard. The older [adaptive boundary selector](../../adaptive_boundary_selector_20261004.md) explains why retaining only one promising alternative can discard an available completed continuation. None of their original geometric or Collatz conclusions is imported here.

The Fibonacci connection below is exact at the constraint interface. [S41's colored-grid packet,§4](../../session41_owned_attachments_fibonacci_20261023/COLORED_GRID_AND_PHASE_QUOTIENT.md) shows that a forced-zero slot must be eliminated existentially, not removed and reconnected. The same distinction is essential in our port table.

## 2. Colored cyclic hierarchies and the missing active coordinate

Build a finite rooted ternary hierarchy by cyclic joins. Leaves carry actual rigid tournament types from a fixed equal-order, equal-full-profile palette, such as the inherited B and D=A▷1. The cyclic three-child quotient is oriented, so the canonical class of an internal node is its ordered triple of child classes **modulo rotations**, not arbitrary permutations. Put

    a(v)=1 if all three child classes at v are isomorphic,
         0 otherwise.

The full automorphism group is recursively a direct product when a(v)=0, and a cyclic wreath extension when a(v)=1. In particular its order is3^E, where E=Σ_v a(v). This is the inherited recursion applied to a supplied colored hierarchy, not a new isomorphism classification.

Choose one named vertex inside each selected rigid leaf. For the selected leaf set S, its pointwise stabilizer has exponent

    e(S)=E−Σ_(v with a descendant in S) a(v).           (1)

Indeed a selected descendant prevents any nontrivial cyclic rotation at its ancestors. Entirely unselected subtrees keep their original colored automorphism groups. Product decomposition over the maximal unselected subtrees proves (1), including heterogeneous shapes. The leaf rigidity hypothesis matters: a fixed point in a nonrigid leaf could leave an extra internal stabilizer.

The valuation chart recovers the tree, but not a(v). Its sidecar must contain the actual child colors/classes. Its leading-unit orientation is also necessary: an unoriented three-child ball has an S3 action, whose one-point stabilizer has order2; a directed cyclic triple has C3, whose one-point stabilizer is trivial. These are different marking problems. The chart is still only a chart for the hierarchy: a D leaf need not admit the full 3-adic determinant orientation, as Session312 explicitly proves.

**Sharp colored hostile to the old raw-distance rule.** At depth two, use the three child color triples

    BBB, DDD, BBD.

The root has a=0; the first two child nodes have a=1 and the third has a=0. The unmarked group has order9. After marking address00, raw common-prefix valuations tie addresses10 and20 as equally far. Marking20 leaves group order3; marking10 gives order1. Thus arbitrary farthest-distance ties are no longer safe after colors change, even though the old address minors and all leaf cover polynomials are unchanged.

## 3. Weighted active-ancestor greedy is exact

The useful extension is broader than (1). Give every internal node of any finite rooted tree a nonnegative real weight w(v), and let L be any specified subset of eligible leaves. For S⊆L define

    C(S)=Σ_(v whose descendant set meets S) w(v).       (2)

For symmetry use w=a. A mark's gain is the weight of its ancestors not already covered. Start with S empty and repeatedly choose an eligible unselected leaf with maximum gain.

**Theorem.** Every tie choice in this greedy rule gives a nested family maximizing C(S) among all k-element subsets of L, for every feasible k. Eligibility is arbitrary. Negative node weights, unequal mark costs, and additional coupled constraints are not included.

**Proof.** For a subtree v let F_v(k) be its optimum with exactly k eligible marks, F_v(0)=0. At a leaf the only available increment is zero, if that leaf is eligible. For an internal node, allocate k among its children and add w(v) once when k>0. Inductively each child's successive optimum increments are nonincreasing. The best allocation is obtained by merging those child increment lists in descending order: a prefix of any list always precedes a later, no-larger entry, and ties can be ordered to respect prefixes. Thus the max-convolution has nonincreasing increments. Adding w(v) only to its first increment preserves that property.

The same induction proves the physical greedy statement. Within a chosen child, a maximum-gain leaf is a valid greedy step for that child. Between children, the next available marginal is their next list entry. The first parent charge is common to every possible first child and therefore does not change that choice. Every tie produces an optimal child-prefix allocation and an optimal global prefix. This proves all k simultaneously. ∎

In particular, no generic submodular approximation loss is needed for this laminar ancestor objective. One can implement the optimum by merged marginal lists; the simple code checks the literal greedy choices instead.

### A decorated maximum-spanning-tree formula

Let D(x) be the sum of weights on all internal ancestors of leaf x. For distinct x,y let c(x,y) be the sum on their common internal ancestors, including their lowest common ancestor. Then for nonempty S,

    C(S)=Σ_(x∈S)D(x)−max_spanning_tree_on_S Σ_edges c(x,y).  (3)

To prove it, let m_v be the number of selected descendants of v. Any tree on S has at most(m_v−1)_+ edges entirely inside that descendant set. A path through the leaves in depth-first order attains every one of these upper bounds simultaneously, since each subtree is an interval. Nonnegative weights therefore give maximum weight Σ_v w(v)(m_v−1)_+. Subtract from Σ_v w(v)m_v to get (2).

For S nonempty, the greedy gain at x is

    D(x)−max_(s∈S)c(x,s).                             (4)

The covered portion of x's ancestor chain ends at its deepest common ancestor with an existing mark. Formula (4) follows even with zero weights. With a homogeneous depth-d tree, D=d and c=1+common-prefix-length, recovering Session312. For colored trees the missing D(x) and active weights cannot be replaced by raw valuation alone. We decorate the original chart; we do not assert that contracting every zero-weight split preserves its signed orientation data.

## 4. Fibonacci admissibility needs a boundary state

Order leaves by their original left-to-right addresses. Impose the actual hard-core constraint that two consecutive original leaves may not both be selected. A subset of eligible leaves is not first compressed into a shorter row. One may also assign an external role color to each leaf and prescribe a quota vector. These role colors are separate from the B/D physical palette colors unless an explicit identification is supplied.

Greedy is no longer a complete selector. Already on three leaves, choosing the middle leaf prevents a second choice, although the target budget two is feasible via the two endpoints. If the middle leaf is forbidden, those endpoints remain jointly legal. Deleting the middle position and treating the two eligible leaves as adjacent gives the wrong answer. This is exactly S41's forced-zero relation A–Z–B, now with named physical ports.

For any subtree interval store states

    (k, first, last, c, q),                           (5)

where k is the number of marks, first/last say whether the **first and last original leaves of the interval** are selected, c is covered internal weight, and q is the role-count vector. For the exact symmetry target, w=a is integral and one requests c=E−e.

At a leaf the unselected state always exists; the selected state exists exactly when the leaf is eligible. Merge consecutive child intervals by adding k,c,q and forbidding only last_left=first_right=1. The merged boundary flags are first_left,last_right, even when a whole interval is unselected. At each internal tree node, add w(v) to c exactly when k>0. Filter the desired budget, quota, and target weight at the root.

**Exact selector theorem.** This table gives all and only feasible physical marking sets, and hence decides every requested attainable residual symmetry exponent under the stated constraints. Counting each state by products of child counts and sums over decompositions gives its exact multiplicity. Backpointers give a rank/unrank bijection with the accepted markings in a declared deterministic table order.

**Proof.** A selected set has unique child restrictions and unique boundary flags. The only new adjacency across a merge is between the two original boundary leaves, which are exactly the two stored flags. Its covered-node weight adds across children and charges the parent iff the set is nonempty. Role quotas add. These observations prove necessity and sufficiency by induction, including witness reconstruction. Products count independent restrictions; distinct decompositions or child lifts produce distinct subsets. Mixed-radix rank/unrank is therefore exact. ∎

For maximization at fixed cardinality and no quotas, retain only the best weight per(k,first,last), giving a standard O(|tree| K²) arithmetic-operation implementation. Exact weight distributions and quotas require larger tables; with integer weights this is pseudo-polynomial in the weight and quota bounds. No claim of polynomial complexity in binary-encoded arbitrary weights or exponentially implicit trees is made. The supplied implementation retains the complete exact-target table and all rank lifts, not one best witness.

### Constraints destroy the old no-gaps conclusion

At depth three, take one top subtree monochromatic and the other two rigid, with all their internal a(v)=0. For example use leaf colors BBBBBBBBB in the first subtree and BBD BDD BBD in each other subtree. Only the first subtree's four internal nodes are active. Permit one leaf in that subtree and one in a rigid subtree, and require exactly one mark. The possible killed weights are2 and0, so the residual exponents are2 and4, never3. Both ends of the interval are attained, but its middle is not. A table retaining only minimum and maximum exponent would invent a false target.

### The actual three-color row can serve as a role assignment

S41 proves that the row-N role at position j is letter N−j+1 of the word generated by C→CA, A→B, B→CA. We use that authenticated row as an input labeling, without identifying it with every no-consecutive decomposition. At N=9 the row is

    C B A C A C B A C.

On the colored depth-two hostile of §2, require three nonadjacent marks, one of each role C,A,B, and residual group order1. The exact table has four accepted markings. Rank zero in the implementation's deterministic order is the zero-based set{1,3,7}. This is a concrete constrained source choice satisfying the joint guards, not three separately successful tests. The original unit colors, physical palette colors, and any mod-two quotient charge remain separately typed.

### Shared dictionary choices must precede the marking table

The sibling [shared-dictionary theorem](../history/SHARED_DICTIONARY.md) supplies an authenticated owner map on the stripe. Every occurrence of one owner must use the same choice; a length-N prefix has only O(log N) visible owners. For a shared-dictionary role labeling, enumerate those legal owner assignments first, construct the resulting whole row, and then run the table above. This is an exact joint selector for a pair consisting of a consistent dictionary and a feasible physical port set. Necessity restricts any joint witness to its owner assignment; sufficiency uses the dictionary theorem's infinite extension and the table's reconstructed marking. Unseen owner bits may be set to their canonical values.

The extra factor is the number of visible dictionaries, polynomial in numerical N under the stated R/D grammars, not2^N independent position choices. It is not a polynomial-time claim in log N or in an implicit description of a huge hierarchy. Different dictionaries can yield the same physical marking; when counting joint witnesses, these pairs must remain distinct. The exact residual target and all quotas are tested on one common owner assignment, not on separately chosen role labels at each occurrence.

`joint_selector.py` performs this composition on row9, retaining the fixed physical B/D hierarchy while allowing either authenticated shared-dictionary grammar to determine the external C/A/B roles. It compares the complete rank-decoded marking sets with literal subset enumeration. Role colors are constraints on supplied ports; no claim identifies them with tournament arcs or Hamilton paths.

## 5. From marks to a genuine induced graph operation

Pointwise marking alone does not change a graph or its unmarked automorphism group. There is, however, a physical lift of every accepted selection.

Assume the initial colored cyclic hierarchy is strong (positive depth suffices). For each selected rigid leaf, choose its declared named vertex and perform one actual root-child birth there:

    r→z→all other current vertices.

Keep every old vertex and arc. Apply the births in a specified order, such as increasing leaf address. The inherited [ANCESTRY A1–A2](../../adaptive_metagraph_reformulation_20261103/ANCESTRY.md) makes each newborn canonically recoverable and gives

    Aut(final tournament)=Aut(initial tournament) fixing all chosen roots. (6)

One can also apply its one-birth isomorphism restriction successively: the legal newborn and its predecessor are fixed, so the automorphism group restricts to the old root stabilizer, and conversely every such stabilizer element extends. Thus the output's **actual global** group has order3^e for precisely the target e selected by the table. The construction adds k vertices and preserves the entire old induced graph. It works at every accepted address, not just a numerical exponent.

For unconstrained budgets, §3 therefore minimizes the residual global symmetry using k births whose roots belong to the supplied eligible original leaf ports. This is not a minimum over arbitrary tournament edits or births at unrestricted fresh ports. It is also not count-neutral: H and the full rooted response generally change. To target a Hamiltonian count one still needs the compatible response data and a separate exact arithmetic entry theorem. The selector resolves the marked symmetry obligation; it does not silently discharge that other obligation.

Incoming [Session313 PHYSICAL_MARKS](../../session313_zeckendorf_target_selection_20261113/PHYSICAL_MARKS.md) already realizes every pointwise mark set on a homogeneous B tower with two exterior probes and a suitable red/blue incidence coloring. That is a different protocol with a smaller added-vertex budget; the k-birth compiler here is not a general minimum. Its useful sidecar is the inherited, canonically recoverable root-child ancestry at each step. Session313 further proves that arbitrary r-vertex extensions retaining a monochromatic depth-d tower have a path-count divisor3^{Σ_(j=0)^(d−1) max(3^j−2r,0)}, even when the final graph is rigid. Thus physical symmetry selection alone cannot remove the path-packet obstruction. Neither our constrained selector nor the neutral sealed-interface repair in selection/ claims otherwise.

## 6. Reproduction and exact scope

Run `python3 -B 05-knowledge/results/zeckendorf_target_selection_20261113/charts/selector.py`.

The controls exhaust all512 binary colorings at depth two and all239344 reachable greedy-prefix sets, including arbitrary ties. A further81 nonnegative weight vectors with restricted eligibility give2511 maximum-spanning-tree checks. An independent enumeration of the81 address rotations checks all512 pointwise marking masks of the9-element color stabilizer fixture. Nine constrained tables are compared with literal subset enumeration, including their complete rank/unrank lifts. Four named birth programs on the actual72-vertex colored tournament check the old induced graph, strongness, and the lifted automorphism action; absence of additional automorphisms is supplied by the all-order ancestry theorem, not an unrestricted75-vertex permutation search. Further controls retain the forced-zero gap and the missing interior exponent. These finite universes support the elementary all-tree proofs; they do not classify unbounded color palettes or infer arbitrary tournament path-count coverage.
