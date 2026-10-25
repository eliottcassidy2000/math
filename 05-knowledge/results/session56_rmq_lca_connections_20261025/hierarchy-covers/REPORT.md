# Minimum same-label hierarchy covers

Status: PROVED by elementary arguments below; FINITE-EXACT controls are supplementary. No theorem ID is claimed or reserved. No novelty/priority claim is made.

## Scope and definitions

Fix a finite nonempty label set U. A rooted hierarchy is a finite rooted tree whose leaves are bijectively labelled by U. Its cluster at a node is the set of labels below that node. Internal nodes need not be binary; unary nodes are unnecessary. Every tree has full-set cluster U and every singleton cluster. A family F consists of nonempty subsets of U. Repeated requests for the same subset can be deduplicated: occurrence multiplicity, node heights, and edge lengths are not constraints in this problem.

A hierarchy cover is a collection of trees, EACH on exactly the same labels U, such that every requested cluster is a node cluster in at least one tree. Additional clusters are allowed. Trees may not replace labels by copies, delete labels, or relabel different requested occurrences independently. Let h(F) be its minimum size, allowing a zero-tree collection when F is empty. If one insists a representation always have at least one tree, its minimum is max(1,h(F)).

The empty subset is not an ordinary node cluster. If it is requested literally, there is no finite cover; either exclude it as here or explicitly declare it a free formal bottom element outside the trees. U empty is excluded because ordinary rooted trees have at least one leaf. For |U|=1, the single-node tree works. Full-set and singleton requests are harmless, but a nonempty family consisting only of them still requires one tree.

Define the crossing graph X(F): vertices are the requested subsets; distinct A,B are adjacent exactly when A∩B is nonempty, A is not a subset of B, and B is not a subset of A.

## Theorem: exact coloring equivalence

For every finite family F of nonempty subsets of a finite nonempty U,

h(F) = chi(X(F)), with chi of the zero-vertex graph equal to zero.

Proof, cover to coloring. Clusters of one rooted tree are laminar: if two nodes have a common descendant leaf, both lie on that leaf's path to the root and hence are comparable by ancestry; their clusters are therefore nested. Otherwise their clusters are disjoint. Given a cover of t trees, assign each requested subset to any one tree containing it. Adjacent crossing subsets cannot be assigned to the same tree. This is a proper t-coloring, so chi ≤ h.

Proof, coloring to cover. A proper t-coloring partitions F into laminar families F_1,...,F_t. For each color, set H_i = F_i ∪ {U} ∪ {{u}: u∈U}, deduplicated. This remains laminar. For every A in H_i other than U, make its parent the minimal strict superset of A in H_i. Such a parent exists because U exists and is unique because supersets of the nonempty A intersect and are therefore a chain. These parent edges produce a rooted tree, since strict inclusion prevents cycles and iteration terminates at U. The minimal sets are precisely the singletons. The maximal proper H_i-subsets below any nonsingleton A partition A: every point has its singleton below A and extends to some maximal proper subset; distinct maximal proper subsets are disjoint by laminarity. Inductively the descendant-leaf set of node A is exactly A. Thus the tree represents every F_i cluster on the SAME leaves U. This proves h ≤ chi. The empty-family case is immediate.

A binary-tree requirement alone does not change the answer: resolve a multifurcation by grouping its disjoint child clusters, adding unions of children. This refinement preserves every existing cluster and introduces no incompatible node. However, restricting the allowed extra clusters, fixing heights, or imposing external tree shapes is a different optimization problem.

## Affine-plane control: four trees, not one

Let U=F_3². The 12 affine lines form four parallel classes, each of three lines: x=b, and y=mx+b for m∈F_3. Different lines in the same class are disjoint. Lines of distinct directions meet in exactly one point: solving their two independent linear equations gives a unique point. Their cardinalities are both three and they are distinct, so neither contains the other. Hence X is exactly K_{3,3,3,3}.

Four colors suffice by the parallel classes. Taking one line from each class gives a clique of size four, so four are necessary. Each representing tree has root U, three line-cluster children in its parallel class, and the nine labelled points as leaves. The crossing graph has 6·3·3=54 edges.

This is a shared-label obstruction. As unlabelled shapes all four trees are identical, but that does not permit their line clusters to coexist in one same-label tree. For example { (0,0),(0,1),(0,2) } and { (0,0),(1,0),(2,0) } cross.

## Orthogonal partitions and CRT

Let U=A×B with |A|=a≥2 and |B|=b≥2. Row blocks {x}×B and column blocks A×{y} are disjoint within a partition and cross between partitions. Their crossing graph is K_{a,b}; exactly two same-label trees are required to represent all these blocks.

For an explicit integer CRT instance, take U={0,1,2,3,4,5}. The mod-2 partition has blocks {0,2,4},{1,3,5}; the mod-3 partition has blocks {0,3},{1,4},{2,5}. Every cross-partition intersection is a singleton, giving K_{2,3}, six crossing edges, and h=2. The pair {0,2,4} and {0,3} is already a hostile against one tree.

More generally, on U=∏_{i=1}^r A_i with all |A_i|≥2, the coordinate-fiber partitions have crossing graph the complete r-partite graph with part sizes |A_i|. Fibers of distinct coordinates intersect, and each has points outside the other because both alphabets have another value. Thus exactly r trees are needed. CRT identifies Z/(∏ p_i)Z with this product for distinct primes, so all r first-level prime congruence partitions require exactly r trees.

## Independent prime levels versus a diagonal tree

On a finite set U of labelled integers, let P_{p,k} be the partition by congruence mod p^k, k≥0. For each fixed prime p this refines as k increases, so all its nonempty blocks form a laminar family. This remains valid when residue representatives are restricted to a subset, and for several labels with the same integer value. Singleton leaves can always be added; if equal integer values have different labels, those singleton refinements are bookkeeping rather than a genuine congruence distinction.

For a finite prime set S, the common refinement

D_k = { ∩_{p∈S} B_p : B_p∈P_{p,k}, intersection nonempty }

also refines with k. Indeed membership together at level k+1 implies membership together in every prime partition at level k. Thus all D_k blocks form a laminar family, representable by ONE rooted tree. Equivalently, for integer labels D_k is congruence modulo ∏_{p∈S} p^k. Repeated partitions or repeated blocks at many k are suppressed in the finite unweighted tree; retaining their levels requires decorations. Since U is finite, there are only finitely many distinct subsets even if k ranges without bound.

This is not a contradiction to the two-tree or r-tree lower bounds. The diagonal tree asks only synchronous common-refinement questions and generally does not have the individual prime blocks as nodes. In the mod-2/mod-3 six-label example, D_0={U}, D_1 is already the singleton partition, and D_k stays singleton thereafter. Its star tree has none of the five nontrivial individual prime blocks as a node. Queries for congruence mod 2 alone or mod 3 alone cannot be answered by a single LCA-cluster query in that star.

Precisely state the loss: it is loss of independent-prime cluster geometry from the bare diagonal hierarchy, not necessarily erasure of numerical data stored as labels. If integer values or the prime partitions are supplied as a sidecar, those queries can be reconstructed externally. Without that sidecar, many distinct independent partitions have the same diagonal star. For example on labels {a,b,c,d}, the pair {{a,b},{c,d}} and {{a,c},{b,d}} has singleton common refinement, and replacing its second partition with {{a,d},{b,c}} gives exactly the same common refinement while changing an independent query. The labelled diagonal tree alone cannot distinguish these inputs.

For a finite prescribed collection of prime-depth blocks, the crossing-graph theorem computes the exact minimum number of trees. One tree per prime always supplies an upper bound; equality need not hold for an arbitrary sparse label sample, since partitions may coincide or fail to cross. The full CRT product first-level example certifies sharpness. Taking all finer prime levels cannot increase the requirement above r and cannot remove the first-level r-clique, so it also has exactly r trees there.

## Arbitrary graphs: the optimization is not an LCA shortcut

For a finite simple graph G=(V,E), make a distinct private label q_v for every v and a distinct label q_e for every e. On U={q_v:v∈V}∪{q_e:e∈E}, define

S_v={q_v}∪{q_e:e is incident to v}.

For distinct vertices v,w, neither S_v contains S_w: q_v is absent from S_w and q_w is absent from S_v. Their intersection is nonempty exactly when vw∈E, in which case it consists of q_{vw}. Thus X({S_v}) is isomorphic to G, not merely a homomorphic image or induced subgraph. The construction uses |V|+|E| labels and |V|+2|E| incidences, polynomial in the graph input. For the graph with no vertices use any dummy nonempty universe and the empty requested family.

Consequently a minimum-hierarchy-cover algorithm, or a decision algorithm for h(F)≤t, gives a graph-coloring algorithm or corresponding decision algorithm via this explicit size-controlled reduction. This is the exact content of “as hard as graph coloring” here; no external complexity theorem is needed for that reduction claim. The reverse direction is also immediate after constructing X(F). Standard LCA routines on a GIVEN tree do not solve the preceding choice/cover problem.

A clique number is only a lower bound, not the general answer. Realize the five-cycle by this construction: its maximum clique has size two but h=3. This hostile prevents replacing chi by maximum pairwise-crossing subfamily size, even though the affine and CRT examples have equality.

## Connection contract and boundaries

Source: a finite family of subsets with fixed leaf labels. Target: its crossing graph. Map: retain precisely pairwise failure of laminar compatibility. Preserved predicate: existence of a cover by at most t unrestricted rooted same-label trees, exactly equivalent to t-colorability. Lost information: actual intersections, inclusions, subset sizes, labels, hierarchy heights, arithmetic origin, and consumer-specific query semantics. Restoration sidecar: the labelled subsets and, for arithmetic applications, prime and depth annotations. A coloring plus this sidecar constructively produces the trees; the graph alone does not identify the labelled clusters.

Closest proved mechanism used: the elementary laminar-family/rooted-tree correspondence, proved here. Canonical hostile: crossing row and column blocks on a shared product. Corrected near miss: multiple prime partitions need not fit one tree, but their synchronous common-refinement filtration always does. Orthogonal control: arbitrary-graph realization and its odd-cycle instance. Concept board: laminarity, graph coloring, CRT product fibers, diagonal filtrations, and label/annotation loss. Their comparison exposes a combinatorial encoding theorem, not an arithmetic cancellation theorem.

No claim is made about generic Harder–Narasimhan filtrations, a tournament map, signed or weighted cancellation, an unrestricted arithmetic hierarchy, or LRC(14). This result is a self-contained finite set-system statement.

## Exact controls and reproduction

Run: python3 verify.py

The script independently compares backtracking graph coloring against direct set-theoretic laminar-subfamily cover dynamic programming; then constructs the inclusion-parent trees and checks that children partition every nonsingleton node, leaves are the singletons, and requested clusters survive.

Universe and filters: all 128 families of nonempty subsets of a three-element label set; all labelled simple graphs on n=0,1,2,3,4,5 vertices (1,1,2,8,64,1024 respectively) under the private-label construction; the 12 affine lines; the five CRT blocks; the diagonal CRT blocks; and the private-label five-cycle. No probabilistic filters or numerical tolerances are used. Expected affine result: 12 vertices, 54 edges, minimum four trees. Expected CRT result: five vertices, six edges, minimum two trees. Expected diagonal result: minimum one tree with five independent nontrivial blocks absent. Expected five-cycle result: minimum three trees.

The proof establishes all finite sizes; these checks are controls, not the source of the general quantifiers. results.json is the matching exact output.
