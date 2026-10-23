# Hamiltonian observer bundles and native port compilation

Status: proved statements on the stated carriers; finite-exact receipts for orders 3–5 and one order-6 hostile. These mechanisms do not classify all tournament metagraphs or establish a compact encoding for arbitrary Hamilton path families.

## 1. The observer is structure, not a canonical choice

Let H(T) be the number of directed Hamiltonian vertex orders of a tournament T. An observer object is (T,P), where P=(v_0,...,v_(n-1)) is one such order. Isomorphisms must transport the ordered path, not merely its unordered edge set.

Every automorphism fixing P fixes every vertex, hence is the identity. Therefore Aut(T) acts freely on its Hamilton paths, and the fibre of observer-isomorphism classes over [T] has size

    q(T) = H(T)/|Aut(T)|.

Index vertices in their path order. Adjacent arcs are then fixed forward, and every other arc is free. Thus observer objects up to isomorphism are in exact bijection with the binary cube of dimension m=(n-1)(n-2)/2, with no further quotient. Consequently

    sum_[T] q(T) = 2^m,
    sum_labelled_T H(T) = n! 2^m.

The familiar existence of a Hamilton path supplies at least one observer, not a preferred observer. In particular, if Aut(T) is nontrivial there cannot be an automorphism-equivariant choice of a single ordered Hamilton path: equivariance would force that path to be fixed by every automorphism, contradicting freeness. A canonical isomorphism-class representative and a natural actual-owner choice are different requests.

Uniform observer masks weight unlabelled tournaments in proportion to q(T). They are not uniform on tournament classes followed by a uniform choice of a path. The distinction matters when restricting the graph/tournament coupling operator.

## 2. Two involutions and exact colour budgets

Use the backward-edge graph in observer order: its consecutive pairs are absent. This is the ordered-graph correspondence restricted to the Hamilton observer cube.

Let F flip all nonpath tournament arcs. It complements precisely the free graph bits. Let R send (T,P) to (T^op,reverse(P)), then normalize in the new observer order. It reflects the free bits through i -> n-1-i. It DOES NOT complement them: full converse at fixed order complements once, and reversing the order complements again.

F and R commute. For n>=3, F is fixed-point-free, and so is FR, because the free pair {0,n-1} is fixed by reflection and would have to equal its complement. Reflection has f=floor((n-1)/2) fixed free pairs, hence

    # Fix(R) = 2^((m+f)/2).

If blue means both endpoints are reflection-fixed and black means both are not, the exact edge budgets are

    blue = 2^((m+f)/2 - 1),
    black = (2^m - 2^((m+f)/2))/2.

The F,R orbit count is (2^m+2^((m+f)/2))/4. These are direct involution counts, not new tournament-class enumeration formulas.

After forgetting the observer and merging converse classes, a self-converse node has q(T) incident observer ports, and a non-self-converse merged node has 2q(T). F pairs ports; loops consume two ports. Thus its graph-theoretic degree equals its port count when loops count twice. A reflection-fixed observer necessarily lies over a self-converse tournament, but self-converseness of the underlying tournament alone does not imply every observer is reflection-fixed. Blue-only nodes require the stronger fibrewise condition.

## 3. Three different native questions

1. **Transport a known actual path through a complete tournament isomorphism.** No new Hamilton test is needed: its image is automatically a genuine path.
2. **Require partial indexed owners on a complete isomorphism transporter.** Each requirement is a point-transporter coset; a nonempty joint solution remains a coset. Mapping an entire specified ordered spanning P to Q determines at most one map.
3. **Assemble candidate actual paths from internal charts not yet known to preserve the whole tournament.** Boundary arc predicates need not be cosets. This is where a branch-aware compiler is needed.

The third problem must not be used to contradict the first two. Conversely, proving a compiler exact for the first two does not prove it complete for arbitrary path assembly.

## 4. Small native union-of-cosets hostile

Take two named directed 3-cycles A and B. Put A_u -> B_v exactly when v-u mod3 lies in {1,2}. The resulting actual tournament is almost regular: A scores are 3 and B scores are 2.

For charts (a,b) in C3 x C3, consider the literal candidate path

    A_a, A_(a+1), A_(a+2), B_b, B_(b+1), B_(b+2).

All internal arcs are genuine. The single boundary arc is genuine exactly when

    b-a belongs to {0,1} mod3.

There are six distinct genuine paths in this family. The valid chart relation is the disjoint union of two diagonal-C3 cosets; it cannot be one coset because 6 does not divide 9. Explicitly, the valid triples (0,1),(0,0),(0,1) have heap x-y+z=(0,2), which is forbidden.

The whole-block stabilizer, in contrast, is the single diagonal subgroup a=b. Its action has three lifts per difference value. Thus whole-block transport, native path assembly, and the retained simultaneous phase are visible separately in the smallest two-triple example. The receipt also enumerates ALL Hamilton paths of this six-vertex tournament, but does not equate that full set with the six-path assembly family.

Converse preserves the family by reversing each literal path in the converse tournament. Reusing the forward local ordering in the converse without reversing it would be incorrect.

## 5. A compact positive all-size native family

Let p be odd. In each of t named cells use a circulant tournament on Z/p containing every arc u -> u+1; its standard cyclic cuts are genuine internal Hamilton paths. Between cells i,j use any actual circulant block

    (i,u) -> (j,v) iff v-u in D_ij,

with the opposite direction its complement. Fix a cell visitation order and concatenate one rotated internal path in each cell. The boundary from cell i to its successor j is valid exactly when

    a_j-a_i in E_ij := D_ij-1.

Assuming all D_ij on this visitation path are nonempty, the number of these literal Hamilton paths is exactly

    p * product_boundary |D_ij|.

Proof: choose the initial phase (p possibilities), then choose each boundary difference independently and reconstruct all later phases. Distinct chart tuples give distinct ordered spanning paths. Empty D gives no such path. This counts a specified, genuinely native family, not all paths of the tournament. It works for arbitrarily many cells and does not require enumerating the host's Hamilton paths.

More generally, take any fixed spanning template word in the named cells, with each actual local vertex appearing once, and shift all occurrences in cell i by one shared a_i. Within-cell arc tests are shift-invariant and must first pass. Every cross-cell consecutive pair then imposes a difference constraint a_j-a_i in an explicit translated D_ij. Parallel constraints must share the same two variables. For any choice of one difference per boundary occurrence, a common phase assignment exists exactly when signed sums vanish on every cycle of the boundary-incidence multigraph (including parallel-edge cycles). A consistent choice has p^c lifts, where c counts its connected components, including isolated cells. This follows from spanning-forest integration; no field division or primality is needed. For a single spanning path template visiting every cell, this incidence graph is connected, so c=1; the more general p^c statement is useful for partial path forests before joining them. Singleton owner marks fix component phases or reject inconsistent lifts.

This is a native version of residual-cycle potentials: differences are candidate arc branches; common phases are retained kernel lifts. Dropping the cycle tests creates spurious paths, while dropping the p^c lifts loses actual owners.

## 6. When is a native boundary relation a coset?

For a nonempty finite subset S of a group, S is a left coset of a subgroup if and only if S is closed under x y^-1 z. Necessity is direct. For sufficiency pick r in S; r^-1 S contains identity and is closed under x y^-1 z, hence inverses and products, so it is a subgroup.

For a cyclic boundary R_E={(a,b):b-a in E}, this says R_E is a coset iff E is a coset of a subgroup of Z/p. For prime p and nonempty proper E it is a coset iff |E|=1. Thus multivalued proper native arc palettes typically leave the incoming Mal'tsev regime even though exact whole-block comparisons remain in it.

A finite branch-aware extension is nevertheless exact. Write each native relation as a DISJOINT union of k_e cosets in the same shared chart variables. For each branch tuple, intersect the selected cosets together with the original complete transporter constraints and markings. Each nonempty intersection is a coset. Different branch tuples have disjoint full assignments because the local branches were disjoint. Therefore there are at most product_e k_e branches, each handled by the original shared-frame compiler, with exact kernel lifts and no double counting. This is an explicit parameterized sufficient criterion, not a claim that branch counts are uniformly small for all tournaments.

For the cyclic difference family, taking branches indexed by E gives k_e=|E| and cycle-potential tests replace full inventory enumeration. For general actual groups an arc predicate is a union of point-transporter cosets; its useful compactness depends on orbit structure and may fail.

## 7. Guards must be checked on the actual lifts

Converse is implemented on the actual tournament and literal reversed path, not by a guessed numerical symmetry. A mod8 invariant is retained as its exact original labelled predicate. It may be pushed to a quotient only after proving it is constant on each discarded kernel fibre. Linear congruence technology does not automatically cover electrical or representative-dependent quantities.

For illustration only, impose the synthetic guard (a-b) mod8 in {0,7}, using representatives a,b in {0,1,2}, on the order-6 example. It leaves five of the six valid chart pairs. The difference class b-a=1 mod3 has three genuine lifts, but this guard keeps two and removes (2,0). Thus a representative-sensitive mod8 guard can destroy constant kernel multiplicity. This is not asserted to be the repository's electrical guard; it is a reproducible reason to keep that guard's actual definition and clocks.

## Provenance and scope

Read at base commit 8b71b746f7a9d90378c702e01608207227b1f83e:
- cyclic_attachment_lattices_20261023/NONABELIAN.md: whole-frame cosets, Mal'tsev feasibility, complete prefix banks, explicit-domain complexity.
- attachment_lattices_square_sum_20261023/KERNELS.md: independent kernel factorization inside allowed projections, lifted owners.
- attachment_lattices_square_sum_20261023/README.md: observer correspondence, shared clocks, retained all-size boundaries.

The new synthesis is the precise separation of the three native problems, a genuine six-vertex branch hostile, and a compact all-size cyclic-cut family with cycle-compatible shared phases. The elementary orbit and heap facts are standard recovered mechanisms, not claimed as novel general group theory. The restricted graph/tournament spectral results belong to the separate Session43 coupling packet; observer existence does not preserve the old 8/9 gap after conditioning.
