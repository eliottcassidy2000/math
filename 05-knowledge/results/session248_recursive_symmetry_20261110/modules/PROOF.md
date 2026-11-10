# Recursive symmetry as a tree code, an odometer action, and a surviving charge

**Session248, 2026-11-10. Status: PROVED statements below; FINITE-EXACT companion probes.**
All automorphism claims explicitly refer to the marked hierarchy unless stated otherwise. No claim about LRC(14), Chan–Pak, or the full unmarked flip metagraph follows.

## 0. Inheritance, portfolio, and the change of object

The anchor is recursive emergence and destruction of tournament symmetry. The niche is the repository's physical odometer/carry wall; the wildcard is error correction on hierarchical block flips. The working board is (i) local rotational freedom, (ii) occupied prefix support, (iii) coherent flip cost, (iv) binary observer types, and (v) arithmetic charge. The surprising outcome is that (iv) can eliminate all of (i), while (v) survives at exponential strength.

Closest mechanisms and provenance:

* `01-canon/theorems/THM-2788-physical-modular-odometer-versus-heisenberg-bockstein-extension.md` distinguishes actual successor from carry-suppressed successor. Their depth-two join is a wreath group. Here both act as actual automorphisms of an explicitly constructed tournament, not just as permutations on an unrelated digit carrier.
* `07-reflections/keller-fixed-linear-full-wreath-all-level-primitivity-codex-20260816.md`, reflecting THM-3535, emphasizes **supported** local permutations rather than just their characters. We use occupied prefix support to measure exactly how a defect removes local rotations. Its symmetric-group/inertia theorem is NOT imported into tournaments: a tournament cannot have a transposition automorphism.
* `07-reflections/the-merged-metagraph-is-a-three-sorted-recursive-incidence-object-codex-S2.md` keeps presentations and incidence instances before passing to merged nodes. Accordingly our metric theorem lives in the labelled physical arc cube; we do not identify it with the old complement-line metagraph or an unmarked isomorphism quotient.
* `05-knowledge/results/session247_metagraph_weave_20261110/marked/PROOF.md`, Theorem 1, supplies the fixed-divisor interface theorem. We connect its charge to tree symmetry by a free action on ordered path covers.

Canonical hostile: a new automorphism of a defect tournament need not belong to the old group. Thus “symmetry lost” below means loss from the inherited marked group, not an unproved formula for the new full automorphism group. Another hostile is deleting the hierarchy: differently marked decompositions may become isomorphic, contracting metric distances. The least-used sidecar is the address of each local rotation and of each flipped rectangle.

## 1. The cyclic tower and its exact supported action

Let `A={0,1,2}`, let the leaves of a depth-d ternary tree be `A^d`, and let `R_d` orient `x -> y` iff, at their first unequal digit, `y_j-x_j=1 mod 3`. Put `R_0` equal to one vertex. Thus

    R_(d+1) = C3[R_d,R_d,R_d].

Every prefix cylinder is a module. Let `G_d` be the tournament automorphisms preserving the entire prefix hierarchy, allowing permutation of the cylinders. Its action at every internal node is a cyclic rotation of the three child cylinders. Conversely independent choices of these rotations define such automorphisms. Consequently

    G_(d+1) = G_d^3 semidirect C3,
    |G_d| = 3^M_d,              M_d=(3^d-1)/2.                 (1)

This is an equality for the marked group, and a subgroup of the unmarked automorphism group; the latter equality is unnecessary here. The marked group is transitive on leaves. There are d successive local freedoms along each leaf address, although there are M_d local freedoms over the whole tree. We call these **symmetry depth** and **total symmetry exponent**, respectively. They are different observables.

### Theorem 1: occupied-prefix law

For a set S of individually marked leaves, let `a(S)` be the number of internal tree nodes with a descendant in S. Then its pointwise stabilizer satisfies

    |(G_d)_(S)| = 3^(M_d-a(S)).                              (2)

**Proof.** At a node containing a fixed leaf, the cyclic action must fix the child containing that leaf, hence is the identity. Recursively the same holds at every occupied internal node. At an unoccupied node every cyclic choice remains free. Encoding a tree automorphism by its local rotations is bijective, so exactly a(S) of the M_d choices are forced to zero. QED.

For one leaf the index is `3^d`. For two distinct leaves whose common prefix has length l, `a({x,y})=2d-l-1`: each path has d internal nodes and their overlap has l+1. This includes the branching node itself.

### Corollary 1: single-arc defect law

Reverse one arc e of `R_d` to obtain `R_d^e`. If its endpoints have common-prefix length l, then

    |G_d intersect Aut(R_d^e)| = 3^(M_d-(2d-l-1)),
    |G_d.e| = 3^(2d-l-1).                                  (3)

Indeed an old automorphism preserves the modified tournament iff it preserves the single disagreement pair. Since it preserves the original orientation, it cannot exchange that pair's endpoints. Apply (2). This is a statement about the common/inherited symmetry group; the modified tournament might have automorphisms outside G_d.

All arcs at fixed l form one orbit: rotations move the branching prefix, the oriented pair of its children, and the two suffix addresses independently. Their number also equals

    3^l * 3 * (3^(d-l-1))^2 = 3^(2d-l-1).

Thus a microscopic one-arc defect distinguishes scales by the logarithmic group loss `2d-l-1`. A bottom defect costs d powers of 3 because it pins one branch at every ancestor, not merely the local bottom C3. A top-level defect costs `2d-1` powers.

For d=2, `|G|=81`: the 27 cross-top-block arcs have surviving group 3; the 9 bottom arcs have surviving group 9. For d=3, the three arc orbits have sizes 243,81,27 and surviving groups `3^8,3^9,3^10`, respectively.

## 2. A weighted repetition code inside the physical flip cube

For every internal node v and unordered pair `{a,b}` of its children, let `E_(v,ab)` consist of all leaf pairs crossing those two child cylinders. These rectangles partition the entire set of leaf pairs: every pair has exactly one first branching node. If v is at depth l, then

    w_v = |E_(v,ab)| = 9^(d-l-1).                           (4)

Allow the three quotient arcs at each node to vary independently. For a bit vector b on the `3M_d` rectangles, let `T_b` be the tournament whose arcs are constant on each rectangle, using the fixed child-order convention to read a bit. Equivalently this permits arbitrary three-vertex quotient tournaments at every node.

### Theorem 2: weighted code isometry and decoding

For any b,c, the labelled arc-flip distance is exactly

    d_flip(T_b,T_c) = sum_(v,ab) w_v |b_(v,ab)-c_(v,ab)|.    (5)

Relative to any fixed codeword, the orientation vectors form a binary linear repetition code of dimension `3M_d` and length `binom(3^d,2)`. Each coordinate is repeated w_v times on disjoint support.

**Proof.** Arc distance in the labelled cube is Hamming distance. The two tournaments disagree on the whole rectangle exactly when the corresponding bits differ. The rectangles are disjoint, so their sizes add without interaction. This proves both assertions. QED.

A weighted edge of cost w lifts to a geodesic of w individual arc reversals. For w>1, intermediate tournaments usually leave the prescribed module locus. The code locus is therefore NOT a graph subcube with unit physical edges at all scales; it is an isometric weighted Hamming vertex set.

Majority on each rectangle is a nearest-codeword retraction from the full labelled arc cube. All w_v are odd, so no ties occur. It exactly recovers a codeword whenever fewer than w_v/2 arcs are corrupted in every rectangle separately. This is an actual error-correction map; the missing sidecar is the chosen hierarchy. The full code has minimum distance 1, because bottom rectangles have size 1. If all levels below depth L are frozen and only levels `0,...,L-1` vary, its minimum distance is `9^(d-L)`. No claim of global robustness is licensed by large top-level blocks alone.

For n=9 the code has dimension12 in the 36-dimensional physical arc cube, with weights nine copies of1 and three copies of9. For general d the weight identity

    sum_(l=0)^(d-1) 3^(l+1)9^(d-l-1) = binom(3^d,2)

is the exact scale decomposition of all physical arcs.

**Scale mismatch.** A one-arc defect in a rectangle of weight w is distance1 from the original coherent state, and distance w-1 from the coherent state with that quotient arc reversed. Completing the proposed coarse flip costs w-1 further reversals; simply undoing the defect costs1. Combining (3),(4), a top-level defect can remove `2d-1` inherited symmetry factors at physical cost1, even though the corresponding opposite coherent quotient costs `9^(d-1)` flips. This is not a claim that restoring *some* symmetry costs w.

### Quotient coarsening and its information invoice

Truncating the tree at depth L replaces each leaf cylinder of size `3^(d-L)` by one supervertex. On the hierarchical code locus, this forgets all deeper bits and produces an order-`3^L` tournament in the same family. Distances contributed by retained bits obey

    d_retained(T_b,T_c) = 9^(d-L) d_flip(coarse_L(b),coarse_L(c)). (6)

The exact discarded distance is the sum of the deeper weighted bit differences. Thus coarsening preserves the retained orientation predicates and marked ancestor relation; it loses internal orientations, their symmetry, and potentially most distance. The full fine-bit sidecar reconstructs both tournament and distance. Majority also defines a noisy version, but no global Lipschitz claim is made: crossing a majority threshold can reverse a coarse bit in one physical step.

## 3. A genuine 3-adic map, with the carry retained

Identify a word `(x_0,...,x_(d-1))` with `sum_j x_j 3^j mod 3^d`, so the FIRST tree digit is the LEAST significant arithmetic digit. Then

    x -> y iff (y-x)/3^v_3(y-x) = 1 mod3.                  (7)

Every translation `x -> x+a` and every multiplication `x -> ux` with `u=1 mod3` preserves (7). These are explicit elements of G_d. The actual odometer `O:x -> x+1` has order `3^d`; suppressing its carry gives a different element of order3.

At d=2, write `x=v+3w`. The transformations

    O:x->x+1,     Y:(v,w)->(v+1,w),     X:(v,w)->(v,w+v)

all preserve R_2. They are exactly the p=3 digit actions studied in THM-2788. Its distinction between the order27 physical metacyclic subgroup `<O,X>` and the order27 carry-suppressed Heisenberg subgroup `<Y,X>` now occurs inside the order81 marked automorphism group of a nine-vertex tournament. One can see the full join directly: `O Y^-1` is a rotation of the second digit supported in one first-digit cylinder; conjugating by Y gives all three supported bottom rotations, and Y supplies the top rotation. Hence `<O,Y>=G_2`. This uses actual supports, not a commutator analogy.

The inverse limit of these hierarchical groups acts on `Z_3`; any two distinct 3-adic integers still have a first differing digit, so (7) defines a tournament on this boundary. Truncation is an honest quotient on distinct residue cylinders, not a tournament map on pairs that collapse to one point. The group filtration by fixing the first L digits is an exact symmetry-resolution filtration. Distances `3^-l` between leaves encode shared ancestry; they do not by themselves encode orientation, which also requires the first nonzero digit's sign.

Source/target/map: ternary addresses -> residues mod `3^d` -> `Z_3`; preserve first-difference orientation and hierarchy; lose none at finite depth, and lose fine digits under reduction; sidecar is digit order and carry. Cheap hostile: omitting carry changes an element of order9 into order3 already at d=2.

## 4. Observer automata and the two rigid binary types

Color every leaf with one of q colors. Consider colored isomorphism under G_d, not arbitrary unmarked tournament isomorphism. A colored depth-(d+1) type is an ordered triple of depth-d types modulo cyclic rotation. It admits a top rotation precisely when all three child types coincide.

### Theorem 3: exact type and stabilizer recursion

Let `P_(d,q)(z)` sum `z^s` over colored G_d-orbits, where the stabilizer order of a representative is `3^s`. Then

    P_(0,q)(z)=q,
    P_(d+1,q)(z)=(P_(d,q)(z)^3-P_(d,q)(z^3))/3
                         + z P_(d,q)(z^3).               (8)

In particular the number of types `b_(d,q)=P_(d,q)(1)` satisfies

    b_(d+1,q)=(b_(d,q)^3+2b_(d,q))/3.                     (9)

The number `r_(d,q)=P_(d,q)(0)` of rigid types satisfies

    r_(0,q)=q,       r_(d+1,q)=(r_(d,q)^3-r_(d,q))/3.       (10)

**Proof.** Nonconstant triples have rotation orbits of size3 and no extra top symmetry. Their stabilizer exponents add. Subtract the diagonal triples from the cube and divide by3, giving the first term of (8). A diagonal triple with child exponent s has exponent `3s+1`, giving the second term. This proves (8), and specializing gives (9),(10). QED.

This is a bottom-up tree-state automaton: compute child types, canonicalize their cyclic triple, and add one stabilizer exponent only on the equal-type diagonal. The state alphabet grows with depth; we do not call this a single fixed finite-state automaton for all depths.

**Binary observer surprise.** With q=2, `r_(d,2)=2` at EVERY depth. At depth0 use the two colors A_0,B_0. Recursively take

    A_(d+1)=[A_d,A_d,B_d],
    B_(d+1)=[A_d,B_d,B_d].                                (11)

These are the two rigid colored types; color complementation exchanges them. All their children are rigid, and their child triple is nonconstant. Thus a binary root neighborhood can kill the entire inherited G_d, even while its size is doubly exponential in depth:

    # rigid binary labelings = 2 |G_d| = 2*3^((3^d-1)/2).

For a uniformly random binary coloring, the exact rigidity probability is

    2*3^((n-1)/2)/2^n,       n=3^d,                       (12)

which tends to zero exponentially in n. Deliberately chosen one-bit observations achieve what a random one usually fails to do.

Realize a coloring by adding a marked root x, with `x -> v` iff v has color1. Its root-fixing, hierarchy-preserving automorphism group is exactly the color stabilizer. Formula(11) therefore supplies one-root extensions with trivial inherited marked symmetry, for every d. There might be an unmarked automorphism moving x or failing to preserve the hierarchy; we do not assert full rigidity of the underlying extension.

The exact n=9 stabilizer-orbit polynomial is

    P_(2,2)(z)=2+10z+8z^2+2z^3+2z^4.

There are24 binary types, only2 rigid. At d=3 there are4624 types, still only2 rigid. Unlike a scalar symmetry count, (8) retains the full distribution and an explicit recursive witness.

## 5. Symmetry can disappear while its arithmetic charge persists

For any nonempty tournament M, let `pc_M(c)` count its partitions into c nonempty directed paths, each path retaining its vertex order but the set of paths unordered. Put

    delta(M)=gcd_(1<=c<=|M|) c! pc_M(c).

### Theorem 4: automorphism charge

For any subgroup K of Aut(M),

    |K| divides delta(M).                                (13)

**Proof.** K acts freely on ordered c-tuples of directed paths covering all vertices. A stabilizer must fix each path in its assigned tuple position and fix each vertex's position along it, hence fixes every vertex. Therefore every orbit has size |K|, which divides the number `c!pc_M(c)` for every c. Taking the gcd proves the assertion. QED.

In particular, `3^M_d | delta(R_d)`, with no need to prove `G_d=Aut(R_d)`. No equality for this charge is asserted.

### Corollary 2: recursive symmetry-charge persistence through arbitrary interfaces

Take k disjoint copies of R_d, orient every intermodule pair uniformly according to any quotient tournament, and add q further vertices with arbitrary arcs among themselves and to individual core vertices. Then

    3^[M_d max(0,k-2q)] divides H(T).                    (14)

**Proof.** Apply Theorem1 of Session247's marked/PROOF.md with each module charge divisible by `3^M_d`. For completeness, partition Hamilton paths into skeletons of maximal core-module runs and individual interface vertices. At most2q modules have a run adjacent to an interface vertex. Every untouched module contributes a number `c!pc(c)` divisible by `3^M_d`; these independent factors give the claimed common divisor in every skeleton contribution. Sum over skeletons. QED.

Apply the SAME binary pattern(11) to every module's adjacency to one added marked root. This kills all of the product subgroup `G_d^k` pointwise on the modules: every factor's coloring stabilizer is trivial. Yet for k>=3,

    3^[M_d(k-2)] divides H(T).

Thus large divisibility is not equivalent to surviving visible symmetry. It can be the arithmetic memory of symmetry that existed before the interface was attached. For d=2,k=3,q=1 this gives a 28-vertex example: no nontrivial inherited internal symmetry survives the chosen root neighborhoods, while `81 | H(T)`. This conclusion is independent of the three-module quotient orientation. The root marked statement is essential; no full-rigidity claim is made.

The mechanism differs from the occupied-prefix law: a pointwise marked leaf constrains only its d ancestors, while an external root broadcasts a binary constraint to every leaf. One external vertex can therefore erase all inherited local rotations, but a Hamilton path sees it through only two boundary incidences. Group-theoretic visibility and path-count divisibility have different interface budgets.

## 6. Reproduction and onward questions

Run from repository root:

    python3 05-knowledge/results/session248_recursive_symmetry_20261110/modules/probe.py

The committed `probe.json` records exact enumeration of G_d, arc stabilizers and all binary colorings for d=1,2; compares these to an independent recursive polynomial calculation; checks orbit-stabilizer reconstruction of all `2^(3^d)` colorings through d=5; and checks the n=9 rectangle partition. The experiments explore the new construction, not an audit of the inherited literature. Positive controls are the uniform colorings with maximal stabilizer; hostile controls include the minimum-weight1 bottom rectangles and the distinct top/bottom arc losses.

Next directions (OPEN, not implications of the proofs):

1. Determine the full automorphism group of one-defect R_d, including new symmetries outside the inherited group. The occupied-prefix theorem already gives the exact common subgroup and a scale parameter for organizing this classification.
2. Find a canonical modular-decomposition marking for larger families, then ask which weighted code patches remain isometric after forgetting labels. Equation(5) alone does not answer this quotient problem.
3. Extend (8) to observers with restricted support or a prescribed number of root outneighbors. A second variable marking color weight would distinguish economical from arbitrary root neighborhoods.
4. Compute whether delta(R_d) has additional arithmetic beyond `3^M_d`; symmetry guarantees divisibility, not equality. Combine that charge with Session247's valuation budget rather than assume the local group exhausts all obstructions.
5. Use the prefix filtration and rectangle code to compare two different persistence profiles: surviving symmetry under sparse physical defects, and robust recoverability of coarse module orientations. They provably measure different things already in(3),(4).
