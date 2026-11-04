# A physical six-port parity collision and a common integer-observability theorem

Session 190, 4 November 2026. **PROVED elementary results; FINITE-EXACT replay.** This packet gives a concrete construction beyond the inherited formal six-port kernel. The integer-lattice language is standard Smith-normal-form algebra, not a priority claim. It provides a common analysis of two response operators, **not** a map from planar graphs to tournaments or even graphs.

## 1. Inheritance and exact scope

The input `session190-planar/PROOF.md`, §§2–6, established the Hamiltonian partial-forest gluing formula, realizable unit probes, unique completion through five ports, and the six-port determinant −2 block. These are specializations of established matching-connectivity / meander methods. The input `session190-bijection/NOTE.md`, §§2–5, records the inherited signed even-graph Walsh basis, its orthogonality, and the square-entry probabilistic coupling. We do not rebrand these as new general theories.

The concrete additions here are:

1. Two **connected simple planar disk graphs** whose entire profile parities differ, while every permitted exterior gives the same Hamiltonian-cycle parity.
2. A two-exterior, divided-count identity that exactly measures a named interior statistic and recovers the missing bit.
3. An explicit comparison of integral response precision: modulo four suffices for the six-port operator; the order-three signed Walsh operator already needs modulo eight.

A planar patch has ordered boundary R=(r0,...,r5), nonempty interior, and a crossing-free disk drawing. Two patches used in a gluing intersect in exactly R, lie in opposite disks on the sphere, and have disjoint edge sets. Every exterior below has nonempty interior. The function h counts **unoriented spanning Hamiltonian cycles**, not arbitrary cycles. In particular, rim cycles that omit mandatory interior vertices do not count.

A partial Hamiltonian forest is a spanning edge subset with degree two at every interior vertex, degrees at most two on R, and no cycles. Its positive components pair the degree-one ports. Write f_X(d,M) for its count. For the all-ones degree pattern use this ordered matching basis:

    M0 = 01|23|45
    M1 = 01|25|34
    M2 = 03|12|45
    M3 = 05|12|34
    M4 = 05|14|23.

Let x_i=f_X((1,1,1,1,1,1),Mi). These are coordinates of X; X may also have states with other degree patterns.

## 2. The two disk gadgets

Start with six interior rim vertices v0,...,v5 forming a cycle in that order. For each i, attach boundary port ri through a private interior vertex ai, using edges ri–ai and ai–vi.

- **Ring patch B:** these 18 vertices and 18 edges.
- **Hub patch A:** add an interior vertex c and the six edges c–vi. It has 19 vertices and 24 edges.

Both graphs are connected and simple. An explicit straight-line drawing is obtained from

    v0,...,v5 = (2,0), (1,2), (−1,2), (−2,0), (−1,−2), (1,−2),
    ai = 2vi,   ri = 3vi,   c = (0,0).

Take the boundary to be the hexagon through the ri. The rim and radial segments have no crossings and all nonboundary vertices are strictly inside. No patch edge joins two boundary vertices, so these gadgets can be placed opposite any permitted exterior without duplicating a boundary edge.

**Proposition 1 (entire profiles).** In the ordered all-ones basis,

    f_A = (3,1,1,3,1),       f_B = (1,0,0,1,0),

and both profiles vanish at every other boundary degree pattern.

**Proof.** Every ai is an interior vertex of graph degree two. Both of its incident edges must be selected. Therefore all six ports have degree one, and every vi already has one selected spoke edge.

For B, each vi must select exactly one rim edge. The two perfect matchings of the rim cycle give M0 and M3, once each.

For A, the mandatory interior hub has selected degree two, so it selects two distinct rim neighbors. Those two vi select no rim edge; the other four vi must be matched by rim edges. This is possible precisely when the two chosen vertices have opposite parity around the rim. If they are adjacent, the remaining four vertices form a path with a unique perfect matching. The six adjacent pairs contribute three copies of M0 and three of M3. If the chosen vertices are opposite, the remaining vertices form two disjoint edges. The three opposite pairs 25, 03, 14 give M1, M2, M4, respectively. There are no other choices. All selected components are paths, because a rim vertex already has a selected spoke and can have only one additional selected edge. ∎

The subdivisions ai are essential. Without them, the leaves ri need not be used in a partial forest; one could not assert that the displayed five coordinates are the full profile.

## 3. Universal parity collision, and a physical information-recovery identity

The all-ones matching-connectivity matrix is

    K = [0 0 0 1 0]
        [0 0 1 0 1]
        [0 1 0 0 1]
        [1 0 0 0 0]
        [0 1 1 0 0].

Indeed, opposite-side matchings form a single alternating cycle exactly for these entries. The pair {M0,M3} forms a two-by-two permutation block; {M1,M2,M4} forms the adjacency matrix of a triangle. Direct multiplication gives

    K(f_A − f_B) = K(2,1,1,2,1)^t = 2(1,1,1,1,1)^t.       (1)

**Theorem 2 (physical parity collision and extraction).** For every permitted six-port exterior X,

    h(X ∪ A) − h(X ∪ B) = 2 S(X),
    S(X) = sum_(i=0)^4 f_X((1,1,1,1,1,1),Mi).               (2)

Consequently A and B have identical Hamiltonian-cycle parity in every such exterior, even though their profile vectors modulo two differ by

    k = (0,1,1,0,1).

**Proof.** A spanning Hamiltonian cycle restricts to a partial forest on each side. A restriction cannot contain a cycle: such a cycle would already use both incident edges of each of its vertices and could not cover the other side's nonempty interior. Degrees must add to two at each port. Since every profile state of A and B has degree one at every port, only the all-ones states of X contribute. Contracting each path component gives the two matchings; the global cycle is connected precisely when their union is connected. Restriction and union are inverse maps, with no orientation factor. Thus

    h(X ∪ A) = x^t K f_A,       h(X ∪ B) = x^t K f_B.

Equation (1) proves (2). The parity assertion and profile difference follow. ∎

This is an **actual pair of connected graph gadgets**, not a formal sum or difference mistakenly asserted to be one graph. The algebraic difference is used only to compare two physical measurements.

For a unit exterior P_Mj consisting of three length-two paths with matching Mj and distinct private interior vertices, the counts are:

| exterior | h(P ∪ A) | h(P ∪ B) |
|---|---:|---:|
| M0 | 3 | 1 |
| M1 | 2 | 0 |
| M2 | 2 | 0 |
| M3 | 3 | 1 |
| M4 | 2 | 0 |

Every pairing is noncrossing, so each unit exterior has a disk drawing. In particular, M1 gives two spanning Hamiltonian cycles after inserting A and none after inserting B. Equality of every exterior **parity** therefore does not preserve Hamiltonian existence. This is stronger than the usual observation that a nonzero count can be even: the entire universal parity-response function agrees for these two patches.

For exterior copies of B and A, respectively, the A/B counts are 6/2 and 24/6. They verify (2) with S(B)=2 and S(A)=9.

**Corollary 3 (sharp port threshold).** Within the stated disk model, six is the smallest number of ports admitting two physical patches with different profile parities but the same Hamiltonian parity in every exterior.

**Proof.** Through five ports the inherited matching-connectivity operator is a permutation matrix, and every live coordinate has a realizable unit exterior. All exterior parities therefore recover all profile parities. Theorem 2 supplies a connected example at six. ∎

### 3.1 A separate nonplanar-exterior corollary

The **parity equality** extends to arbitrary finite simple graph exteriors X, without a disk embedding, provided X intersects A or B exactly in the six named ports, has a nonempty interior, and has a disjoint edge set. The clean identity (2) is retained only for planar exteriors.

To prove the extension, an unrestricted all-ones exterior state is one of the 15 perfect matchings N on the six ports. The even part 2(e_M0+e_M3) of f_A−f_B makes no difference modulo two. It remains to show that N has a connected union with an even number of M1,M2,M4.

Partition the ports into E={0,2,4} and O={1,3,5}, in those orders. The matchings M1,M2,M4 are precisely the three odd-permutation bijections E->O. If N is also bipartite, a union with such a matching is connected precisely when their relative permutation is a three-cycle. There are two such odd matching partners when N is odd, and zero when N is even. If N is not bipartite, it consists of one E–E edge, one O–O edge, and one E–O edge. Its union with a bipartite matching is disconnected precisely when that single E–O edge is shared. Exactly one of M1,M2,M4 shares that edge, so exactly two give a connected union. In every case the count is even.

The full integer count differences against the 15 unit matchings are 2 for five matchings, 4 for seven, and 6 for three. They are all positive and even, but not constant. This proves universal parity equality for arbitrary exteriors and shows why the planar extraction identity must not be carried over unchanged. The nonplanar unit exteriors are still ordinary simple graphs consisting of three internally disjoint two-edge paths; only their prescribed cyclic disk drawing fails.

## 4. What exactly is recovered, and at what precision?

Let y=Kx be the five counts against the unit matching exteriors. Over the integers,

    x0 = y3,      x3 = y0,
    x1 = (y2+y4−y1)/2,
    x2 = (y1+y4−y2)/2,
    x4 = (y1+y2−y4)/2.                                    (3)

Hence the exact image lattice is

    K Z^5 = { y in Z^5 : y1+y2+y4 is even }.                (4)

The Smith invariants are (1,1,1,1,2), and the cokernel is Z/2Z. For nonnegative formal profiles, (4) is supplemented by y_i>=0 and the three triangle inequalities among y1,y2,y4. These are conditions for a nonnegative **formal** profile, not a claim that every such profile is realized by one graph.

Over F2 the response has one-dimensional kernel spanned by k. If

    s = S(X) mod 2,

then the missing coordinates are recovered from ordinary response parities by

    t = s + y0 + y3,
    x1 = t + y1,     x2 = t + y2,     x4 = t + y4,          (5)

with all operations in F2. Equation (2) provides this one extra bit physically:

    s = (h(X ∪ A) − h(X ∪ B))/2 mod 2.                     (6)

Knowing the two physical counts modulo four suffices to compute (6), because their difference is guaranteed even. This does **not** say that upgrading one arbitrary unit-probe count to modulo four suffices. The added observable is the specified divided difference.

All other six-port degree sectors have at most four degree-one ports; their connectivity blocks are permutation blocks. Thus their profile parities are already recovered by unit-probe parities. Equations (5)–(6) recover the only missing bit of the full 365-coordinate six-port profile. One additional bit is necessary because the kernel is one-dimensional, and the construction supplies one sufficient bit. Equivalently, the full ordinary response modulo four recovers the full profile modulo two by (3).

## 5. The common integer-observability theorem

Let R: Z^s -> Z^t be any integer response operator, and let p be prime. Write R_p for its reduction modulo p and C=coker R.

**Theorem 4.** There is an exact sequence

    0 -> (ker R)/p(ker R) -> ker R_p -> C[p] -> 0,          (7)

where C[p]={c in C : pc=0}. The last map sends x modulo p, with Rx divisible by p, to the class of Rx/p in C.

**Proof.** The map is well-defined: replacing x by x+pz changes Rx/p by Rz. Its value is killed by p because p(Rx/p)=Rx. Any p-torsion class represented by y has py=Rx for some x, so the map is onto. Its kernel consists precisely of classes with Rx/p=Rz, or x−pz in ker R. Reduction of the integral kernel is injective modulo its own p-multiple: if a kernel element is pz in Z^s, then pRz=0 in the torsion-free group Z^t, so Rz=0. ∎

Thus modular blindness is **not** always just cokernel torsion. An existing rational/integral kernel contributes the first term in (7). When R is injective over Q, that first term vanishes and the modular kernel is canonically C[p]. For a nonsingular square matrix, its dimension equals the number of Smith invariants divisible by p.

**Precision consequence.** Suppose R has full column rank, with positive Smith invariants d1,...,ds. Let a=max_i v_p(di). Then the least uniform exponent m>=1 for which R x modulo p^m determines x modulo p, for every integer x, is

    m = a+1.                                              (8)

To prove sufficiency, use unimodular Smith coordinates and divide each observed d_i x_i by its p-power before reducing modulo p. Unimodular coordinate changes preserve all relevant congruences. For necessity, choose a Smith coordinate with valuation a; changing that domain coordinate by one is invisible modulo p^m for every m<=a, but changes the domain modulo p.

Equation (8) concerns known genuine responses to arbitrary integer vectors. It is not an assertion that each Smith-coordinate test vector is a graph profile. The explicit construction in §2 supplies that missing physical step for the planar parity collision.

## 6. Applying the same theorem to the signed tournament Walsh response

Let W_(T,G)=Phi_G(T) be the integer matrix from the inherited signed even-graph Walsh basis, with tournament types T as rows and automorphism-even graph types G as columns. Here "even" means that every graph automorphism reverses an even number of oriented edges, as defined in `session190-bijection/NOTE.md`. It does not mean Eulerian or even-order automorphism group.

The coordinates of a vector a in this section are **integer coefficients of signed graph Walsh sums**; W a is the resulting integer-valued invariant function evaluated on tournament types. They are not planar path-count coordinates. The two settings share an integer observation mechanism, not their combinatorial objects.

Two elementary facts make the comparison precise.

1. W is square and invertible over Q by the inherited Fourier-basis theorem. Its inverse has denominator dividing 2^N, N=binom(n,2). Indeed, Fourier-expand the indicator of a labelled tournament orbit. The coefficient at a labelled graph support is an integer sum of signs divided by 2^N. Grouping these coefficients into signed graph orbits gives an integer matrix C with W C=2^N I. Consequently every Smith invariant of W is a power of two, at most 2^N.
2. Modulo two, every signed monomial is one. Therefore Phi_G(T) is congruent to the number r_G of labelled copies of G, independently of T. Every row of W has the same reduction. Its empty-graph column is one, so rank_F2 W=1 exactly. If there are q tournament types, the modulo-two kernel has dimension q−1.

Thus the same exact sequence detects dyadic loss in both systems, but it does not identify their kernels. The first planar exceptional block has only one factor two. Even the first nontrivial tournament example has a larger precision obstruction.

At n=3, in row order (transitive, cyclic) and column order (empty, P3),

    W = [1  1]
        [1 −3].

One direct derivation is Phi_P3=x01*x02−x01*x12+x02*x12: its value is 1 on a transitive tournament and −3 on a cyclic tournament. Hence

    Smith(W) = (1,4),
    W Z^2 = { y : y_transitive − y_cyclic is divisible by four },
    a_P3 = (y_transitive − y_cyclic)/4.

The signed Walsh response must be known modulo eight to recover arbitrary coefficient vectors modulo two. Modulo four is insufficient: the coefficient vectors (0,0) and (1,−1) have the same response modulo four, but different parities. By contrast, modulo four suffices for the six-port planar response. This is a concrete obstruction to identifying the integral operators by unimodular coordinate changes, even after appending identity blocks: their nontrivial Smith factors are 2 and 4.

Finally, this is distinct from the information lost by the **nonlinear squaring** Q_(T,G) -> Q_(T,G)^2 in the probabilistic coupling B. Squaring removes signs even over the real numbers and, for instance, identifies converse rows. The linear cokernel theorem applies to W, not directly to that squaring map. Neither the wheel gadget nor this theorem constructs a natural tournament/even-graph bijection.

## 7. Reproduction and independently inspectable certificates

Run from the workspace root:

    python3 session190-interface/verify.py > session190-interface/RESULTS.json
    python3 -O session190-interface/verify.py > session190-interface/RESULTS.optimized.json
    cmp session190-interface/RESULTS.json session190-interface/RESULTS.optimized.json

The script uses only the standard library and imports no repository producers. Runtime checks are explicit exceptions, not assertions disabled by optimization. It performs:

- Exact degree-constrained edge-subset enumeration of the **full** profiles of A and B, not only a chosen degree sector.
- Direct degree-two edge-subset and connectivity enumeration of spanning Hamiltonian cycles in 377 physical exterior cases: all 365 live six-port planar unit probes, ten additional crossing-matching unit probes, plus exterior copies of A and B. This enumeration does not invoke the matching-connectivity formula to decide the cycles.
- Exact noncrossing checks for the integer-coordinate gadget drawings.
- Smith invariants independently obtained from gcds of all minors, not guessed from rank.
- All 32 binary domain vectors for the parity fibres and all 1,024 vectors modulo four for profile-parity recovery.
- Independent evaluation of the order-three signed Walsh polynomial on all eight labelled tournaments, its Smith invariants, and exhaustive precision checks modulo four and eight.

`RESULTS.json` lists every graph edge, all drawing coordinates, every probe gluing count, and the two complete Hamiltonian-cycle witnesses for the M1 existence hostile. It records zero for the corresponding ring-patch gluing. These finite checks supplement the proofs; they do not stand in for the universal gluing argument.
