# Blue fixed-frame folding: paths and odd-cycle packings

Let T be a tournament with involutory anti-automorphism r, with m swapped pairs and optional fixed vertex c. Write a signed vertex for a choice of one specified endpoint of a swapped pair.

## Exact folded blue path count

For every subset S of pair indices and signed endpoint v over some i in S, let D(S,v) count directed paths selecting one vertex from each pair of S and ending at v. Initialise D({i},v)=1. The recurrence is

    D(S,v) = sum_{u over S minus i, u -> v} D(S minus i,u).

The number of r-palindromic Hamiltonian paths is precisely the sum of D(all,v) over v with v -> r(v) in even order or v -> c in odd order. Indeed the first half selects a transversal, the second half is its image under reversal and r, and the seam test is exactly the displayed arc; at odd order c -> r(v) is forced by v -> c. For n=1 the answer is 1. This computes the exact fixed-path count in O(m² 2^m) integer operations and O(m 2^m) storage. It is not a count of all Hamiltonian paths.

After fixing r, divide this answer by |C_Aut(T)(r)| to recover the number B(T) of blue tilings above the isomorphism class, including in nonrigid tournaments. Indeed blue presentations arise from fixed paths for various anti-involutions, and all involutory anti-automorphisms are conjugate by the automorphism group (Sylow, or THM-024 with its conjugacy refinement). Once r is fixed, two such paths give identical labelled tilings exactly when transported by a tournament automorphism, which must commute with r because r is uniquely determined by either path and position reversal. This action on paths is free. THM-643 supplies the path-colour dictionary; the DP retains only half as many independent pair labels as an ordinary subset DP.

## Folded fixed-point polynomial for odd-cycle packings

The odd-cycle formula writes H(T)=P_T(2) where P_T(z)=sum z^k over unordered packings of k vertex-disjoint directed odd cycles of length at least three. Reversal and r send directed cycles to directed cycles and preserve disjointness. Write P_T=F+2G, where F counts fixed packings and G counts two-element packing orbits, by grade. A fixed cycle, being odd, must contain a vertex fixed by r and has exactly one. Thus no such cycle exists in even order; in odd order there is at most one per packing.

Every nonfixed cycle bundled into a fixed packing occurs together with its image, disjoint from it; it therefore chooses at most one signed endpoint per pair. A bundle consumes both endpoints of all its selected pairs. For an odd-cardinality set S of at least three pair indices, let q(S) be the number of directed cycles selecting exactly one endpoint over each pair in S, with the selected endpoint at min S prescribed to be a chosen first sign. Exactly one member of each reversal bundle has that sign, so no quotient division is needed. Set q(S)=0 otherwise. A fixed central cycle using pair set S has the form

    c, v_1,...,v_j, r(v_j),...,r(v_1), c.

For each nonempty S let a(S) count such cycles. It is obtained by running the same signed endpoint DP with singleton initialization [c -> v], and summing over endings v with v -> r(v). No rotation or orientation ambiguity remains because c is fixed and the directed cyclic order specifies the first half. Set a(empty)=0.

For allowed pair subset U, put L(empty)=1 and

    L(U) = L(U minus i) + z² sum_{S subset U, i in S} q(S) L(U minus S),
    i = min U.

This is the fixed-point generating polynomial for packings of only paired cycles. The fixed packing generating polynomial is L(all) at even order, and

    F(z) = L(all) + z sum_{S subset all} a(S)L(all minus S)

at odd order. Every term has a unique combinatorial decomposition, so this folds **all** fixed packings, not only the triangle terms. All q(S) can be computed with at most m rooted signed DPs (anchored at min S and closed by its final arc), and the last convolution uses at most 3^m subset incidences, up to polynomial factors in m.

## What the fold does and does not determine

The exact dissipated information is 2G(2); it is even but not generally zero modulo higher powers of two. A fixed-packing-only oracle cannot recover H: at even order four, every r-anti tournament has F=1 because a bundle needs at least three pairs; the transitive tournament and the strong order-four tournament nonetheless have H=1 and 5. Both have blue fixed frames, e.g. natural path 0,1,2,3 with all edges increasing in the first, and just edge {0,3} reversed in the second. Their fixed Hamilton path counts are 1 and 3, respectively. Grades make the precision transparent:

    H = 1 + 2 F_1 + 4 (G_1+F_2) (mod 8).

Here G_1 counts all nonfixed odd-cycle orbits, not merely triangles. Since even order has no fixed odd cycle, one recovers the known H=1 (mod 4) constraint for every even blue node (THM-791). More sharply, at even order F_2 exactly counts disjoint cycle pairs C,r(C); hence

    H = 1 + 4 (# odd-cycle reversal orbits with C intersect r(C)) (mod 8).

Overlap is equivalent to the cycle containing both endpoints of at least one swapped pair. This gives a concrete obstruction at precision eight and identifies *which* nonfixed information is needed. It does not assert that the overlapping orbit count is constant or that triangle counts suffice. At odd order the mod-4 digit is exactly parity of the central fixed odd-cycle count, available via the signed DP in O(m² 2^m) arithmetic operations.

## Full exact count: a separate meet-in-the-middle reflection fold

Write E(S,a) for the full Hamilton path count on S ending at a, computable by the ordinary subset recurrence. At even order n=2m, cutting any Hamilton path at its middle arc gives the exact formula

    H(T) = sum_{|S|=m} sum_{a in S, b not in S} [a -> b] E(S,a) E(r(S^c),r(b)).

At odd order n=2m+1, cut at its middle vertex x:

    H(T) = sum_x sum_{S subset V minus x, |S|=m}
              sum_{a in S, b in V minus (x union S)}
              [a -> x][x -> b] E(S,a) E(r(V minus (x union S)),r(b)).

The second factor is a start count transformed to an end count by r and path reversal. In either case the outer sum's split data are involuted by S -> r(complement), with x -> r(x) in the odd formula. Pair off unequal splits to evaluate one copy with weight 2, and retain fixed splits with weight 1. This yields an exact ungraded computation retaining nonfixed contributions rather than mistaking the faster signed half-path DP for all of H. It avoids subsets larger than m and duplication of endpoint calculations, but no superior exponential base for arbitrary exact H is claimed.

## Blue degree and the score sequence still lose count information

There is a small exact obstruction to recovering H from the local
blue degree even if the order, full sorted score multiset, strongness,
and automorphism order are supplied. The two five-vertex outgoing
mask rows

    (18,4,9,19,6)     (2,12,9,17,7)

are both strong and anti-preserved by index reflection. They are
rigid, each has three reflection-fixed Hamilton paths (hence blue
degree three), and both have scores (1,2,2,2,3). Their total H
counts are 11 and 13 respectively. The companion replay checks both
directly. This is not a counterexample for an oracle that sees the
whole blue graph with additional boundary data.
