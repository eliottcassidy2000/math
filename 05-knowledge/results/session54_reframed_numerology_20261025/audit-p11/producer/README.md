# P11 stress test: sharp Hall theorem and failed field extrapolation

## Status and outcome

PROVED: in every doubly regular tournament of order n>=11, the only size-(n-1) Hall blockers of Hamiltonian paths are the two-source and two-sink deletions, numbering n(n-1). FINITE-EXACT: at P11 this is110 blockers; eleven non-Hall stars give121 distinct ten-arc HP blockers, each independently verified inclusion-minimal.

REFUTED: the natural equivariant extrapolation of the P7 multiplicative carrier to F1024×. The order-five multiplier x->4x on P11 fixes exactly ONE ten-arc blocker, while the corresponding Frobenius z->z^4 fixes THREE nonzero elements of F1024. The full ten-arc HP-blocker family cannot be equivariantly identified with that field carrier, regardless of its unknown cardinality. The same obstruction excludes the full minimum-blocker family even without knowing its minimum size, as detailed below.

OPEN HERE: the global HP-blocking number beta(P11), beyond beta<=10, and the number of all ten-arc HP blockers. In particular, this work DOES NOT refute the bare cardinality1023. It refutes a specified equivariant field model.

## General Hall theorem

Let a doubly regular tournament have n=4lambda+3=2r+1 vertices, r=2lambda+1, indegree/outdegree r, and exactly lambda common in-neighbors per vertex pair. A Hall certificate consists of S,T with |S|-|T|=2; delete all arcs from S into V\T. This makes maximum bipartite matching size at most n-2, whereas a Hamiltonian path supplies n-1 distinct tail/head matches.

The deletion cost is

    r|S| - e(S,T) >= r|S| - r|T| = 2r=n-1.

Equality holds iff every in-neighbor of every vertex in T lies in S. Write t=|T| and U=union_{v in T} N^-(v). If d_x counts arcs from x into T, then

    sum d_x = rt,
    sum d_x^2 = rt + lambda*t*(t-1),
    |U| >= r^2*t/(lambda*t+r-lambda).

For 1<=t<=2r-2=n-3, the numerator of the strict comparison with t+2 is

    D(t)=r^2*t-(t+2)*(lambda*t+r-lambda)
        =r*(r-3)+lambda*(t-1)*(2r-2-t).

For n>=11, r>=5 and D(t)>0. Hence |U|>t+2=|S|, so equality is impossible for every interior t. Only t=0 and t=n-2 remain. They give, respectively, deletion of all out-arcs at two vertices (two sinks), and all in-arcs at two vertices (two sources). Each family has binomial(n,2) distinct members, and they are disjoint.

This also exhausts arbitrary size-(n-1) deletion sets that destroy every one-path-cycle factor, not only certificates chosen a priori. By the bipartite Hall deficiency theorem, such a deletion has some neighbor set T with |T|<=|S|-2. Enlarge T to exact deficiency2. Its required deleted arcs are a subset of the original deletion and already cost at least n-1, forcing equality and the classification above.

The stars are not Hall blockers: their residual digraph has exactly one isolated vertex and no other source or sink, whereas the complete Hall classification requires two sources or two sinks. They do block every Hamiltonian path. This proves n^2 distinct size-(n-1) HP blockers in any such tournament, but not global minimality or exhaustive HP-blocker classification. Inclusion-minimality is checked here only for P7 and P11, not asserted in general.

## Why seven is exceptional

At n=7, r=3 and lambda=1. The same identity becomes

    D(t)=(t-1)*(4-t).

The interior endpoints t=1 and4 now allow equality. These are precisely the common-neighbor Hall families omitted by a source/sink-only description. Enumeration recovers21+21+7+7=56 Hall blockers; the seven stars bring the known total63. At n=11 the positive term r*(r-3)=10 closes those extra Hall branches. The lost operation is a sharp neighborhood-expansion equality, not a loss of the valid P7 carrier.

## Exact P11 census

All sum_k binomial(11,k)*binomial(11,k-2)=binomial(22,9)=497420 certificates are tested. Their minimum costs for |S|=2,...,11 are

    10,12,14,15,15,15,15,14,12,10.

There are exactly110 distinct cost10 Hall masks, each with an explicit source-labelled arc list and Hall certificate in result.json. The independent equality-case checker enumerates unions of in-neighborhoods instead of the producer's edge deletions and recovers exactly the same masks. It verifies all121 examples with endpoint-subset HP dynamic programming and verifies that restoring any one deleted arc restores some HP. Explicit directed cycles in the complement of each star also witness its surviving one-path-cycle factor.

No enumeration of binomial(55,10) deletion masks is attempted. The unclassified case is a deletion retaining a one-path-cycle factor but destroying every Hamiltonian path. Such exotic blockers may exist below size10 or in additional size10 families.

## Decisive fixed-point hostile

The vertex multiplier h:x->4x mod11 has order5. It fixes no arc, since its only fixed vertex is0. Its55 arcs therefore split into eleven orbits of size5. Every h-invariant ten-arc deletion set is the union of exactly two such orbits: only binomial(11,2)=55 possibilities.

`fixed_point_stress.py` enumerates all55. Independent endpoint DP and direct depth-first path search agree: only the star at0 destroys every HP. For the other54, the output records explicit eleven-vertex surviving paths and verifies every consecutive arc against the deletion mask. It also checks all eleven invariant five-arc deletion sets and records surviving paths for each.

On F1024× the map z->z^4 fixes the solutions z^3=1, exactly three. Any equivariant bijection preserves fixed-point counts. Hence the full ten-arc blocker set is not this G-set. In fact the minimum-blocker set cannot be that G-set either: beta<=10 from stars; h-invariant arc sets have sizes divisible by5; no invariant five-set blocks; therefore the number of h-fixed minimum blockers is zero if beta<10 and one if beta=10, never three.

The mismatch persists for any nontrivial power of the order-five Frobenius, since gcd(2j,10)=2 for j=1,2,3,4. It also persists for affine-semilinear variants z->a*z^(4^j) that genuinely have order5: their order-five condition makes the fixed-point equation solvable, with exactly three nonzero solutions. No arbitrary action or bare1023 cardinality is ruled out.

## Reproduction and scope controls

Scripts use exceptions rather than assertions. Run hall.py to produce result.json, check_independent.py to produce independent.json, and fixed_point_stress.py to produce fixed-points.json. Their normal and -O outputs match byte-for-byte. The fixed-point test and Hall theorem are independent mechanisms and neither silently substitutes a partial census for a full one.

All files are scoped to this separate packet; the frozen P7 bridge packet is unchanged.
