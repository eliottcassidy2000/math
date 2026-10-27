# Exact carrier, lower bounds, and complete geodesic cones

Named vertices are retained throughout. A directed Hamilton word P is blue iff
its positional reversal is an anti-automorphism: for every i<j,
T(p_i,p_j)=T(p_(n−1−j),p_(n−1−i)). All other words are black. R_P reverses
exactly the directed path arcs. K_P complements exactly the nonpath pairs.
No Hamilton count, endpoint-closing condition, or scalar score test is used
as a replacement for this full predicate in any certificate replay.

For a nonempty proper shore A, c_T(A) counts its outgoing arcs. The path's
signed crossings telescope, giving

    c_(R_P T)(A)−c_T(A)=1_A(last P)−1_A(first P).

The singleton case gives the degree current: the start loses one outdegree,
the end gains one, and every other named score is unchanged. Define λ as the
minimum outgoing cut and δ as minimum semidegree.

**Lower bounds.** Every nonstrong tournament has a zero outgoing cut, so its
R-distance from strong T is at least λ. Every singleton-plus-strong-core target
has a source or sink; making a specified source v needs at least its initial
indegree moves, and making a specified sink needs at least its initial
outdegree. Therefore distance to that target is at least δ.

**All-colour λ equality.** Fix a minimum shore A. Before λ moves, every cut
is still positive. In the current strong tournament choose a Hamilton cycle
and cut at an arc entering A; the resulting Hamilton word starts in A and
ends outside. Reversing its path reduces c(A) by one. Repeat λ times. This is
an all-colour argument; it does NOT prove those words black.

**Complete λ cone.** If a lengthλ route reaches a zero cut A, that shore had
initial value at leastλ, and each move lowers it by at mostone. Equality forces
initial c(A)=λ and every path to start in A and end outside. Thus all initial
minimum shores, with those endpoint restrictions at every step, form a complete
lower-bound-depth search universe. A route may not switch its witness shore
without independently retaining the corresponding complete cones.

**Complete δ cone.** If a lengthδ route reaches a source v, its initial
indegree must equalδ and every step must end at that same named v. A sink case
is the converse. Hence enumeration of all initial minimum-indegree and
minimum-outdegree ports, with the fixed-endpoint condition, is complete for
attainment of the δ lower bound. A failure at one port is not a global lower
bound. Intermediate strongness is checked, not assumed from scores.

**Literal K conversion.** For an R-route U_j→U_(j+1) with word P_j, set
T_j=U_j for evenj and T_j=U_j^op for oddj. Use word P_j for evenj and its
reverse for oddj. This is a literal K-route, since K_P(T)=(R_P(T))^op.
Both target predicates are converse-invariant. The independent replayer
constructs K by complementing nonpath unordered pairs, rather than merely
assuming this identity; it checks actual directedness and full colour at both
ends. Consequently all distances in this packet are valid in both gauges.
