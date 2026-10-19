# Signed affine potentials: square-sum paths, graceful trees, and signed Collatz words

**2026-10-19. Status: elementary exact reformulations and a proposed cross-domain proof framework. The named general conjectures remain open. No novelty claim for standard signed-graph switching or tree potentials.**

## 1. A precise version of the proposed analogy

There is a useful common language behind square-sum arrangements, graceful tree labelings, and the two signed Collatz maps: assign an affine transition to each edge, propagate a root coordinate, then enforce the admissible integer labels and any cycle closure. The edge transitions differ, so the admissibility data must stay attached to the maps.

For a rooted tree T, root coordinate t, and an oriented edge u→v, write

    x_v=a_e x_u+b_e.

There is a unique propagated expression x_v=A_v t+B_v, with root coefficients (A,B)=(1,0) and recursion

    A_v=a_e A_u,     B_v=a_e B_u+b_e.

A closed directed walk has affine monodromy x→Ax+B. Closure is exactly (1-A)x=B. If A=1, either the whole coordinate line closes (B=0) or none does; if A≠1, precisely one candidate root closes. Integer, positivity, distinctness, and branch conditions are additional exact requirements.

The significant distinction is between a whole-family identity and an isolated root. The same distinction appears in the new [integer Collatz stabilizer atlas](integer_stabilizer_20261019.md), where a literal integer merge can live on a nonidentity affine stabilizer.

## 2. Square-sum edges are reflections; bipartite switching is exact

Prescribe a square q_e² on each edge of a rooted tree. The sum constraint is

    x_u+x_v=q_e²,    x_v=q_e²-x_u.

Thus a_e=-1. At vertex v,

    x_v=s_v t+B_v,    s_v=(-1)^(depth v),
    B_root=0,         B_v=q_e²-B_u.                         (1)

Multiply coordinates by their depth signs: y_v=s_v x_v. Then

    y_v-y_u=s_v q_e².                                      (2)

This is a difference-potential problem with signed square increments. It is a bijection on all integer label assignments. Positive x labels become positive y on one color class and negative y on the other; that sign restriction is part of the transformed problem.

For a general connected graph, a vertex-sign assignment converting every sum edge to a difference edge exists exactly when the graph is bipartite: the signs must be opposite across every edge. Odd cycles obstruct this switch. Along an even square-sum cycle, closure requires the alternating sum of edge squares to vanish. Along an odd cycle it fixes twice the root label to an alternating sum of squares. This explains algebraically why triangles pin a scale/center while a tree leaves one root parameter free.

## 3. An exact interval-and-collision criterion for consecutive coverage

Let T have N vertices, and fix its edge squares. Compute B_v and s_v by (1). There is a square-sum labeling that uses every integer 1,...,N exactly once if and only if the following finite conditions hold:

1. B_u≠B_v whenever distinct u,v have the same sign s.
2. The integer interval

       I = intersection_(s_v=+1) [1-B_v,N-B_v]
           intersect intersection_(s_v=-1) [B_v-N,B_v-1]

   contains an integer t avoiding every value

       t=(B_v-B_u)/2 for s_u=+1, s_v=-1.                   (3)

**Proof.** The intervals are precisely the N separate inequalities 1≤s_vt+B_v≤N. Same-sign collisions are independent of t and occur exactly when B_u=B_v. Opposite-sign collisions occur exactly at the roots in (3). Once all N labels lie in an N-element interval and are distinct, coverage is automatic. Conversely a valid labeling provides such a root and satisfies the conditions. ∎

For the square-sum arrangement problem, take T to be the N-vertex path. Its edge squares must lie between 3 and 2N-1, since adjacent labels are distinct positive integers in the interval. The theorem converts the existence question into choosing those edge squares so that an interval survives and its integer points are not exhausted by the collision roots. Square values may repeat on different edges.

For a concrete illustration, the square-sum path

    8,1,15,10,6,3,13,12,4,5,11,14,2,7,9

has edge-square word `9,16,25,16,9,16,25,16,9,16,25,16,9,16`. Propagation gives

    B=(0,9,7,18,-2,11,5,20,-4,13,3,22,-6,15,1).

The interval test yields I=[7,8]. The root t=7 makes the last two labels both 8 and is excluded; t=8 yields exactly the displayed permutation. This is a worked illustration of the criterion, not a new existence result for N=15.

This is also a useful algorithmic certificate: once edge squares are supplied, the root interval and all forbidden values can be checked exactly in polynomial time. The hard search over the edge-square word has not disappeared.

For unconstrained positive labels, replace the upper bounds by infinity. Choosing edge squares successively large enough can always build a positive distinct tree labeling. The difficult consecutive-coverage requirement is exactly the simultaneous bounded interval in (3), not the existence of square-sum tree embeddings with unlimited labels.

## 4. Graceful labelings are translation potentials with an exact span condition

Let T have m edges and m+1 vertices. Choose a bijection d:E(T)→{1,...,m} and choose a sign ε_e∈{±1} for each rooted edge. Set

    B_root=0,    B_v=B_u+ε_e d_e.

**Proposition.** T is graceful if and only if some such edge assignment has all B_v distinct and

    max_v B_v-min_v B_v=m.                                 (4)

**Proof.** A graceful labeling supplies the signed differences and its vertex labels are exactly 0,...,m, so (4) holds. Conversely translate the potentials by -min B. They are m+1 distinct integers in an interval of m+1 integers, so they use precisely 0,...,m. Their absolute edge differences are the prescribed distinct values 1,...,m. ∎

Equivalently one may ask for span at most m, because distinct integral potentials already force span at least m. This gives the exact parallel with Section 3: both problems require affine root propagation, exclusion of vertex collisions, and tight interval packing. Their palettes differ: graceful edge differences are a permutation of every value 1,...,m, while square-sum edges belong to a sparse square set and may repeat. Switching square sums to signed differences preserves neither palette nor the positive-label interval automatically.

The inherited [graceful exchange packet](graceful_square_exchange_20261019.md) and [frozen-tail obstruction](graceful_fiber_obstructions_20261019.md) already show why a fixed short-color buffer does not suffice; those finite tests are not rerun here.

A concrete joint target is therefore a **palette-preserving local surgery** that improves interval feasibility without introducing potential collisions. For graceful trees it must keep the edge-difference permutation; for square-sum paths it must keep every selected edge value square. A generic affine conjugacy does not ensure either requirement.

## 5. Signed Collatz has the same closure equation, with different guards

For ε=±1, use the parity-accelerated map

    T_ε(x)=x/2 if x is even,       (3x+ε)/2 if x is odd.

On the signed integer line, reflection R(x)=-x gives the exact conjugacy

    T_- = R T_+ R.

It exchanges the positive and negative half-lines. It therefore explains why studying positive 3n-1 dynamics probes negative 3n+1 dynamics, without implying the same cycle count on the positive half-line. The familiar cycles do not establish that there are no others.

A fixed lawful parity word of length n with k odd steps composes to

    x → (3^k/2^n)x + ε B,

where B>0 whenever k>0 and B=0 for the all-even word. For a nonzero positive cycle with k>0, closure gives

    x=ε B/(1-3^k/2^n).

Consequently the plus map requires 3^k<2^n, while the minus map requires 3^k>2^n. These are necessary sign conditions, not sufficient conditions for an integer cycle. Integrality and every parity branch must also hold. The sign flip reverses which side of the multiplicative resonance can contain a positive candidate.

This is different from the square-sum/difference switch: square-sum edges have slope -1 and graceful edges slope +1, while both signed Collatz maps have positive branch slopes. The common mechanism is affine propagation plus preserved admissibility constraints, rather than an identification of the signs in the names.

For comparisons with period statements elsewhere, this section uses the Terras clock, not the standard unaccelerated Collatz clock. The finite eleven-point golden dictionary uses the standard clock and should not be read as a five-step assertion for T_ε.

## 6. A concrete next obstruction to attack

The new contact constructions solve increasingly flexible rational/integer square-sum embeddings, but their common-center labels are sparse and can be enormous. They do not yet solve consecutive coverage. The exact feasibility theorem above specifies what a coverage-preserving construction would need to control:

- the intersection interval for the root must contain an integer;
- same-color potential differences must remain nonzero;
- opposite-color collision roots must leave an integer uncovered;
- edge squares must remain within the allowed consecutive-label palette.

The graceful analogue asks for a distinct signed potential set of minimum possible span while using every edge magnitude once. This suggests transporting **local moves and defect functions**, rather than asserting an isomorphism of the conjectures. A candidate defect is excess span plus a vertex-collision count, with the palette treated as a hard constraint. It becomes a proof only after proving a terminating improvement rule or a complete obstruction classification; neither is supplied here.

The useful synthesis is exact: tree structure removes cycle consistency but leaves global packing; odd square-sum cycles impose isolated-root consistency; Collatz merges can be whole-cylinder identities or isolated stabilizer hits. Each setting exposes a different part of the same affine feasibility framework.
