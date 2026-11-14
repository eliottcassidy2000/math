# An exact base-four inverse for prescribed H and Q(2)

This is an arithmetic relaxation of the physical tournament inverse. It rejects impossible responses and describes every slot polynomial allowed by the relaxation; feasibility of a relaxed polynomial need not imply a tournament realizes it.

At order n≥3 put δ=1 for odd n and δ=2 for even n, and r=(n−δ)/2−1. The positive odd-cycle packing formula forces

    Q(z)=z^n+2 z^δ Σ_(j=0)^r c_j z^(2j),   c_j nonnegative integers.

The leading term uses no cycles; a selected odd cycle consumes an even number of powers relative to its singleton vertices and adds an even multiplicity. The minimal exponent is positive and has parity n. For a prescribed H,Q2, set

    S=(H−1)/2,    W=(Q2−2^n)/2^(δ+1).

The lower coefficients are exactly representations W=Σ c_j 4^j with coin count S=Σ c_j. Require nonnegative integral S and W. The bounded greedy representation a of W uses floor(W/4^r) top coins and the ordinary lower base-four digits; let G=Σ a_j. Then the response passes this arithmetic inverse **if and only if**

    G≤S≤W,    S≡W (mod 3).

Indeed merging four equal coins into the next denomination reduces the count by three and eventually reaches the unique bounded greedy representation (normalize from bottom to top). Conversely split a coin to increase the count by three until every coin is a unit.

## Exact carry variables and higher slots

Let t=(S−G)/3. Every permitted polynomial is obtained by precisely t splits from the greedy one. If u_j is the number of splits of denomination 4^j, then

    c_j = a_j+4u_(j+1)−u_j,   0≤u_j≤a_j+4u_(j+1)   (1≤j≤r),
    c_0 = a_0+4u_1,   u_(r+1)=0,   Σ u_j=t.

These are necessary and sufficient: do the splits from top to bottom. At fixed slack t there are at most binomial(r+t−1,t) vectors (with the elementary empty-index boundary understood), not a free search over all coefficients. In particular Q3 has the exact range obtained by evaluating

    Q3 = 3^n + 2·3^δ [Σ a_j 9^j − 5 Σ u_j 9^(j−1)]

on those integer vectors. This gives sharp extrema in the relaxation. It is also a direct NP-style arithmetic certificate for every proposed higher slot value, independent of how the tournament is presented. To require at least τ directed triangles, append c_r≥τ: c_r counts exactly the triangle terms. A strong tournament needs at least ceil(n/3) triangles, since each vertex lies on one (otherwise its in-neighbors dominate its out-neighbors, making a one-way cut). Higher physical constraints are not represented by these carry inequalities.

At zero slack the entire Q polynomial is forced. Thus a same-order, same-H, same-Q2 fibre with zero slack has no hidden unpointed continuation information at any higher slot.

## A full order theorem at the first nonzero slack

When t=1, each allowed polynomial is one split at an occupied positive denomination. If j<k are two possible split levels, subtracting their polynomials gives

    Q_j(z)−Q_k(z)
      =2 z^δ (z²−4)(z²−1) z^(2j−2) Σ_(i=0)^(k−j−1) z^(2i).

This is z(z−1)(z−2)g(z−2) with g having nonnegative coefficients and positive leading coefficient. If F_m is the ordered nonempty m-cover count, inclusion-exclusion yields

    F_m(j)−F_m(k) = m(m−1)(m−2) Δ^(m−3)g(1)   (m≥3),

where Δ is the unit forward difference. Every such difference is nonnegative: unit differencing preserves nonnegative polynomial coefficients on positive arguments. It is strictly positive from m=3 through the degree of the difference, and zero thereafter; slots one and two are tied. Thus **every physical fibre at slack one is totally ordered in its full ordered-cover profile**, even if many distinct graphs realize the same profile.

In a common homogeneous context, each quotient itinerary reads a nonnegative product of nonempty cover counts. Any feasible visit count inside the stated strict range separates the two profiles with a fixed sign; no opposing higher-cover responses can compensate there. A context that only uses equal coordinates still ties. This is a first-carry obstruction to designing a signed return in the same-(H,Q2) fibre, complementary to the [literal ordered-cone obstruction](../session321_selective_returns_and_fano_wiring_20261113/FOUR_VISIT.md). It is not a theorem that every relaxed fibre is physically populated.

## Two elementary companion bounds

The leading term and its next allowable exponent give 2^n≤Q2≤(H+3)2^(n−2). Equality at the upper edge requires every nonempty packing to consist of just one directed triangle; equivalently there are no directed odd cycles longer than three and no disjoint directed triangles. This already sharply separates the small-H and large-order axes before invoking a strong-order theorem.

There is also an exact path/cycle identity without polynomials. Group ordered nonempty two-covers into unordered pairs of paths, and let B0 count those for which neither directed concatenation is legal. Let C be the number of directed spanning cycles modulo rotation. Both concatenations legal means exactly a spanning cycle with two deleted arcs, so there are binomial(n,2)C such unordered covers. There are (n−1)H legal ordered concatenations by splitting Hamilton paths once. Hence

    Q2 = 2nH − n(n−1)C + 2B0.

The earlier injection bound Q2≥(n+1)H is precisely the consequence of nC≤H. Any lower bound on nonclosing Hamilton paths H−nC, or on doubly blocked seams B0, strengthens it. Unlike a bare moment estimate, this identity states exactly what lost endpoint information is needed for a sharper physical bound.

Replay `replay_slot_carry.py` independently exhausts small integer coefficient boxes, checks the one-carry order on nonempty-slot profiles, and verifies the seam identity by literal permutations on example tournaments.
