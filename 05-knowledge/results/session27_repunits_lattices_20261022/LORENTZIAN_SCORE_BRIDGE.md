# The 25-vertex score lift: intransitivity is a Lorentzian norm

Status: PROVED elementary score identity; CITED classical lattice construction. The same norm identity was proposed concurrently by the S52 lane in message MSG-1913 at commit a08f25ecbf8abd97959305f0836c1c25890edf65. It is independently verified here, without treating its reserved theorem IDs as proved dependencies or claiming priority.

For a tournament T on25 vertices, let s=(s0,...,s24) be its outdegree vector and c3(T) its number of cyclic triples. Set

    v(T)=(s0,...,s24;70),
    Q(x;t)=sum_i x_i^2-t^2.

Since sum s_i=300 and c3=binom(25,3)-sum_i binom(s_i,2),

    Q(v(T))=sum_i s_i^2-4900=-2c3(T).

This proves that the lift is always causal, in the signature convention with negative time norm. It is null precisely when T is transitive, because a tournament is transitive exactly when it has no cyclic triangle. Otherwise it is strictly timelike. The identity is all-order-in-objects at the fixed order25, not a numerical fit.

## Actual even Lorentzian lattice membership

Use the standard coordinate model of II_(25,1): coordinates are all integers or all half-integers, and have integral Lorentzian pairing with (1/2,...,1/2;1/2). The integral score lift satisfies this condition because (300-70)/2=115 is an integer. Thus it is an actual vector of the even unimodular Lorentzian lattice, not merely a point in a real Minkowski space.

For a transitive tournament the scores are a permutation of0,1,...,24. Hence its lift is a coordinate permutation of

    w=(0,1,...,24;70),

whose null norm is the cannonball identity sum_(j=1)^24 j^2=70^2. The classical construction identifies the roots r with r^2=2 and(r,w)=-1, under their difference metric, with the Leech lattice. The primary coordinate model and this statement are given on printed p76 of Borcherds–Conway–Queen–Sloane, *A Monster Lie Algebra?*, Advances in Mathematics53(1984),75–79, [author-hosted paper](https://neilsloane.com/doc/Me99.pdf), citing Conway–Sloane's1982 Lorentzian construction.

Equivalently the positive null quotient w-perp/Zw is Leech. To see the affine/root-to-quotient correspondence, choose r0=e0-e1, with r0^2=2 and(r0,w)=-1. Differences of such roots lie in w-perp. For any y in w-perp, a unique integral multiple a of w makes r0+y+aw have norm2, since a=(y^2+2(r0,y))/2 is integral in an even lattice. This identifies the difference metric with the quotient metric, modulo the choice of root origin.

This is a precise25→26→24 bridge: tournament order25, Lorentzian rank26, and the positive quotient rank24. It does not arise by relabelling a forbidden Hamiltonian value.

## The map's strict domain boundary

For a cyclic tournament v(T)^2<0, so v(T) is not contained in its own orthogonal complement. The expression v(T)-perp/Zv(T) is then not a defined lattice quotient of that form. One cannot assign the same Niemeier/null-quotient construction to every tournament by copying the transitive case.

Even the exact Lorentzian point is not a faithful Hamiltonian invariant. Take the two eight-vertex examples in SCORE_ROOT_COLLISION.md and prepend a transitive chain of17 vertices dominating the original tournament. Their common score vector becomes

    (4,4,4,3,3,3,3,4,24,23,...,8).

Both lifts have norm-40 and c3=20, yet their Hamiltonian counts remain657 and643: the source chain has exactly one possible prefix and contributes no cycles. Thus the Lorentzian norm measures a real intransitivity statistic, but not all intransitivity or complete metagraph structure.

The preserved predicate is transitivity/nullness and the preserved scalar is c3. The lost data include higher directed cycles, Hamiltonian count and attachment response. The finite-context reconstruction theorem supplies a different, faithful observable rather than asserting that the score lift secretly retained these coordinates.

## Checks

The companion check_lorentzian_score.py tests the score identity on explicit examples and deterministic random25-vertex tournaments, checks the lattice parity condition and transitive null vector, and verifies the padded hostile algebraically. Finite controls corroborate the elementary proof; no25-vertex Hamiltonian enumeration is claimed.
