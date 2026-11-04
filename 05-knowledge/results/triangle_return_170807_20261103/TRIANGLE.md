# A rotating triangle is a precise return obstruction

2026-11-03. **PROVED elementary consequences of inherited run and cut
theorems.** Finite controls are separate. The unrestricted {7,21} spectrum
classification remains OPEN; no historical-priority claim is made.

## 1. The hidden module in the actual 3-ready source

Use the literal named family T_m and the second feedback vertex q from
[TRANSPORT](../return56941_two_feedback_20261103/TRANSPORT.md), m>=1.
Its transitive base consists of m clones followed by
32,22,2,1,6,8,4,13,5,14,7, and apex21 has the inherited1243 incidence word.
Write U_m=T_m+q for the exact saved birth, not an isomorphic replacement.

The three named vertices {32,22,q} form a homogeneous directed triangle:
32->22->q->32. All clones beat all three; all three beat every later base
vertex and apex21. Thus the second feedback operation literally replaces
the TT2 module {32,22} by C3, with every other old name retained.

This identifies a physical symmetry that the scalar identity
H(U_m)=3H(T_m) alone does not specify.

## 2. Where the first missing response occurs

For a fixed exterior, let b_j be the number of valid exterior skeletons
with j visits to the substituted module. Visits are in ordered path bins;
inside each bin adjacent module pieces are combined into one maximal run.
All prescribed endpoints in this statement lie outside the module.
The inherited maximal-run bijection gives response=sum_j b_j C_j(A),
where C_j(A) counts ordered nonempty j-path covers of the module A.
See [run-content T1](../uniform_prime_transport_20261103/THEORY.md#t1-arithmetic-content-of-a-homogeneous-module).

The two actual modules have

    (C1,C2,C3)(TT2)=(1,2,0),
    (C1,C2,C3)(C3) =(3,6,6).

Consequently, for every such fixed exterior and outside-port response,

    response(C3)-3 response(TT2)=6 b3.                 (TR1)

There are no higher terms, since the larger module has only three vertices.
In particular exact tripling holds iff the context has no valid three-run
skeleton. This is a nonnegative, realization-preserving criterion, not a
formal cancellation among signed coordinates. In the slot-polynomial basis
the same difference is X^3+2X-3X^2=X(X-1)(X-2).

For the old root31, the inherited all-size identities imply b3=0 for H,
E1,S1,deletion-H,E2; for S2 they give the exact positive skeleton count

    b3=22546(3^(m-1)-2^(m-1)).                          (TR2)

Thus the higher-run correction already recorded in TRANSPORT is the first
actual way the exterior can use the extra triangle vertex. TR2 explains
that correction; it does not claim the earlier formula as new.

## 3. A mark-stabilizer divisibility rule

If G is a subgroup of Aut(T) fixing every named endpoint, then |G| divides
each unweighted ordered spanning path-cover response with those endpoints,
and divides H(T). The action on ordered vertex words is free: fixing a
cover fixes every vertex at its position in an ordered bin. This applies
to joint responses with individually fixed ports; merely preserving the
set of port names is insufficient. Arbitrary non-invariant weights need
not satisfy the divisibility statement.

For a strong tournament put b_r=E_r(2)-H. The inherited rooted-defect law
gives b_r>=d_r(d_r-1), with d_r the indegree of r. In U_m the triangle
rotation fixes every exterior root, so b_r is a multiple of3 there.
All exterior roots except the first clone have indegree at least2. The
first clone has b_r=3 for m>=2 by the inherited transport theorem. Every
triangle root has indegree m+1>=3. Therefore

    min_r (E_r(2)-H)(U_m)=3, for every m>=2.             (TR3)

In particular **no re-marking anywhere in U_m makes it1-ready or2-ready**.
The m=1 boundary is excluded: its first-clone defect is zero. Naming a
triangle vertex breaks the stabilizer but incurs the indegree obstruction;
this is why simply retaining a symmetry frame does not produce a reset.

## 4. No number of births can fill the original two-unit deficit

More generally let T be strong and contain a proper homogeneous C3 module.
Then no induced extension of T, at any order, has H(T)+2 or H(T)+4 paths. (TR4)

Proof. If the new incidence is constant on the triangle, its rotation
remains an automorphism. Both counts are divisible by3, excluding increments
2 and4. Otherwise its out-neighbor cut splits the triangle and crosses at
least one directed triangle arc u->v. The positive cut law and its stronger
weight bound give an increment at least

    2 w_uv >= 2 H(T-{u,v}) >= 6.

Here T-{u,v} is exactly the quotient contracting the proper triangle to its
remaining vertex. That quotient is strong, has at least three vertices,
and so has at least three Hamilton paths. The cut formula and bound are
inherited from [Session124](../session124_tournament_cut_response_20261031/README.md),
not new claims about generic weighted graphs. This proves the one-birth
case and supplies the uniform first-split cost.

For any longer birth sequence, consider the first birth whose incidence
splits the original triangle. Before it, the triangle rotation is still an
automorphism. Every prefix count is an odd multiple of3, so any positive
increase from H(T) is already at least6. If all those prefixes instead
have count H(T), the first split costs at least2H(T-{u,v})>=6: its current
pair-deleted graph contains the original T-{u,v} induced, and H cannot
decrease under births. All later counts are no smaller. If no split ever
occurs, the final count remains divisible by3. This proves TR4 for every
finite extension, including nonstrong intermediate and final graphs.

Applied to U_m, TR4 rules out **every birth-only continuation** producing
3H(T_m)+2 or3H(T_m)+4, including arbitrary nonuniform incidences and any
number of feedback vertices. A repair must change an old incidence or
delete an old vertex. This is stronger for this particular source than
excluding a fixed uniform two-feedback template. It does not exclude
other sources or the larger targets3H(U_m)+2/+4.

## 5. Ports, phases, and the hostile control

Keep the full relation, physical names31,6,21,32,22,q, and the old clock m.
An exterior-only summary sees the tripling but misses TR2. A summary that
also forgets the original m=1 boundary cannot infer readiness from field
phase: the inherited m=1/m=7 collision still applies.

The existing current F49/F343 codecs, all five Q4 presentations, and their
transport laws are reused unchanged. These are observations, not permission
to relabel or reset a source. The12-vertex81/27 pair is replayed with its
literal relations and distinguished birth; matching unrooted profiles do
not identify their joint response carriers. None of TR1--TR4 asserts that
the saved two-run responses form a carrier closed under all complements.

## 6. A general symmetry barrier for every induced extension

Let G be a nontrivial automorphism subgroup of a tournament T, of odd order
g. Let mu_G be the minimum H(T-{u,v}) over distinct vertices in the same
G-orbit. Every finite induced extension V of T satisfies

    H(V)=H(T), or H(V)-H(T)>=2 min(g,mu_G).              (TR5)

Proof. Order the new vertices in any way. As long as each birth is constant
on every original G-orbit, G extends by fixing the newborn vertices. Every
prefix H is divisible by g, and tournament parity makes each increase even.
Thus a positive cumulative increase already costs at least2g. At the first
birth that splits an orbit, that induced orbit is a vertex-transitive
tournament and hence strong: its uniquely ordered strong components would
otherwise be preserved separately, contradicting transitivity. The cut
therefore crosses an arc u->v within the orbit. Its cost is at least
2H(current-{u,v})>=2H(T-{u,v})>=2mu_G. Future births cannot decrease H.
These two cases prove TR5. No transitivity of G on all vertices is required.

All tournament automorphism groups have odd order: a nontrivial involution
would swap a pair and reverse its edge, and Cauchy's theorem excludes even
group order. Alternatively odd order can simply be retained as a hypothesis
of TR5. The group must be actual; divisibility of H by g is insufficient.

For the inherited81/27 pair, independent counts of all66 two-vertex
deletions of each graph have the same minimum7821. Using their actual
automorphism groups, TR5 excludes positive increments below162 for the
order81 source and below54 for the order27 source, at every extension
order. These are **guaranteed windows**, not claimed sharp gaps or unequal
true minimum increments. Both scalar H values are divisible by81; this
does not authorize the81-group argument for the second source. The full
named-cut hostile still distinguishes them. The pair-deletion receipts
are in hostile_pair_deletions.json, in lexicographic pair order.
