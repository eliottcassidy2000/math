# Suspension needs a source-incidence record, not just the old target class

Status: PROVED by an explicit native example, with FINITE-EXACT controls.

## Inheritance

Incoming `sixslot570_20261024/PORTS.md` proves the general suspension identities
and correctly says its construction is on source/path records, not bare target
classes. This note supplies a small, explicit obstruction to removing that
qualification, including when all new outer ports are intrinsically recognizable.
It is a new witness to an already identified boundary, not a new proof of the
incoming suspension theorem.

For a tournament T and Hamiltonian order P, let K_P(T) reverse all arcs not
on P. Write S(T)=a -> T -> b, with a->b. The suspended path is (a,P,b), and
write H(T,P)=K_(a,P,b)(S(T)). The question is whether

    [K_P(T)] -> [H(T,P)]

can be a well-defined map of ordinary or converse-merged tournament classes
independent of the source/path record. In general it cannot.

## Explicit order-five witness

Use owners0..4 and P=(0,1,2,3,4). All path arcs point forward. List off-path
pairs in the order

    (0,2),(0,3),(0,4),(1,3),(1,4),(2,4).

For b=2 or3, orient a listed pair i<j as i->j precisely when its corresponding
little-endian bit is1; otherwise use j->i. Call the two sources T2,T3.
Their labelled degrees are respectively

    (2,1,2,2,3),   (3,1,1,2,3).

Neither source has a source vertex or a sink vertex. Their old target
matrices K_P(T2) and K_P(T3) are isomorphic: the actual owner map

    (0,1,2,3,4) -> (1,2,3,0,4)

preserves every old target arc.

Nevertheless their enlarged seven-owner target score multisets are

    H(T2,P): (1,2,3,3,3,4,5),
    H(T3,P): (1,2,2,3,4,4,5).

Their degree-square sums are73 and75. For a tournament on N owners,
converse replaces d by N-1-d; because sum d=N(N-1)/2, the sum of squares
is unchanged. Thus these targets are neither isomorphic nor converses.
This disproves the proposed map even after converse merging.

In BOTH enlarged targets a is the unique degree-one owner and b is the
unique degree-five owner. Their outgoing/incoming singleton incidences
recover the old endpoints0 and4. The failure therefore cannot be attributed
to ambiguous outer markers. The missing data is the source/path incidence
used to construct the enlargement, not just the ability to recognize its ports.

## What survives, and what this does not refute

The map from full named records (T,P) to the enlarged full record is still
well-defined and faithful in the stated sense. On the distinguished enlarged
targets, deleting the identified ports recovers the old inner target. Our
example does NOT refute that reverse readout. A many-to-one readout from new
targets to old ones is compatible with the failure of a forward single-valued
map. A harmonic pullback through a reverse readout could still exist if its
full Laplacian equations hold.

Nor does the example disprove the existence of a more elaborate incidence
carrier, correspondence, or compensated transport. It requires retaining
those extra choices explicitly. Counting matching classes or recording an
inner restriction is not a proof of harmonicity on the enlarged metagraph.
The all-level matrix/readout test in the carrier-transfer packet addresses
that different obligation; no order-ten harmonic-bank certificate is claimed.

## Exact boundary-cycle conservation at every order

There is a useful positive invariant in this enlargement. Write c(T) for the
number of cyclic directed triangles, and let f,l be the first and last owners
of P. Then for every T,P on n>=2 owners,

    c(K_P(T)) = c(T) + d_T(f) - d_T(l) - 1,
    c(H(T,P)) = c(T) + n.

Consequently the number of cyclic triangles of H that touch a new outer owner
is exactly

    c(H) - c(K_P(T)) = n+1-d_T(f)+d_T(l).

Proof. For every tournament c(T)=choose(n,3)-sum_v choose(d_T(v),2), since
each transitive triangle has a unique source. The inner complement degree
vector is n-1-d_T(v)+1[v=f]-1[v=l]. Its squared-degree sum differs from that
of T by 2(d_T(l)-d_T(f)+1), proving the first identity. The enlarged degrees
are n-d_T(v),1,n. Since sum d_T(v)=choose(n,2), their transitive-triangle
count is the old transitive-triangle count plus n(n-1). Also
choose(n+2,3)-choose(n,3)=n^2, proving the second identity. Deleting the two
outer vertices leaves exactly K_P(T), so subtraction counts precisely the
boundary triangles.

For a FIXED source, moving path endpoints may change the inner target's
cyclic-triangle count, but the new boundary triangles compensate exactly.
The enlargement converts endpoint-dependent degree corrections into literal
boundary incidence while keeping its total c(T)+n fixed. This is an all-size
combinatorial conservation law. It neither transports an electrical harmonic
character nor establishes that the native packet route respects its currents.
Keeping these two conservation problems distinct is the point of the example.

## Controls

`check.py` constructs both actual source and target matrices, checks the
explicit old isomorphism, authenticates the Hamiltonian paths, verifies the
outer-port recovery and both score obstructions, and independently compares
full canonical forms under every owner permutation and converse.
It also verifies the suspension identities on all1,099 forward-path tournaments
of orders2 through6, including direct triangle counts. This is not a full order-seven metagraph enumeration.

Run `python3 check.py --output normal.json` and repeat with `python3 -O`.
The normal/optimized receipts agree. The example is explicit; no external
classification or numerical approximation is needed.
