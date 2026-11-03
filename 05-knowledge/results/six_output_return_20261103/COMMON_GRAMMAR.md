# A matching-based grammar for actual module returns

**PROVED, scoped:** G1--G4. The grammar has an exact acceptance test for a
supplied return program. Existence of both ternary programs at every grammar
object remains **OPEN**. This is not a classification of all strong atoms.

## G1. The exact run bound is a rooted triangle matching

Let T be a tournament and q a distinguished vertex. Form a bipartite graph
L_q with left vertices N+(q), right vertices N-(q), and edge a--b precisely
when a->b. Its edges are exactly the directed triangles q->a->b->q. Let
nu_q be its maximum matching size. This is an intrinsic binary graph
relation; no ties are being forced into a tournament.

Collapse a substituted module at q in a spanning Hamilton path. A *visit*
is one maximal run inside that module. For a nonempty child of order m,

    largest possible number of visits = min(m,nu_q+1).    (G1)

In particular the transitive-substitution formula is

    H(T[q <- TT_m]) = sum_(r=1..nu_q+1) c_r r! S(m,r),     (G2)

with nonnegative integer coefficients and c_(nu_q+1)>0.

**Proof, upper bound.** Between successive visits there is an exterior
directed path starting in N+(q) and ending in N-(q). Along that path, some
consecutive pair crosses from N+(q) to N-(q). Choose that pair. The exterior
paths are vertex-disjoint, so their chosen edges form a matching of size one
less than the number of visits.

**Proof, attainment.** A matching of size nu_q gives vertex-disjoint pairs
(a_i,b_i) and a legal collapsed walk

    q,a_1,b_1,q,a_2,b_2,q,...,a_nu,b_nu,q.

Insert every remaining exterior vertex without deleting any existing visit.
The standard tournament insertion argument still applies: prepend if the
new vertex beats the first vertex, append if the last beats it, and otherwise
use the first adjacent transition from beating the new vertex to losing to
it. Every comparison with a repeated q has the same orientation, so the
argument is well-defined. This produces a spanning collapsed walk with
nu_q+1 visits. A smaller matching gives any smaller positive visit count;
any child of order m can be split along a Hamilton path into r<=m nonempty
paths. Assigning m ordered clone labels to r nonempty increasing
runs has exactly r! S(m,r) choices, proving (G2). ∎

For k ordered path slots the same argument gives at most **nu_q+k** visits.
Every path visiting the module contributes one possible extra visit beyond
its exterior excursions. The bound is attained for the unconstrained query
allowing empty slots and sufficiently many child vertices: use the construction above in one slot and put one
singleton module visit in each of the remaining k-1 slots. Endpoint guards
at exterior names can lower the maximum, but cannot increase it.

Consequently nu_q+k coefficients suffice for any fixed ordered k-slot query
whose endpoint guards are exterior to the module. This is a bound on a fixed
query's run expansion, not a claim that joint endpoint incidence can be
forgotten. The dominant scalar exponential base is exactly nu_q+1; lower
exponential modes can have cancellations, so a full minimal recurrence order
is not inferred merely from that dominant base.

**Dual certificate.** The number nu_q is also the minimum size of an
exterior vertex set meeting every directed cycle through q. A set meets
all such cycles iff it meets every triangle through q, by the transition
argument above. This is exactly a vertex cover of L_q. A maximum matching
gives a cover of the same size: start alternating paths at unmatched left
vertices, let Z be all reached vertices, and take (left minus Z) together
with (right intersect Z). Every edge is covered and exactly one endpoint
of each matched edge is selected. A matching is also a lower bound on
every cover. Thus the grammar's run capacity has both a disjoint-cycle
certificate and a minimum deletion certificate.

**Capacity cannot change in a uniform affine return.** Suppose two strong
rooted cores satisfy H(W[q<-TT_m])=a H(T[q<-TT_m])+c for all m>=1, with
a>0 and fixed c. The functions r!S(m,r) are linearly independent on the
positive integers (their evaluation matrix at m=1,...,d is triangular
with diagonal 1!,...,d!). The constant function is the r=1 term. Hence
the largest nonzero run index is unchanged, and nu_q(W)=nu_q(T). For
the current return multiplier a=3 this proves that capacity preservation
is necessary, not just a convenient search restriction. More generally,
for two supplied cores, d=max(nu_q(T),nu_q(W))+1 literal evaluations
decide the uniform scalar identity. This statement does not assert the
existence of a legal complement word between those cores.

## G2. Capacity one has a physical cycle guardian

For a strong T with at least three vertices, nu_q>=1. The following are equivalent:

1. nu_q=1.
2. Every edge of L_q shares one vertex w.
3. Some w!=q meets every directed cycle through q.
4. In T-w, the strong component containing q is the singleton {q}.

To see 1=>2, two disjoint edges would be a matching of size 2. In a bipartite
graph, if two intersecting edges have different non-shared ends, any edge
avoiding their common end is disjoint from one of them. The other implications
follow because every directed q-cycle contains a triangle through q: take
an N+(q)-to-N-(q) transition along its exterior segment. The equivalence with
4 is the definition of a cyclic strong component.

This gives a common structural input grammar:

    finite strong core T;
    module slot q with nu_q=1 and a named guardian w;
    a positive transitive block size m and its ordered vertex addresses;
    retained exterior observer ports, including r=2 in the current bank;
    full core relation and selected-path program.

It strictly enlarges the previous one-in/out-neighbour rule. For example,

    T=(16,1,3,7,14), q=2, w=0

is strong; q has two in-neighbours and two out-neighbours, but nu_q=1 and
its run coefficients are c_1=9,c_2=4. Thus H(T[q<-TT_m])=4*2^m+1.
`port_rank.py` checks G1 independently on **all 5405 rooted labelled
tournaments of orders 1 through 5**, using matching, permutation/gap
enumeration, and exhaustive exterior deletions with reachability. This
finite check supports the proof; it does not replace it.

## G3. A common input type, with an exact return guard

The earlier common-frame theorem in
`adaptive_four_returns_20261103/LIFTING.md`, L1, works whenever two successive
actual path complements use the same ordered predecessor and successor at
q. A uniform birth preserves the module. Any finite word of such natural
operations, and explicit deletions of designated exterior helper vertices,
therefore commutes with transitive substitution at q. A deletion is an
additional operation and is never treated as a free class reset.

If the starting and ending cores are strong and have nu_q=1, both are objects
of the grammar in G2. Write their two scalar run coefficients as (c_1,c_2)
and (d_1,d_2). The supplied word implements 3H+c at **every m>=1** iff

    d_1=3c_1+c,             d_2=3c_2.                     (G3)

Thus two literal size evaluations certify the entire scalar return. Three
evaluations certify a fixed ordered two-slot response, since nu_q+2=3.
All evaluation divisions retain their factorial integrality guards.

This is a decidable certificate for an actual word, not an existence theorem
for every core. The grammar is closed under composition of accepted words;
we have not proved that both required offset words are available everywhere.
The full core and chosen paths remain necessary for deciding legal commands.

The previous six output families all belong to this grammar. The newly saved
continuations remain at the same module slot and old root, and satisfy G3.

Here are explicit common-grammar input certificates. The displayed guardian
is one choice; it need not be unique. All six retain observer root 2.

| Head | Module q | Guardian w | H(m), m>=1 |
|---|---:|---:|---|
| 4373 | 0 | 2 | 3564*2^(m-1)+809 |
| 4375 | 0 | 2 | 3564*2^(m-1)+811 |
| 4379 | 0 | 3 | 5538*2^(m-1)-1159 |
| 4381 | 0 | 3 | 5538*2^(m-1)-1157 |
| 4393 | 7 | 5 | 3720*2^(m-1)+673 |
| 4399 | 6 | 7 | 3228*2^(m-1)+1171 |

`port_rank.json` records the matching and guardian witnesses on the actual
labelled cores, not just these scalar rays.

For searches that widen the one-neighbour rule, `search.cpp` accepts the
opposite newborn direction only if all rooted triangle edges still share
an explicitly recorded guardian. The guardian may change; the module slot
and observer root do not silently change names.

## G4. Arbitrary-child naturality identifies the missing constant

Let S be any nonempty tournament with a selected Hamilton path R. In each
of the six established programs, replace q by S and expand both selected
paths using the **same R**. They remain legal. The first complement changes
S to tau_R(S), and the second changes it back to S; common entry and exit
hinges cancel exactly as in the transitive case. Uniform newborn relations
also commute with the substitution. This proves an actual relation-level
commuting square for every S, not just an equality of counts.

Let Q_r(S) count covers by r ordered **nonempty** directed paths. For each
fixed exterior query with collapsed coefficients c_r, the response is

    sum_r c_r Q_r(S).                                    (G4)

This is the inherited ordered-run substitution mechanism, now applied to
the actual six complement programs. The complete inherited theorem is
**THM-3121**, `01-canon/theorems/THM-3121-path-cover-walk-content-substitution-kernel.md`.
The scalar coefficient identities from the previous packet give, for all S,

    H(W[S]) = 3 H(T[S]) + c H(S).                        (G5)

Therefore the desired constant c is retained **iff H(S)=1**, equivalently
iff S is transitive. The latter follows because any nontransitive tournament
contains a directed triangle and induced extension cannot decrease H, while
a transitive tournament has exactly one Hamilton path. In particular a C3
child gives errors 4 and8 beyond the desired +2 and+4. A scalar repair would
have to remove exactly c(H(S)-1), with an actual realization guard.

The 1457 observer law extends to

    (H,X,Y) -> (3H+cH(S), 3X, 3Y+cH(S)).                 (G6)

For the 1459/+2 program it extends to

    (H,X,Y) -> (3H+2H(S), 3X+F_2(S), 3Y+F_2(S)),
    F_2(S)=Q_2(S)+2H(S),                                 (G7)

where F_2 allows empty slots. These are universal homogeneous transports.
Replacing the H(S) factor by the weight of **one selected child path** is
a precise new selector problem. A selected path name is already present,
but no physical operation extracting just that path's contribution has
been proved. We do not subtract positive counts formally and call the
result a tournament.

## The hostile that even the complete commutative response cannot repair

The older **Session36** result,
`session36_information_curl_contexts_20261023/context/CURL-AND-CLAN-COMPLETION.md`,
already proves that equality of labelled cyclic triples is equivalent to
equality of every vertexwise tournament-substitution Hamilton count. Its
structural classification input is explicitly cited there. This stronger
old result is inherited, not rediscovered here.

For the actual 81/27 pair we directly recover

    (Q_1,Q_2,Q_3)=(462915,9011898,80614764)

on both sides, and matching capacities `(4,5,5,5)` repeated three times.
G1 explains why their largest one-port run orders are respectively5 or6
at the corresponding labels. All six current root/Q4 observations agree
after inserting either child, by (G4); the exact 24-vertex target values are
obtained symbolically, not by a claimed direct 24-vertex enumeration.
Their inherited two-run matrix ranks 24/22 and three-run ranks 430/435 modulo 3
still differ. Full joint incidence and actual source transport remain outside
this commutative count grammar. Retaining them is essential to the return
problem, even after all scalar clone probes have been supplied.
