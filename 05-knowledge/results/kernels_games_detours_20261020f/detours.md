# Endpoint-changing black detours at cyclic singleton extensions

Research continuation: October 20, 2026.

**PROVED:** an explicit three-black detour for the prescribed blue word
`(a,V,z)` in every cyclic singleton extension `a -> V -> z -> a`, when
`|V|>=4`; the same holds for an open blue middle of order three. The
construction retains the entire middle tournament and reaches the literal
labelled complement target. It is a family of strong-to-strong detours.

The closed three-vertex middle is a sharp exception: its target is the
pure-blue regular five-vertex tournament. Arbitrary strong mixed-to-mixed
blue edges and all-order black-core connectivity remain **OPEN**. In
particular, this theorem concerns the specified aligned word, not every
Hamiltonian word of the same unmarked cyclic extension.

## 1. Inheritance, current target, and retained coordinates

The audited [outer-SCC theorem](../spectrum_basins_20261020e/outer.md),
[its independent audit](../spectrum_basins_20261020e/outer_audit.md), and
[the pure-hub reduction](../spectrum_basins_20261020e/pure_hubs.md) already
settle every blue line incident to a nonstrong mixed vertex. The remaining
connectivity obligation is exactly the strong mixed-to-mixed blue case,
relative to the cited combined-metagraph connectivity theorem.

The closest mechanism is the common-cut identity in
[defect_transitions_20261020/blue.md, Section 5](../defect_transitions_20261020/blue.md).
Its pinned-endpoint assumption fails even on genuine five-vertex strong
blue edges. The corrected near miss is that a legal relocation is an
auxiliary move inside one Hamiltonian-path fibre, not a complement edge;
the [relocation theorem](../hamilton_path_relocation_connectivity_20261020/THEOREM.txt)
does not itself transport black components. The least-used sidecar is a
middle module together with its two actual boundary vertices.

The cyclic object `C3[1,V,1]` itself is inherited: see
[the symmetry-support note, Theorem 6](../continuation_width_symmetry_20261020/README.md)
and [the singleton-padded cyclic theorem, moves Section 7](../scc_basins_20261020d/moves.md).
The latter repairs blue lines at `1 -> C3[1,V,1] -> 1`; it does not already
give the strong-source detours proved here. The current board consists of
marked paths, endpoint charges, cyclic modules, actual arc guards, and
the black-component quotient.

Write `K_P(T)=T^op xor E(P)` for nonpath complementation along a directed
Hamiltonian path P. Its exact outdegree change is

    d_(K_P(T)) = (n-1)1-d_T+e_first(P)-e_last(P).       (1)

The path is blue precisely when its positional reversal is an
anti-automorphism. In particular, a blue path has complementary endpoint
scores summing to n-1. This is a necessary condition, not an iff test.
The path is blue at T iff it is blue at K_P(T).

All words below are evaluated after the preceding K moves. Internal
middle path arcs remain directed forward at every visit. An odd-length
ledger reaches K_P(T) literally whenever its path-edge XOR is E(P).
This preserves labelled arcs and the precise blue edge being bypassed;
no free change of marked path or tournament isomorphism is used.

## 2. A guarded three-move endpoint release

Let a blue closed word be `P=(a,b,W,z)`, where W has at least two vertices.
Write x and y for the first and last vertices of W. Suppose

    z->a,  b->y,  x->z.                              (2)

The first arc is the closing condition. Consider

    Q1=(b,W,z,a),
    Q2=(z,W,b,a),
    Q3=(a,z,W,b).                                    (3)

**Lemma (PROVED).** All three are directed at their successive sources,
and their composite is K_P(T). Q2 is always black. Q3 is black whenever
`d_T(a)+d_T(b) != n`. Q1 is black whenever the two-step cyclic shift of P
is not an automorphism; in particular `d_T(a)+d_T(b) != n-1` suffices.

For validity, Q1 is a cyclic cut of P. Its K move reverses the nonpath
arcs x->z, b->y, a->b, giving the three new joins z->x, y->b, b->a of
Q2. Here x!=y is essential. Q2 retains z->x and y->b while reversing
z->a, giving Q3. Internal W edges occur three times. The new boundary
edges az, zx and yb occur twice; the surviving edges are ab, bx and yz,
exactly the original boundary. This proves the literal target.

Put s=d_T(a)+d_T(b). Formula (1), and the original blue equality
`d_T(a)+d_T(z)=n-1`, give the endpoint-score excesses over n-1:

    Q1: s-(n-1),   Q2: -1,   Q3: s-n.                (4)

Thus `s not in {n-1,n}` is a simple complete sufficient guard for this
three-black ledger. At equality one must test the actual anti-symmetry;
a zero excess is not a proof of blueness.

### The whole-cycle symmetry guard has a parity boundary

If rho is the positional reversal of P and pi its cyclic successor,
the positional reversal of Q1 is `pi rho pi^(-1)`. Since rho is anti,
Q1 is blue iff `pi rho pi^(-1) rho=pi^2` is an automorphism. This is
**shift by two**, not necessarily shift by one. For odd n these generate
the same cyclic group. For even n the two parity orbits remain separate.

The smallest even hostile has outgoing rows

    (26,36,41,18,38,9),

and score vector `(3,2,3,2,3,2)`. Its standard path is closed and blue,
as is its one-step cyclic cut. Shift by one is not an automorphism (it
changes scores), while shift by two is. Thus a nonautomorphic one-step
shift alone would be an invalid whole-path blackness guard. The endcap
proofs which keep an entire opposite endcap fixed use a different
restriction argument and are unaffected.

## 3. The strong cyclic family and the open-middle branch

Let V have r>=3 vertices with prescribed blue path M. On n=r+2 vertices
define U by

    a -> V -> z,   z -> a,   P=(a,M,z).              (5)

Both U and B=K_P(U) are strong. U has the displayed Hamiltonian cycle.
For B, the cycle `(a,z,M,a)` is directed when r>=2: a->z and z->first(M)
are complemented cross-arcs, the middle path survives, and last(M)->a
is complemented. The source P is blue because M is blue.

Suppose M is open: its first vertex b beats its last vertex y. Write
M=(b,W), with first(W)=x. The guards (2) hold: b->y is the open-middle
arc, and every middle vertex beats z. Because r>=3, b beats its distinct
path successor and endpoint, so `d_V(b)>=2`.

Use (3). Since `d_U(a)=r` and `d_U(b)=d_V(b)+1`, the three score
excesses (4) are

    d_V(b),  -1,  d_V(b)-1.                         (6)

All are nonzero. This proves a three-black detour for every such middle,
including arbitrary nonstrong blue middles and open paths of strong ones.
For V=TT_r it repairs `TT_(r+2)` with its first-to-last arc reversed,
uniformly at every total order at least five.

## 4. Closed middles: two complementary local ledgers

Now r>=4 and M=(b,W,c) closes: c->b. Thus V is strong, W has at least
two vertices, and every V score is between 1 and r-2. Write x,y for
the first and last vertices of W, and let epsilon=[b->y]. Blue reversal
of M sends b to c and x to y, so epsilon=[x->c] as well.

If epsilon=0, apply

    Q1=(c,W,b,z,a),
    Q2=(z,c,W,b,a),
    Q3=(a,z,b,W,c).                                 (7)

If epsilon=1, apply

    Q1=(c,b,z,a,W),
    Q2=(z,c,b,a,W),
    Q3=(a,z,b,W,c).                                 (8)

**Direction and XOR.** In (7), epsilon=0 supplies c->x and y->b;
the other first-word joins follow from (5). K along Q1 complements
c->z and a->b, supplying z->c and b->a for Q2. The remaining second
word joins are retained. In (8), the first word uses c->b, b->z, z->a,
a->x; the second uses the complemented c->z and a->b plus retained
c->b and a->x. In both cases direct boundary cancellation gives

    E(Q1) xor E(Q2) xor E(Q3) = E(a,b,W,c,z).         (9)

For (7), the only added edges are cx, yb, bz, za, zc and ba; their
even cancellations leave ab,bx,yc,cz. For (8), cb and bz occur twice,
za occurs twice, ax occurs twice, and zc occurs twice; ab,bx,yc,cz
remain odd. Internal W edges occur three times in each case.
Equivalently, after Q1 and Q2 the tournament is K_Q3(B). Q3 is directed
in B along `(a,z,M)`, so it is directed there as well and reaches B
exactly. This also verifies the final joins without guessing orientations.

**Blackness of (7).** At U, Q1 has endpoint sum
`d_U(c)+d_U(a)=r+1+d_V(c)>n-1`. After Q1, the Q2 endpoints z,a have
scores r,0, so their sum is n-2. Hence Q1 and Q2 are black.

**Blackness of the first two words in (8).** Since V is strong, a is
the unique score-r vertex of U and z the unique score-one vertex. Any
anti-automorphism exchanges them. Their positions in Q1 are 3 and 2
(zero based), so a blue Q1 would force `n-1=5`, hence r=4. In that
boundary W=(x,y), but Q1 blue would require c->x iff b->y. The former
is false and the latter true because epsilon=1. Q1 is therefore black
at every r>=4.

After Q1, the endpoint-score excess for Q2=(z,c,b,a,W) is

    r-2-d_V(y) >= 1.                                (10)

Indeed y has the two distinct inneighbors b and its preceding W vertex,
so `d_V(y)<=r-3`. Thus Q2 is black as well.

**The common final word is black.** Test Q3 at its other endpoint B,
where it is `(a,z,b,W,c)`. In B the outneighbors of a are exactly b,z.
The vertex c has both b and y as inneighbors: c->b was a nonpath arc in
U, while y->c was a retained path arc. If Q3 were blue, its reversal
would send a to c and z to y; hence it would have to fix b, in order
to send a's two outneighbors onto c's two inneighbors. But b has
position 2 in Q3. It is fixed by reversal only when n=5. Here n>=6,
a contradiction. This proves Q3 black and completes both detours.

No property of V beyond the actual blue path and its indicated boundary
arcs was discarded. Its size and all interior arcs remain unrestricted.

## 5. Sharp small boundaries and the surviving obstruction

**Intrinsic endpoint corollary (PROVED).** At total order n>=6, every
closed blue word whose first vertex has score n-2 has the three-black
detour above. Its last vertex has score one by blue reversal. The closing
arc makes these vertices each other's sole exceptional neighbor: the
first loses only to the last, and the last beats only the first. Thus
the remaining vertices are exactly the middle V in (5), in their actual
marked order. Conversely, that cyclic extension has these endpoint
scores. Reading the same edge from its open complement partner gives
the equivalent sufficient endpoint pattern `(2,n-3)`, by (1).

The next precise endpoint layer has closed scores `(n-3,2)`. It retains
one additional incoming arc at the first vertex and its mirrored outgoing
arc at the last. Those two exceptional incidences are new retained ports;
the uniform middle-module proof above does not erase them. This is a
specific next construction target, not a claim that every remaining
strong blue edge lies in that layer.

If r=1, U=C3 and its target is TT3; both vertices are pure blue. If
r=2, U is the strong four-vertex tournament and its blue target is in
the same merged class: the blue edge is a loop, so the empty merged
walk suffices. Neither boundary is included in the literal three-move
theorem.

At r=3 the open middle TT3 is covered by (6). The closed middle is C3.
Then U=C3[1,C3,1], and its marked target is the regular five-vertex
tournament R5. This endpoint is pure blue, so no black detour exists.
This is the precise point at which the final-word proof above allows
blueness: b occupies the reversal's fixed middle position. The exact
checker independently verifies that every target Hamiltonian path is
blue. The exceptional object is also the inherited five-vertex bypass
hostile; it is not a new tournament classification.

For arbitrary strong sources the common-cut method is still only
sufficient. Given two open source paths with common endpoints and a
shared directed edge u->v, their partner cycles can be cut at that edge
to meet at `T xor {az,uv}`. This identity needs openness, not nonstrongness.
If both cut words are black and one source word is black, it gives a
three-black detour. Fixed-endpoint black alternatives can nevertheless
be absent. One five-vertex witness has outgoing rows

    (18,4,9,19,6).

The standard blue path is open, both endpoints of its blue line are
strong and mixed, and every path with standard endpoints is blue. Its
closed partner is the near-transitive family covered by Section 3.
Thus the new endpoint-changing ledger repairs an actual failure of
the old pinned-endpoint test.

The remaining general obligation is to find legitimate endpoint-changing
ledgers outside these guarded cyclic configurations. Relocation
connectivity inside a fibre does not supply their blackness or targets.
No conclusion about arbitrary strong blue lines, or about black-core
connectivity, is inferred from the finite success counts.

## 6. Exact replay scope

The standard-library checker independently constructs every fixed-path
blue middle of orders 3 through 8. Of 4,694 presentations, precisely the
closed order-three middle is the retained pure-target hostile. All other
4,693 receive their appropriate ledger: 2,347 open-middle cases and 1,173
in each closed-middle branch. It verifies every successive directed path,
literal positional blackness, both strong endpoints, and the complete
labelled target matrix. No Python assertion is used.
These same exhaustive blue-middle controls cover the intrinsic endpoint
corollary, whose coordinate characterization is proved above.

Its separate local probe exhausts open blue presentations of total orders
4 through 8 whose endpoints are both strong and mixed. Strong pure-blue
exclusion uses the inherited classification at orders at least six;
orders four and five test all Hamiltonian paths directly. The probe
retains exact hostiles rather than interpreting a failed certificate as
a disconnected pair. Loop detection uses an independent degree-pruned
labelled isomorphism search.

Reproduce from the repository root:

    python 05-knowledge/results/kernels_games_detours_20261020f/detours_check.py
    python -O 05-knowledge/results/kernels_games_detours_20261020f/detours_check.py

The adjacent `detours_check.out` records the finite universe and controls.
Normal and optimized runs both passed with byte-identical stdout (2,188
bytes). The checker also verifies the symbolic endpoint charges in all
three branches and the shift-by-two equivalence throughout its small blue
presentation universe.
