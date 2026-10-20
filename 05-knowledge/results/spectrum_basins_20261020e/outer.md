# Outer-SCC detours and the endpoint-rigidity reduction

**PROVED:** every blue line incident to a nonstrong mixed tournament has
a black-only detour of length at most five. The new outer-SCC argument
uses the proved rigid-strong alphabet, itself a deduction from the CITED
Benhocine–Wojda theorem. The whole black core is not claimed connected:
blue lines with two strong endpoints remain outside this result.

The new part is the arbitrary-order outer-SCC theorem in Sections 2–6.
Its finite controls are separate from the universal proof.

## 1. Inheritance and the obstruction retained

The closest mechanisms are the fixed-endpoint orbit and common-hub lemmas
in [the previous moves note, Sections 2–4](../scc_basins_20261020d/moves.md),
and [the generalized terminal ear](../deficiency_three_20261020/ear_probe.md).
The latter needs both a closing arc and a terminal chord. They cannot be
silently conflated: the strong R4 presentation `(6,12,8,1)` has the closing
arc `3->0` but lacks `3->1`, and has only one Hamiltonian path ending at 3.
The regular cyclic R5 is the different hostile to omitting the ear's
rotation-automorphism guard.

The least-used sidecar is now the ordered SCC decomposition after deleting
one endpoint. Concept board: endpoint path orbits; near-source deletion;
rigid SCC words; consecutive cyclic triangles; two-port path ledgers;
black bypasses of blue lines.

Write `K_P(T)=T^op triangle E(P)`, where P is a directed Hamiltonian path.
A path is blue exactly when its positional reversal is an anti-automorphism.
For a K move, the outscore vector changes by

    d'=(n-1)1-d+e_first(P)-e_last(P).                    (1)

Thus an endpoint-score sum different from `n-1` proves a path black. All
path words below describe actual directed paths, not arbitrary edge-XOR
identities. Their targets are literal labelled K targets.

Fix a blue presentation of

    U=S^op -> Z -> S,

where S is strong of order r>=4. Read S along the prescribed right-hand
path Q; the left presentation is its reversed converse. Z may be empty
and its chosen path is blue. Put n=2r+|Z|. Every vertex in the right S
has the same score in U as in S, hence at most r-1.

The source-target map is complementing the edges outside this prescribed
path. The preserved predicate is existence of a black path between these
two metagraph vertices. Full adjacency, path order, and the positional
blue predicate are retained; replacing S by only its score sequence would
lose the required information.

## 2. Endpoint rigidity and its recursive structure

Call `(S,z)` endpoint-rigid if all Hamiltonian paths of S ending at z are
in one orbit under `Aut(S)_z`. Equivalently, reading S along any such path
gives the same adjacency matrix.

If `(S,z)` is not endpoint-rigid, choose a second terminal path in a
different orbit. Change only that right-hand path in the blue presentation
of U. The new presentation is black and retains the global first, second,
and last vertices. The inherited nonstrong common-hub lemma gives a
three-black detour around the original blue line. Hence assume endpoint
rigidity from now on.

**Lemma 1 (closure and nonrotation).** Every HP of S ending at z closes
into a Hamiltonian cycle. If S has order at least four, cyclically rotating
any such path by one position is not an automorphism of S.

**Proof.** Camion's theorem gives an HC; cut it just after z. Endpoint
rigidity transfers its closing arc to every HP ending at z. If full cyclic
rotation were an automorphism, it would be transitive on vertices. Move
the endpoint of any HP to z, then use endpoint rigidity: all HPs of S
would belong to one Aut orbit. The [all-order rigid-strong classification,
Section 3](../merged_metagraph_finite_pure_atoms_20261020/FINITE_ATOMS_AND_RIGID_CLOSURE.md)
implies that the only nontrivial strong tournament with one HP orbit is
C3. That contradicts r>=4. QED.

That classification uses the CITED Benhocine–Wojda Hamiltonian-bypass
theorem for order at least seven and explicit smaller cases. It does not
depend on a conjectural pure-blue strong alphabet. Camion is the same
CITED classical input already used in the descent proof.

**Lemma 2 (missing-chord decomposition).** Let Q=(a,b,...,z) be an
endpoint-rigid strong presentation with `b->z`, and assume Q contains no
consecutive cyclic triangle after its first vertex. Then

    S-a = W -> D,                                       (2)

where W is a nonempty ordered word of singleton and C3 SCCs, D is the last
strong component containing z, and `(D,z)` is endpoint-rigid. The vertex a
beats every vertex of S except z. If D is nontrivial, its restricted path
ends at z and closes. In particular D has smaller order than S.

**Proof.** Write Q=(v0,...,v_(r-1)), with a=v0,b=v1,z=v_(r-1).
By Lemma 1, `z->a`; by assumption `b->z`. A legal cyclic rotation of a
proper prefix containing a and b would, by endpoint rigidity, be an
automorphism fixing z. Its cycle puts a and b in the same orbit, which
is impossible because z distinguishes them.

First there is no index i with `v_i->a->v_(i+1)`. Such an index lies
between 2 and r-3, and left rotation of the prefix through v_i would give
the directed HP

    (b,v2,...,v_i,a,v_(i+1),...,z),

contradicting the preceding orbit observation. Consequently the
outneighbors of a in Q occur in a single initial interval
`v1,...,v_t`, followed by inneighbors. Here t>=1 and t<=r-2.

If t<=r-3, the absence of a consecutive cyclic triangle after a gives
`v_t->v_(t+2)`. Right rotation of the prefix through v_(t+1) then gives
the directed HP

    (v_(t+1),a,b,...,v_t,v_(t+2),...,z).

Its prefix cycle again contains a,b and fixes z, a contradiction.
Thus t=r-2: a beats every vertex other than z.

Deleting the forced initial a gives a bijection
between HPs of S ending at z and HPs of S-a ending at z. Every automorphism
fixing z fixes a (it is the only possible initial vertex of such a path).
Restriction and extension by fixing a identify the two endpoint
automorphism groups. Thus S-a inherits endpoint rigidity. It cannot be
strong: Lemma 1 applied there would imply `z->b`. Consequently its last
SCC D is proper and the preceding SCC word W is nonempty.

Vary the HP in one preceding SCC while fixing every other component path.
Endpoint rigidity supplies an automorphism of S-a relating the two full
paths; ordered SCCs must be preserved. Thus that SCC has one HP orbit and
is 1 or C3. Vary only the path in D to obtain endpoint rigidity of (D,z).
QED.

The no-consecutive-triangle hypothesis is used only in the second prefix
rotation. It is precisely the branch left after Lemma 3 below. An earlier
attempt tried to prescribe both endpoints of a Hamiltonian path in the
first SCC of S-z by cutting an HC. That is invalid when S-z is strong;
the two explicit prefix rotations above are the repair.

## 3. A consecutive-triangle three-move ledger

This lemma is unconditional and will also handle C3 terminal components.
Suppose the prescribed HP of S is closed and its one-step rotation is
not an automorphism. Suppose moreover that it contains a consecutive
cyclic triangle after its first vertex. Write

    Q=(a,L,b,c,d,R),        b->c->d->b,

where L,R are possibly empty lists. Let z be the last vertex of Q, so
z=d if R is empty. Let A be the concatenated left-S and middle path.
Use the three K paths

    P1=(A,L,b,c,d,R,a),
    P2=(c,A,L,b,d,R,a),
    P3=(c,A,a,L,b,d,R).                              (3)

**Lemma 3 (PROVED).** They are directed black paths and their composite
equals `K_(A,Q)(U)` on the original labelled vertices.

**Proof.** P1 rotates the closed right path. For P2, K along P1 turns
`A_first->c` into `c->A_first` and turns the nonpath arc `d->b` into `b->d`.
Its other links are retained path arcs. For P3, `c->A_first` and `b->d`
are retained, while the links `A_last->a` and `a->first(L,b)` have each
been complemented twice and return to their original directions.

Every internal edge of A,L,R is traversed three times. The added links
`c A_first`, `b d`, `z a`, and `A_last first(L,b)` each occur twice; the
remaining links are exactly the original path. Thus the XOR of the three
path-edge sets is E(A,Q), proving the target identity.

P1 is black by the rotation guard: relative to the original blue
reversal, rotating only the right path is blue iff that rotation is an
automorphism. Put `s_v=d_S(v)`. In the second and third intermediate
tournaments, the endpoint-score excesses over n-1 are respectively

    n-2-s_c-s_a >= n-2r+1 > 0,
    s_c+s_z+2-n <= 2r-1-n < 0.                       (4)

Here two distinct r-vertex scores sum to at most 2r-3. Both moves are
black. QED.

## 4. A guarded ear inside a terminal strong component

Let D be a terminal endpoint-rigid strong component arising after one or
more applications of (2), with chosen path `(a,b,C,z)`, C nonempty. If
`z->b`, set F to the entire path before this D, including the original
left S and middle. Use

    (F,b,C,z,a),    (C,F,b,z,a),    (C,F,a,b,z).        (5)

These are the generalized-ear R paths converted to K paths by reversing
the second path. Their directedness needs only the original path, the
arcs `z->a,z->b`, `F_last->b`, and `F_first->last(C)`; all hold. In
particular it does not require every vertex of F to beat z. The original
leftmost vertex still beats every right-S vertex.

The edge-XOR identity is unchanged by placing the earlier right-hand
components into F. The first move is black because a nonautomorphic
rotation of D cannot be the restriction of an automorphism of S; all
other right-hand vertices are fixed. The remaining two moves have the
same strict score bounds (4), using their vertices' full S scores.
Thus (5) is a three-black detour to the original labelled K target.

## 5. The five-move terminal near-transitive ledger

Suppose a right-hand suffix is

    B=a -> V -> c -> z,

with the arc between a and z reversed, so `z->a`, and with V any nonempty
tournament on k vertices with its chosen HP. All earlier right-hand
vertices beat B except that they may lose to z. Let F be everything
before B. The first vertex of F is the original leftmost vertex.
Use

    (F,V,c,z,a),
    (z,F,V,c,a),
    (a,z,F,V,c),
    (a,c,z,F,V),
    (c,z,F,a,V).                                     (6)

**Lemma 4 (PROVED).** These are five black K moves whose exact composite
is K along the original path `(F,a,V,c,z)`.

**Proof.** Every word retains F and V as forward contiguous paths. The
needed cross-arcs are `F_first->z`, `F_last->a`, `F_last->first(V)`,
`a->V`, `V->c`, `c->z`, and `z->a`. The successive complements supply
the other displayed links. The path-edge XOR is the original path:
internal F and V edges occur five times, unwanted boundary links occur
evenly, and exactly the four original links between F,a,V,c,z occur oddly.

The first move rotates only B. This induced B is nonregular: c and z
have score one inside B, while a has score k+1>=2. A full cyclic rotation
cannot be its automorphism, hence cannot extend by the identity to an
automorphism of S. The original presentation is blue, so this move is
black.

Write q=d_S(z), and let t be the internal V outscore of its final vertex.
Then `q<=r-k-2`, `0<=t<=k-1`. The S scores of a,c,and that V vertex are
respectively k+1,1,t+2. Applying (1), the endpoint-score excesses of moves
2–5 are

    n-k-q-3,   k+3-n,   n-k-t-3,   t+4-n.             (7)

The first is at least n-r-1>0; the second and fourth are negative; and
the third is positive. Thus all five moves are black. QED.

## 6. All-order closure at nonstrong sources

**Theorem (PROVED).** For every r>=4 and every blue presentation of
`S^op->Z->S`, where S is strong of order r, the corresponding blue line
has a black detour of length at most five. The algorithm is:

1. At the original S, use an inequivalent same-terminal path if available.
2. In the endpoint-rigid case, a consecutive cyclic triple after the
   original first vertex supplies (3). Henceforth assume none exists.
3. Use an ear whenever the current terminal strong component has its
   terminal chord. Otherwise apply Lemma 2 and recurse on D.
4. If D=C3, its consecutive triple occurs after the original first vertex;
   use (3). If recursion ends at singleton z, inspect the final rigid W.
   A C3 atom supplies (3). If W ends in a singleton c, write W=V->c.
   When V is nonempty use (6). When V is empty, `(a,c,z)` is a consecutive
   C3 after the original first vertex, unless S itself is C3, excluded.

The no-consecutive-triangle condition is inherited by every restricted
terminal component. In fact it excludes all C3 blocks encountered after
the original first vertex; those cases were already solved at step 2.
The induction strictly reduces the terminal strong component. It produces
at most five moves, not five per recursion level: recursion only selects
a final boundary ledger; it performs no intermediate metagraph moves.

Every nonstrong self-converse tournament has converse-paired outer SCCs.
For outer order at least four use this theorem. For outer order one use
the [all-order singleton-padding theorem, Sections 9–11](../scc_basins_20261020d/moves.md).
For outer C3, mixedness forces its middle to be mixed because C3 has one
path orbit; the inherited mixed-center padding hub gives three black
moves. No strong SCC has order two. Thus every blue line at a nonstrong
mixed source has a black detour of length at most five. In particular
none is a mixed-to-mixed blue bridge. Blue loops are also covered as
walks; they were never bridges to begin with. The construction preserves
the actual two endpoints in the converse-merged metagraph.

## 7. Reproducible finite information

`outer_probe.py` enumerates all fixed-path presentations through order 7,
tests strong connectivity, and compares every path ending at the prescribed
terminal vertex by its complete induced matrix. Counts are:

| order | strong presentations | endpoint-rigid | rigid with missing chord |
|---:|---:|---:|---:|
| 4 | 5 | 3 | 1 |
| 5 | 50 | 10 | 3 |
| 6 | 903 | 29 | 7 |
| 7 | 30,773 | 83 | 19 |

All missing-chord rigid presentations in this universe have their first
vertex dominating every vertex other than the terminal one, even without
the no-consecutive-triangle condition. That stronger unrestricted
statement is only FINITE-EXACT here; the proved restricted lemma suffices.

Reproduce from the repository root:

    python 05-knowledge/results/spectrum_basins_20261020e/outer_probe.py
    python -O 05-knowledge/results/spectrum_basins_20261020e/outer_probe.py

The preserved sidecar is the complete list of chord-failing matrices in
`outer_probe.json`, including endpoint path counts and score vectors.

`outer_verify.py` constructs the complete black path rather than just
testing its selection criterion. It checks every intermediate directed
path, every positional blue predicate, and the literal final adjacency
matrix. All 31,731 strong fixed-path presentations of orders 4–7 were
replayed with three middle presentations (empty, singleton, C3), giving
**95,193** exact controls. The constructor branches were:

| order | alternative endpoint orbit | consecutive triangle | terminal ear | five-move suffix |
|---:|---:|---:|---:|---:|
| 4 | 2 | 2 | 0 | 1 |
| 5 | 40 | 8 | 1 | 1 |
| 6 | 874 | 24 | 3 | 2 |
| 7 | 30,690 | 73 | 7 | 3 |

The five-move lemma was also checked on every arbitrary fixed-path V of
orders 1–6, with zero, one, or three earlier near-source layers and an
empty or singleton middle: **6,600** further labelled controls. A separate
family checks recursive depths 0–4, reaching outer order 12. A deeper
naive endpoint-path enumeration was stopped; no coverage beyond those
five depths is claimed. The recursion has no theoretical depth cutoff.
The complete small replay and the additional family controls passed with
normal Python and with `-O`; no check uses an `assert` statement.

    python 05-knowledge/results/spectrum_basins_20261020e/outer_verify.py
    python -O 05-knowledge/results/spectrum_basins_20261020e/outer_verify.py

The independent `outer_inventory_audit.py` replays the parent's
cross-feedback four-path-reversal identity on all 1,100 fixed-path A
presentations of orders 1–6. It verifies the general labelled target and
middle score charges for every A; all four moves are additionally checked
black for the 77 presentations whose first vertex is a source. Its
separate JSON receipt does not turn that sufficient guard into a claim
about arbitrary A. Normal and optimized Python both passed.

Independent mathematical audit: the parent agent and `snark_reductions`
checked the repaired prefix-rotation argument and SCC induction. The
initial invalid endpoint-prescription shortcut was retracted before
promotion. The current proof uses no such shortcut and no unproved
pure-blue strong alphabet. See the separate `outer_audit` receipt for
the independently implemented local-ledger audit.
