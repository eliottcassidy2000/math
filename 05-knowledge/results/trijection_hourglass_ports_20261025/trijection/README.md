# Movable-cut repairs as affine graph certificates

**PROVED:** a lossless affine compiler for literal path-complement ledgers;
an exact finite-menu covering criterion by inconsistent signed cycles;
and an overlapping-prefix extension of the endpoint-four odd-cut repair.
**FINITE-EXACT:** the bounded independent controls in `receipt.json`.
The unrestricted endpoint-four layer and black-core connectivity remain
open. A failed menu is not a disconnected black component.

## Inheritance and information ledger

The closest mechanism is the ordered-graph observer and path-agreement
trijection in
[the tree-gauge packet](../../tournament_graph_trijections_20261020/TRIJECTIONS_AND_TREE_GAUGES.md),
especially equations (1), (14), (19) and its Theorem 2. The inherited
odd-prefix ledger and its four honest final-blue hostiles are in
[the endpoint-four frontier](../../endpoint_four_frontier_20261021/README.md).
The [second-port/growing-prefix repair](../../ports_owners_f1_20261025/endpoint.md)
closes a further branch. The all-order endpoint-three theorem and audited
order-ten head are already inherited; this note does not rediscover them.

Canonical hostile: path-edge XOR alone permits an illegal intermediate
join. Corrected near miss: complementing in two different tree charts is
not the same operation. Least-used sidecar: a signed constraint graph on
the *unordered pairs*, retaining the orientation correction at each pair.
No cut normalization, clone contraction or isomorphism quotient occurs.

The active objects are literal paths, observer bits, exceptional ports,
anti-reflection equations, signed constraint cycles and overlapping
permutation supports. Source and target are labelled tournaments; the map
to pair bits is invertible. It preserves every arc and hence the actual
target. Forgetting the observer or any constant parity destroys this
property. Its cheapest decisive hostile is a reversed path join.

## 1. An exact affine compiler

Fix the observer `0<...<n-1` and let `x_{uv}=1` mean `u->v` for `u<v`.
Write `E` for all unordered pairs and `p(Q)` for a word's path-edge set.
Each `K_Q` complements precisely `E\p(Q)`. Before the i-th word of a
specified ledger, put

    d_i = (i-1 mod 2) 1_E + sum_{j<i} 1_{p(Q_j)}        (over F2).

The actual intermediate orientation bits are `x+d_i`. Therefore:

1. A join `u->v` in `Q_i` is legal exactly when

       x_{uv} = 1 + [u>v] + d_i({u,v}),

   where the subscript is sorted. All legality conditions are affine.
2. If `r_i` is positional reversal of `Q_i`, that move is blue exactly
   when, for every sorted pair `e={u<v}`,

       x_e + x_{r_i(e)}
       = 1 + [r_i(u)>r_i(v)] + d_i(e)+d_i(r_i(e)).      (1)

3. For an odd number m of moves, the literal target is `K_P(T)` for
   every source iff `sum_i p(Q_i)=p(P)`. For arbitrary m the exact
   criterion is

       (m mod 2)1_E + sum_i p(Q_i) = 1_E + p(P).

These statements follow by substituting the bit flips; they are iff
claims, not sufficient numerical proxies. In particular (1) tests all
arcs, including fixed pairs, rather than a degree sequence.

Each equation has at most two variables. Adjoin a constant vertex `*`
with `x_*=0` and make an edge `e--f` labelled c for every equation
`x_e+x_f=c`; a forced bit is an edge to `*`. The system is consistent
iff every cycle has XOR label zero (parallel edges and loops included).
To prove this, propagate a value along a spanning forest. Path
independence is exactly the cycle condition. Each component disjoint
from `*` contributes one free bit, giving `2^d` solutions. Thus a failed
blue test has a small, checkable contradiction certificate: a cycle of
equations with odd total constant. This is the graphic-matroid cycle
criterion from the trijection, now applied to *arc variables*, not to
the original tournament vertices.

This compiler decides a specified affine family without enumerating its
`2^d` tournaments. Its input must include every source filter and every
actual intermediate legality condition. It does not find a finite menu
covering all orders by itself.

## 2. Exact finite-menu coverage

Let F be a nonempty affine source family. Suppose r candidate ledgers
have already been proved legal on every member of F and to have the
same exact target. For ledger j let `B_{j,i}` be its i-th blue system
intersected with F. Then these ledgers cover F by all-black detours iff

    intersection_j B_{j,s(j)} is empty
    for every selection of one index s(j) from each ledger.           (2)

Proof: an uncovered source has at least one blue move in each otherwise
legal ledger; choose those moves. Conversely a point in one intersection
fails every ledger. This is distributivity of finite unions/intersections.
In the pair-constraint setting, each required empty intersection is
certified by an odd signed cycle. Three-move menus have at most `3^r`
such combinations, independent of the number of free source bits.
Redundant branches may be eliminated by conflict propagation.

When legality varies, first partition by the distinct boundary-join
bits. On each resulting affine cell keep only its legal ledgers and use
(2). This partition retains exact ports; replacing it by mere path-edge
XOR is invalid. The number of cells can grow without bound. Accordingly
(2) is an exact structural proof target for connectivity, not an asserted
small unrestricted recurrence or a solved connectivity theorem.

The same mechanism explains why several numerical lookalikes should
not be identified: a code quotient avoiding three crossing subspaces
and a repair menu avoiding blue subspaces share affine incidence, but
need an explicit map preserving their forbidden loci before one proof
can transfer to the other.

### Sharp boundary of the pure reflection test

For two positional reversals r,s on n>=2 labels, there exists a
tournament in which both are anti-automorphisms **iff rs has odd order**.
This is an elementary existence statement without prescribed path arcs.

Necessity is the usual odd-order automorphism obstruction. For
sufficiency put g=rs and suppose its order m is odd. The generated group
consists of the automorphisms `g^k` and the intended anti-automorphisms
`g^k r`. All the latter are conjugate to r, since multiplication by2
is invertible modulo m; hence each has at most one fixed vertex.
On any unordered-pair orbit choose an orientation and propagate it,
preserving orientation for g and reversing it for r. This is
well-defined: an element `g^k` stabilizing a pair cannot exchange its
ends because its order is odd; an element `g^k r` stabilizing a pair
must exchange its ends because it cannot fix two vertices. These are
exactly the compatible stabilizer actions. Each pair orbit supplies
one independent orientation choice.

Thus an odd-order product is a real limitation of a *permutation-only*
blackness proof, not evidence of a blue path or a failed detour.
Prescribed joins and ports can still make its affine blue system empty.
For example P=(0,1,2,3) and Q=(1,0,3,2) have the same reversal, so
their product is the identity and simultaneous anti-automorphisms exist.
They cannot both be directed paths: P demands0->1 and Q demands1->0.
The compiler retains exactly this missing sidecar.

## 3. Overlapping odd prefixes: three positions of extra reach

Use the inherited endpoint-four setting:

    P=(a,b,x,3,...,y,c,z),  a=0,b=1,x=2,y=n-3,c=n-2,z=n-1,
    N^-_T(a)={z,h,y,c},  N^+_T(z)={a,1,2,n-1-h}.

P is closed and blue; all original labels remain distinct. Let ell be
positive odd, `L=ell+3`, `D=(3,...,ell+2)`, `U=(ell+3,...,c)`. Assume

    n >= max(11,2ell+3),       (ell+3)->x,
    h not in {2,ell+2},
    and, if ell>=3, h not in {4,ell-1}.                 (3)

Then the inherited words

    Q1=(a,x,D,U,z,b),
    Q2=(D,a,x,U,z,b),
    Q3=(D,a,b,x,U,z)                                  (4)

are a three-black detour to the same literal `K_P(T)`.

The predecessor required `n>=max(15,2ell+6)`. The new theorem permits the
two rotated endpoint supports to overlap in up to three positions and
uses the sharper cyclic-cut bound from the sibling endpoint proof.

**Legality and target.** Unchanged from the predecessor: the new joins
are controlled by `h!=2,ell+2` and `(ell+3)->x`; each added edge occurs
twice in the XOR. Internal D and U edges occur three times. No new
source condition concerning x versus y or c is imposed.

**Q1.** Its reversal is the one-step cyclic-cut reversal. Use the
sharpened endpoint argument in
[the sibling cyclic-prefix theorem](../endpoint/THEOREM.md): if this
reversal and r_P were both anti, shift by two modulo n would be an
automorphism. For odd n the tournament would be regular, forcing n9.
For even n=2m, the only incoming neighbours of a in its even-label
orbit are c and possibly h, because y,z are odd. Regularity of that
orbit gives `(m-1)/2<=2`, forcing n<=10. Thus Q1 is black at n>=11.

**Q2.** At `S=K_Q1(T)`, `N^-_S(z)={a,c,x,n-1-h}`. For ell=1 the
inherited degree mismatch applies. For ell>=3 reversal sends 4 to z
and a to `n-ell`. The weakened bound still puts `n-ell` in U, since
`n-ell>=ell+3`. The h exclusions give `S(4,a)=1` and
`S(n-ell,z)=0`, a literal anti contradiction.

**Q3.** Work at `B=K_P(T)`, where P is blue and
`N^+_B(a)={b,h,y,c,z}`. If Q3 were blue, `g=r_Q3 r_P` would be an
automorphism. Let A be the first L positions and Z the last L.
Directly from Q3, g is the product of the rotations by +3 on these
two ordered intervals (first apply the Z rotation). If they are
disjoint, the original even-cycle proof applies because L is even.
Otherwise put `j=2L-n`, so `j=1,2,3`.

The following cycle calculation finishes the proof. The calculation can
be obtained by cutting each +3 cycle at the j shared positions and
splicing its outgoing strands; this works also when 3 divides L.

| Case | Cycle information needed |
|---|---|
| `3` does not divide L, j=1 or3 | g is a single n-cycle |
| `3` does not divide L, j=2 | the cycle of a has length L-1 and meets `{b,y,c,z}` in at most two points |
| `3` divides L, j=1 or2 | at least one unchanged cycle of length L/3 remains |
| `3` divides L, j=3 | three cycles of length `2L/3-1`; the cycle of a meets `{b,y,c,z}` exactly in `{y}` |

Here is a detailed strand justification for the first two rows. When
3 is invertible mod L, each interval has one cycle. For j=1 their
single common strand joins them into one cycle. For j=2, when
`L=2 mod6` the cycle of zero is
`(0,3,...,2L-4,2,5,...,L-3)` and its named endpoint is c; when
`L=4 mod6` it is `(0,3,...,2L-5,1,4,...,L-3)` and its named endpoints
are b,y. Both have length L-1. For j=3, the three shared vertices
appear in the same cyclic order in the two +3 cycles, because their
coordinate lists differ by a translation modulo L. The first-return
permutation on the three cut strands is the square of this 3-cycle,
which is still a 3-cycle; hence the splice is one circuit. The same
tracing by residues modulo3 proves the
last two rows: at j<=2 an untouched residue class persists in each
interval; at j=3 matching residue classes splice once, giving length
`2(L/3)-1`. The cycle containing zero then contains y and no other
listed endpoint.

An automorphism of a tournament has odd order: an even-order power
would be an involution swapping a pair and reversing its arc. Every
orbit of a cyclic automorphism induces a regular tournament, because
the generator is transitive on that orbit.

In the first row regularity would give `d_B(a)=(n-1)/2`, whereas
`d_B(a)=5`. The only candidate n11 cannot occur in this row: j1 would
give L6, divisible by3, while j3 would give odd L7. Thus n>=13 here.
In the second row a has at most three outneighbors
inside its orbit (at most two named endpoints and h), while regularity requires `(L-2)/2`.
For L>=10 the requirement is at least4. The only additional case is
L8,n14; its a-orbit `(0,3,6,9,12,2,5)` contains only c among the four
named endpoints, hence at most two outports including h, whereas
regularity requires three.
In the third row L/3 is even. In the last row a has at most two
outneighbors (y and h), while regularity requires `L/3-1`; here
L is a multiple of6 and n=2L-3>=11 forces L>=12. Every case contradicts
the assumed automorphism. This proves (3)--(4).

## 4. Exact scope and next question

The signed-graph compiler exhaustively quantifies a family through
linear constraints; it does not turn a finite range of n into an
all-order proof. Its longer-cut probe finds no additional blue failures
at n19..40, but that is only finite evidence and is not included in
Theorem (3). The exceptional long-cut cells at n15 and18, ell9, h3 or6
retain genuine blue final words, so dropping all guards is false.

The next structural target is to classify the prefix/suffix rotation
product for overlap beyond three and use exact port occupancy, or to
cover the remaining affine cells by a second ledger using criterion
(2). The latter can advance connectivity without treating every
tournament completion as a separate object.

Run `python check.py` and `python -O check.py` in this directory. The
checker uses explicit exceptions, independently compares the compiler
with literal tournaments in small complete universes, tests the cycle
table and the full guarded theorem on bounded affine families, and
retains positive and hostile controls. It does not invoke a SAT solver
or a probabilistic source sample.
