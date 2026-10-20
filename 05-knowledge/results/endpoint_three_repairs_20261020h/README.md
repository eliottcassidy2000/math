# Three-black detours for the endpoint-three layer

Research continuation, October20,2026. **PROVED structural theorem for
every n>=11**, with separate finite controls and explicit smaller hostiles.
This extends the [endpoint-two theorem](../endpoint_two_repairs_20261020g/README.md).
No claim of literature novelty is made.

## The theorem and its exact scope

Let P=(a,M,z) be a closed blue Hamiltonian word in a tournament T on n
vertices. Write K_Q(T)=T^op xor E(Q); a word is blue when its positional
reversal is an anti-automorphism, and black otherwise. Suppose

    n>=11,       d_T(a)=n-4,       d_T(z)=3.

**Then there exist three successive directed black K moves whose
literal labelled target is K_P(T).** Thus this marked blue line has its
two endpoints in the same black component and cannot be a bridge. There
is no additional assumption on the interior, the target's SCCs, or its
mixedness. Reading the line from its open partner gives the equivalent
endpoint score pattern (4,n-5).

The proof handles all placements and coincidences of the two mirrored
exceptional ports. It is a case partition by actual arcs, not by score
sequences alone. The general black-core connectivity conjecture and the
remaining endpoint layers are **OPEN**. This packet does not classify
every possible strong-to-strong blue line. The range6..10 is not claimed
by this structural theorem.

## The retained coordinates and exhaustive partition

Write P=(a,b,x,W,y,c,z), with W arbitrary, and let r_P be its reversal.
There is a two-element incoming set I in the middle such that

    N^-_T(a)={z} union I,       N^+_T(z)={a} union r_P(I).

The first middle vertex b cannot lie in I. The two sets I and r_P(I)
may overlap. An open middle means b->c; a closed middle means c->b.
These words describe the boundary chord, not a strong-component test.

The following rows are disjoint and exhaustive for n>=11. Every row
supplies three black words, with arbitrary orientations inside W.

| Middle and ports | Further split | Proof |
| --- | --- | --- |
| Open, y not in I | none | [R2](root.md#r2-open-middle-detours-away-from-the-penultimate-port) |
| Open, y in I, c not in I | none | [R4](root.md#r4-complete-penultimate-port-block-repair) |
| Open, I={y,c} | none | [Open end-pair ledger](probe_new_ledgers.md) |
| Closed, c not in I | c->x | [Theory section3](theory.md) |
| Closed, c not in I | x->c, x not in I | [Theory section4](theory.md) |
| Closed, c not in I | x->c, x in I, 3 not in I | [Theory section5](theory.md) |
| Closed, c not in I | x->c, I={2,3} | [Theory section6](theory.md) |
| Closed, c in I, y not in I | none | [Movable-cut theorem, section9](theory.md) |
| Closed, I={y,c} | x->c | [R3 two-block repair](root.md#r3-last-port-two-block-repair-with-a-second-exceptional-port) |
| Closed, I={y,c}, c->x | y->x | [Adjacent-port repair](lift.md) |
| Closed, I={y,c}, c->x | x->y | [Adjacent-port repair](lift.md) |

Here numeric vertices3,4,... mean positions in P, starting a=0. In each
ledger the three undirected path-edge sets XOR to E(P), so the final
object is the exact marked target, not merely an isomorphic tournament.
Blackness is unchanged across a move's own K pair; this permits checking
the third word at B=K_P(T). The proofs separately check directedness,
path-edge cancellation, and actual anti-reversal failure.

## The main new move: put the cut at the actual port

For the closed-middle case I={h,c}, h!=y, let k=r_P(h) and split
M=(L,R,c), with R starting at k. Then use

    Q1=(a,L,z,R,c),
    Q2=(L,z,R,a,c),
    Q3=(L,R,a,c,z).

L has at least two vertices. The crucial boundary join is last(L)->z
followed by z->k. Moving the cut with k discharges the guard that defeats
a fixed cut. Q1 cannot be blue because its reversal equates the true
c->b with the false y->a. Q2's endpoint score sum is n-2, not the n-1
required by an anti-reversal. In the third source, Q3 reversal equates
the false b->c with the true x->z. Thus all three words are black.

This mechanism generalizes beyond endpoint three: under the same closed
boundary, require z->b, a->y, and a further middle outneighbor of z.
Choose the first such outneighbor after b; its predecessor beats z.
The same proof works with any number of additional exceptional ports.
See the exact hypotheses in [theory.md](theory.md).

The adjacent end-pair needs a different repair. When I={y,c}, c->x and
x->y, a particularly short uniform ledger is

    Q1=(a,x,W,y,c,z,b),
    Q2=(z,W,y,x,b,a,c),
    Q3=(b,z,W,y,x,a,c).

For n>=9 its three blackness witnesses are literal arc contradictions:
first(W)->z versus x->c in T; a->P5 versus x->P3 after Q1; and
a->P3 versus x->z in B. The opposite y->x case uses a whole-block cut,
with the two possible z-inneighborhoods retained explicitly. These two
cases remove the last closed end-pair obstruction.

## Sharp boundaries and what seventeen means here

For a general closed blue endpoint pair (n-1-k,k), k>=2, the one-step
cyclic cut is necessarily black when n>4k-2. This bound is sharp for
every k. The proof composes two path reflections: their product shifts
positions by two, so its vertex orbits must induce regular tournaments.
The endpoint's limited incoming capacity bounds those orbit sizes.
[R1 gives the infinite extremal family](root.md).

At k=3 the sharp order is10. The displayed order10 source has a blue
cyclic cut but nevertheless has another three-black detour, recorded in
[probe_search.json](probe_search.json). Thus the sharpness of a method
is distinct from an actual metagraph bridge.

A second fixed-word obstruction occurs at order17. The proposed final
word has reflection-product orbit lengths3,3,5,3,3; a port may balance
the three-element orbit of a. [A3 in lift.md](lift.md) proves that17
is the only remaining order above10 for that word and gives an exact
blue example. The movable-cut theorem bypasses it. This is an explicit
permutation explanation for17 in this problem, not a numerical transfer
from the earlier seventeen-port or105-matching constructions.

There is also a necessary low-order warning. The six-vertex source rows
`(6,4,8,51,39,7)` have I={3,4}, but the blue partner is C3->C3 and has
no black Hamiltonian words. Hence the present unfiltered theorem cannot
simply be extended to every n>=6. Small strong/mixed filters and any
remaining orders must be stated explicitly.

## The requested arithmetic and foundations connection

The separately requested session completed and pushed a
[marked-graph and arithmetic packet](../trijection_reset_20261020/README.md).
We read its main proof and recovered its exact maps here.

For a tournament with a marked path, nonadjacent position pairs i<j
map to graph pairs {i,j-1}. Choosing backward arcs gives a simple graph
on ordered [n-1], or the same bits give a tournament on ordered [n-1].
This is an actual trijection. The present move K_P becomes graph
complement in these coordinates, while blue means invariance under
order reflection. Changing the Hamiltonian word changes the coordinate
system; that sidecar is exactly what the detour ledgers retain.

The ten proper nontrivial divisors of p^2qr admit the faithful map

    p^a q^b r^c -> (a-b-c, b-c).

Rows have sizes3,4,3; the negative, zero and positive first-coordinate
regions have sizes4,2,4. Divisor complement negates both coordinates,
and swapping q,r negates the second. Thus both requested decompositions
live on one object, with explicit symmetries. The three prime shapes
solving F=U+S are p,p^3,p^2qr with p,q,r distinct in the last shape.
No divisor-to-tournament detour map is asserted.

The square calculation distinguishes two actual operations: removing
the diagonal from a25-square leaves600 cells, whereas removing one
distinguished row and column leaves576. The balanced path-deletion
fibres4,2,1 come from keeping the exceptional pointed state. Reset-debt
composition similarly requires the multiplier gauge and legal source
domain. The useful common mechanism is retaining the data needed to
transport an operation; scalar equalities alone do not do that.

The incoming [fixed-anchor and endpoint-response work](../zero_heavy_ports_20261020g/README.md)
provides a useful hostile: equal complete path-cover profiles, scores,
feedback minima and automorphism cycle indices still permit different
maximum one-vertex H responses,159 versus147. Its missing data are actual
component endpoints. Likewise the incoming
[unit-gap matching boundary](../charge_palettes_boundary_20261020g/zero_heavy.md)
preserves full weight words and the roles of unit endpoints. Our two
mirrored ports are the analogous retained incidence data for K moves.
These observations connect the mechanisms without identifying the
different state spaces or importing their scalar105 into this proof.

## Validation

The three collaborating proof lanes were cross-read. The final partition
was independently implemented with directed arc sets in [root_check.py](root_check.py),
without importing any producer or previous checker. It replays every
marked source at n9 (7680), plus10976 deterministic sources at n11..60,
with no isomorphism, strong-partner, or mixedness filter. All larger
replays are three-black detours. The small run records blue template
steps rather than pretending the larger theorem applies there.
Normal and optimized Python give the same SHA256 receipt.

Separately, [probe.md](probe.md) certifies all8427 strong-to-strong mixed
marked presentations at n6..9:11,96,640,7680 respectively. Only eleven
presentations need a first-success search beyond the structural bank,
each resolved within82 candidate pairs. Their full selected certificates
are replayed by an independent set-arcs implementation. The unique
excluded source is the pure-target order-six boundary displayed above.
This **FINITE-EXACT** statement does not fill the unclassified order10.

The companion theory, lift, and probe checkers separately verify their
literal constructions, affine constraints, and hostile examples. These
finite controls support the structural proofs; they do not replace the
all-order arguments or extend their stated quantifiers.

Reproduction, from the repository root:

    python3 05-knowledge/results/endpoint_three_repairs_20261020h/root_check.py
    python3 -O 05-knowledge/results/endpoint_three_repairs_20261020h/root_check.py
    python3 05-knowledge/results/endpoint_three_repairs_20261020h/theory_verify.py
    python3 05-knowledge/results/endpoint_three_repairs_20261020h/lift_check.py
    python3 05-knowledge/results/endpoint_three_repairs_20261020h/probe_complete.py

The inherited theorem remains the endpoint-two result. This packet is
the current endpoint-three proof; exploratory fixed-word limitations
inside the sidecars must not be mistaken for the final frontier.
