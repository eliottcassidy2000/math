# Cyclic attachment lattices, shared alignment, and exact elimination

23 October 2026. **PROVED scoped**, with **FINITE-EXACT** tests in
`check.py` and `receipt.json`. The carrier is a supplied named-cell partition,
actual internal charts, and actual binary cross blocks. A chosen cyclic group
must act by genuine internal automorphisms. It need not be the full internal
group; if it is smaller, the result classifies that restricted comparison
problem only. A cell name, a generator, and its actual action are retained.

## Inheritance and the new question

The nearest proved mechanisms are the shared-root criterion in
[Session41, SHARED_GAUGE_CRITERION.md](../../session41_owned_attachments_fibonacci_20261023/attachments/SHARED_GAUGE_CRITERION.md)
and the fixed-local-type quotient in
[odd-cell completion, NETWORKS.md](../../odd_cell_completion_20261023/NETWORKS.md).
The latter used full internal groups `C_p` and a vector-space quotient over
`F_p`. The canonical hostile is the same file's actual regular nine-cell
`R9=C3[C3,C3,C3]`: a nonuniform cross block leaves 27 target charts after the
entire source chart is fixed. A corrected near miss is independent edge
canonicalization: even a tree with two rigid five-cell attachments has five
different shared alignments. The least-used sidecar here is the *number* of
hidden middle-chart extensions, rather than only their existence.

Our live concepts are actual cyclic actions, one-sided kernels, integer
lattice quotients, owned boundary variables, and elimination multiplicity.
The extension below admits unequal and composite cyclic orders, nonfunctional
relations, and exact elimination. It does not turn a nonabelian chart group
into a cyclic group by numbering its elements.

## 1. Exactly when two congruences suffice

Let the row and column chart groups be `C_m=<g>` and `C_n=<h>`. On a cross
block act by `(a,b) B = g^a B h^(-b)` in permutation-matrix notation. The
same convention must be used for all incident blocks. Given blocks `B,B'`,
the relation

    R(B,B') = {(a,b): (a,b) B = B'}

is empty or a coset `o+H`, where `H=Stab(B)` is a subgroup of
`C_m x C_n`. This follows directly from the group-action identity; the
product is abelian, so left/right coset distinctions disappear here.

Take the inverse image `L` of `H` in `Z^2`. It contains `(m,0),(0,n)`.
There is a unique basis

    (h0,0), (s,v),   h0>0, v>0, 0<=s<h0.

Here `h0 | m`, `v | n`, and `h0 | s(n/v)`. To prove the basis form without
assuming a normal-form algorithm, take the least positive horizontal
element `(h0,0)`, the least positive second coordinate `v`, and reduce its
first coordinate modulo `h0`. Subtracting multiples of these two vectors
reduces every lattice point to zero. Containment of the two periods gives
the divisibilities. Conversely these conditions make a lattice containing
the periods and hence define a subgroup. Thus the normal forms exhaust the
abstract cyclic-pair subgroup cosets. This does not claim that every such
coset occurs for every fixed pair of actual tournament actions and blocks.

For integer representatives `x=a-o_a`, `y=b-o_b`, membership is exactly

    y == 0 (mod v),
    v*x-s*y == 0 (mod h0*v).                         (1)

Equivalently, after the first test, `x-s(y/v)==0 (mod h0)`. The divided
quantity is integral; no division modulo a nonunit is being performed.
In particular there are precisely

    |R| = mn/(h0*v)                                 (2)

allowed chart pairs. The number of target charts above each *admissible*
source chart is

    k_right = n*gcd(h0,s)/(h0*v),                    (3)

while the number of source charts over an admissible target chart is
`k_left=m/h0`. Formula (3) is the size of `{b:(0,b) in H}`: write
`b=v*t`, then `s*t==0 (mod h0)`. It may be greater than one, so (1) does
not silently assume functional propagation. The admissible source charts
form a coset of the subgroup of multiples of `gcd(h0,s)` in `C_m`.

## 2. The integer replacement for the prime-field alignment quotient

Fix a finite support `E` on named cells `V`, their cyclic orders `m_i`, and
one actual local block orbit on each edge. Choose an orientation and a
reference block for each edge. Let `L_e` be the preceding rank-two lattice.
All possible edge configurations form a torsor for

    W = direct-sum_e Z^2/L_e.

Use one shared integer variable `a_i` per actual cell, not one per incident
edge. The incidence map `D:Z^V -> Z^(2E)` sends `a` to the ordered endpoint
pair `(a_i,a_j)` on every edge. Then joint attachment classes are exactly

    Q = Z^(2E) / (direct-sum_e L_e + image(D)).       (4)

This is a finite abelian group. An offset vector `o` is globally realizable
by chart changes iff `o` belongs to the denominator lattice in (4). If it
is realizable, the number of global chart witnesses is exactly

    K = (product_i m_i) * |Q| / product_e(h0_e*v_e). (5)

Proof: each local choice is a coset of its stabilizer, hence precisely one
element of `W`. The same cell action on all incidents is `D`, so the orbit
set is its cokernel. Cell periods map into the local lattices, making this
an action of `product C_(m_i)`. Every nonempty witness set is a coset of
the kernel of this homomorphism. The first isomorphism theorem gives (5).
No regularity assumption on cross blocks is needed, though all tournament
controls below use actual biregular blocks.

For a checkable integer certificate, make `M` with the two embedded basis
columns of each `L_e`, followed by the cell-incidence columns of `D`.
`M` has full row rank. Its Smith invariant factors describe `Q`; their
product is `|Q|`, also the gcd of the maximal minors. A proposed offset
belongs to its column lattice iff adjoining that offset as one further
column leaves this gcd unchanged. An explicit integer solution `Mz=o`
is a positive witness. A Smith/Hermite residue is a negative witness.
The equality-of-index criterion is sufficient because the old lattice is
contained in the augmented lattice, with finite index. Rational row
reduction alone cannot test this membership: it loses torsion.

This specializes exactly to `W/image(D)` over `F_p` in the prior packet.
Unlike a dimension, the Smith factors preserve composite torsion and
unequal periods. The exact information that was formerly `delta` field
coordinates becomes the finite group `Q`; its cardinality is the minimum
alphabet size of a lossless label above the supplied local types.

## 3. Eliminate a hidden middle chart without forgetting how many

Let `R12=o12+H12` be a relation on `C_a x C_b`, and
`R23=o23+H23` on `C_b x C_c`. Their composition consists of `(x,z)` for
which one shared `y` satisfies both. It is either empty or a coset of a
subgroup of `C_a x C_c`, so it again has the two-congruence form (1).
Indeed the feasible triples form an affine coset of a subgroup of
`C_a x C_b x C_c`, and composition is its endpoint projection.

Every feasible endpoint pair has the SAME number of middle extensions:

    multiplicity = gcd(k_right(H12), k_left(H23)).  (6)

The difference of two middle witnesses lies in the intersection of the
two one-sided kernels inside `C_b`. Two subgroups of a cyclic group have
intersection size the gcd of their orders. This proves (6), even with
offsets. Emptiness still must be checked: the projected middle cosets can
be disjoint. Dropping this multiplicity preserves existence but corrupts
counting and weighted sums. Degree-two boundary contraction can therefore
replace two factors by their composed relation and this scalar, provided
the eliminated variable has no other incident constraint or marked query.

### A stronger cyclic-only result: exact fill-in at every degree

Suppose a middle owner `y in C_b` has incident relations
`R_i subset C_b x C_(a_i)`, all nonempty affine subgroup cosets. For any
allowed value of `x_i`, its possible middle values are a coset of the
fixed one-sided kernel `K_i={y:(y,0) in H_i}`. Any subgroup of `C_b` is
the multiples of a divisor `d_i` of `b`, so these are congruences
`y==r_i (mod d_i)`. Generalized CRT says they have a common solution iff
every pair does: `r_i==r_j (mod gcd(d_i,d_j))`. When nonempty their
intersection has `b/lcm_i(d_i)=gcd_i |K_i|` elements.

Consequently, at ANY degree, elimination of this cyclic owner is exactly:

1. retain the projected unary domain on every neighbor;
2. insert the composed relation `R_i^op composed-with R_j` between EVERY
   pair of neighbors, intersecting with any existing relation;
3. multiply the counting weight by `gcd_i |K_i|`.

Every new pair relation is again empty or a coset with the two-congruence
form (1). Intersections preserve the same class. First intersect any
unary restriction on `y` into its incident relations; if it is isolated,
its elimination factor is just its unary-domain cardinality. Empty
relations/domains terminate with zero. Repeating the process solves
compatibility and counts all actual shared cyclic chart choices exactly.
There are at most quadratically many binary relations at any stage;
straight variable elimination makes at most cubically many relation
compositions. These relation operations use exact gcd/CRT or integer
lattice elimination, not enumeration of all chart assignments.

More explicitly, encode the two incident relations as four linear integer
congruences in `(x,y,z)`, with integer slack variables. A fixed-size
Smith/Hermite solve first decides whether this affine lattice is empty;
otherwise project its integer basis and offset onto `(x,z)` and put that
rank-two lattice into the normal form of section 1. Intersection uses the
same construction with shared endpoint variables. This gives a compact
arithmetic implementation of each fill-in operation. The reference
checker deliberately uses explicit finite sets instead, providing an
independent small-universe test of this compact-operation theorem. The
cubic operation count assumes the initial authenticated relations have
already been supplied; it does not bound the cost of discovering cells or
unknown internal groups from an arbitrary tournament.

The proof is the preceding CRT intersection theorem at each step. It also
gives back-substitution: retain the old incident relations and choose an
actual middle witness in their common coset. The product of the recorded
factors counts all such choices. Counting here starts with weight one on
each compatible assignment, or a common scalar; arbitrary nonconstant
edge weights need not survive contraction in this form. This is a special
extra property of the
cyclic coset carrier, not a theorem for arbitrary attachment constraints.
Dropping some fill-in pairs, retaining just satisfiable individual
relations, or replacing noncyclic charts by integer labels destroys the
argument. Finite abelian groups need more coordinates and can have
cosets that are pairwise intersecting but not jointly intersecting.

This is an arithmetic specialization of the incoming
[owner-indexed enumerator rule](../../session41_owned_attachments_fibonacci_20261023/residual-cycles/OWNER_INDEXED_ENUMERATORS.md).
An eliminated actual owner that is later queried must be recovered by its
stored relations, or left on the boundary. The scalar multiplicity is not
its identity, nor does it replace the Smith coordinates of (4).

## 4. Actual unequal/composite tournament controls

For odd `r`, use the regular tournament `C3[R_r,R_r,R_r]` on `3r` owners,
where each `R_r` is the cyclic interval tournament on `r` vertices. Label
owners `t=3*x+a`, with outer coordinate `a=t mod3`. The permutation
`t -> t+1 mod(3r)` is an actual automorphism, giving a declared cyclic
subgroup of order `3r`. Take `r=1,3,5` for orders 3,9,15.

Between any two such cells use the literal block `B(t,u)=1[t==u mod3]`.
Its row/column degrees are the other/own cell order divided by three.
The cyclic comparison stabilizer is `a==b mod3`, with normal form
`(h0,s,v)=(3,1,1)`. The 9-by-9 block is exactly the prior repeated-profile
hostile, now restricted to the authenticated `C9 x C9` subcarrier:
27 compatible pairs and THREE target shifts above each fixed source.
The FULL internal groups are larger; they give 2187 compatible pairs and
TWENTY-SEVEN target charts over each fixed source, as inherited.

For a 3--9--15 path, composing equality-mod-three relations leaves 15
endpoint pairs, each with exactly three middle witnesses. For a triangle
of these three cell sizes, the local orbit space has 27 configurations;
there are three shared-rotation classes, each with nine configurations,
and every nonempty global comparison has 45 witnesses. The Smith group
is `C3`. A single incompatible cycle offset is rejected even though every
individual attachment relation is nonempty. All values are obtained both
from the actual blocks and by the lattice formulas.

## 5. The minimal regular nonabelian boundary and the repair

The Paley seven-tournament has actual automorphisms

    t(x)=x+1,  d(x)=2x  on F7, with d*t*d^(-1)=t^2.

This is the nonabelian group of order 21. It is already present in the
regular-seven carrier, the first odd regular order with a nonabelian
internal automorphism group (the complete regular-seven census and the
three/five classifications are in the prior packet). The matching block
`I7` has diagonal comparison group consisting of all 21 pairs `(g,g)`.
Its composition law cannot be represented faithfully by addition in
`C21`: the displayed commutator is nontrivial. Labeling the 21 elements
with integers does not authenticate an additive chart action. The checker
also shows an actual translated matching relation failing affine-coset
closure under the natural labels `7*k+b` for `x -> 2^k*x+b`.

The full nine-cell group is an additional larger-kernel boundary, not a
cyclic group of order 81. The repaired general carrier is a finite group
comparison relation (a coset of a subgroup of the product), with actual
shared group elements and all one-sided kernels retained; or a recursive
module decomposition with those kernels explicitly carried. For finite
ABELIAN groups, one coordinate per invariant cyclic factor gives the same
integer-cokernel construction in higher dimension. Two coordinates per
attachment are justified specifically by cyclic endpoint domains.

The [non-normal kernel control](NON_NORMAL_KERNEL.md) gives a second
literal nine-cell block: its kernel is normal in the actual projected
chart domain but not in the full group. It also distinguishes source
kernels acting on the right from transported target kernels acting on the
left, and exhibits a nine-lift quotient whose actual owner image still has
three possibilities. This authenticates the restrictions in the root
[kernel-factorization theorem](../KERNELS.md).

Actual marked owners restrict chart sets further and must be encoded as
constraints on the SAME variables. Actual Hamiltonian ports and their
integer positions are transported by an accepted global chart. Native
mod-eight electrical data, converse, and cell names remain separate
authenticated inputs. This packet neither replaces those inputs by a
lattice index nor claims unrestricted source-only tournament recovery.

An extra path or electrical guard need not be a binary affine-coset relation.
In that case keep its queried variables on the boundary or retain its full
joint factor. The cubic relation-operation bound does not absorb arbitrary
additional guards by dropping them or replacing them with pair projections.

## Reproduction

From the repository root:

    python -B 05-knowledge/results/attachment_lattices_square_sum_20261023/lattice/check.py
    python -O -B 05-knowledge/results/attachment_lattices_square_sum_20261023/lattice/check.py

No external package or network is used. `--write` deliberately regenerates
the receipt. Tests enumerate every subgroup normal form for cyclic orders
1 through 9, every coset for each such subgroup, actual 3/9/15 tournament
blocks, a composite Smith example, compatible and incompatible cycle
offsets, and all two-relation compositions for middle/endpoint orders at
most 3. The any-degree contraction checks all 37,975 triples of cyclic
cosets through order 12, 600 seeded full networks with exact direct
enumeration, and a four-port star with a three-choice hidden middle.
The `C3^2` lines `x=0,y=0,x+y=1` are pairwise compatible but have no common
point, providing an odd-order noncyclic boundary. Exceptions, rather than
removable assertions, enforce checks.
