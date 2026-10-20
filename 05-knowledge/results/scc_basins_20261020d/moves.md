# Releasing endpoints at four- and five-vertex SCC boundaries

Status: PROVED for the bounded black detour theorem below; FINITE-EXACT
controls are labelled separately. All-order black-core connectivity is OPEN.
Scope: the converse-merged complement-line metagraph, not single-arc flips.

## 1. Inheritance and target

The [all-order descent theorem](../black_descent_20261020c/cycles.md) puts
every black component in contact with a nonstrong tournament. The
[terminal-basin reduction](../black_descent_20261020c/scc.md) therefore
leaves the problem of joining those SCC boundary states. Its four-step
interchange has a two-singleton prefix guard.

The closest complementary mechanism is the exact common target in
[the nonstrong rotation lemma](../defect_transitions_20261020/blue.md),
Section 2. Its fixed-boundary implementation fails on

    R4 -> Z -> R4^op,

with Z pure blue: the standard first and last vertices force every path
between them to be blue. That is the canonical hostile, proved in Section
5 of that note. The least-used coordinate is the choice of endpoint inside
the terminal R4. The present construction releases it explicitly.

Concept board: triangle basins; ordered SCCs; a retained middle path;
endpoint degrees; path-edge XOR; and black detours around blue edges.
The new result joins basins across a previously obstructed family of
blue lines. It does not infer connectivity from degree inventories.

Write K_P(T) for complementing every arc outside a directed Hamiltonian
path P. Thus

    K_P(T) = T^op triangle E(P),
    d_(K_P(T)) = (n-1)1-d_T+e_first(P)-e_last(P).       (1)

For an odd number q of successive path presentations P_1,...,P_q,

    K_(P_q)...K_(P_1)(T)
       = T^op triangle (E(P_1) xor ... xor E(P_q)).    (2)

Every line is unchanged by replacing this operation with full path
reversal and then taking converse. A presentation is blue precisely when
its positional reversal is an anti-automorphism.

## 2. Uniform detours for every four-vertex outer SCC

**Theorem (PROVED).** Let U be a nonstrong self-converse tournament whose
first and last strong components have order four. Every blue line
incident to U has a black-only detour of length at most five between its
endpoints. In particular every such blue line is nonbridging, and its
endpoints belong to the same black component.

There is a unique strong four-vertex tournament R4. In path coordinates
one representative has outgoing-neighbor bitmasks

    R4=(2,4,9,3).

Thus U has the form R4^op -> Z -> R4, where Z may be empty. We do not
require Z to be pure blue. Fix the blue path defining the line under
consideration. Its restrictions to the two outer SCCs are paired by
positional reversal, and its middle restriction is blue.

The five Hamiltonian paths in this representative of R4 are

    0123, 1230, 2301, 3012, 3120.

Reading R4 along these paths gives the following five right-hand path
presentations B. The left R4 presentation is its reversed converse.
Put A=(the left R4 path, the middle Z path). In every word below A is
kept contiguous and in the same forward order. Write 0,1,2,3 for the
four right-hand vertices; these labels are disjoint from those in A.
The blue path whose complement is the target is

    P_0=(A,0,1,2,3).

Apply K successively along the words in the indicated row:

| B outgoing bitmasks | successive black path words |
|---|---|
| (2,4,9,3) | (A,1,2,3,0); (2,A,1,3,0); (2,A,0,1,3) |
| (2,12,9,1) | (A,2,0,1,3); (1,3,A,2,0); (1,2,3,A,0) |
| (6,12,8,1) | (A,1,2,3,0); (3,A,1,2,0); (0,3,A,1,2); (0,2,3,A,1); (2,3,A,0,1) |
| (6,4,8,3) | (A,1,2,3,0); (2,A,1,3,0); (2,A,0,1,3) |
| (10,4,9,2) | (A,1,2,0,3); (2,0,3,A,1); (2,3,A,0,1) |

### Path validity and exact target

These are finite boundary identities that extend to an arbitrary A path.
Its consecutive internal path arcs are retained by every K operation,
so the forward A word is always available. Its other internal arcs are
complemented at every step; the number of steps is odd, exactly as in
K_(P_0). Every arc from an interior A vertex to B that is never a path
boundary likewise has the correct odd complementation parity.

The only remaining information is the first and last A vertices and
the four vertices of B. Replace A temporarily by a two-vertex path
(L,R), initially dominating B. For each of the five rows, direct
six-vertex substitution checks every displayed path arc and

    E(P_1) xor ... xor E(P_q) = E(P_0).                (3)

The checker below performs these five finite identities with no
isomorphism search. Since every visit to A has the same direction and
the only outside attachments are at L or R, the same boundary identities
hold for A of any length at least two. Here |A|>=4. Equations (2) and
(3) give the literal labelled target K_(P_0)(U).

### Blackness

First consider P_1. It changes only the path through the final R4,
leaving the initial R4 and Z words fixed. If both P_0 and P_1 were blue,
the composition of their positional reversals would induce a nonidentity
automorphism on the final R4. The restriction is nonidentity because
the two terminal words differ. Every automorphism preserves each SCC,
and R4 has trivial automorphism group. For completeness, its score
classes have sizes two and two, so an automorphism has order dividing two;
a tournament has no nontrivial involutory automorphism, since it would
reverse the arc between a swapped pair. This contradiction proves that
P_1 is black in all five rows.

For the subsequent steps it is enough to test endpoint degree sums.
A blue presentation would require that sum to be n-1. Formula (1)
gives the following differences from n-1, in chronological order after
the first step:

| B | endpoint-sum differences for P_2,...,P_q |
|---|---|
| (2,4,9,3) | n-5, 6-n |
| (2,12,9,1) | n-4, 4-n |
| (6,12,8,1) | n-5, 4-n, n-4, 4-n |
| (6,4,8,3) | n-5, 5-n |
| (10,4,9,2) | n-4, 4-n |

All are nonzero because n>=8. These computations use only the four
terminal scores and the initial score of the first A vertex. The latter
is n-1-d_B(3), by the original blue reversal; no hidden hypothesis on Z
or its internal degrees is used.

Thus every displayed move is black, and the endpoint is the original
blue neighbor. Removing repeated vertices from the walk gives a black
path of length at most three or five. This proves the theorem.

## 3. Regular five-vertex outer SCCs need only three black moves

**Theorem (PROVED).** Let the first and last strong components of a
nonstrong self-converse U both be the regular five-vertex tournament R5.
Every blue line incident to U has a black-only detour of length at most
three. The middle Z is arbitrary; it need not be pure blue.

The mechanism is the same, with A=(the first R5 path, the Z path).
R5 on residues modulo five with positive increments {1,2} has 15
Hamiltonian paths in three automorphism orbits. Representatives are
01234, 02341, and 02413. Reading along these words gives the three
right-hand presentations in the table. The target blue word is
(A,0,1,2,3,4).

| B outgoing bitmasks | successive black path words |
|---|---|
| (6,12,24,17,3) | (A,3,0,1,2,4); (1,2,4,A,3,0); (1,2,3,4,A,0) |
| (18,12,9,17,6) | (A,0,1,3,4,2); (4,A,0,3,2,1); (2,4,A,0,3,1) |
| (10,20,9,18,5) | (A,0,1,4,2,3); (4,3,A,0,2,1); (3,A,0,2,4,1) |

The completeness of these three path presentations is a finite
five-vertex fact: direct enumeration gives precisely these three
induced matrices, each five times. The accompanying checker verifies
this over all 120 permutations. Equality of two induced matrices means
exactly that the map between their ordered words is an automorphism.

The seven-vertex boundary ledger, using a two-vertex A proxy, proves
path validity and the odd edge-XOR identity exactly as in Section 2.
The first move changes the right-hand Hamiltonian path to a different
automorphism orbit: in each row, reading B along that new terminal word
gives a matrix different from B. If this move and the original path
were both blue, their reversals would induce an automorphism between
those two paths in the last SCC, which is impossible. Thus the first
move is black even though its endpoint degree sum is balanced.

Every right-hand vertex initially has score 2, and the first A vertex
has score n-3. Formula (1) gives endpoint-sum differences from n-1
equal to n-5 and 5-n at the second and third moves. Both are nonzero
since n>=10. Hence those moves are black too, and their labelled target
is precisely the original blue partner.

## 4. What this adds to SCC normalization

The theorem repairs the entire infinite hostile family
R4 -> Z -> R4^op, including every blue path presentation, not merely
the standard one. Its moves change global endpoints while preserving
the chosen middle path as a contiguous block. A fixed-endpoint search
cannot discover this route.

The middle may be mixed. Consequently the theorem also handles examples
outside the incoming mixed-center detour theorem's padding argument:
here blackness can come solely from choosing incompatible outer R4 paths,
even when the middle has no black path at all.

The theorems identify the two endpoints of every such blue line in the
black-component quotient. It does not prove that every SCC basin can be
transformed to one with an outer R4 or R5. Source/sink components of other
types, and strong mixed blue sources, still need additional moves.

## 5. Reproduction and exact scope

Run from any directory:

    python 05-knowledge/results/scc_basins_20261020d/moves.py
    python -O 05-knowledge/results/scc_basins_20261020d/moves.py

The self-contained checker verifies five six-vertex R4 boundary
identities and three seven-vertex R5 identities, the complete path
presentation lists, R4's trivial automorphism group, and the affine
endpoint-degree expressions. It then tests every fixed-path
blue middle presentation through order seven, with all five outer path
presentations for R4 and all three for R5, checking each actual path,
its blackness, and literal
equality to the original blue partner. The resulting counts and receipts
are in [moves.json](moves.json). These finite tests audit the symbolic
argument; the unbounded conclusion follows from the boundary identities
and the degree formulas above.

There are 601 blue middle presentations in the stated range. They give
3,005 complete R4 detours and 1,803 complete R5 detours, comprising 15,626
individually checked black steps. No isomorphism search is used to verify
their final targets.
