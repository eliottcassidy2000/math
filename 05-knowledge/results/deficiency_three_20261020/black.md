# Endpoint-changing black detours through a terminal R4 ear

Status: **PROVED** for the explicit three-move constructions below;
**FINITE-EXACT** for the stated audits. All-order black connectivity remains
**OPEN**. This note uses the repaired all-order descent and terminal-basin
normalization in [black_descent_20261020c/scc.md](../black_descent_20261020c/scc.md)
and [cycles.md](../black_descent_20261020c/cycles.md). It does not replace them
by an earlier finite descent observation.

## Inheritance and the new target

The nearest proved mechanism is the two-complement common hub from
[defect_transitions_20261020/blue.md](../defect_transitions_20261020/blue.md).
Its sharp hostile is R4 -> Z -> R4^op: when Z is pure blue, fixing both
endpoints of the standard source path can force every replacement path to
remain blue. The corrected near miss is the n=6 shared-internal-cut repair;
it still cannot handle that two-endpoint hostile. The less-used sidecar is
the exact endpoint score charge of a path reversal, from the new SCC note.

The present construction lets endpoints change. It repairs every one of
the five R4 endcap path types, carries the complete labelled middle, and
then weakens the endcap hypothesis to a guarded induced four-vertex ear.
The live concepts are: exact path-edge parity, endpoint charge, terminal
SCC size, positional blue symmetry, and contraction to descent basins.

Write R4 for the strong tournament on a,b,c,d with arcs

    ab, bc, cd, ca, da, db.

This is the same fixed-path code-zero R4 as the earlier hostile. Its
outdegrees are 1,1,2,2. It is rigid: the two score classes and the arc
inside each class fix all four vertices. Its five Hamiltonian paths are

    q0=abcd, q1=bcda, q2=cdab, q3=dabc, q4=dbca.             (1)

Reversing labels a<->d,b<->c gives its converse. In particular R4^op is
isomorphic to R4, with this specified identification, not an implicit
arbitrary relabelling.

## 1. Exact path reversal and the five templates

For an actual directed Hamiltonian path P of T, let R_P(T) reverse just
its n-1 path arcs. This is the converse of the nonpath complement K_P(T),
so represents the same merged-metagraph edge. Its colour is black exactly
when reversal of the positions of P is not an anti-automorphism of T.

Let left vertices 0,1,2,3 induce another labelled R4, let right vertices
a,b,c,d induce R4, and let M be an ordered list of the remaining vertices.
For each q_i in (1), define L_i by reversing q_i and replacing
a,b,c,d by 3,2,1,0 respectively. Thus

    L0=0123, L1=3012, L2=2301, L3=1230, L4=3120.

Put A_i=(L_i,M), J=(0,1,2,3,M), and P=(A_i,q_i). Concatenation of lists is
written without commas below; rev(A) reverses the entire list A.

| i | first path P1 | second path P2 | third path P3 |
|---|---|---|---|
| 0 | A0 b c d a | a d b rev(A0) c | c A0 a b d |
| 1 | A1 a b c d | d b a rev(A1) c | c A1 b d a |
| 2 | J c d a b | a d b c rev(J) | A2 c b d a |
| 3 | A3 a b c d | d c a rev(A3) b | b A3 d a c |
| 4 | A4 b c d a | b rev(A4) a d c | c a A4 d b |

In every row, as sets of unordered edges,

    E(P1) triangle E(P2) triangle E(P3) = E(P).             (2)

Here triangle denotes symmetric difference, not triangle count. Each
internal M path edge occurs three times, and no other internal M edge
occurs. The finite boundary ledger cancels all unintended connecting
edges. Consequently, whenever the successive paths are directed,

    R_P3(R_P2(R_P1(T))) = R_P(T)                            (3)

as literal adjacency matrices on the original vertex labels. No
isomorphism search, converse ambiguity, or unrecorded Hamiltonian path
change enters (3). The checker independently verifies each edge ledger.

**Theorem 1 (all middle tournaments).** Let

    T = R4 -> Z -> R4,

where Z is any tournament, possibly empty, and M is any Hamiltonian path
of Z. Every row of the table is a sequence of three black moves satisfying
(3). The original path P is blue if and only if its middle presentation
M is blue. Hence every blue line at every R4 -> Z -> R4^op endcap source
has an explicit three-black replacement, even when both-endpoint-preserving
replacement fails.

Proof of directedness: internal M path arcs alternate in direction, so
its successive forward/reverse/forward traversals are valid. The labelled
R4 edge lists verify the remaining arcs; initial cross-arcs all point
left to right. Equation (2) proves the final labelled equality.

For blackness, the first move changes the right R4 path while retaining
the left one, except in row 2 where it changes the left while retaining
the right. Both cannot remain positional reversals: that would give a
nonidentity automorphism of the rigid opposite R4. For rows 0,1,3,4 the
second and third moves have both endpoints in the right endcap. Their
score sums in their current tournaments are respectively

    row 0: 4,5;  row 1: 5,4;  row 3: 4,4;  row 4: 3,3.

These are strictly below n-1, since n>=8. A blue path requires every
paired score sum, including its endpoint sum, to equal n-1. Row 2 is
covered by the stronger structural argument in Theorem 2 below. This
proof uses no property of the internal arcs of Z beyond the chosen path.

For the original standard hostile, row 0 reads especially simply. With
A=(0,1,2,3,M),

    (A,b,c,d,a), (a,d,b,rev A,c), (c,A,a,b,d).               (4)

The first and last global endpoints visibly change. For the more general
source A->R4 with |A|=m>=3, (4) still gives the exact three reversals and
all are black whenever d_A(first(A)) != m-2. The first endpoint sum is
d_A(first(A))+5; the next two are 4 and 5, compared with n-1=m+3. This is
a sufficient guard, not a necessary condition for a black detour.

## 2. A larger family: one guarded terminal R4 ear

**Theorem 2.** Let T be nonstrong and let P be a blue Hamiltonian path.
Suppose the last four vertices of P induce R4. Label them a,b,c,d so that
their order is q_i from (1), and write w for the vertex immediately
preceding this four-vertex suffix. Require the following one extra arc:

| i | extra boundary arc |
|---|---|
| 0 | w -> b |
| 1 | w -> a |
| 2 | w -> a |
| 3 | w -> a |
| 4 | w -> b |

Then the corresponding three-black template gives exactly R_P(T).

Proof. Since the terminal four-vertex subtournament is strong, it lies
inside the terminal SCC S of T. Blue reversal identifies it with the
first four vertices inside the initial SCC, with precisely the left
labels used above. Thus n>=8. The entire terminal SCC has size s<=n/2,
because its reversed mate is a distinct initial SCC of the same size.

For rows other than 2, let A be the entire prefix preceding the last four
vertices. The only new prefix-to-suffix arc needed by P1 is the extra arc
in the table. The arc from first(A) to the appropriate right vertex is
present because distinct SCCs are ordered. P2 and P3 use reversals of
these same two boundary arcs. All other boundary arcs are original path
arcs. In row 2 the first four vertices are changed instead: reflection
of w->a supplies the needed new exit from the initial four-vertex block.
The same edge ledger proves directedness and (3), without requiring any
uniform relation between the rest of A and the terminal four vertices.

P1 is black by the rigid-R4 argument already given. For distinct u,v in
S, their original scores satisfy

    d_T(u)+d_T(v) <= 2s-3,                                 (5)

because S has no outgoing arcs to earlier SCCs, and one of u,v loses
their mutual arc. Reversing a directed spanning path subtracts 1 from
the score of its first vertex and adds 1 to its last vertex. In terms
of original scores, the endpoint sums for P2 and P3 are

| row | P2 endpoint sum | P3 endpoint sum |
|---|---|---|
| 0 | d(a)+d(c)+1 | d(c)+d(d)+1 |
| 1 | d(d)+d(c)+1 | d(c)+d(a)+1 |
| 3 | d(d)+d(b)+1 | d(b)+d(c)+1 |
| 4 | d(b)+d(c) | d(b)+d(c) |

By (5), all are at most 2s-2<=n-2. These paths are black.

For row 2 a score inequality is unnecessary. Before P2, the left four
vertices induce the transitive tournament: P1 reversed their standard
R4 path. The right four still induce a strong tournament, with the
cycle a->d->b->c->a. They are respectively the last and first four
positions of P2, so no positional anti-automorphism is possible.
Before P3, the left R4 has been restored; the right four have a vertex
a of internal outdegree zero. They occupy the first and last four
positions of P3. Again a strong and a nonstrong induced tournament
cannot be paired by an anti-automorphism. Thus both moves are black.

Theorem 1 is recovered because a terminal R4 endcap makes every listed
extra boundary arc automatic; the initial and terminal four-vertex
patterns already agree even if the middle is not blue.

## 3. Sharp limitation and connection to terminal basins

The extra boundary arc is substantive. Let F5=T*_5 be the fixed-path
code-zero tournament with rows (2,4,9,19,7), and take

    T = F5^op -> F5,

using the converse-reversed labels so that the standard Hamiltonian path
is blue. Its last four vertices form the q0 presentation of R4, but the
preceding vertex beats a and loses to b. The first proposed move requires
the absent arc w->b. Exact rows of this ten-vertex hostile are

    (994,996,1001,1011,999,64,128,288,608,224).

The invalid proposed P1 is (0,1,2,3,4,5,7,8,9,6). Its missing arc is 5->7.
The first failed implication is "an induced terminal R4 suffices";
the strongest survivor retains the explicit one-arc boundary guard.
This is a failure of this template, not a disconnected pair.

Under the proved all-order descent, every strong black-active node lies
in a basin attached to a nonstrong terminal. Theorems 1 and 2 say that
the strong neighbour reached by any of the covered blue terminal edges
already lies in the SAME black component as that terminal. Equivalently,
these blue incidences collapse to loops after black-component contraction.
They cannot be a missing link between two distinct black components.

This removes a genuine infinite obstruction family from the remaining
blue quotient, but it neither connects all terminal basins nor covers
every blue path whose source has a larger initial/terminal strong SCC.
The boundary-failing F5 test is now resolved in the
[enlarged-ear supplement](ear_probe.md), which proves a uniform detour for
every standard full-tiling endcap of order at least four. The remaining
question is whether arbitrary terminal SCC basins admit a suitable
endpoint-changing ear and blackness certificate. F5 is not the regular
five-vertex tournament, denoted R5 elsewhere in the repo.

## 4. Exact controls and scope

Reproduce with

    python 05-knowledge/results/deficiency_three_20261020/black_check.py
    python -O 05-knowledge/results/deficiency_three_20261020/black_check.py

The checker imports only the existing elementary tournament/path helpers.
It verifies path validity, every positional anti-automorphism predicate,
the independent unordered-edge XOR, and literal final matrix equality.
No canonicalization or removed Python assertion is used.

* Endcaps: all fixed-path middle presentations of orders 0 through 6.
  They have 1,1,1,2,8,64,1024 presentations and 5,505 three-black
  certificates across all five boundary path types.
* Guarded ears: every strong fixed-path presentation of orders 4,5,6,
  as the terminal SCC, paired with its converse and a middle of order
  0,1,2. Strong presentation counts are 5,50,903. Tail-R4 counts are
  5,35,600; guarded counts are 5,15,300; all 960 applicable three-move
  certificates pass. These are presentations, not isomorphism counts.
* The complete inherited nonstrong-blue presentation generator gives
  70 sources at n=8 and 519 at n=9; the five R4-endcap source presentations
  at each order are all covered. Other source types are outside this
  specific ear theorem, not failed certificates.
* Hostile control: dropping the one-arc boundary guard produces the
  explicit invalid F5^op->F5 move above.

The [receipt](black_check.out) is finite evidence for these exact programs.
The all-order conclusion comes from the symbolic path ledger and SCC
score/rigidity arguments, with the stated quantifiers and boundary guard.
