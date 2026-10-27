# SCC minima, terminal basins and a uniform black interchange

Status: PROVED for the statements explicitly labelled below.
The basin theorem uses the new all-order descent theorem in [cycles.md](cycles.md).
Complete black connectivity and unrestricted SCC normalization remain OPEN.

## 1. Inheritance and the coordinate that a local rewrite loses

The closest proved mechanisms are
[THM-4729, closed/open metagraph lines](../../../01-canon/theorems/THM-4729-closed-and-open-tilings-hamiltonian-cycle-bipartition-of-the-metagraph-lines-and-the-forced-arrangement.md),
the [endpoint transport formula](../merged_metagraph_endpoint_transport_20261020/endpoint-transport-and-complement-leakage.md),
and [ordered SCC converse](../tournament_converse_context_20261011d.md).
The hostile is a proposed interior SCC interchange whose vertices never
become global Hamiltonian-path endpoints. The corrected target retains
the actual paths, their blackness and the unchanged exterior.
The least-used sidecar is the labelled endpoint-charge vector.

Our board has five objects: triangle potential; SCC words; endpoint
charges; actual black path sequences; and the quotient by descending
trees. The first gives termination, but the fourth is needed to connect
different terminal words. A terminating descent alone cannot prove
irreducibility.

Write T^P for reversal of every arc of the directed Hamiltonian path P.
This is the same converse-merged complement line as the inherited
nonpath-complement operation. A path is blue exactly when positional
reversal is an anti-automorphism; otherwise it is black. Write s(v) for
outdegree and c3(T) for the number of directed triangles.

For P=(a,...,b), the exact labelled score change is

    s(T^P)-s(T) = e_b-e_a.                                  (1)

Consequently, in any sequence of such moves, the total score change is
the incidence boundary of the directed multigraph of global endpoint
pairs a->b. A vertex never used as an endpoint cannot change score.
For example, replacing an interior C3->1 by 1->C3 changes its internal
score multiset from (2,2,2,0) to (3,1,1,1). A fixed exterior adds the same
constant to every score, and cannot repair this discrepancy. Thus a
contextual implementation must transport endpoint charges into that block.
This is an obstruction to the proposed local implementation, not to the
existence of a longer global black route.

## 2. All local minima of the full triangle potential

**Theorem (PROVED, all n>=3).** In the full converse-merged metagraph, the
strict local minima of c3 are exactly the nonstrong tournaments.
Every line incident to a nonstrong tournament strictly increases c3.

Let the first and last SCCs of a nonstrong T have sizes s and t. A
Hamiltonian path begins at a in the first and ends at b in the last.
When s>1, the first SCC is strong and every vertex has an internal
outneighbor; when t>1, every last-SCC vertex has an internal inneighbor.
Therefore

    s(a) >= n-s + 1_{s>1},
    s(b) <= t-1 - 1_{t>1},
    c3(T^P)-c3(T)
        = s(a)-s(b)-1
        >= n-s-t+1_{s>1}+1_{t>1} > 0.                       (2)

The final inequality holds because s+t<=n; if both are singletons it
is n-2>0, and otherwise an indicator is positive.

Conversely, a strong tournament has a directed Hamiltonian cycle.
The endpoint score differences of all its cuts sum to zero, so one cut
has s(last)>=s(first) and hence triangle change at most -1. This argument
allows blue moves; [cycles.md](cycles.md) supplies the additional
blackness at every order when a black line exists.

The theorem is consistent with the inherited stronger incidence fact:
the nonstrong periphery is independent and every incident complement
line goes to a strong node. See
[the exact strong-core cut counts](../merged_metagraph_purity_and_strong_core_20261020/STRONG_CORE_COUNT_RECURSIONS.md).
Strict increase is the new potential formulation needed for descent.

**Corollary (PROVED, all n>=3).** On the black-active subgraph, the
strict triangle minima are exactly its nonstrong nodes. Every black
component contains a nonstrong node. Starting from any strong
black-active class, repeatedly choose the descending move from
[cycles.md](cycles.md); the integer potential decreases at every step,
so the route must terminate at a nonstrong class.

In particular, the terminal set is much larger than the rigid words
over {1,C3}. Treating just those rigid leaves as all terminal states
would be false: every nonstrong black-active class is a strict minimum.

## 3. The terminal-basin quotient

For n>=3, let M be the set of nonstrong black-active merged classes.
For every strong black-active vertex choose one incident line to a
strictly smaller c3 value. These chosen lines form a forest oriented
towards M: there are no directed cycles, and each connected tree has
exactly one root in M. Let r(v) be the root reached from v.

Define the **terminal-basin multigraph** by retaining one vertex for
each member of M and sending every black line uv to r(u)r(v), preserving
multiplicity and allowing loops.

**Theorem (PROVED).** This quotient preserves black connected
components. In particular,

    black-active graph connected
      iff terminal-basin multigraph connected.             (3)

The proof contracts the chosen forest. Equivalently, a quotient route
lifts by joining each endpoint to its root along its tree, and every
original route projects. The quotient depends on the descending choices;
its component partition does not. This is an exact reduction of the
remaining obligation to SCC boundary states, not an irreducibility proof.

There is also a useful linear form. Give black lines positive
conductances and let L be the graph Laplacian, ignoring loops. Partition
vertices into strong active vertices S and the terminals M. Every
component meets M, so L_SS is positive definite: if x vanishes on M,
the quadratic form is a sum of positive multiples of (x_u-x_v)^2 and
can vanish only for x=0. Thus the Schur complement

    K = L_MM - L_MS (L_SS)^(-1) L_SM                         (4)

is defined. It is a positive semidefinite Laplacian on M and preserves
the component partition. Indeed for given boundary data y its energy
is the minimum of the original energy over all extensions; zero energy
means constancy on each original component. Hence

    rank K = |M| - number of black components.              (5)

When S is empty, take K=L_MM. Empty M means there are no active black
vertices, by the corollary. Formula (4) is a standard Laplacian
elimination; the new descent theorem is what guarantees its grounding
here at every order. Positivity does not imply full rank |M|-1.

## 4. Exact all-order size of the terminal set

These counts are **isomorphism types**, not the fixed-path tiling counts
with similar letters in the strong-core cut note.
Let a_n,b_n,s_n,r_n count, respectively, all ordinary tournament types,
self-converse types, strong types, and strong self-converse types.
Let p_n count pure-blue merged types. They are self-converse.
For every n>=3 the number of nonstrong black-active merged classes is

    m_n = (a_n-s_n+b_n-r_n)/2 - p_n + 1_{n=3 or n=5}.        (6)

Converse pairs the nonstrong ordinary types, giving the first term by
Burnside. Subtract the nonstrong pure-blue types. The only strong
pure-blue types are C3 and R5, which explains the correction.
This is a counting identity at all orders; its interpretation as all
black descent terminals has been proved here at every order.

The inherited unique SCC decomposition gives the formal identities

    A(z)=sum a_n z^n = 1/(1-S(z)),
    B(z)=sum b_n z^n = A(z^2)(1+R(z)),

where a_0=b_0=1, S(z)=sum_{n>=1}s_n z^n and
R(z)=sum_{n>=1}r_n z^n. The second is the converse-palindrome formula in
[the SCC context note, Section 3](../tournament_converse_context_20261011d.md).
The pure-blue classification gives

    sum_{n>=0} p_n z^n = (1+z+z^3+z^5)/(1-z^2-z^6).

Thus (6) is an exact coefficient prescription, with the stated
small-order convention. It does not turn the ordinary tournament type
counts into a simple closed form.

## 5. A uniform four-black interchange with arbitrary prefix data

**Theorem (PROVED).** Let A=1->1->B, where B is any tournament, possibly
empty. There is a four-black-move route

    A -> 1 -> C3    to    A -> C3 -> 1,                     (7)

which restores every arc of A literally. All arcs between A and the
four remaining vertices are restored as well. The orientation of the
triangle uses a different labelled triple at the end.

Write m=|A|>=2 and choose a Hamiltonian path
A=(u_0,u_1,...,u_(m-1)). Its first two vertices are the two dominating
singletons. Put x,a,b,c after A, initially with x beating the other
three and triangle a->b->c->a. The successive paths are

    P1 = ( A,     x,c,a,b ),
    P2 = ( rev A, c,x,b,a ),
    P3 = ( x,a,b,c, A     ),
    P4 = ( c,b,a,x, rev A ).                               (8)

Every path arc is present at its step. One can verify this by the
following small ledger: internal edges of the chosen A path reverse at
every step; nonpath edges of A never reverse. Edges joining A to the
four vertices occur in cancelling pairs:
u_(m-1)x in P1,P4; u_0c in P2,P3.
Among the last four vertices, the edges ac and xb reverse once;
xa,ab,xc,bc reverse an even number of times or never. The final triangle
is x->a->b->x, all three beat c, and the displayed target follows.
In detail the four internal words are respectively
x,c,a,b; c,x,b,a; x,a,b,c; c,b,a,x, so their validity is checked in the
six-edge four-vertex tournament independently of B.

Blackness requires a separate check. A blue positional reversal would
pair all scores to n-1=m+3.

* For P1 the endpoint scores are m+3 and 1, whose sum is m+4.
* For P2 the endpoint sum is d_A(u_(m-1))+5.
  For m>=3 the last vertex has d_A<=m-3, so this is below m+3.
  For m=2 the second and penultimate scores are 4 and 2, sum 6!=5.
* For P3 the endpoint sum is d_A(u_(m-1))+6.
  For m>=4, the last vertex lies in B and has an inneighbor on its
  Hamiltonian path, hence d_A<=m-4; the sum is below m+3.
  For m=2 it is 6!=5. For m=3 the second and penultimate scores are
  2 and 5, sum 7!=6.
* For P4 the second vertex b has score 2 and its reversed partner u_1
  has score m+2, sum m+4!=m+3.

All four paths are therefore black. Taking converses gives the dual
interchange with two trailing dominated singletons.

The guard is substantive. If A consists of just one singleton, the
same four paths are valid but P4 is blue. The target 1->C3->1 is pure
blue. The theorem therefore cannot be silently promoted to arbitrary
contexts or used to cross the pure-blue boundary.

This is a local rule on **SCC words** realized by global paths. It is
not a local arc operation. Formula (8) explicitly transports the charges
needed in Section 1 and retains arbitrary prefix data B.

## 6. Finite audits and remaining normalization problem

Reproduce the independent macro and inherited-census audits with

    python 05-knowledge/results/black_descent_20261020c/scc_normalization.py
    python -O 05-knowledge/results/black_descent_20261020c/scc_normalization.py

The [receipt](scc_normalization.json) covers every labelled B of order
0 through 5 and every Hamiltonian path of B: 7,888 actual four-move
certificates. It checks all four positional anti-automorphism predicates,
all path arcs, and literal equality of the final adjacency matrix.
It also verifies the one-singleton hostile.

The complete inherited class/line tables give:

| n | merged nonstrong types | black-active types | terminal types m_n | terminal quotient component sizes |
|---|---:|---:|---:|---|
| 3 | 1 | 0 | 0 | empty |
| 4 | 2 | 2 | 1 | 1 |
| 5 | 4 | 7 | 2 | 2 |
| 6 | 12 | 32 | 10 | 10 |
| 7 | 56 | 268 | 52 | 52 |

These validate the count formula, the full strict-minimum classification,
and one explicit descent forest at each order. These are FINITE-EXACT audits of the separate all-order proof.

The earlier [rigid-word distance probe](scc_probe.py) and
[four-step discovery probe](scc_splice_probe.py) retain their actual
paths and receipts. They are provenance, not substitutes for the proof.

What remains is a complete set of black normalization moves connecting
arbitrary members of M, or a pair in distinct components. Rule (7)
handles a genuine infinite family of changes and arbitrary B, but does
not eliminate arbitrary strong SCC letters or treat every prefix/suffix.
The all-order descent theorem and the count formula make this outstanding
obligation precise.
