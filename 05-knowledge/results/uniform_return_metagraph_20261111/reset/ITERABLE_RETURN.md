# An indefinitely repeatable ready-return chart

**PROVED actual graph operation, all-k rooted response, and exact iteration.
OPEN:** preserving the ternary target during the return; Chan–Pak.

Inheritance: the [older nonzero preparation packet]
(../../nonzero_preparation_20261103/README.md) already gives actual
1429→4289/4291 repairs and a two-round nonzero tree, while
[multiport defect grammar]
(../../multiport_defect_grammar_20261103/README.md) gives an indefinitely
reusable strong tripling channel and unbounded port transport. The result
here is a different, fully explicit dyadic ready chart with one-old-vertex
loss and the exact return4H−33. It is not the first reusable operation
in the repository, and it does not supply unbounded nonzero-digit closure.

## 1. Recovering a dyadic ready family from the 425 continuation

Let R_k be the1-ready source in [READY_REBASE](READY_REBASE.md). Perform
its standard two successive births, delete the old port a, and add a
vertex w beaten only by the retained old apex z. Call the result Q_k.
For k>=1 it is the following fixed ten-vertex template with vertex0
replaced by TT_k:

    Q_1=(54,20,48,55,32,578,159,319,127,479).

The ordered labels of Q_k are

    P_0,...,P_(k−1), B,A,r,y,z,v,u,v2,w.

Here P is the transitive module. The symbol v2 is one vertex name.
This template description and the preceding literal surgery agree on
every comparison: they agree on the fixed k1 core and treat the whole
prefix uniformly. In particular Q_2 has H851 and is a source obtainable
from the actual425→799→1277 continuation by one deletion and one birth.

Writing D=2^k, the exact marked-root response at w is

    H=210D+11,       E=3,              S=42D+14,
    U=966·3^k+720D−45,                 V=210D+12,
    d=42D+15,       K=966·3^k+768D−30.                (1)

The notation(H,E,S,U,V,d,K) is the seven-coordinate interface in the
incoming [root amalgamation theorem]
(../../h13_paired_descent_20261111/modular_cycles/ADDITIVITY.md).
Thus V−H=1: every Q_k is1-ready. Its total count is dyadic even though
its full continuation response also has a3^k channel.

## 2. A physical return that can be repeated without limit

From Q_k, perform these four operations:

1. Add x beaten only by w.
2. Add t beaten only by x. This attains the exact ternary target3H+2.
3. Delete old vertex v (label k+5), preserving all other arcs.
4. Add w' beaten only by old z.

The final rooted graph is EXACTLY Q_(k+2), up to the following explicit
block relabelling. The old root w remains as a vertex, but is not the
new marked root.

|Retained/new vertex|Label in Q_(k+2)|
|---|---|
|old P_i|P_(i+2)|
|old v2, old r|P_0,P_1 respectively|
|old B,A,y,z|B,A,y,z respectively|
|old u|r|
|old w|v|
|new x,t,w'|u,v2,w respectively|

Every comparison in this table follows from the template rows and the
three specified birth cuts. It proves the graph identity for all k;
the equality is not inferred from equal Hamilton counts. Exactly one
old vertex is lost, all other old arcs persist, and the order rises by2.

The complete count trace is

    210D+11 → 378D+9 → 630D+35 → 168D+15 → 840D+11.

If h=H(Q_k), its last three entries are

    3h+2 → (4h+31)/5 → 4h−33.                       (2)

The two parity families k odd and k even are therefore two infinite
orbits of the same physical return program. Starting at Q_k and
repeating j times gives Q_(k+2j), with

    H_j=210·2^k·4^j+11.

Its complete sidecars follow from(1) by replacing k with k+2j; their
only exponential modes are4^j and9^j, with constants. For example
the rooted deletion counts satisfy d_(j+1)=4d_j−45. All vertices and
root changes in this iteration are explicitly known.

Every ternary excursion in this family is composite:

    3H(Q_k)+2=35(18·2^k+1).

Thus this particular closed return mechanism does not supply the prime
or seven-times-prime part of the outstanding spectrum frontier.

This supplies an indefinitely composable ready-return family, but the
return does not preserve the ternary target. Its final count exceeds
that intermediate target by h−35. In particular it does not prove that
an arbitrary ternary address can be followed while maintaining an exact
prescribed target. The lost coordinate is the absolute count, despite
recovery of readiness and a known root state.

## 3. Proof of readiness and the full response

Put P=Q_k−w. Starting at z, its unique Hamilton path is

    z,v,u,v2,r,P_0,...,P_(k−1),B,A,y.

Vertex y can only go to z, already used, so y must be last. Choosing B
immediately after z would force B,A,y too early, hence z must choose v.
Only v can enter u and only u can enter v2, forcing the next two steps.
The remaining vertices form the displayed transitive chain. The last
edge y→z closes a spanning directed cycle, and inserting w between z
and v gives one in Q_k. Both graphs are strong.

The exact positive module-run coefficients are

    H(P): (99,42,0,...),       H(Q): (431,210,0,...).

Hence H(P)=42D+15 and H(Q)=210D+11. The indegree-one birth identity
then gives d=H(P), S=d−1, V=H(Q)+1. Counting paths in P ending at z
has coefficients(3,0,...), giving E=3.

For U and K, insert a temporary separator into the path-counting
representation: it has an arc in each direction with every vertex of P.
A spanning path through this separator is exactly an ordered two-path
cover of P, allowing either bin empty. This auxiliary directed graph is
a representation of covers, not a claimed tournament. Its complete
module-run coefficients are

    K: (4404,3666,966,0,...).

Forcing the whole spanning word to start at z by removing incoming arcs
to z gives coefficients(111,48,0,...). Thus the weighted root-start
cover sum at z in P is48D+15. Removing newborn w from a cover counted by
U(Q,w) gives every two-cover of P except those whose first nonempty
path starts at z. Consequently

    U(Q,w)=K(P)−U(P,z).

Substituting j!S(k,j) for each run coefficient proves both remaining
formulas in(1). This supplies exact all-k sidecars, not sampled fits.

As a small sequence consequence, every H(Q_k) is3 modulo4, since
H−1=10(21·2^k+1) has2-adic valuation1. No member is one more than a
square. Thus this family is also disjoint from A002496, independent
of any assumption about which of its values are prime.

## 4. All ready roots, including the converse

The only1-/2-ready root of Q_k is its displayed1-ready w. Neither
orientation has a2-ready root, for every k>=1. In the forward graph the
only vertices of indegree at most2 are u,v2,w. For u and v2 the two-port
residual still contains the cycle z→B→A→z. The rooted defect bound
delta>=2H(residual) therefore gives gap>=6. Vertex w has gap1.

In the converse the only degree-eligible vertices are B,A,y. Deleting
the root and its two in-neighbours leaves z→v→u→z for B and
v→u→v2→v for A, again giving gap>=6. The converse root y has indegree1
and gap equal to the number of paths ending at z in Q_k−y. Such paths
must have w first and B,A,z last. Between them lie the directed triangle
v,u,v2, then r and the transitive prefix. Exactly three paths occur,
so this converse gap is3. All other roots have indegree>=3 and hence
gap>=6 by the same inherited defect bound. These arguments cover every
root, without an extrapolation from the finite controls.

Reproduce with `python -B iterable_return.py`. It independently counts
literal graphs for k1..5, verifies every retained arc and the actual
block isomorphism, proves the fixed response by complete marked run
tables, and checks ten consecutive physical iterations at each starting
k without fitting a recurrence.
