# Leading digits become response determinants

**PROVED IN NOTES:** an exact marked-response difference under a neutral
cyclic reversal, its arbitrary-module specialization, and an all-depth
boundary for old ports. **FINITE-EXACT:** the declared controls in
chirality.json. No all-target spectrum statement follows.

## Inheritance and the actual operation

The [Session312 chart](../session312_cluster_pythagorean_symmetry_20261113/CLUSTER_CHART.md)
retains three different kinds of information: ancestry valuations, signed
leading units, and actual leaf modules. Its unmarked count does not see
all of this information. The underlying neutrality is inherited from
[THM-4803 / labelled cyclic triples](../../../01-canon/theorems/THM-4803-intransitivity-is-the-relation-modulo-its-modular-tree-and-a000568-from-the-prime-cycle-index-census.md)
and the [THM-3121 substitution kernel](../../../01-canon/theorems/THM-3121-path-cover-walk-content-substitution-kernel.md).
[Session188 covered-quotient §2](../session188_structural_bridges_20261104/covered_quotient/PROOF.md)
already gives the adjacent-minor identity for aggregate starts in a block.
Equation(1) below is its named-root refinement, with L_A replacing the
unconditioned cover sequence. The new four-coordinate specialization and
readiness boundary make the signed orientation an exact selector
for a marked response, rather than assuming it selects a new H.

Let A,B,C be nonempty tournaments, with a named root r in A. Start with
T+=C3[A,B,C], oriented A→B→C→A. Reverse every arc BETWEEN the three
blocks, retaining every internal arc and every vertex name, to obtain T−.
The labelled cyclic-triangle hypergraph is unchanged: triangles inside
one block are untouched, triangles meeting two blocks remain acyclic,
and those meeting all three remain cyclic. Thus full Q is unchanged,
including after any fixed vertex deletions. In particular H and the
root-deletion count q are the same. The root star generally changes.

The operation applies independently at any node of a recursive C3 tree:
each pair of leaves has only one least common ancestor, so these edge sets
are disjoint. In the address chart, interchange two child digits at that
node. Common-prefix valuations remain fixed and the leading orientation
units reverse at exactly that split. For fixed original lexicographic
labels i, the chart x_i=3^d i+r_new(i) still has positive minors, since
0≤r_new(i)<3^d. This is a physical orientation change on the quotient,
not a relabeling of the protected root.

## The exact root-start determinant

Let c_X(k) count ordered covers of X by k nonempty directed paths, and put
c_X(0)=0 for nonempty X. Let L_A(k) count the same ordered covers of A
with the FIRST path starting at r. Write s± for the numbers of Hamilton
paths starting at r in T±. Then

    s+ − s− = Σ_(k≥0) L_A(k+1)
                    [c_B(k+1)c_C(k) − c_B(k)c_C(k+1)].       (1)

This is a weighted sum of adjacent two-by-two minors of the two response
sequences. It is not an identification of the response minors with the
numeric minors of the address chart.

Proof: compress a Hamilton path into its maximal block runs. A word
starting in A has the cyclic pattern A,B,C,A,B,C,... . Its possible visit
counts are

    (k+1,k,k), (k+1,k+1,k), (k+1,k+1,k+1).

Each pattern contributes the product of its ordered cover counts, with
L_A in the first coordinate. Reversal interchanges B and C. The first
and third products cancel in the difference, leaving (1). All sums are
finite, bounded by the respective module orders. This is the exact
maximal-run decomposition, not a numerical fit.

For A=TT2=(r,a), with r→a, one has L_A(1)=L_A(2)=1 and no higher term.
Writing h_X=H(X) and p_X=Q_X(2), c_X(2)=p_X−2h_X gives

    s+ − s− = p_B h_C − h_B p_C.                           (2)

Consequently the direction minimizing s is selected by comparing the
exact ratios p_B/h_B and p_C/h_C. The difference is even because each
p_X is even; the normalized determinant using p_X/2 is integral. A
nonzero minimum difference is two, not one. This parity boundary matters
when using the move to change a readiness defect by two.

## An arbitrary-module preparation with all four responses

Use the root-state convention (H,s,q,F), where s counts starts at r,
q=H(T−r), and F counts ordered covers by two paths with the first nonempty
and ending at r and the second possibly empty. Put

    h=H(M), p=Q_M(2), w=Q_M(3),
    T+=C3[TT2,M,1], T−=C3[TT2,1,M], root r first in TT2.

The second notation describes the reversed orientation, with the same
internal modules and the corresponding fixed names. The exact states are

    T+: (3p−h, p,  p+h, 3p),
    T−: (3p−h, 2h, p+h, 2w+p−2h).                        (3)

The complete maximal-run coefficients, in the basis c_M(1),c_M(2),...,
are respectively

| Observer | T+ | T− |
|---|---|---|
| H | (5,3) | (5,3) |
| s | (2,1) | (2,0) |
| q | (3,1) | (3,1) |
| F | (6,3,0) | (6,7,2) |

To certify completeness, use the three exterior vertices r,a,b. One
path has at most four maximal M-runs; two paths have at most five.
Enumerate every word in those exterior vertices and identical M-run
symbols, excluding adjacent M symbols within a path. For F split into
two ordered paths, with the first ending r. This finite exterior-word
enumeration gives the table and zero coefficients at every higher
possible run. It is reproduced in chirality.py; substituting
c1=h,c2=p−2h,c3=w−3p+3h proves (3) for every nonempty M.

T+ has exactly one incoming edge at r. Its readiness defect is
q−s=h, whereas the other marked presentation has q−s=p−h and generally
has more than one incoming edge. Thus chirality can prepare a unit port
when h=1 (M transitive), but does not make an arbitrary module unit-ready.
For M=TT2 the same H11 and q5 have s4 versus s2: the defects are1 and3.
The latter root is not an indegree-one root, so this does not contradict
the previous [indegree-one cyclic-triple congruence](../symmetry_port_return_20261113/history/MARKED_HYPERGRAPH.md).

The inherited universal NR birth identity reads H_NR=H+2F. The two
prepared outputs therefore have exact targets

    Y+=9p−h,              Y−=4w+5p−5h,                    (4)
    Y−−Y+=4(w−p−h).

The selector now has an explicit arithmetic obligation: supply an actual
module with the requested response, then select the appropriate neutral
orientation. It does not supply that module for every integer Y.

## A boundary that optimal marking cannot remove

For an old vertex r in a leaf of a C3 hierarchy, its indegree is

    indeg_leaf(r) + Σ_(ancestor v) |incoming sibling block at v|.  (5)

Every summand is positive. Hence depth≥2 forbids indegree one at every
old leaf vertex, regardless of local chirality choices or optimal marking.
At one fixed root path, the complete degree-choice enumerator is

    z^indeg_leaf(r) ∏_v (z^a_v+z^b_v),                    (6)

where a_v,b_v are the two sibling orders. Coefficients count orientation
assignments, not distinct graph isomorphism classes.

There is also a response obstruction to one fresh active-root birth. If
r lies inside a nontrivial STRONG block A (|A|>1) at a C3 node, then L_A(1)≥1
and L_A(2)≥H(A−r)≥1. The first is a consequence of a Hamilton cycle in
a strong tournament; for the second take r as the first singleton path.
The k=0 and k=1 terms in the start formula give

    s_T(r) ≥ [L_A(1)+L_A(2)] H(B)H(C) ≥2.              (7)

At depth≥2 the child containing r at the top split is itself strong, so
this applies even if the terminal leaf was a singleton. The new active
root after one N birth has q_new−s_new=s_old, by

    N(H,s,q,F)=(H+F−s,H−s,H,H+F).

Thus one such birth does not produce unit readiness from an old port in
this deep hierarchy. A newly selected port may still break global
symmetry and be primitive; neither fact implies unit readiness. The
boundary disappears only after an operation leaves this fixed hierarchy
or changes the interface more substantially. The new neutral module
repair in selection/ is an independent physical operation of that kind,
not a consequence of marking alone.

## Source, map, losses, and test

Source: the signed child digit of an ancestry chart. Map: reverse the
three between-block orientations. Preserved: labelled cyclic triples,
full Q, every fixed deletion count, names and internal arcs. Changed:
the root star, start count, and two-path endpoint response. Required
extra data: L_A and the ordered cover sequences of B,C. Decisive test:
(1), specialized to the literal H11 pair and independently enumerated
permutations. The main consequence is a lawful neutral response selector;
the missing arithmetic source supply is kept explicit.
