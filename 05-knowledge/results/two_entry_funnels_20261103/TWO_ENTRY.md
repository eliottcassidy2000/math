# Two starting paths: a complete rooted grammar

**PROVED all sizes, using the cited complete unique-cycle grammar.** A
tournament has one directed arc per unordered vertex pair. Paths are
ordered, directed and spanning. Converse is not quotiented here. Let
`s_T(r)` count paths starting at r, and count isomorphism classes with r
distinguished. Write `D=1-3x+x^2`. All generating functions below are ordinary
class generating functions, not labelled exponential generating functions.

## 1 The inherited SCC inverse and its two-entry boundary

Delete r, and order the strong components `U1 -> ... -> Uk`. Every Hamilton
path of the deletion passes through these components in that order. Hence

    s_T(r) = (product_(i>=2) H(Ui))
             * sum_(v in U1, r->v) s_U1(v).             (1)

A nontrivial strong component has a directed Hamilton cycle, and its
rotations give at least three paths. Every vertex in such a component
starts at least one. Thus, when s is1 or2, every component after U1 is a
singleton. For s=2, the nonempty first-entry set J is exactly one of:

* one vertex of U1 with two starting paths;
* two vertices of U1, each with one starting path.

This is an iff decomposition. Build T by taking a strong core U, a
transitive tail `v1->...->vt`, with U dominating the tail, and a new r
which beats exactly J in U. Encode the tail by `epsilon_i=1` iff `vi->r`.
For t>0 the source is strong iff the last bit is1 (J is already nonempty).
For t=0, J must be a nonempty proper subset of U. The two-entry cases
automatically have proper J. The inverse is unique because the SCC order,
all orientations at r, and the distinguished root are intrinsic.

## 2 Exactly which unique-cycle vertices are unique starts?

If a strong T has `s_T(r)=1`, it has exactly one Hamilton cycle: cutting
any cycle immediately before r injects its cycles into these starting
paths. The converse at a *specified* root is not automatic.

Use the complete classical maximum-outdegree insertion grammar in
[the inherited synthesis](../operation_closed_synthesis_20261029/THEORY.md#ii-a-classical-grammar-with-a-new-finite-observer-experiment).
At n>=4 a unique-cycle tournament has one or two vertices of indegree one.
In type A there is one maximum-outdegree vertex a; legal parents of an
indegree-one newborn are `pred(a),a`. In type B the two maximum vertices
are `a0->a1`; legal parents are `pred(a0),a0,a1`. The new vertex is the
canonical newest vertex of the resulting unique-cycle object. Conversely
these choices are exactly the grammar's children.

Adding a vertex x beaten only by p creates exactly `e_T(p)` Hamilton
cycles, where e counts paths ending at p: each cycle uses `p->x`, and
deleting x leaves precisely such a path. Therefore the listed parent
sets are exactly the vertices with e=1. Apply converse to obtain:

**Entry-set theorem.** On the unique cycle, let M be the outdegree-one
vertices of T. For n>=4, M has size one or two; in the latter case its
vertices are consecutive. Then

    {v : s_T(v)=1} = M union succ_C(M).                 (2)

Its size is two or three. For C3 all three vertices are unique starts.
At n>=4 these tournaments are rigid: every automorphism rotates the
unique directed cycle, and preserves the characteristic set M of size
one or two. A nonidentity rotation acting on M would have order two in
the latter case. No tournament admits an involution (a swapped pair's
arc would reverse). Hence the rotation is the identity.

The classical type counts obey `A'=A+B, B'=A+2B`, with `(A4,B4)=(0,1)`.
Converse gives the same counts when types are defined by outdegree-one
vertices. Their generating functions are `A=x^5/D`, `B=x^4(1-x)/D`.
This is an inherited classical mechanism, not a new ordinary UHC census.

## 3 Entry-pair inventories and the new two-start recurrence

Let P_n count unique-cycle strong cores marked with an unordered pair of
unique-start vertices. An automorphism cannot interchange the vertices
of such a pair, since it must preserve their mutual arc. At n>=4 rigidity
gives `P_n=A_n+3B_n`; at n=3 the three pairs form one orbit. Thus

    P(x) = x^3(1-x^2)/D.                              (3)

Let Q count only pairs adjacent on the unique cycle. A type A entry set
has one adjacent pair; a type B entry set has two. C3 contributes one
pair orbit. Consequently

    Q(x) = x^3(1-x)/D.                                (4)

Every legal core and entry set can have an empty tail or a nonempty binary
tail ending in1. Including the new root, their attachment series is

    K(x) = x + x^2/(1-2x) = x(1-x)/(1-2x).

Let V count strong rooted two-start classes. The unique inverse in (1)
and the two disjoint entry cases give `V=K(V+P)`, so

    V(x) = x^4(1-x)(1-x^2)/D^2.                        (5)

At orders4 through12 its coefficients are
`1,5,18,60,191,589,1775,5257,15360`. After the numerator range,

    v_n = 6v_(n-1) - 11v_(n-2) + 6v_(n-3) - v_(n-4).

All four denominator roots with multiplicity are necessary over Q:
the numerator has no common factor with D. The squared Fibonacci
denominator comes from feeding the classical pair inventory into the
same recursive funnel operator. It is a compositional convolution,
not a numerical coincidence or a finite fit.

A rooted two-start tournament is rigid. Its stabilizer acts freely on
its two starting paths, so its order divides2, and a tournament's
automorphism group has odd order. Thus each class contributes exactly
two fixed-spine charts whose first vertex is the root. The full fixed-order
faithful *class-address* set has V_n elements. This cardinality statement
does not claim a minimum response carrier under arbitrary complements.

For arbitrary, possibly weak tournaments with s(r)=2, r's component is a
strong two-start source and all later components are singleton. Thus the
all-source rooted class series is `V/(1-x)=x^4(1-x^2)/D^2`.

## 4 The Hamilton-cycle refinement

In either one-entry case, and in every nonempty-tail case, both paths
starting at r return to r. Thus there are two Hamilton cycles. With an
empty tail and two unique entries, the path from entry u ends at its
predecessor on the unique core cycle. It fails to return to r exactly
when that predecessor is the other entry. Exactly one failure occurs for
adjacent entries; none occurs otherwise. Consequently

    V_HC=1 = x Q = x^4(1-x)/D,
    V_HC=2 = V-xQ = x^5(1-x)(3-2x)/D^2.                (6)

This gives a separate geometric interpretation of the repeated denominator:
one-cycle objects form the ordinary Fibonacci layer; the two-cycle layer
adds the iterated entry choice.

## 5 Recognition and canonical labels

Equation(1) recognizes and recovers all starting paths when their number
is at most two, without general Hamilton-path counting. Delete the root;
an empty entry set means zero paths. If a later SCC is nontrivial or there
are more than two entries, the count exceeds two. For one entry recurse;
for two entries call the inherited unique-start recognizer twice. These
are complete tests by(1). The only branching is the terminal pair of
unique-start calls, so the bound is O(n^3) adjacency inspections. The
reference implementation returns the literal path words.

For a two-start source compare its two adjacency encodings labelled along
these paths; the lesser is a canonical rooted code. The encodings cannot
be identical because the rooted stabilizer is trivial. A relabelling is
retained as physical port transport. Ordinary tournaments possessing a
root of start budget at most two can be canonically labelled by doing this
at every eligible root, in O(n^4). This is a restricted-class algorithm;
the 81/27 hostile pair is outside its domain, not collapsed by it.
