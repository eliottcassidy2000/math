# A golden unit switch that cannot be flattened to one five-core birth

**PROVED IN NOTE:** an all-address unit-switch formula, its full rooted
sidecars, an exact deletion floor, and the exact scope of common O/Z
continuation. **FINITE-EXACT:** the bounded image and graph controls below.
No novelty claim for the inherited O/Z maps or the known H13/69/71 source.
Five-core value coverage and Chan–Pak remain **OPEN**.

The starting object is the old July spectral hostile
`Q=(18,8,3,5,14)`, H13, whose directed-triangle arc-overlap is the golden
P4. The previous [history note](../INCIDENCE_AND_HISTORY.md) connects it to
THM-1780 and THM-5117. Here the new question is constructive: does its
69/71 unit switch persist at arbitrary address length, and can the actual
graphs be flattened into the five-core one-birth grammar?

## 1. A portable switch at every right address

Use the exact fixed-core compiler of
[Session251](../../../session251_creative_completion_20261111/ternary/PROOF.md),
with left address A, right address B, parity bits α,β, and

    g(0)=g(1)=1,
    g(2q+1)=g(q),
    g(2q)=g(q)+q+(q mod2).

For the two apex out-cuts `O+={0,1}` and `O−={1,2}`, its nine signatures
`(d,k+,k−,i+,i−,o+,o−,j+,j−)` are

    O+: (13,91,−13,50,−8,32,−6,18,2),
    O−: (13,91,−13,47,−9,32,−6,15,3).

Their difference in the exact two-flank output is therefore

    H+(A,B)−H−(A,B)=3A−α+2.                      (1)

In particular the difference is independent of B. At A=0 put

    L(B)=13g(B)+32B+6β+18.

Then **every B≥0 supplies actual tournaments with counts L(B),L(B)+2**.
The first examples are31/33 at B0,69/71 at B1,121/123 at B2,
133/135 at B3,211/213 at B4. This is an unbounded word family, not just
a fixed numerical unit gap. Formula(1) also identifies exactly what fails
when a nonempty left address is introduced: the extra source response is
`3A−α`, so the gap stays2 precisely when A=0.

The signatures are finite integer sums of weighted simple paths in Q.
Their derivation follows directly from the compiler definition; the local
script computes them from the literal arcs. Once these constants are fixed,
equation(1) holds for all address lengths by the inherited compiler theorem.

## 2. The whole rooted sidecar agrees

Write a rooted state as `(H,E,U,d)`, where E counts full paths ending at
the root, U counts ordered spanning two-path covers whose distinguished
path starts at the root, and d is the root-deletion path count.

At B0 the states are

    (31,7,45,13),    (33,7,45,13).

For B>0, let b be its bit length and `rev_b(B)` its bit reversal in b bits.
The states are exactly

    (L(B),  13, 45·2^b−13rev_b(B),13),
    (L(B)+2,13, 45·2^b−13rev_b(B),13).             (2)

**Proof.** Both cuts have full-path ending selected sum6 and weighted
path-start sum32 in Q. Thus before a tail is attached, E=13−6=7 and
U=13+32=45. The deletion graph is always `Q→TT_b`, so d=13.
The canonical right word ends at an unselected vertex (B's leading bit1),
so every deletion Hamilton path ends at a vertex that beats the root,
giving E=13. Its selected-word integer is
`2^b−1−rev_b(B)`. Weighted starts in the ordinal product are
`32·2^b+13(2^b−1−rev_b(B))`; adding the root-alone contribution13 gives U
in(2). This proves all sidecars without assuming equality of H.

This specializes the inherited pure-H common-continuation mechanism; it
does not establish that any desired absolute baseline can be prescribed.

## 3. Sharp all-length obstruction to flattening the actual graphs

**Theorem.** For every B≥2, both displayed tournaments have minimum
one-vertex-deletion count exactly13. The lower family already has this
property for B≥1. At B0 both minima are5; at B1 the upper minimum is9.

**Proof.** Deleting the outer apex leaves `Q→TT_b` and hence gives13.
Deleting a tail vertex merely deletes one binary digit in B; after removal
of harmless leading zeroes, the output is the same family at an address B'.
Its count is at least31 in the lower family and33 in the upper family.

For deletion of a core vertex v, the count has the form

    d_v g(B)+o_v B+t_v β+c_v.

The exact coefficient table, with rows v=0,...,4, is:

| v | Lower `(d,o,t,c)` | Upper `(d,o,t,c)` |
|---:|---|---|
| 0 | `(3,6,2,2)` | `(3,3,1,2)` |
| 1 | `(5,6,2,4)` | `(5,5,1,4)` |
| 2 | `(5,6,2,4)` | `(5,13,3,8)` |
| 3 | `(3,9,3,2)` | `(3,9,3,2)` |
| 4 | `(5,14,2,10)` | `(5,11,3,4)` |

These follow by deleting v in Q and applying the same compiler. All
coefficients are nonnegative and g(B)≥1. Every lower entry is≥13 for
B≥1. Every upper entry is≥13 for B≥2, except that the crude bound on
row0 at B2 is only11; the exact value g(2)=3 gives17. For B≥3 even that
row is at least14, hence at least15 since a tournament count is odd.
The displayed boundary values follow by direct substitution. QED.

Every literal one-birth five-core construction has some vertex whose
deletion count is1,3,5 or9: its construction apex. The theorem therefore
rules out **every isomorphic presentation** of our B≥2 graphs in that
grammar, regardless of which vertex is chosen as the purported apex.
This is stronger than failure of one preferred source representation.

Yet two births suffice. Delete original vertex1 of Q; the remaining
four-vertex tournament has rows `(8,1,3,6)`, H5, so it is S4 up to
relabelling. Reinsert vertex1, letting it beat the same tail uniformly,
to obtain `Q→TT_b`; then insert the outer root with the selected cut and
arbitrary right address. Thus these graphs form a strict depth-two family
over the same small-parent ingredients. The scalar counts may still have
different one-birth five-core witnesses; this theorem concerns actual
graphs and does not refute five-core **value** coverage.

The coordinate that prevents flattening is now precise: the minimum
deletion count of the current source, not H alone or the golden overlap
spectrum. Any target-preserving flattening must change the tournament,
rather than only re-root or relabel it.

## 4. Exact common O/Z continuation: genuine local intervals, sparse global image

The inherited physical operations are

    O(H,E,U,d)=(H,0,2U,d),
    Z(H,E,U,d)=(H+U−E,d,2U−d,d).

O appends a sink; Z appends a vertex beating just the retained root.
Read B from least to most significant bit, with0=O and1=Z. Applied to
the B0 pair, the resulting graphs are exactly the family above. Conversely,
every O/Z word has this form after deleting terminal O's, which preserve
H. Thus the **entire** common O/Z continuation image of these twins is

    S={L(B),L(B)+2 : B≥0}.                           (3)

This is substantially larger than one repeated-letter ray, but it is not
tail coverage. Since `L(B)≥32B+31`,

    |S∩[1,N]| ≤ 2 floor((N−31)/32)+2                 (N≥31).

Its upper density is at most1/16 among all positive integers, or1/8 among
the odd integers. This is an elementary all-height bound even if the image
has collisions or a density limit fails to exist.

Named subfamilies are particularly transparent:

    Z^k, k≥1:       32·2^k+{5,7},
    O^k Z, k≥1:     45·2^k+{31,33}.

Branches do interlace locally. Complete inversion of(3) through target
1,000,000 gives49,861 distinct values, about9.97% of positive odds. The
longest consecutive-odd run in that bounded image has nine members,

    473005,473007,...,473021.

This is a useful **FINITE-EXACT** interval signal; it does not bound or
prove growth of the maximal run at larger heights. The exact density
ceiling proves that common continuation alone cannot establish the full
spectrum. A proof using these local intervals still needs a way to place
them near an arbitrary target, or a target-preserving change of source
outside common O/Z continuation.

## 5. Scope and reproduction

`python3 -B 05-knowledge/results/five_core_web_20261111/history/continuation/probe.py`

The script pins the inherited compiler by hash, computes both signatures
and every coefficient in the deletion table, checks literal full graphs
for B0–15 by subset DP, and independently follows the four-coordinate
O/Z maps. The bounded image through1,000,000 is complete because the
proved lower bound limits B to31,249. Its finite run length is never used
to prove an all-height claim.

The new structural contribution is the strict actual-graph nesting
hierarchy and its explicit source coordinate, together with the golden
H13 lineage. The inherited pure-H continuation laws retain their credit.
