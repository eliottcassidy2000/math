# From the ordinal triangle to the two-root normal form

**PROVED elementary normal form; INHERITED two-birth count law, rederived by
forced blocks; FINITE-EXACT literal controls.** No all-target construction is
claimed. This isolates the smallest way to retain two indegree-one vertices
while breaking the ordinal C3 factor required by the short mutation return.

## 1. Three roots are intrinsically ternary

The following degree argument is inherited from
[ANCESTRY A1](../../adaptive_metagraph_reformulation_20261103/ANCESTRY.md).
A tournament has at most three vertices of indegree one: four would have
six internal edges but total indegree only four. If there are three, their
three internal edges exhaust the total indegree, so they form a directed
triangle and every exterior edge points out from that triangle. Conversely
the first SCC of C3 -> R supplies exactly three indegree-one vertices.
Thus three such roots force H=3H(R) under any construction, not just mutation.

## 2. Exactly two roots have a canonical exposed-triangle form

Let a,b be two indegree-one vertices, orient their edge a->b, and let c be
the unique in-neighbour of a. Then necessarily

    a -> b -> c -> a,
    a,b -> R,             R=V minus {a,b,c}.

Indeed b is beaten only by a, and a only by c. Hence every remaining
vertex loses to both a and b. The vertex c may have arbitrary contacts
with R. Put I={r in R : r->c}, and W=T minus {a,b}=R+c.

Conversely any such choice of R and I gives indegree one at a,b. Every
vertex of R has indegree at least two, and c has indegree1+|I|. Therefore
I is empty exactly when there are three roots and an ordinal C3 factor;
I nonempty gives exactly the two roots a,b. Since the oriented edge a->b
distinguishes the two roots, W,c and I are canonically recovered from an
unmarked graph with exactly two indegree-one vertices.

This normal form is exactly two successive inherited adaptive births from
the rooted source (W,c): append a beaten only by c, then b beaten only by a.
It is therefore a structural reading of an existing constructor, not a
new unrestricted operation. See [READY](../../two_entry_funnels_20261103/READY.md)
and [ADAPTIVE](../../adaptive_metagraph_reformulation_20261103/ADAPTIVE.md).
The literal graph and ordered names a,b,c remain necessary sidecars.

## 3. The count exposes the missing attachment statistic

Write h=H(R), d=H(W), and let U_I count ordered spanning two-path covers
of R whose first path is nonempty and ends in I; the second path may be
empty. Both paths use disjoint vertex sets covering R. Put U_empty=0 and
H(empty)=1. Then

    H(T)=d+2(h+U_I).                                    (1)

Proof by first vertex. Paths beginning at a must next visit b, because b
can only be entered from a; their suffix is an arbitrary Hamilton path
of W, giving d. Paths beginning at b must contain c,a as a consecutive
block; contraction of that block leaves an arbitrary path of R before it
ending in I (or empty), and an arbitrary complementary path after it.
There are h+U_I choices. Paths beginning at neither a nor b must contain
the consecutive block c,a,b; they have exactly the same h+U_I choices.
These three classes are disjoint and exhaustive. QED.

Equivalently let F_W(c) count ordered spanning two-path covers of W with
nonempty first path ending at c and possibly empty second path. Removing c
from the first path gives F_W(c)=h+U_I. Hence (1) is the inherited identity

    H(T)=H(W)+2F_W(c)=3H(W)+2(F_W(c)-H(W)).               (2)

This displays the actual ternary digit: it is the attachment defect
F_W(c)-H(W), not an arithmetic choice independent of the rooted source.
The exposed cut I is the datum destroyed by retaining only H(R) or the
fact that a directed triangle is present.

For I empty, c is a source of W, d=h, and (1) reduces to3h. For R a single
vertex and I=R, d=h=1 and U_I=1, so H(T)=5 with exactly two roots. This
is the unique strong four-vertex tournament and is the smallest positive
escape from the factor3 constraint. Nonempty I does not guarantee that H
is coprime to3; the formula retains the exact integer data instead.

The [return lane](../returns/PROOF.md) studies which of these roots has
the next required readiness defect. The count identity itself supplies no
repeatability theorem. In particular degree one and defect one are
different predicates.

## Connection audit

Source: tournament with two indegree-one roots. Target: rooted deletion
graph (W,c), or exposed triangle (R,I). The normal form is a bijection on
the stated exact-two-root class when I is nonempty. It preserves the full
graph if the cut and names are retained. The scalar equation loses the
individual path-cover attachments and subsequent root choices. The minimal
hostile to an ordinal reading is H5 above; the minimal positive ordinal
control is C3 at R empty. Formula (1) extends to arbitrary digraph R, with
the two-root classification itself remaining tournament-specific.

The checker tests every tournament R through order3 and every cut I,
using literal Hamilton permutations and independently counted path covers.
It also recovers the normal form from every labelled tournament through
order5 with at least two indegree-one vertices. These finite checks verify
the implementation; the degree and forced-block proofs supply all orders.
