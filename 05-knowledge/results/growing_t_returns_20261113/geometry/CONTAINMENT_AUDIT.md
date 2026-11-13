# Independent audit of the ready-target containment obstruction

**PROVED:** the weak-source reduction and operation implication below.
**FINITE-EXACT:** the independent all-subset recount in `containment_audit.json`.
**INHERITED:** all-order completeness of the unique-start parent grammar
and its strong-order bound. This audit does not independently re-enumerate
all parent classes. No global spectrum hole is claimed.

## Weak unit-ready reduction

Let T have H(T)>1 and a unit-ready root r, meaning F_T(r)-H(T)=1.
Write its ordered SCC decomposition as A->C->B, where r belongs to C;
either exterior factor may be empty. Put a=H(A), h=H(C), b=H(B), with
H(empty)=1. Let P2(A) count two ordered spanning path bins, either allowed
empty, and P2(empty)=1. Concatenating the restrictions of the two paths
across the SCC order gives exact identities

    H(T)=ahb,    F_T(r)=P2(A) F_C(r) b.

The path ending at r cannot visit B, so every vertex of B lies in the
other path. The two paths can independently use their A restrictions;
these are exactly the P2(A) choices. This proves the factorization with
empty-bin boundaries included.

Thus b(P2(A)F_C-ah)=1, forcing b=1. If A is nonempty, P2(A)>=2a
because either bin can carry a Hamilton path. Also F_C>=h: cut every
Hamilton path of C immediately after its occurrence of r. Consequently
P2(A)F_C-ah>=ah=H(T)>1, a contradiction. Hence A is empty, C is strong
and unit-ready at the same H, and B has H1, so B is transitive. The
singleton C possibility would force H(T)=1 and is excluded.

The exact hereditary-product law in [PROOF.md](PROOF.md) now shows that
the induced count spectrum of T equals that of C. This proves the
all-order reduction, rather than assuming weak padding harmless.

## Census and operation quantifiers

The inherited complete target counts are1691,3311,6551,13031,25991, with
respectively3,3,4,2,2 strong rooted classes. Their predecessor counts are
881,1691,3311,6551,13031. The source
[next_target.py](../../multicore_metagraph_20261113/next_target.py) enumerates
every rank in the unique-start grammar through its stated order18. For
each target N, the parent H is at most(N+1)/2; the cited increasing
strong-order minima exclude all orders above its bound. The receipt
records that each required bound lies inside the enumerated range.
This is the completeness dependency; the new subset test does not prove
completeness by itself.

The root's [count_spectrum.py](../containment/count_spectrum.py) takes
every subset of all fourteen literal target graphs, including the empty
set. It compares endpoint-deletion and endpoint-extension recurrences,
checks literal permutation controls, verifies the whole-target count and
unit readiness, and rejects no subset before testing its Hamilton count.
Its integer arrays are safe: the largest graph has17 vertices, so even
the crude17! bound is below signed64-bit capacity.

The independent checker here recomputes every subset count using Python
unbounded integers and compares the complete sorted count multiplicities
against the root receipt. It also checks the missing predecessor and
all-order completeness flags. Total universe:321536 subsets of fourteen
literal graphs, of orders11 through17.

Therefore **every** target in each complete unit-ready fibre lacks an
induced predecessor-count tournament. No birth-only program from **any**
source of that predecessor count reaches the fibre. This allows arbitrary
initial same-count preparation, arbitrarily many births, arbitrary cuts,
and any final root. Initial preparation does not remove the obstruction:
the prepared source itself must survive induced. The conclusion concerns
five adjacent count pairs; it is not an induction over the entire405 ray.

Old-arc changes after a birth, deletion of vertices, or replacing the
whole source can escape the invariant. A proved endpoint repair of one
of these kinds does not contradict the birth-only exclusion.

## Terminal-newborn bound

There is also a universal bound when the final ready root is the last
newborn. For H>1, the ready classification gives root indegree1. Its
sole predecessor x lies in the deletion parent U, and the newborn formula
gives F_new-H_new=s_U(x)=1 and

    H_new=H(U)+F_U(x)-1 >= 2H(U)-1.

If U was obtained from a source of count h by vertex births alone, it
contains that source induced and H(U)>=h. Hence H_new>=2h-1. The desired
next405-ray count2h-71 violates this bound. This argument needs no bounded
source census, but only excludes programs whose final root is the last
newborn. The five finite-head hereditary exclusions remove that root
restriction for their stated count pairs.

Reproduce with `python -B .../geometry/containment_audit.py`, and again
with `-O`. The source hashes are recorded with the audit receipt.
