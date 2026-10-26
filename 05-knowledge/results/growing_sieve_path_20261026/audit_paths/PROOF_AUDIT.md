# Independent mathematical audit of used-arc reconstruction

**PASS.** Reviewed the frozen [path proof](../paths/RESULT.md), SHA256
`67e1b911741f09d29fbf072e9075ff06a438060732849405759e6b7dea1c53dd`.
No mathematical or scope correction is requested. The initially reviewed
candidate had hash `4332eae549913d0575ed668c579ec45a480f4e9e663ece3b705512eebde4be1f`.
Its promotion changed only the status line and added this audit link. Reversing
those two edits reproduces the candidate hash exactly, confirming that the
mathematical body is byte-identical. Substantive changes require renewed review.

The reviewer is the field/tournament lane, separate from the path producer.
This is a proof audit, not an inference from finite counts. The root lane's
independent dynamic-programming verification is packaged separately in this
directory and is not claimed as this reviewer's computation.

## 1. The exact primary dependency is sufficient

I independently located Arthur H. Busch's 2005 thesis *Arc-traceable
Tournaments* through the primary university repository, then read the pinned
local PDF's text on printed pages 45--47 and inspected the rendered statement
page. [Primary PDF](https://digital.auraria.edu/downloads/ma7yx-pcr25/AA00005831_00001.pdf).
The PDF hash matches the producer's pin:
`76c57cfb33b8ad32e7e830e811bb3a3a38d4a2eeea72416d03b45412f51cdf19`.

Theorem 4.12 applies to a **strong tournament** and an arc x->y unused by
all Hamilton paths. It supplies a vertex whose deletion separates the reverse
reachability, placing x and y in the first and last of at least four ordered
strong components. The theorem gives further dominance conditions, which the
producer also tests. Its proof begins with the reverse-path separator supplied
by Theorem 4.7 and Menger's theorem. The source explicitly warns that the full
displayed structure is necessary and not sufficient for nontraceability.

The reconstruction uses only this necessary implication. It does not apply
the statement to arbitrary digraphs, omit strongness, infer nontraceability
from a guessed cut, or reverse the direction x->y. The endpoint-component
clause also ensures that the separator differs from both endpoints. The thesis
number is correctly distinguished from the later 2006 journal article.

## 2. P1: component restriction is exact

A Hamilton path visits the tournament's strong components in their unique
linear order and cannot return to one already left. Arbitrary internal
Hamilton paths concatenate, so the restriction of the global used-arc union
to a component is exactly its own used-arc union. This is an equality, not
only one containment.

For a nontrivial strong component, breaking each possible edge of a Hamilton
cycle shows that all cycle arcs survive in the union; hence the component
remains strong there. Singletons are handled without pretending they have
cycles. Arbitrary prescribed starts and ends are available in a strong
component, so every cross arc between successive components can be used.
Non-successive components cannot be consecutive in a Hamilton path. Thus
both the component sets and their order are recoverable from the union.

## 3. P2: the deletion lemma survives the hostile cases

The claim is `R(U(T)-z)=R(T-z)` for a **strong** tournament T. The easy
containment follows from U(T) being a spanning subgraph. The difficult
containment was checked in both cases.

When T-z has at least two components, strongness of T forces a vertex a in
the first component with z->a and a vertex b in the last with b->z. Two
different extension families are needed. Prepending z to a path starting at
a preserves arbitrary chosen cycle arcs in every component except possibly
the first. Appending z to a path ending at b preserves arbitrary cycle arcs
in the first component, because first and last differ. A nontrivial strong
component has at least three vertices, so an edge different from the one to
retain can always be broken. Singleton components need no internal edges.
One concatenated path supplies every consecutive component cut. These facts
give precisely the same reachability; no backward cut can appear.

When T-z is strong, fix an arc e of one Hamilton cycle. If every break at an
edge f different from e failed both prepending and appending z, each such f
would force its head to dominate z and its tail to be dominated by z. Among
the remaining cycle edges are two consecutive edges, even when the cycle
has length three. Their common vertex would have both orientations with z.
Thus every cycle edge survives in the global union after deleting z.

Boundary orders are sound: a two-vertex tournament is not strong, a strong
three-vertex tournament falls in the multiple-component deletion case, and
the singleton is treated separately. The transitive three-vertex hostile
correctly demonstrates that the strong-source hypothesis cannot be dropped.
The proof does not silently identify U(T)-z with U(T-z), which is a different
and generally stronger statement.

## 4. P3: existence and consistency are separate, and both hold

For a missing pair inside a recovered component C, the classical theorem
gives a separator in the actual strong tournament T[C]. P2 transfers its
one-way reachability to the observable U[C]-z. Hence the procedure can always
find a qualifying z using only the supplied union.

For **any** qualifying z, a reverse actual edge would itself supply the
forbidden reverse reachability in T[C]-z. P2 rules this out. Consequently
all qualifying separators give the same orientation, independently of the
order in which missing pairs are processed. The implementation correctly
uses the original U for its witnesses, rather than using previously guessed
arcs and creating circular certificates.

Outside a component, its recovered order gives every omitted orientation.
Ordinary transitive closure alone is not enough inside a strong component;
the deletion relations supply the extra information. The stated O(n^4)
time and O(n^3) storage bounds are conservative for precomputed deletion
reachability. The claim remains a promised inverse of a complete used-arc
union; it does not claim arbitrary-digraph recognition or polynomial listing
of all Hamilton words.

## 5. Elementary existence facts and the all-order continuation consequence

The insertion proof of a Hamilton path is correct. The maximum-cycle proof
of the strong-tournament Hamilton-cycle fact is also correct: every outside
vertex must either dominate the whole cycle or be dominated by it. A path
from the latter class to the former must contain an arc between them before
returning to the cycle. Inserting those two vertices yields the contradiction.
This covers the prescribed-start/end uses in P1 and P2.

I reread the original strong-atom C5 state definition. Its operations are
named permutations p of all positions, with legality exactly membership in
the current labelled Hamilton language and a distinguished failure state.
Equality of all guard bits therefore is equality of those languages. The
proved inverse recovers the labelled tournament, and hence the literal
presentation in that fixed-path gauge. This does close the formerly open
all-order named-guard injectivity statement. It does not address re-marking
after forgetting operation names or aggregating successor isomorphism
classes; the producer correctly leaves HYP-9416's corresponding questions
open.

Relabelling commutes with the union and with reconstruction. The claimed
equivalence of tournament isomorphism and union-digraph isomorphism follows.
No conclusion about black-core connectivity, an unrestricted recurrence,
Hamilton counts, or minimal unnamed response compression is imported.

## 6. Code and hostile-certificate review

I reviewed `check_paths.py` at its frozen hash in the manifest. Its two
enumerations are distinct: filling constraint subcubes for each permutation
versus a direct vertex-by-vertex DFS in each tournament. Every labelled
tournament of orders one through six is included. Reachability is restricted
to the requested vertex set, component sorting uses the strict downstream
order, and internal certificates exclude the endpoints. The code tests the
strong deletion lemma componentwise and all four classical separator clauses.
It contains explicit failure gates that survive optimization.

The two frozen producer receipts are byte-identical and report 321,237
gates. I checked those bytes and source hashes; this paragraph does not claim
an additional full producer rerun. Separately, I directly enumerated all 5!
vertex orders of the displayed first strong hostile `(20,1,2,7,14)`. It has
exactly nine Hamilton paths and exactly one unused arc, 4->1. Deleting 0
orders the components as `[4],[3],[2],[1]`, giving the correct orientation.
The legal path `(1,0,4,3,2)` also verifies the warning that nonuse does not
force a fixed relative order of its endpoints.

The five-vertex lower boundary agrees with the classical requirement of at
least four components after deleting one vertex. The incomplete directed
path versus transitive tournament control correctly marks the completeness
assumption. All reviewed claims pass at their stated quantifiers.
