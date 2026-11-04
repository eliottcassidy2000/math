# Two actual next returns and a direct native guard compiler

**PROVED:** the pair-block guard construction below, composing with the
inherited one-helper arc-flip theorem. **FINITE-EXACT:** the two next named
endpoints and all saved guards/counts. No arbitrary-depth seed policy or
coverage theorem.

The actual repaired parent has H13835447 and23 named vertices. Its two
chosen reader targets have H41506343 and41506345, respectively orders26
and27. They retain every parent physical name, add3 and4 permanent vertices,
and use40 and46 native complements, with20 and23 temporary helper births
and deletions. The outer root39, active mark4, old mark31 and spare6 remain
distinct names. Prefix50→51 remains a prefix; the former additional prefix
vertices1,2 are now part of the actual component. Native intermediates need
not preserve H, strongness or the complete port stars.

## A direct guard avoids subset search

For any tournament T and distinct vertices a,b, construct a Hamilton word Q
on `T-{a,b}` by inserting each new vertex immediately before its first
out-neighbor in the current word, or appending if none exists. Before an
insertion point every old vertex beats the inserted vertex, so this elementary
procedure preserves legality.

Find the first vertex v of Q beaten by a or b. Insert the consecutive block
(a,b) or(b,a) immediately before v, choosing its last vertex to beat v.
Every preceding Q vertex dominates both a and b, so the incoming join is
also legal. If no such v exists, append either block: every Q vertex then
dominates both endpoints. Empty Q is allowed.

Every arc of the resulting full word except possibly the a,b arc is legal
in T. Choose that pair's orientation to agree with the word. The word is
therefore a Hamilton path containing the designated pair, in either T or
its single-pair reversal. Constructing it uses O(n²) adjacency tests and no
Hamilton-path search.

Apply the existing traced-arc gadget from **METAGRAPH section5,
`adaptive_metagraph_reformulation_20261103/METAGRAPH.md`**. It adds one helper,
uses two legal off-path complements, and deletes the helper to reverse
exactly the chosen old pair. If the constructed word was legal in the
reversed tournament, traverse the two involutions backwards. This is the
inherited theorem and gadget; the direct block guard supplies a simple
polynomial implementation of its one-pair separation step. No new claim
about general reconstruction from many changed arcs is needed.

The compiler keeps actual labels throughout; it performs no physical
relabelling. Its complexity statement concerns realizing a supplied target
matrix, not finding an arithmetic target witness. Repeating it over the
explicit discrepancy list is polynomial in the vertex order and list size.

## Exact graph and interface checks

`followup.py` builds every permanent vertex from its recorded cut and retains
every pair guard, both native words and every intermediate endpoint digest.
It replays saved programs on every run. The same script independently counts
the actual endpoint graphs using the previously audited small-feedback gap
engine, with structural feedback sets of size3. It checks the full marked
component, the exact wrapper, strongness and the seed-reader recurrence.

Each target retains the four ordered pairs `(39,6),(31,6),(4,6),(50,51)`.
For every pair, the receipt contains H, P2, both endpoint/deletion profiles,
ordered W2 and all17 universal static future-H queries. It also records
rooted golden and current F49/F343 codecs for all six named vertices.
These codecs are observers, not tournament automorphisms.

`followupaudit.py` imports only the old elementary independent native
interpreter, not this compiler or its helpers. It replays all86 guarded
complements,43 helper deletions and7 permanent births from the exact same
repaired parent. It verifies every final matrix, the surviving names and
roles, and recomputes all named field codecs. It explicitly does not claim
a second H recount in that native interpreter; the independently implemented
gap counts are bound to those exact matrices.

The direct compiler also supplies410 literal witnesses covering every
labelled tournament and pair at orders2–4. The independent interpreter
reconstructs each tournament, checks the pair word in one orientation,
checks both native guards and verifies the exact single-pair endpoint.
Illegal directed and repeated-vertex guards are rejected separately.

Run both scripts normally and with `python3 -O`. `followup.py` writes its
receipt and labels partial checkpoints explicitly; `followupaudit.py`
prints JSON. Final normal/optimized receipts are byte-identical. All sources
are pinned, and every gate remains active under optimization.
