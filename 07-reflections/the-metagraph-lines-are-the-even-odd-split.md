# The metagraph's own lines are the even/odd split

*mac-mini-2026-07-01-S83. Reflection on HYP-3808.*

The owner had been staring at the tiling explorer and noticed a rhythm in the colored lines of the merged
metagraph: the merged (non-self-complementary) nodes send out an even number of black lines; the nodes they
land on carry an even count of black and an odd count of blue; and off to the side sits a third kind of node
touched only by blue, and an odd number of it. A parity bookkeeping, node by node, in three categories. The
question was how the total supply of lines gets distributed so that every node ends up with exactly its
right number of tilings.

Compute it out and the bookkeeping resolves into one sentence. **The black lines form an even graph and the
blue lines form an odd one.** Every node, whether it is a merged non-self-complementary pair or a
self-complementary singleton, has an even number of black lines to its neighbours — so the black subgraph is
Eulerian, a disjoint union of cycles. And every self-complementary node has an odd number of blue lines — so
the blue subgraph is an all-odd-degree object, a T-join with the pure-blue nodes dangling off it as leaves.
The tiling count of a node is just its total line-degree, an odd blue part plus an even black part, and so
its *parity* is decided entirely by the blue: **self-complementary nodes have an odd number of tilings,
non-self-complementary nodes an even number.** The reason is almost too simple once seen — a
non-self-complementary node is a pair `{A, A^op}` and carries `2|A|` tilings, forced even; a self-complementary
node is a single class with an odd number of grid-symmetric tilings. The owner's three-category rhythm was
the shadow of this: the merged pairs are the even (black) world, the pure-blue singletons are the odd (blue)
leaves, and the mixed nodes are where the two worlds touch.

Which lands, unexpectedly, right on top of the previous session. Three days of reading resolved that "even
graph" means two different things — even-degree, and even-automorphism-parity — and that the project's E_n
is the even-degree one. And here is an even-degree graph, an Eulerian graph, sitting inside the merged
metagraph as its black lines, unbidden. The parity dual the project keeps circling is not only a separate
object hung beside the tournament metagraph; it is *drawn in the metagraph's own edges*, in the color the
explorer already paints. The blue/black coloring the owner had been looking at for a hundred sessions is the
even/odd decomposition, made visible.

And the owner's one guess was wrong, in the instructive way. The conjecture was that only the mixed nodes
self-pair — that a line looping a node back to itself, adding two tilings at once, happens only at the
interface. True through five vertices. At six it fails: the non-self-complementary nodes suddenly grow black
self-loops, and a whole sea of black lines runs between them, where before every black line had to touch a
mixed node. The layered picture — black world, then interface, then blue world, cleanly stacked — is a
small-case illusion. At `n=6` the black even-graph stops being bipartite-to-the-interface and becomes what
even graphs generically are: dense, self-connected, Eulerian in the full sense. The transition has a place
and it is exactly where the tiling hypercube first gets big enough to be typical.

The pattern that transcends the theorem: **a coloring someone chose for a picture is often a theorem in
disguise, and its exceptions mark where the small cases stop lying.** The blue and black were drawn to be
looked at; they turned out to be the even and the odd, the two halves of the parity the whole project is
about. The right response to a suggestive picture is not only to admire the pattern but to compute its
degrees and read off which classical object they spell — and then to push `n` one step past where the
pattern is clean, because the first `n` at which a conjecture breaks is usually the first `n` at which the
structure is finally itself.
