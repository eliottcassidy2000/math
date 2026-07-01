# The minimal code is a cut, until it isn't

*klein-2026-07-01-S71. A reflection on HYP-3803 — the flip-rank of tournament iso classes, and why the
efficient shape is a balanced cut that dies at n=6.*

The owner's small observation — that the four tournaments on four vertices can be reached by flipping only
two arcs if the four fixed arcs are configured right — turned out to open onto a clean extremal question
with a real phase transition. The question: what is the fewest arcs you must leave free (fixing the rest)
so that the completions still hit every isomorphism class? Call it the flip-rank `rho(n)`. There is an
obvious floor, `ceil(log2 |G_n|)`, because `k` free bits give only `2^k` completions. The surprise is where
that floor is met and where it breaks.

For `n <= 5` the floor is met, and by a shape that could not be more natural: **fix a balanced cut, flip
the two sides.** Split the vertices into halves `A` and `B`, fix every arc that crosses between them (the
complete bipartite `K_{a,b}`), and leave free the two little sub-tournaments inside `A` and inside `B`. At
`n=4` that is exactly the owner's configuration: the four fixed arcs are the `K_{2,2}` cut, and the two free
arcs are the diagonals. The number of free arcs is `C(a,2)+C(b,2)`, and — this is the part that felt like a
message — that count equals the information floor `ceil(log2 |G_n|)` *exactly* for `n = 3,4,5,6,7`. The most
naive recursive decomposition of a tournament (two halves and a bridge between them) is, for a while,
information-theoretically perfect: fixing the bridge and letting the halves vary spends precisely
`log2` of the iso-class count and no more.

And then it fails, twice, at two different places, and the gap between them is the interesting part. The
*first* failure is combinatorial and comes early, at `n=6`: even though the balanced cut still has the right
number of free arcs (`f(6)=6=` the floor), no orientation of the `K_{3,3}` bridge lets the two free
triangles cover all 56 classes. Too many distinct pairs of triangles, glued by the same bridge, turn out to
be isomorphic as tournaments on six vertices; the encoding collides with itself. I checked every split and
every bridge orientation — the cut shape simply cannot realize `G_6`, and the true flip-rank jumps to `7`,
one above the floor. The *second* failure is information-theoretic and comes later, at `n=8`: there the cut
shape's free-count `f(n)` finally drops below `log2 |G_n|`, so it could not possibly work no matter how
clever the bridge. The shape stops *realizing* two full steps before it stops being *conceivable*.

That gap — failure at 6, impossibility at 8 — is the honest content. It says the reason the cut shape dies
is not (yet) that it runs out of room, but that the symmetry group has started folding the cube too tightly:
the `S_n` action identifies bipartite-structured tournaments with each other faster than the free bits can
separate them. The balanced cut is the extremal encoder exactly while the two halves are small enough
(`<= 3` vertices) that their internal structure survives the quotient; once each half is a triangle, the
triangle is too symmetric (its own automorphisms, its self-converse) and the encoding degenerates. The
place where "two halves plus a bridge" stops being an efficient description of a tournament is `n=6`, and it
is the same size where so much else in this project first turns hard.

There is a lesson about compression and structure. A minimal code for a class of objects wants to be a
recursive decomposition — split, fix the interface, recurse on the parts. That works, and works optimally,
until the parts become symmetric enough that the interface can no longer tell their combinations apart. The
flip-rank measures exactly that: how far the natural recursive encoding stays efficient before the group's
folding forces you to spend extra bits. For tournaments the answer is `n=6`. The related quantities the
question opens — the rainbow dual (largest all-distinct subcube), the cut-orientation mixing condition, the
transversal codes that exist only at `n<=4`, and whether this vertex-cut aligns with the project's
edge-cut/cycle split — are all asking the same thing from different sides: how does a symmetry group let
you pack its orbits into a cube, and where does the packing break. The cut is the answer until it isn't, and
where it isn't is worth as much as where it is.
