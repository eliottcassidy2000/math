# Self-complementary classes cannot be compressed

*mac-mini-2026-07-01-S85. Reflection on HYP-3810.*

The owner asked a sharp question: does the T-join parity — the all-odd-degree blue subgraph on the self-
complementary nodes — obstruct low-dimensional covers of those classes? The answer came back clean and a
little surprising: yes, and by a lot. To realize every self-complementary tournament class inside an axis-
aligned subcube of the arc-hypercube, you need `floor((n-1)^2/4)` free coordinates — the entire dimension of
the complement fold — while information theory says `ceil(log2 of the count)` bits would suffice. At six
vertices there are twelve self-complementary classes; four bits could name them; but no four-dimensional
subcube contains all twelve, nor five — you need six. The gap grows like `n^2/4`. The self-complementary
classes are, in the strongest sense, incompressible.

The reason is a parity, and it is the same parity that has been surfacing all week. Every self-complementary
class holds an odd number of grid-symmetric tilings — the fixed points of the complement mirror. Odd, always
odd. In the blue subgraph that is odd degree at every self-complementary node, which is the definition of a
T-join with the whole vertex set as its boundary. And an odd count is exactly what a low-dimensional subcube
cannot promise to catch. The grid-symmetric tilings live in the fold, a subspace cut out by the mirror's
pairing constraints, and they are the essential representatives of the self-complementary classes; strip away
one dimension of the fold and some class, carrying its odd handful of fixed tilings, falls out of reach. The
upper bound is easy — the fold itself is a cover, the grid-symmetric subcube realizes exactly the self-
complementary classes and nothing else. The content is the lower bound: you cannot do better, and the
obstruction is the odd parity of the fixed points.

There is a moral here about what self-complementarity costs. These classes are defined by a symmetry — they
are their own complements. One expects a symmetry to *help*: fewer degrees of freedom, a smaller
fundamental domain, a cheaper description. And it does help in one direction — the fold halves the tiling
count, the invariants compute twice as fast. But in the covering direction the same symmetry *hurts*: the
self-complementary classes are pinned to the fold, they fill it, and they inherit its full quadratic
dimension. The symmetry that makes them special is the symmetry that makes them incompressible. You cannot
name them with fewer bits than the fold has cells, because to be self-complementary is to live everywhere in
the fold at once.

And the fold keeps its schedule. The blue T-join is bipartite through five vertices and turns non-bipartite
at six — the same six where the black sea first connects to itself, where the pure-black nodes grow self-
loops, where the minimal-flip gauge over-identifies. Four faces of one threshold, seen from four sessions.
The T-join gaining an odd cycle, the sea gaining an edge to itself, the covering gauge gaining a redundancy —
these are not four coincidences at `n=6`. They are the single moment when the tiling hypercube becomes large
enough that its structure stops being a small-case accident and starts being generic. The bipartite blue
graph was a courtesy of smallness; at six the courtesy ends.

The pattern that transcends the theorem: **a defining symmetry is a double agent — it compresses the object
in one coordinate and forbids compression in another, and which face you meet depends on the direction you
push.** Fold the triangle and you save half your work; try to *cover* the folded classes and the fold's full
dimension is the floor you cannot get under. The information-theoretic bound counts what you would need to
*distinguish* the classes; the parity obstruction counts what you need to *contain* them; and for a self-
symmetric family these two numbers diverge quadratically. When a set is closed under a symmetry, ask not only
what the symmetry lets you forget but what it forces you to keep — the answer to the second is usually the
whole invariant subspace, and that is exactly as large as the symmetry is deep.
