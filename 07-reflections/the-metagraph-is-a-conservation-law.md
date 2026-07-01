# The merged metagraph is a conservation law

*klein-2026-07-01-S75. A reflection on HYP-3808 — formalizing the owner's blue/black line-pairing process,
proving its checksum, and refuting one of its conjectured rules.*

The owner described the merged metagraph as a *process*: start with every node holding zero tilings, then
assign lines — each blue or black line raising the tiling count of two nodes by one (or one node by two, a
self-loop) — until every node reaches its correct count, subject to eligibility rules between three kinds of
node. Asked to make those rules "as precise and explanatory as possible," I computed the whole structure for
n up to 6, and it resolved into something cleaner than a list of rules: a system of conservation laws.

The nodes carry a conserved bit. Every self-complementary merged node has an *odd* tiling count; every
non-self-complementary one has an *even* count. This is not a coincidence to be catalogued — it is Rédei's
theorem in disguise. Each unmerged isomorphism class holds an odd number of tilings, because tilings are
Hamiltonian paths wearing the base path as a costume, and Rédei says a tournament has an odd number of
Hamiltonian paths. A self-complementary node is one class, hence odd; a non-self-complementary node is two
complement classes of equal size, hence even. So the parity of a node's tiling count is a Z/2 *charge*, and
Rédei is the law that assigns it. The owner's three categories are, at bottom, two charges and an interface.

The lines respect a color law that makes the structure tripartite. Blue lines — the grid-symmetric ones —
touch only self-complementary nodes; black lines touch only the non-self-complementary side; and the mixed
nodes are the sole place both colors meet. Pure blue is a blue-only world, pure black a black-only world,
mixed the membrane between them. That immediately yields a second conservation law by the oldest argument in
graph theory: every self-complementary node has odd blue-degree, blue lines pair them up, so there must be
an even number of them. The count of self-complementary merged nodes is even, and it is even *because* the
blue handshake demands it. Two conservation laws — a per-node parity from Rédei, a global parity from the
handshake — plus an exact count of blue and black lines from the grid-symmetry formula, and the "process"
is no longer a mystery to be simulated but a degree sequence to be realized under gauge constraints.

And then the part I most needed to get right: one of the conjectured rules is false, and false in an
informative way. The owner conjectured that only mixed nodes can carry self-loops — the self-contained
`+2` moves. For n = 4 and n = 5 this is exactly true; every self-loop sits on a mixed node. At n = 6 it
collapses. The pure-black nodes begin to self-loop, and not timidly — twenty-four of the twenty-six
self-loop lines are theirs, the mixed nodes keeping only two. Pure blue never self-loops at all. So the true
rule is "self-loops live on mixed and pure-black, never on pure-blue," and the pure-black self-loops switch
on at n = 6. I have seen that number before this week: the flip-rank's efficient balanced-cut shape also
dies at n = 6, and the packing law breaks there too. Something about six vertices — the first size where the
two triangles are symmetric enough, where the group folds the cube too tightly — keeps being the threshold.
A self-loop here is a tiling whose all-tiles-flip is isomorphic to itself or its complement, a flip
symmetry; those symmetries are rare until the tournaments grow structured enough to admit them, and six is
where they arrive on the black side.

The lesson is one about how to honor a conjecture. The owner asked for precise rules and offered a
conjecture as scaffolding. The respectful thing was not to confirm it but to test it to the point where it
breaks and report exactly where and how — n <= 5 true, n = 6 false, pure-blue never, pure-black dominant —
because the break is the discovery. The conjecture pointed at the right object (self-loops as a distinct
move) and the right suspicion (that eligibility differs by category); it was wrong only about which
category, and being wrong at n = 6 rather than always is itself a signal, tying this structure to the same
transition the flip-rank showed. A process becomes understood not when its rules are listed but when its
conserved quantities are named and its rules are pushed until they fail; here the conserved quantities are
two parities, the failure is at six, and the parity that runs through all of it is Rédei's.
