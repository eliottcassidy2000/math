# Everything is a coloring (S633)

The instruction was to see everything as a graph coloring — some nodes, some edges, some both — and
to reframe as much as I could. I expected a loose analogy session. Instead the colorings kept turning
out to be *the same object* I had been computing under other names, and by the end the word "coloring"
had merged with the word "partition function" the cluster has been circling for weeks.

Start with the cleanest one, because it is exact. Last session the tight Lonely Runner tie-graph turned
out to be the cycle `C_n`. A sieve is a proper coloring of that cycle, and a corrector is one color
class — an independent set. So the sieve arity is `χ(C_n)` and the corrector size is `α(C_n)`, and these
are not metaphors: `χ(C_n)` is `2` for even `n` and `3` for odd `n`, `α(C_n) = ⌊n/2⌋`, and the
chromatic polynomial `P(C_n,2)` is `0` exactly when `n` is odd — the odd cycle that refuses to
two-color. The single-versus-multi-sieve dichotomy the cluster spent a dozen sessions on is the
statement that an odd cycle is not bipartite. The 2-adic seam is `χ` of a cycle. I have rarely seen a
year of structure collapse into a sentence that cleanly.

The chromatic polynomial is where coloring stops being a lens and becomes the master object. `P(G,k)`
counts proper `k`-colorings, and it is exactly the zero-temperature Potts partition function. The
covering-depth `Z(z)` I built two sessions ago is the same animal at a different temperature. So
"partition functions everywhere" and "colorings everywhere" are one claim: every problem carries a
conflict graph, and its chromatic / Potts partition function is the thing whose special values are the
answers. Loneliness is a free color in the ground state; the chromatic number is the smallest palette
that avoids every conflict; an independent color class is a corrector or a lonely set. The deformation
parameter that was fugacity is now temperature; the questions are unchanged.

Then the edge side, which gave the surprise of the session. The pair-sum sieve — the reason the modulus
`2n−1` runs through everything — labels each pair of runners by `v_i + v_j` mod `2n−1`. I checked
whether two pairs sharing a runner ever get the same label, and they never do: it is a *proper edge
coloring* of the complete graph. Which is to say it is the round-robin tournament schedule, the
circle-method one-factorization that has organized chess tournaments for a century. The pair-sum sieve
was a scheduling algorithm in disguise. The `2n−1` modulus is the number of rounds. That is a genuinely
new way to hold THM-401, and it hands the whole sieve over to the theory of edge colorings and Latin
squares — and, on the engineering side, to scheduling and code design, where one-factorizations are
bread and butter.

And then unit distance, which needed no reframing at all, because it was already the most famous coloring
problem there is: the chromatic number of the plane, Hadwiger–Nelson, somewhere between five and seven.
The grid disproof, read this way, is a construction that forces chromatic complexity by packing
unit-distance edges — and the tie-graph of the optimum is a Cayley graph of a complex-multiplication
lattice, whose chromatic number is governed by its symmetry. The same sentence as the Lonely Runner: the
conflict graph is a Cayley graph, and its symmetry sets its coloring invariants. One problem on the
discrete circle `ℤ/m`, the other on the plane; the same family of graphs, one dimension apart.

What I take from the session is the shape, finally explicit. Every problem here is: lay down a *conflict
graph* whose edges are the resonances — binding pairs for runners, unit distances for points, equal
endpoint-types for the metagraph — and *color around the conflicts*. The kind of coloring sorts the
problem: color the nodes and you are sieving or asking Hadwiger–Nelson; color the edges and you are
scheduling or running the pair-sum sieve; color both and you are in the trienement. The chromatic
polynomial of that graph is the partition function; its chromatic number is the sieve arity; its
independence number is the corrector; its symmetry group is the perspective key. Trienement, symmetry,
partition function, coloring — I had been calling one structure four names. The session's real result is
that I can now say which name to use for which question, and that the answer to "how hard" is always a
chromatic invariant of a Cayley graph.
