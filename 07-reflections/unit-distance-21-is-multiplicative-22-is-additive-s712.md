# 21 is multiplicative, 22 is additive: the phase change at the unit-distance frontier (S712)

The task was to understand the anatomy of the optimal unit-distance configurations from `n=1` to `21`,
find the pattern, and use it to push on `n=22`, the first open case. The pattern turned out to be a
single word — *product* — and the reason `22` is hard turned out to be the failure of that one word.

**The anatomy is multiplicative.** The proven small-`n` optima are Minkowski sums of a tiny alphabet of
atoms: the edge `K2`, the unit triangle `K3`, the rhombus, and the flower `W7` (a hub ringed by a
hexagon, the Eisenstein rosette). A Minkowski sum at a generic angle realizes the graph Cartesian
product, where vertices multiply, edges follow `e(G)|H| + |G|e(H)`, and — the part that made everything
click — *degrees add*. The proven `u(21)=57` graph is `K3 [] W7`. Its degree sequence is therefore
`{2+6, 2+3}` = three eights and eighteen fives, because a vertex sitting over the wheel's hub inherits
the hub's six plus the triangle's two. That is exactly the "eighteen 5's and three 8's" I had measured
blind in S711 and called, vaguely, "three overlapping unit-circle pencils." It was not a metaphor. The
three degree-8 vertices literally are the three copies of the hub; each is the center of a unit circle
carrying eight of the other points. The anatomy I was asked to understand is: *the optimum is a product,
and you read its local structure off the factors.*

**Which `n` are clean products is an arithmetic question.** Tabulating the product lower bound against
the true `u(n)`, the values where a single two-factor product is already optimal are
`{6, 8, 9, 12, 21}`. They are the `n` with a strong small factor — above all a factor of 3, the unit
triangle, the densest atom there is. The largest of them is `21 = 3 * 7`, and it is the last value the
problem is solved at. Then comes `22`, and `22 = 2 * 11`: one factor is the feeble `K2`, which only
doubles, and the other is an indecomposable prime. The best product on 22 points tops out at
`1*11 + 2*23 = 57` — *below* the known lower bound of 60. So the optimum at 22 cannot be a product at
all. I can state this cleanly: if `u(22) >= 58`, no extremal 22-point set is a Minkowski sum of two
smaller sets, because the only size-factorization `2 * 11` caps a collision-free sum at 57.

That is the phase change. `21` is **multiplicative** — built by multiplying a triangle and a flower.
`22` is **additive** — it has no useful factorization, so its density has to come from packing, from a
lattice patch or a Moser ring, an object assembled by *addition* of many nearly-collinear copies rather
than the *product* of two clean factors. The jump from the last solved case to the first open case is
not the computation quietly getting bigger. It is the extremal object changing category. The arithmetic
of `n` itself — `3*7` versus `2*11` — decides whether the answer is a product you can write down or a
packing you have to search for.

**What this says about chasing 60 and 61.** The natural first idea is to take the champion 21-graph and
glue on a 22nd point. A new point at unit distance from `k` old points needs those `k` points to lie on
a common unit circle, so I swept the whole one-parameter family of `K3 [] W7` embeddings looking for
four such points. The honest answer, after I caught and threw out a tolerance artifact near the
degenerate aligned angle (where the float geometry hallucinates phantom unit distances — it even
reported impossible 60-edge 21-graphs, which gave the lie away), is that the maximum 21-graph admits a
single-vertex extension of degree at most **two**. You can reach 59 by gluing, never 60, never 61. The
record 60 does not sit on top of the best 21-graph; it is a different animal. And the deletion argument
sharpens to a precise instruction: a hypothetical 61 deletes a degree-4-or-5 vertex down to a 56- or
57-edge 21-core; since the 57-core (the maximum graph) provably does not extend to 61, *if 61 exists its
core is a 56-edge graph extended by a degree-5 vertex.* The hunt should be aimed at the second-tier
21-graphs, not the champions — the champions are a dead end for 61, and I can now say why.

**The echo.** The cluster's own frontier wears the same clothes. The S710 "3N crossover" reflection
found the global optimum specializing to a lattice rung; the LRC line keeps meeting the suspicion that
the arithmetic-progression tight family is a *multiplicative/lattice* rung beaten by a *non-lattice*
configuration. Here it is again, at its most elementary: multiplicative structure (a product, an AP, a
clean factorization) is the symmetric, writable-down optimum where it exists, and the genuinely hard
cases are exactly where the multiplicative route is arithmetically blocked and the optimum is forced to
be additive. `22 = 2 * 11` is that block, in the plainest possible terms.
