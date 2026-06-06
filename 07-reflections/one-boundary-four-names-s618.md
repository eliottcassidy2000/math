# One boundary, four names (S618)

The user handed me four phrases and a claim that they are siblings: Helly is "how many orders of overlap you must
keep," the Vitali wall is "finite moments cannot decide p₀," Collatz's two-block is "the correlated residue where
density is blind," and OCF/partition functions are "the sibling world where independence polynomials already play
the depth-distribution role." I spent the session checking whether they are four facts or one, and they are one.

The thing they are about is the overlap hierarchy `S_k` — the inclusion-exclusion terms of the forbidden arcs. The
lonely measure is the full alternating sum `p₀ = Σ(−1)^k S_k`; I verified it closes exactly. Once you have that
identity, the four phrases are four things you can ask about *the same series*.

Helly asks: where can I truncate? The Bonferroni partial sums bracket `p₀`, and the order-3 lower bound `1 − S₁ +
S₂ − S₃` is a genuine certificate — when it's positive, loneliness is *proved*. For arcs on a circle that's the
Helly number, 3, and the data obeys it: the generic configs are decided at order 3. What I didn't expect was the
honesty the computer forced on me — `{1,5,8,11,13}` is lonely (`p₀ ≈ 0.09`) but order 3 says `−0.03`. So three is
not a wall, it's a *floor*: order 3 decides everything except a thin layer hugging the collapse set. The exceptions
are exactly the configs creeping toward resonance.

And that thin layer is the Vitali wall. At an actual collapse config the truncations refuse to close — `+1, −.6,
+.281, −.057, +.000` — reaching zero only at the very last term. No finite number of moments can tell you `p₀ = 0`;
worse, even knowing `p₀ = 0` exactly, the measure can't tell you whether the lonely set is a single tight point (LRC
holds) or empty (LRC fails). That is the whole reason LRC is not a measure-theory problem at the boundary: measure
is a continuous functional and the answer is discontinuous there. You have to *construct* the witness. Helly's
truncation and Vitali's non-closure are the same series read at the generic point and at the boundary.

Collatz is that boundary wearing a multiplicative mask. "Almost all orbits" sees a density-1 set; the cycles and the
2-adic two-block are density zero and density cannot see them. That is the Vitali wall again — the structured
residue invisible to the first moment — transported across the additive↔multiplicative dictionary I built two
sessions ago. The blind spot of density on the multiplicative side and the blind spot of measure on the additive
side are literally the same blind spot.

The fourth lens is the one that tells you what tool to reach for. `P(z) = ∫ z^{depth}` is a partition function;
`p₀ = P(0)`; the independence polynomial of the hard-core gas is its combinatorial twin. I tried the crude version —
the pairwise dependency graph — and it collapsed to the complete graph, because at the LRC gap the arcs are so wide
that every pair is correlated. So the pairwise independence polynomial is as vacuous as the first moment, for the
same reason. The lesson is not that the sibling world is wrong; it is that the only honest object on either side is
the *full* one — the whole depth generating function, the whole overlap hierarchy. Which is exactly what Helly and
Vitali already said. Four lenses, one series, one boundary: the collapse family, which is the additive chains, which
are the resonances, which empty the apex sheaf's `H⁰`. I formalized the partition function and its factorization and
its value at zero, and the small surprise is how little the formal core needs to know — one product identity carries
all four readings.
