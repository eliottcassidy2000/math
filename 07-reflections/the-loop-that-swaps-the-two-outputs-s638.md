# The loop that swaps the two outputs

*S638 reflection. On reading "a loop of the input causes a swap of the two outputs" as the
self-converseness of the universal tournament, and finding that the smallest avatar orients its arcs by
the cube root of unity.*

The user handed me two fragments. "Consider the Rado graph as a tournament." And, separately, "a loop
of the input causes a swap of the two outputs." I spent the session convinced they were one sentence,
and I think they are.

The Rado graph is the generic countable graph — flip a coin for every pair and you almost surely build
it; it contains everything, and it is so symmetric that any local picture extends to a global one. The
tournament version is the same idea with arrows: orient every pair by a coin and you almost surely get
*the* countable universal homogeneous tournament. It is the place where every finite tournament lives
as a sub-pattern, and where the automorphism group is so large that its orbits on `k`-subsets are
exactly the isomorphism types of `k`-vertex tournaments. That last fact is the one I have been circling
for twenty sessions under the name "perspectives." On finite tournaments the perspective count and the
iso-type count *almost* agree and then break apart at `n = 5`; I always read the break as the content.
But on the generic object they agree exactly, forever — because agreeing exactly is the *definition* of
the generic object (Ryll–Nardzewski). The finite breakage was the shadow of finiteness, not of the
phenomenon. The phenomenon, cleanly stated, lives upstairs.

Then the second fragment. Every vertex has two outputs — the set it beats and the set it loses to, its
out-set and its in-set. The operation that swaps those two sets at every vertex at once is reversal:
turn every arrow around. And the universal tournament is *self-converse* — reverse all its arrows and
you get an isomorphic tournament, because reversal is a symmetry of the class of finite tournaments and
the limit inherits it. So "a loop of the input causes a swap of the two outputs" is the
self-converseness of the Rado tournament, and the "loop" is the order-two involution that has been the
spine of this whole arc: `σ : v ↦ −v`, the antipodal map, complex conjugation, the complement, the deck
transformation of the antipodal double cover. Apply it once and the two outputs swap; apply it twice
and you are home. A loop.

What made the session worth formalizing was the avatar. I wanted the smallest concrete tournament that
is vertex-transitive and self-converse — the finite stand-in for the generic one — and it is the Paley
tournament on the residues mod a prime `q ≡ 3 (mod 4)`. The smallest is `q = 7`. And `7` is the number
this entire arc keeps returning to: `7 = Φ₃(2)`, the forbidden value; `7 = N(3+ω)`, the Eisenstein
prime and the plane's chromatic bound (last session); `7 = M₃`. So I computed the connection set — the
differences that count as a win — expecting nothing in particular, and it was `{1, 2, 4}`. The nonzero
squares mod 7. Which are *exactly the cube roots of unity*. The Paley-7 tournament orients an arc
precisely when the difference between the two vertices is a cube root of unity. The `ω` that builds the
hexagon, that splits the 7, that is the eigenvalue of the 3-cycle, that the whole arc converges on —
that very `ω` is the rule that says who beats whom in the smallest random-tournament avatar.

And the two fragments meet there, mechanically. The thing that makes Paley-7 a *tournament* and not a
graph is that `−1` is a non-residue mod 7 (because `7 ≡ 3 (mod 4)`), so the residues and their negatives
are disjoint — every nonzero difference is a win in exactly one direction. That same `−1` is the swap:
`σ : x ↦ −x` carries each arc to its reverse, because negating a difference moves it from the residues
to the non-residues. So the condition that *creates* the tournament (`−1` non-residue) and the loop that
*swaps its outputs* (`x ↦ −x`) are the same fact about `−1`. The user's two sentences are one sentence,
and the prime where it is cleanest to say is `7`, and the arcs are the cube roots of unity.

I formalized the avatar — six sorry-free lemmas: the arc set is `μ₃`, `−1` is a non-residue, the three
tournament axioms, the antipode reverses every arc, the antipode is an involution. The infinite
statement — that reversal preserves the extension property so the universal tournament is self-converse
— I left as prose and a handoff; it needs the Fraïssé machinery, and the finite avatar is the
down-payment. But the shape is settled: the generic tournament is where perspectives become exact and
the swap becomes a symmetry, and its smallest face wears the cube root of unity for arcs.
