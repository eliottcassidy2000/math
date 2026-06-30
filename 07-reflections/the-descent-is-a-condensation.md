# The descent is a condensation

*mac-mini-2026-06-30-S35. The owner asked to study more irreducible paradoxes, work the descent, and come to know the finite families. The three turn out to be one thing seen at two levels.*

## What an irreducible paradox is

If a tournament is the intransitivity among `n` things, an *irreducible* paradox is an intransitivity you cannot break apart — no way to split the things into a top group that beats a bottom group. That is exactly strong connectivity: from every thing you can reach every other through the dominance relation. Two things are never irreducible (one beats the other), so the smallest irreducible paradox is the 3-cycle, the Condorcet atom; after that they are `1, 6, 35, 353, …` up to relabeling. Few, and classifiable. These are the tournaments where the odd-cycle structure is genuinely tangled — where `H > 1`, where there is a Hamiltonian cycle, where no ranking gets any traction at all.

## Every tournament is an order of paradoxes

The thing I had not appreciated until I computed it: a tournament is not "orderable or not." It is *always* a total order — of its irreducible blocks. The condensation collapses each strongly-connected component to a point, and the points form a transitive tournament, a clean ranking. So every preference among `n` things decomposes, uniquely, into a ranking of blocks plus one irreducible paradox inside each block. Order on the outside, irreducible intransitivity on the inside. The fully orderable tournament is the one all of whose blocks are singletons (zero paradox, the unique transitive class); the fully tangled one is a single block (the SC types, 6 of the 12 at `n=5`, 35 of the 56 at `n=6`).

That decomposition is a *descent*. You peel the orderable skeleton — the ranking of blocks — and what is left, what carries all the content, is the finite family of irreducible cores.

## The same move, one level down

And that is precisely what the 2-adic descent does to the Lonely Runner. THM-580 peels the even part of a speed set — divides it by two, the doubling, the orderable 2-adic skeleton — and exposes the odd cores, which mod the apex prime are subsets of `Z_7`. klein-S17 just showed those cores realize all 127 nonempty subsets of `Z_7`: the apex finite family is complete. So the two descents line up exactly. The condensation peels the transitive ranking and exposes the strongly-connected blocks; the 2-adic descent peels the even doubling and exposes the odd `Z_7`-cores. On both sides, a large object — a whole tournament, an infinite covering family — finitizes into a finite family of irreducible paradoxes.

This is why the proof can live there at all. klein-S16 made the point sharp: an infimum is only provably attained over a *finite* family; over the infinite covering family the lonely measure sinks to zero. The descent is what manufactures the finite family. It does it by throwing away the orderable part — which was never the obstruction, because an orderable relation is a coboundary, a ranking, and rankings carry no paradox. Everything that can obstruct lives in the irreducible cores, and there are finitely many, so a minimum exists and is attained at the smallest one: the 3-cycle on the tournament side, the doublet `C_7` with gap `4cos²(3π/7)` on the apex side (THM-590).

## Coming to know the finite families

So "the finite families" are the irreducible-paradox atoms, and now they have two faces. The tournament face is the strongly-connected tournaments — `0, 0, 1, 1, 6, 35` — the catalogue of ways `n` things can be irreducibly tangled. The apex face is the 127 `Z_7`-cores, the catalogue of irreducible resonances the runners can descend to. The descent is the bridge: it takes an unbounded problem and hands you one of finitely many irreducible paradoxes, each of which you can simply check. The smallest atom binds, and it does not vanish.

What I would chase next is the matching: which strongly-connected atom is the tournament-side image of the binding doublet `C_7`, and whether the descent provably lands the worst covering on the minimal-gap core. If it does, the whole infinite problem reduces to inspecting one irreducible paradox — the smallest one — and confirming, finitely, that it is intransitive enough.

See [[tournaments-are-intransitivity-among-n-things]] (HYP-3599), klein-S17 (HYP-3598, the `Z_7`-core finite family), klein-S16 (HYP-3597, finite vs infinite), THM-580 (the 2-adic descent), THM-590 (the apex gap). New: HYP-3603.
