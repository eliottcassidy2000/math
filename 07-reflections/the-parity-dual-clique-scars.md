# The parity-dual clique scars: where tournaments and even graphs disagree at K_3

*kind-pasteur-2026-06-22-S31h. Chasing the owner's "does the even-graph spectrum have a unique
forbidden class mirroring K_3?" with the recursion lens. Companion to
`seven-is-the-unique-forbidden-clique-value.md`.*

## The shared spine

Both spectra are organized by the clique values `I(K_r, 2) = 2r+1` (a clique has independent
sets only of sizes 0 and 1). The clique spine is `3, 5, 7, 9, 11, 13, 15, ...`. Each world is a
selection rule on this spine, dictated by what its objects ARE.

## The two selection rules (computed, <=7 vertices)

**Tournament H-spectrum** = `{ I(Omega, 2) : Omega a REALIZABLE conflict graph }`.
Every clique `K_r` is realizable as a conflict graph EXCEPT `K_3` (THM-200: three pairwise-
conflicting triangles force a `C_5`). So the tournament spectrum is **the spine minus the single
odd clique `K_3`**: `3,5,9,11,13,15,...` -- the **unique permanent gap at 7**.

**Even-graph I(G,2)-spectrum** = `{ I(G, 2) : G an EVEN graph (all degrees even) }`.
A clique `K_r` is an even graph iff `r` is ODD (degree `r-1` even). So the even-graph spectrum
realizes ONLY the **odd-clique** values `7, 11, 15, ...` and MISSES the **even-clique** values
`5 = I(K_2,2), 9 = I(K_4,2), 13 = I(K_6,2)`. Computed gaps `<=45`: `{5, 9, 13, ...}`.

## The duality

> The tournament forbids the **one odd clique `K_3`** -> a UNIQUE gap at `7`.
> The even graph forbids **all even cliques `K_2, K_4, K_6, ...`** -> a PARITY FAMILY of gaps at
> `5, 9, 13, ...`. And it **HEALS the tournament's 7-scar**: `K_3` (the triangle) is a perfectly
> valid even graph, so `7` IS in the even-graph spectrum.

So `7` is exactly the value where the two object-constraints DISAGREE: tournament-*realizability*
kills it, even-graph-*parity* keeps it. The apex prime is the unique point where orientation and
parity diverge on the smallest odd clique.

## Why this is the deep answer to "the even-graph analogue"

The owner asked whether `E_n` has a unique forbidden class mirroring `K_3`. The answer is sharper
than yes/no: **the even-graph spectrum is the PARITY-DUAL of the tournament spectrum on the clique
spine.** The tournament's scar is unique (one odd clique) because forbiddenness is a *realizability*
property of orientations; the even-graph's scars are a *parity family* (every even clique) because
they are a *membership* property of un-oriented even-degree graphs. Forgetting orientation does not
reproduce the 7-scar -- it heals it and opens the complementary even-clique holes.

This sits exactly on the `K_3 <-> C_5` SPGT duality (mac-mini S37 verified the LITERAL `C_5 = H=7
K_3` match via the cycle space, HYP-2880): the odd clique `K_3` is the tournament-side imperfection,
the odd hole `C_5` is the even-graph-side imperfection, and `n=7` is where BOTH first bite. The
clique spine carries the whole story: tournaments puncture it once (at the forbidden odd clique),
even graphs puncture it by parity (at every even clique), and the Lonely Runner inherits the scar
through `14 = 2*7`. -> HYP-2878, HYP-2879, HYP-2880, THM-200,
`seven-is-the-unique-forbidden-clique-value.md`.
