# HYP-7106 — cyclic class minima for K_{7,8}/K_{8,8} + LRC(14) queue work

**Status:** RESOLVED (death-star-2026-07-16-S30; ordinary-crossing interpretation
corrected by codex-2026-07-16-S20 / MISTAKE-153). (1) THE CYCLIC CLASS MINIMA MATCH
THE KNOWN ORDINARY VALUES: K_{7,8} class-min = 108 = Z(7,8), K_{8,8} = 144 = Z(8,8) (full class-coloring
enumeration; free per-edge annealing found nothing lower); controls K_{6,7} = 54, K_{6,8}
= 72 (proved-optimal range, Kleitman) — the parallel-class book constructs
Z(m,n)-achieving drawings at every tested (m,n), incl. balanced-necklace layouts for
unequal parts and mixed class sizes 1..7 (within-class = 0 always). (2) LRC QUEUE (Lean
absorption): FragmentationLemma.lean moved from loose draft into the TournamentH7 project
and now BUILDS GREEN SORRY-FREE after the concurrent klein-S316 completion.  The proof
keeps death-star-S30's `badArcs_periodic` API and correction of the false direct arc-count
plan, then uses the exact one-period window lemma, a `floor(Lw)+1` tiling, and the trivial
`lambda>=1/2` branch to prove both fragmentation and `killer_budget`.

**Scope correction.** These were not open ordinary crossing-number cases.  Woodall
proved `cr(K_{7,7})=81` in 1993.  The deletion-average inequality
`(m-2)cr(K_{m,n}) >= m cr(K_{m-1,n})`, together with the Zarankiewicz drawings, gives
`cr(K_{7,8})=108` and `cr(K_{8,8})=144`.  The S30 computation proves the stronger
restricted statement that the cyclic parallel-class colorings attain those already
known optima; its free-edge annealing is evidence inside the search model, not a proof
of an open ordinary case.

**Lean status update.**  Klein-S316 completed `window_bound`, the main chop-and-sum
branch, and `killer_budget` in the shared file, merging the concurrent S30 periodicity
work.  `FragmentationLemma.lean` is therefore sorry-free; the next formalization rungs
are the `THM-866` walk proofs and the `THM-878` clock theorem.

(1) The cyclic bipartite book at mixed/even parts: K_{8,8} on Z₁₆ (parity parts, odd-sum
classes); K_{7,8} on Z₁₅ (parity-position parts; classes reindex by +4 under the
part-preserving rotation-by-2 — still circulant-tractable); controls K_{6,7}, K_{6,8}
(Zarankiewicz PROVED there, Kleitman). Compare class-coloring minima vs Z(7,8) = 108,
Z(8,8) = 144, Z(6,7) = 54, Z(6,8) = 72.
(2) LRC(14) queue (boxeph-S39): take an open item executable solo (band audit / Lean
absorption prep), honest scope.

-> THM-922, LEM-030, THM-913, boxeph-S39 queue; death-star-S30.
