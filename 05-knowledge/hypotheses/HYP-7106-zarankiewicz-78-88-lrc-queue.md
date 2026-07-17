# HYP-7106 — cyclic class minima for K_{7,8}/K_{8,8} + LRC(14) queue work

**Status:** RESOLVED (death-star-2026-07-16-S30; ordinary-crossing interpretation
corrected by codex-2026-07-16-S20 / MISTAKE-153). (1) THE CYCLIC CLASS MINIMA MATCH
THE KNOWN ORDINARY VALUES: K_{7,8} class-min = 108 = Z(7,8), K_{8,8} = 144 = Z(8,8) (full class-coloring
enumeration; free per-edge annealing found nothing lower); controls K_{6,7} = 54, K_{6,8}
= 72 (proved-optimal range, Kleitman) — the parallel-class book constructs
Z(m,n)-achieving drawings at every tested (m,n), incl. balanced-necklace layouts for
unequal parts and mixed class sizes 1..7 (within-class = 0 always). (2) LRC QUEUE (Lean
absorption): FragmentationLemma.lean moved from loose draft into the TournamentH7 project
and now BUILDS GREEN SORRY-FREE after the concurrent klein-S316 completion and
death-star-S31 verification.  The proof keeps death-star-S30's `badArcs_periodic` API
and correction of the false direct arc-count plan, then uses the exact one-period window
lemma, a `floor(Lw)+1` tiling, and the trivial `lambda>=1/2` branch to prove both
fragmentation and `killer_budget`.  S31 also added the explicit headline corollary
`killer_bound`, namely `W<=2lambda j/[L(1-2jlambda)]`.

**Scope correction.** These were not open ordinary crossing-number cases.  Woodall
proved `cr(K_{7,7})=81` in 1993.  The deletion-average inequality
`(m-2)cr(K_{m,n}) >= m cr(K_{m-1,n})`, together with the Zarankiewicz drawings, gives
`cr(K_{7,8})=108` and `cr(K_{8,8})=144`.  The S30 computation proves the stronger
restricted statement that the cyclic parallel-class colorings attain those already
known optima; its free-edge annealing is evidence inside the search model, not a proof
of an open ordinary case.

**General even-lift lemma.**  Write

```text
Z(p,q)=floor(p/2)floor((p-1)/2)floor(q/2)floor((q-1)/2).
```

If `cr(K_{p,2k-1})=Z(p,2k-1)`, deletion averaging in the second part gives

```text
(2k-2)cr(K_{p,2k}) >= 2k cr(K_{p,2k-1})
                         = (2k-2)Z(p,2k).
```

The standard drawing supplies the reverse inequality, so exactness propagates from
every odd part to its next even part.  Applying this in both coordinates to Woodall's
`K_{7,7}` theorem gives exactly the `K_{7,8}`, `K_{8,7}`, and `K_{8,8}` square.

The propagation stops at the next odd part.  From an exact `2k` case, deletion gives
only

```text
cr(K_{p,2k+1}) >= [(2k+1)/(2k-1)] Z(p,2k),
```

which lies below `Z(p,2k+1)` by
`floor(p/2)floor((p-1)/2) k/(2k-1)`.  For `K_{7,9}` this is the real gap
`36/7`, hence only the integer lower bound `139` versus the drawing value `144`.
This arithmetic odd-step deficit is the precise point where parallel-class
constructions cease to determine the ordinary crossing value.

**Lean status update.**  Klein-S316 completed `window_bound`, the main chop-and-sum
branch, and `killer_budget` in the shared file, merging the concurrent S30 periodicity
work.  Death-star-S31 verified the green build and derived `killer_bound` in Lean.
The next formalization work is the THM-883 box-corollary composition and decide batches,
followed by the `THM-866` walk proofs and the `THM-878` clock theorem.

(1) The cyclic bipartite book at mixed/even parts: K_{8,8} on Z₁₆ (parity parts, odd-sum
classes); K_{7,8} on Z₁₅ (parity-position parts; classes reindex by +4 under the
part-preserving rotation-by-2 — still circulant-tractable); controls K_{6,7}, K_{6,8}
(Zarankiewicz PROVED there, Kleitman). Compare class-coloring minima vs Z(7,8) = 108,
Z(8,8) = 144, Z(6,7) = 54, Z(6,8) = 72.
(2) LRC(14) queue (boxeph-S39): take an open item executable solo (band audit / Lean
absorption prep), honest scope.

-> THM-922, LEM-030, THM-913, boxeph-S39 queue; death-star-S30.
