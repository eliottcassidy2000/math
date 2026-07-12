# Two cases, a constant, and a Freiman lead — the combinatorial simplifications the endgame was hiding

*klein-2026-07-12-S265. Owner directive: "spend another similar session searching for combinatorial
connections and simplifications." Four Explore agents swept the ≤6-coprime crux, the finish-map
structure, the tournament/tiling bridges, and the additive-combinatorics toolbox. The theme that
came back was not a new weapon but a smaller problem: the endgame is carrying more cases and heavier
machinery than the mathematics needs.*

---

## Four simplifications, honesty-graded

**1. The case split is TWO cases, not five (rigorous).** The finish map reads as
[non-covering / covering-DC / bounded-diameter / large-diameter / AP-wall]. Two of those are not
intrinsic. The **AP-wall collapses into non-covering**: the tight locus `M=1/14` is exactly
`{AP, GW, dilates}` (THM-612 + klein-S206, exhaustive `n=4..7`), and it is *entirely non-covering*
(`{1..13}` has no multiple of 14). So `#(tight ∩ covering) = 0` — the `1/14` wall lives wholly in
the branch the `t=1/q` sieve already dispatches, and every divisor-complete family sits in a strict
cushion `M > 1/14` (band-edge, opus-S235). The bounded/large-diameter split is a *technique* seam
(finite check below, multi-scale above), not a case. Honest minimal split: **non-covering** (holds
all tight families) **+ divisor-complete** (strict `>1/14`).

**2. Both routes are one moment ladder (rigorous reduction).** `bandCount` (Route B, clearing) is
the same empty-count `N` as the seven-sector residue (Route A, density); `B5 = Σ(−1)^d S_d` is the
alternating factorial-moment majorant both build, and both live on pair-sum rulers (THM-671/668/707).
The much-cited "coverage-clearing duality" as a *quantitative identity* is heuristic (its own
correlation is +0.398); what is rigorous is the moment-ladder identity plus one shared, still-open AP
inverse theorem. At most one route's bookkeeping needs to reach Lean.

**3. The target is a CONSTANT, and I had been chasing a mirage.** This is the one that cost me. Last
session (S264) I built a wider-band Parseval floor and reported, with verification, that it "reaches
the true M, growing with diameter." The growth was real in my sample — and the sample was the
artifact. death-star's THM-721 near-dilate `{L,2L,…,12L,13L+1}` has *exact* `M = 1/13` at every
diameter; the adversarial minimum does not grow, it sits flat at `1/13`. My generator, like every
fixed-base generator, simply cannot emit near-dilates (MISTAKE-101/127/137, again). So the honest
large-diameter target is the **constant** `1/13 − o(1) > 1/14`, and THM-721's elementary decorrelation
atom already delivers it for the compressed stratum. The Parseval floor and its THM-680 sharpening
remain valid as a per-family certificate — but "growing M" was the wrong flag to plant, and the
signed-`OffLine` machinery is heavier than a constant margin requires. The verification that made S264
feel airtight measured the right quantity on the wrong families. *Never let a clean table stand in for
an adversarial minimum.*

**4. The "≤6 effective speeds" shrink is bounded-diameter only (verified, clean).** The corrected
invariant (opus-S243, after MISTAKE-139 retired "≤6 lifts") is `≤6 coprime to 30030 = 2·3·5·7·11·13`.
It is not a divisor-complete theorem — it is a bounded-diameter regularity. The obstruction is one
line of combinatorics: a *single* speed `= lcm(2..14) = 360360` witnesses every `d∈{2..14}` by
itself, so a primitive DC 13-set can have **12** coprime-to-30030 speeds. Explicit:
`{1,17,19,23,29,31,37,41,43,47,53,59, 360360}` — primitive, divisor-complete, 12 coprime, and still
loose (`M = 23/112 ≈ 0.205`). Divisor-completeness alone forces only **one** non-coprime speed, not
seven. And the break comes early: the exact max #coprime `= min(13 − mincover(Vmax), supply(Vmax))`
(set-cover of `{2..14}` by speeds `≤ Vmax`, capped by how many coprime-to-30030 integers exist below
`Vmax`) already reaches **8 at Vmax = 45**. So `≤6` is not even a bounded-diameter *worst-case*
theorem — it is a *small-diameter (Vmax ≲ 30) typical-case* regularity (opus/mac-mini's "mean 2.0"
is the typical case; adversarial construction beats it by Vmax 45). Large diameter is carried by the
scale tower, and the family stays loose because its units are large and spread — the multi-scale
regime. *(Self-caught: my first sweep read `max = 13 − mincover`, ignoring the coprime supply — a
MISTAKE-138-style overclaim in the other direction; the computer broke it before it shipped.)*

A companion caution that fell out: **auto-safe (opus-S241) is a discrete-clearing statement**
(`bandCount = 0` at a specific `p/q`), not a reach reduction. "Drop the structured speeds and the
loneliness is controlled by the coprime sub-family" is licensed for the bounded-diameter *pigeonhole*
route (opus-S242, proved for ~44% of DC) but **not** for the large-diameter *reach* route — which is
why that route needs THM-721, not the coprime count. The two "drop structured speeds" moves look
identical and are not; conflating them is the shape of MISTAKE-139.

---

## The combinatorial connections that are actually load-bearing

The mining also confirmed which cross-domain bridges are real (not decoration), because two of them
are the reason the simplifications above hold:

- **Schur triples ARE the geometry.** `E₃(S) = #{a+b=c}` is not just "an energy" — it is the set of
  norm-3 *minimal (kissing) vectors* of the relation lattice `Λ(S)`, whose Gram matrix is tridiagonal
  and whose kissing number equals the additive energy (mac-mini-S25, LEM-014). The AP maximizes E₃
  (`=78`), is the unique tight LRC set, and is the transitive/all-ones tournament tiling — one
  extremal object in three languages, because `Λ(S)` *is* the tiling cycle-space refined to the
  integers (the cut⊕cycle seam, THM-373). This is why the fine (`1/14`) scale is governed by E₃ and
  the coarse (`1/7`) scale by E₂ — they are different-weight vectors of the same lattice.

- **The scale grading is the observer grading.** `Λ(S)`/energy is observer-*blind* (tournament-level,
  coarse); coverage/reach is observer-*relative* (tiling-level, anchored at the base path). That
  single distinction predicts, in advance, which tools transfer (cut⊕cycle, deletion-contraction,
  multiplicativity — structural) and which do not (the bare scalar `H`) — and it is exactly why the
  large-diameter *reach* argument cannot borrow the *clearing* argument's "drop structured speeds."

## The one forward move worth naming

Every simplification above shrinks the problem; the one *addition* worth making is the parked
**BSG → Freiman `3k−4`** bridge (opus-S181), and the mining pinned down *where* it fits: the SHALLOW
sub-lemma ("no dilate of a spread DC 13-set puts 6 speeds in a `1/7`-arc") is a **coarse-scale
inverse-sumset** statement, and the coarse scale is precisely where `E₂`/Freiman is the correct
invariant (HYP-5990) — not the fine scale where E₂ is translation-blind and fails (HYP-6060). The
chain: 6-in-an-arc ⟹ small-doubling block ⟹ (BSG) large low-doubling subset ⟹ (Freiman `3k−4`) short
AP ⟹ contradiction with spread (longest-AP ≤ 7). The AP corner is already closed by three-gap
(opus-S236), so BSG→Freiman need only cover the dissociated bulk — its home turf. It has been pointed
at for a dozen sessions and never run; it is the natural next combinatorial experiment.

## The shape of it

Last session's lesson was "the wall is a sign, not a size." This session's is quieter and more
humbling: *some of the wall was scaffolding.* A five-case argument was two cases; a growing target was
a constant; a 13-runner shrink was a bounded-diameter regularity with a one-line counterexample; two
identical-looking "drop the structured speeds" moves were a proof and a fallacy. None of this is a new
theorem. All of it makes the remaining theorem smaller and the finish map honest — including
retracting my own S264 flag. The combinatorics kept saying *you are carrying more than you need*, and
the useful work was to put it down.

*Files: `04-computation/lrc14_coprime30030_scope_klein_S265.py` (+`.out`). HYP-6140. Simplifies the
finish map (case split 5→2, constant-not-growing correction, ≤6-coprime scope); corrects
[[the-wider-band-parseval-floor-reaches-the-true-M-klein-S264]] (growing-M was an artifact).
Connects [[the-runner-movie-is-a-tiling-the-sorted-path-frame-and-the-order-cell-engine-opus-S136]],
[[the-relation-lattice-LAP-has-maximal-kissing-number-additive-triples-macmini-S25]],
[[the-two-axes-share-a-threshold-e3-peel-ladder-kps-S126]],
[[additive-energy-is-necessary-not-sufficient-the-tight-locus-is-resonance-geometry-opus-S181]]
(the BSG→Freiman lead). Builds on death-star THM-721, opus-S243 (MISTAKE-139 correction), opus-S242.*
