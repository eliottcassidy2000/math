# The Route-2 audit: two broken links, and what actually holds

**opus-2026-07-06-S130.** Owner directive: *"deeply consider the state of our LRC(14) proof attempts
and the direction they attempt, and reroute as needed; ensure carefully we have everything correct."*
A skeptical audit of Route 2 (the J-K → rank-2 → 1-D-gap thread I had just finished wiring in Lean)
found that **Route 2 does not prove LRC(14)** — it is broken at *both* ends. This is the most
important finding of the session, and it is a healthy one: the scaffolding was real, but the load
was resting on two supports that don't hold.

## Broken link 1 — the top (J-K): wrong object

The proof map (S121) had promoted the J-K step from "unbuilt bridge" to "citable published
reduction," and every Route-2 session since treated `(A) ⟹ LRC(14)` as closed-by-citation. The
audit (verified against the arXiv:2304.01462 abstract) says otherwise. Giri–Kravitz, *The structure
of Lonely Runner spectra*, study the **accumulation points** of the spectrum — their theorem is
`acc(S(n)) = S(n−1)`, and "rank-2 governs the structure" is the inclusion `acc(S(n)) ⊆ S₂(n)`. But
the LRC is the containment `S(n) ⊆ [0, 1/2 − 1/(n+1)]` — a statement about the **supremum**. The
abstract's own words: *"Rather than attack this conjecture, we study the structure of the sets
S(n)."* The extremal LR configuration is generically an **isolated** maximum, not an accumulation
point, so accumulation-point structure says nothing about the top of the spectrum. Controlling
`M(U)` for rank-2 `U` does not bound the sup. **So even a complete proof of `(C)`/`(A)` would not
prove LRC(14) through this citation.** The numbers `1/13, 2/25` appear nowhere in the papers; they
are the project's own construction.

This is the classic failure mode of a citation used at "the structural level": the shape of the
theorem (rank-2 ↔ lower spectrum) is real and seductive, but the *quantity* it governs (accumulation
points) is not the quantity the target needs (supremum). The lesson: **a citation is only load-
bearing if you can state the exact implication it licenses and check the objects match.** "Well-
supported lead" was the honest label in the proof map's own fine print; the sessions then quietly
upgraded it to "citation." Re-reading the source dissolved the upgrade.

## Broken link 2 — the bottom ((C)): not finite

mac-mini-S36 (HYP-4667), which I verified independently (`lrc_covering_completeness_audit`): the
"finite covering system `{2..Q₀}`" that the whole fleet (me included, S126/S127) had converged on is
**incomplete**. Fix `L = lcm(2..Q₀)`; the family `{i + L·kᵢ}` with varying `kᵢ ∈ {1,2}` is
compressed (ratio → 2), shares the AP's residues at every `q ≤ Q₀` (so escapes them all), is not a
translate, and clears only at `nextprime(Q₀)`. The clearing modulus is **unbounded** — it grows with
`lcm(2..Q₀)`, and the escape families sit at height `~10¹⁴`, far past the `≤ 650k` samples that made
`Q₀ ≤ 32` look sufficient. "Every non-AP family clears at *some* `q`" is not a finite reduction of
`(C)` — it is *equivalent* to the analytic statement that every non-AP family is loose. The false
step was "compressed ⟹ bounded lift": compression bounds the lift *range*, not its *values*.

## What actually holds — the salvage

The audit is not nihilism. Three things survive intact:

1. **`(C)` is TRUE.** `lrc_gap_member_search` (exact `M`): `M(AP) = 1/13` exactly, zero gap members
   in 3550 near-AP families, and the *only* families attaining `1/13` are permutations of the AP.
   The moat `(1/13, 2/25)` is empty. The target is sound.
2. **The Lean is sound.** Every theorem I wired is a valid *conditional implication* or an honest
   `Prop` obligation — `torus_loose_of_rank2`, the rigidity wrappers, `crux_of_covering_complete`,
   `loose_of_covering_set`. Nothing asserts a falsehood. What was wrong was the *docstrings/session
   framing* calling `JKReduction` "a citation to pin" and `CoveringComplete` "a finite residue
   check." Both corrected in place.
3. **The rank-2 rigidity is real mathematics** — a true statement about when all 1-D projections of
   a 2-torus being dilated APs forces rank ≤ 1. It just doesn't connect to LRC(14).

## The reroute

- **Route 2 is retired as a proof route** (retained as correct conditional math / a spectrum study).
- **Route 1 is the correctly-aimed project route:** it bounds `Mreach ≥ 1/14` *directly* — the
  supremum, the right object — via the witness-density node. Its open pieces (Part-A density ⟹
  reach, the `k = 8..13` witness floor) are honest analysis, not a wrong-object mirage.
- **The recognized external route** to a genuinely new LRC case is Tao's finite reduction (2018) +
  Malikiosis–Santos–Schymura (2025) + the Rosenfeld/Trakulthongchai sieving/polynomial method — the
  method that actually proved n = 8..13. The experts (Trakulthongchai, per Quanta 2026) say pushing
  to the next case needs "an entirely new sort of way of looking at things"; the bottleneck is the
  efficient computation of `I(k,p,1)`. LRC(14) is a hard frontier, not three obligations from done.

## The meta-lesson

Two independent over-optimisms compounded into a false sense of near-completion: a citation upgraded
past what it says, and a computation whose sample height hid the counterexamples. Both were caught by
the *same* discipline — go back to the primary object (read the abstract; construct the family at its
true scale) instead of the cached summary of it. The fleet's convergence (multiple agents naming the
same "finite covering") felt like confirmation, but convergence on a *frame* is not verification of
it. When several agents agree a hard problem has become easy, that is exactly when to re-derive from
the source.
