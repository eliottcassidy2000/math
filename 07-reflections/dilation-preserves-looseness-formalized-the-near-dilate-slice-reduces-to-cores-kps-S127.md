# Dilation preserves looseness, formalized — the near-dilate slice reduces to bounded cores

*kind-pasteur-2026-07-11-S127 cont.50. Owner: "keep working the remaining open math." boxeph's MISTAKE-140
just corrected the whole large-diameter picture — "min M grows with diameter" is a sampling artifact, because
`M` is dilation-invariant. That correction rests on a transport fact that was only cited (THM-531) or measured.
This session machine-checks it.*

---

## The correction it underpins

boxeph-S19 (MISTAKE-140, the third recurrence of the "size-indexed extremal over a dilation-closed class"
genus) refuted the claim that `min M` over spread DC families grows with diameter: the spread 12-core
`H* = {1,2,3,4,8,9,10,11,12,13,14,16}` has `M = 1/11`, and its near-dilates `v_c = 2c·H* ∪ {δ_c}` are
primitive, divisor-complete, spread, of diameter `30c → ∞`, and **`M(v_c) = 1/11` at every scale**. So the DC
class stratifies by **structure**, not diameter. This is the same lesson as my own cont.47 (the band-clearing
window is bounded only under bounded diameter) and cont.49 (the ≤6 coprime-core is a diameter/availability
artifact): diameter is a proxy; dilation carries structure to every scale. Three "bounded-at-bounded-diameter,
unbounded in general" findings, one cause.

## What was formalized

The transport under dilation — `reach(c·v) ≥ reach(v)` — is the load-bearing fact, and its *looseness-carrying*
direction is elementary. `LRCReachTransport.lean` (kernel-pure `[propext, Classical.choice, Quot.sound]`):

- **`reach_dilate_ge`** — for `c ≥ 1`, `reach(v) ≤ reach(c·v)`. The proof is a single scaled witness: if `t₀`
  attains `reach v`, then `t₀/c ∈ [0,1]` attains the *same* margin for `c·v`, because
  `(c·vᵢ)·(t₀/c) = vᵢ·t₀`. No periodicity, no equidistribution — the reach can only go up under dilation.
- **`loose_dilate`** — the consequence: if a family `v` is loose (`reach v > 1/14`) then every integer dilate
  `c·v` is loose, at every scale, with no recomputation.

## Why the `≥` direction is the right half

The full dilation-invariance is `reach(c·v) = reach(v)` (the `≤` half needs 1-periodicity of the margin over
a longer interval — real work). But for the endgame only the `≥` direction is load-bearing: it **propagates
looseness upward**. A bounded structural core that clears `1/14` makes its entire unbounded dilation orbit
clear `1/14` — so the near-dilate slice (`v_c = 2c·H* ∪ {δ}`, all `c`) collapses to the single bounded core,
which lives in the finite base. This is exactly the reduction MISTAKE-140 identified ("re-anchor to the flat
`[1/13,1/11]` floor; the class stratifies by structure") — now machine-checked as the arrow from core to
orbit. The finite check no longer has to see infinitely many dilates; it checks the cores, and `loose_dilate`
carries the rest.

## Honest scope

This formalizes the transport, not the endgame. Two things stay outside: (a) the core's own looseness — that
`H*` and the other structural cores have `reach > 1/14` — which is the detuned-dispatch / bounded-diameter
finite check (the genuine remaining work); and (b) the `≤` half of dilation-invariance (not needed for
looseness, so deliberately skipped). What is banked: the "reduce the unbounded near-dilate slice to bounded
cores" step is now a kernel-checked lemma, so no rigor plan has to re-verify the orbit of a loose core — a
recurring source of the size-indexed-extremal mistakes (MISTAKE-137/139/140). The transport that keeps
tripping the fleet up as an *informal* aside is now a theorem.

*Files: LRCReachTransport.lean (kernel-pure; `reach_dilate_ge`, `loose_dilate`), root-wired. HYP-6160.
Formal underpinning of MISTAKE-140 (boxeph-S19) / THM-531 (dilation-invariance); reuses the reach framework of
[[the-13-runner-decorrelation-atom-formalized-the-large-diameter-half-in-lean-kps-S127]]; same "structure not
diameter" lesson as
[[the-clearing-window-is-unbounded-the-finite-check-is-only-bounded-diameter-kps-S127]] and
[[divisor-completeness-is-one-integer-so-the-coprime-core-size-is-a-diameter-artifact-kps-S127]].*
