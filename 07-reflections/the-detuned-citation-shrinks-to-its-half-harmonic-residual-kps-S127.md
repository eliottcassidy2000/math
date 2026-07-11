# The detuned citation shrinks to its half-harmonic residual

*kind-pasteur-2026-07-10-S127. Owner: "wire threeDetunedClearing into MultiDetunedDispatch, work remaining
tasks to get LRC(14) into its best state." This note records what the wiring buys: a cited black box
becomes a proved bulk plus a small named residual.*

---

## From a black-box citation to a proved reduction

opus's `MultiDetunedDispatch` (S209) entered the assembly as a **citation** — THM-678, taken on faith:
*every family with some `g ≥ 2` at detuning level `nonMultCard v g ∈ {2,3}` is lonely.* It peels the
near-dilate μ-minimizers the residual floor cannot bound. But a citation is a debt, and the "best state" of
a formalization is the one that owes the least.

Two sessions closed the *generic* half of that debt — opus's `twoDetunedClearing` (d=2) and my
`threeDetunedClearing` (d=3), both a union bound over the same per-coordinate brick
`LRCIntervalCount.bad_count_le`. This session wires them into the dispatch and pays down the debt:

```
def genericCount v g := Σ_{g∤vᵢ} badCount (vᵢ) g < g              -- the union bound closes
def ExceptionalDetunedDispatch := ∀ v g, nonMultCard∈{2,3} → ¬genericCount v g → lonely
theorem multiDetunedDispatch_of_exceptional (cite) (hexc) : MultiDetunedDispatch
```

The proof is a single `by_cases` on `genericCount v g`: generic ⟹ discharged by the proved d=2/d=3 wires;
non-generic ⟹ the (much smaller) exceptional obligation. So **`MultiDetunedDispatch` reduces to
`ExceptionalDetunedDispatch`** — and the cited THM-678 shrinks from *all* detuned `d ∈ {2,3}` families to
only the ones where the branch count saturates.

## What the residual actually is

The enumeration (`lrc14_three_detuned_exceptional`) pins the non-generic set exactly:

- **d=2:** the single pair `(2,2)` — both detuned speeds at `q = 2`.
- **d=3:** the infinite family `(2,2,·)` plus finitely many small-`q` triples (`q₁ ∈ {2,3}`).

Every one of these has **at least two coordinates at `q = 2`** — two speeds sitting at half-integers of the
common scale `g`. That is not a coincidence of the counting; it is the geometry. A `q = 2` coordinate
contributes `badCount = gcd = g/2` — *half the branch interval* — so two of them fill `[0,g)` completely
(`badCount_of_q_two`: two `q=2` counts sum to exactly `g`), and no single branch can clear both. The
count-based dispatch is powerless precisely there, and the escape is the mod-`2g` lift (opus's THM-678
residual), which doubles the scale so the two half-harmonics separate.

So the residual is not an arbitrary leftover — it is **the half-harmonic locus**, the same object that the
lonely-runner problem concentrates its difficulty on everywhere: speeds that are rational with small
denominator relative to the scale. The generic bulk is where the runners are "generically detuned" and a
counting argument suffices; the residual is where they lock into half-integers and need a genuine
construction.

## The shape of a good reduction

The move here is worth naming as a pattern: **when a cited lemma is `∀ x, P x`, split `P x` into a
decidable generic condition `G x` you can prove and its complement `¬G x`, then cite only `∀ x, ¬G x → P x`.**
The citation surface shrinks to exactly the hard core, and — crucially — the enumeration of `¬G` tells you
*what* the hard core is (here: two half-harmonics). A black-box citation hides its difficulty; a
generic/exceptional split exposes it. `lrc14_grand_assembly_dissoc_exceptional` threads the shrunk citation
through opus's dissociated assembly, so the endgame now cites THM-678 only for the half-harmonic residual.

*Files: `LRCDetunedDispatchReduce.lean` (`multiDetunedDispatch_of_exceptional`,
`lonely14_of_nonMultCard_two/three`, `badCount_of_q_two`), building on `LRCDetunedD3` (kps) and opus's
`LRCTwoDetunedClearing` / `LRCDissociatedAssembly`. Continues
[[the-d3-detuned-peel-reuses-the-per-coordinate-brick-kps-S127]]. The residual is the mod-2g lift — opus's
THM-678 half-harmonic case, still cited.*
