# The d=3 detuned peel is one union-bound wider — the counting brick is per-coordinate

*kind-pasteur-2026-07-10-S127. Owner: "take the d=3 generalization next, but pull from other agents and
synthesize." This note records the synthesis that made the d=3 case a short extension of opus's d=2, and
the clean structure of its exceptional set.*

---

## The synthesis that shrank the work

opus built the d=2 detuned dispatch across two sessions: the construction core + reduction
(`LRCDetunedD2`, S210) and the counting `twoDetunedClearing` (`LRCTwoDetunedClearing`, S211). The counting
rests on one lemma — `LRCIntervalCount.bad_count_le` — which bounds, **for a single detuned speed**, the
number of bad branches `c ∈ [0,g)`:

```
|{c : δ(u+c)/g fails to clear 1/14}|  ≤  gcd(δ,g) · (⌊q/7⌋ + 1),   q = g/gcd(δ,g).
```

The moment I read that this brick is stated **per coordinate**, the d=3 generalization collapsed to a
mechanical exercise: call it a *third* time, and replace opus's two-set union bound with a three-set one.
No new analysis — `bad_count_le` already contains all the geometry (the de-circled `ψ`-injection into a
width-`q/7` interval). The only genuinely new lines are:

* `lonely14_of_three_detuned_good` — the construction core, one `rcases` case wider (three detuned
  coordinates cleared at the branch, the rest harmonic);
* the three-set cover bound `(bad₁ ∪ bad₂ ∪ bad₃).card ≤ Σ badCount δⱼ g` via two applications of
  `Finset.card_union_le`;
* the assembly through LRC(**10**) instead of LRC(11) — ten harmonic speeds, clearance `1/11 ≥ 1/14`.

That is the whole file. `threeDetunedClearing` is proved kernel-pure, and `lonely14_of_three_detuned'`
makes the generic d=3 case **unconditional from the LRC(≤13) citation**, exactly as opus's
`lonely14_of_two_detuned'` did for d=2. This is what "pull from other agents and synthesize" looks like when
it works: the reused brick was designed per-coordinate, so the dimension count was never the hard part.

## The generic condition, and where it stops being generic

The three bad sets fail to cover `[0,g)` iff their sizes sum below `g`:

```
Σⱼ gcd(δⱼ,g)·(⌊qⱼ/7⌋+1) < g   ⟺   Σⱼ (⌊qⱼ/7⌋+1)/qⱼ < 1
```

(the `dⱼ = g/qⱼ` cancel against `g = dⱼqⱼ`, leaving a condition on the `qⱼ` alone). I stated
`ThreeDetunedClearing` with this sum as an explicit hypothesis — the exact shape `bad_count_le` produces —
so the proof is a direct union bound with nothing hidden.

The **exceptional set** (where the sum is `≥ 1`, so the generic bound fails) has a clean shape
(`lrc14_three_detuned_exceptional_kps_S127`):

- **`min qⱼ ≥ 4` ⟹ always generic.** Three terms each `≤ 1/4` sum to `≤ 3/4 < 1`.
- **`(2,2,·)` is the one infinite family.** `term(2)+term(2) = ½+½ = 1` already, so *any* third `q₃` fails.
  This is the **double-half-harmonic**: two detuned speeds sitting at `q = 2` (half-integer relative to
  `g`). It is the exact d=3 analogue of opus's d=2 `(2,2)` residual — and it fails for the same reason, so
  it needs the same escape: the **mod-`2g` lift**, not the branch count.
- Everything else (`q₁ ∈ {2,3}`, e.g. `(2,3,q₃≤42)`, `(2,4,·)`, `(3,3,3)`) is a **finite** small-`q` set.

So the counting closes the d=3 case down to `(2,2,·)` plus finitely many small triples — the same residual
shape opus faces at d=2, one dimension up. The lonely-runner "danger" concentrates, as always, where speeds
are half-integers of the common scale.

## The pattern for d ≥ 4

The recursion is now visible and cheap: at `d = k` strays, call `bad_count_le` `k` times, union-bound `k`
sets, assemble through LRC(`13−k`) (clearance `1/(14−k) ≥ 1/14`). The generic bound `Σⱼ Nⱼ/qⱼ < 1` holds
once `min qⱼ ≥ ⌈k/(1 − (k−1)/…)⌉`-ish — concretely, for `d` strays, `min qⱼ ≥ d+1` already gives
`Σ ≤ d/(d+1) < 1` in the range-free part. Each new dimension is one more `card_union_le` and one lower
citation index. The strays-at-one-scale ladder is uniform; only the residual half-harmonic locus (`several
qⱼ = 2`) ever needs the mod-lift.

*Files: `LRCDetunedD3.lean` (the d=3 dispatch), `lrc14_three_detuned_exceptional_kps_S127.py`/`.out`
(the exceptional set). Builds on opus's `LRCDetunedD2` / `LRCTwoDetunedClearing` / `LRCIntervalCount`.
The dispatch consumer is opus's `MultiDetunedDispatch` (S209), which cites d=2 **and** d=3 — the generic
half of the d=3 citation is now discharged.*
