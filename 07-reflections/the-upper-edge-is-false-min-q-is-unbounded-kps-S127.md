# The upper edge is false — the minimal strictly-live modulus is unbounded

*kind-pasteur-2026-07-10-S127. Owner: "prove the upper edge q ≤ 27 on the residual class." I could not
prove it because it is false. This note is the refutation, a correction of my own S127cont9 conjecture,
and the recalibration it forces.*

---

## What I was asked to prove, and why I didn't

Last turn I measured `min q ∈ [15,27]` on the residual class and conjectured `BoundedStrictlyLiveSupply B`
(a fixed `B` bounding the minimal strictly-live modulus), flagging it explicitly as *measured, not proved*.
Asked to prove the upper edge, I tested it first. It is **false**, and so is any fixed bound.

## The refutation

The strictly-live condition at modulus `q` depends only on `v mod q`. So a family whose residues are
"tight" (admit no live multiplier) at every `q ∈ [15,B]` fails all of them — its first strictly-live ruler
is above `B`. The question was only whether such families are *residual* or excluded by the hypotheses.

The naive adversary — the `μ = 0` locus, the dilated APs `c·{1,…,13}` — **is** excluded, by
difference-primitivity: `c·{1,…,13}` has difference-gcd `c > 1`, so it is dispatched by the common-residue
branch. That is genuine structure, and it is why my earlier sampling (small `Vmax`, generic families) saw
nothing above 27.

But a broader search found a **fully residual** counterexample:

```
v = {210, 1378, 1379, 2106, 2222, 2247, 3650, 3773, 4123, 5083, 5561, 5680, 6000}
    covering ✓   ratio 28.57 > 13 ✓   compressed ✓   distinct ✓
    not detuned ✓   difference-primitive ✓   not near-AP ✓   some |v| ≥ 23 ✓
    min strictly-live q = 28
```

It satisfies every residual hypothesis, and its first strictly-live ruler is at `q = 28`. And maximizing
`min q` over the residual class reaches **43** at `Vmax ≈ 10⁵`. The minimal strictly-live modulus is
**unbounded**.

Certified in Lean (`LRCSmallRuler`, kernel `decide`, no `native_decide`):

- `cex_covering` — covering (13 explicit witnesses),
- `cex_gapFamily` — ratio `> 13` (reaches the residual branch),
- `cex_no_ruler_below_28` — no strictly-live multiplier at any `q ∈ [15,27]`.

Axioms `[propext, Classical.choice, Quot.sound]`. The full residuality (compressed, not-detuned,
difference-primitive, not-near-AP) is verified computationally in `lrc14_upper_edge_refuted_kps_S127`.

## What survives, and what breaks

- **Survives (proved):** the *lower* edge. `strictlyLive_modulus_ge_15` — covering forces `q ≥ 15` — is a
  theorem and stands. The witness window opens at 15.
- **Breaks (refuted):** the *upper* edge. `BoundedStrictlyLiveSupply B` is false for every fixed `B`. The
  implication `lrc14_of_boundedStrictlyLiveSupply` remains a true theorem, but its hypothesis is
  unachievable, so it is a dead route. The docstring now says so.

## The recalibration — the transfer is load-bearing

My S127cont9 headline — *"the witness lives at the bottom of the window, not the transfer tail"* — is
**true for generic families and false for the hard ones.** Near-tight residual families (residue-close to
the dilated AP, but with coprime differences and broken AP structure) have all their small-`q` rulers dead
and only become live at large `q`. For them, klein's THM-685 transfer — live rulers at `q > Σv/μ` — is
not a superfluous safety net that the small-ruler law makes redundant. **It is the actual mechanism.**

So the endgame's shape is the reverse of what I sketched last turn:

- The measure floor `μ(S) > 0` is *not* replaceable by a bounded small-modulus search. It is the real
  content, precisely because `min q` is unbounded.
- The correct remaining obligation is `StrictlyLiveSupply` (the `B = ∞` form) or, equivalently, the
  measure floor. Both are the open case of LRC(14); neither collapses to a finite check.
- The two edges are honest: `15 ≤ min q < ∞`, lower edge proved, upper edge open and *provably not
  uniform*.

## The lesson, again

I keep relearning it and it keeps paying: **before proving a crisp bound, test it — especially the
adversarial limit you did not sample.** My "witness at the bottom" was a small-`Vmax` artifact. The
dilated-AP adversary looked like the worst case and was excluded by difference-primitivity, which almost
sold me the bound; the real counterexample was a *coprime-difference near-tight* family, one filter deeper.
Two turns ago I nearly published "covering ⟹ μ > 0" (false); this turn I nearly tried to prove
"`min q ≤ 27`" (false). The wall does not have a small door. It is the conjecture, and the transfer is the
only way through.

*Files: `LRCSmallRuler.lean` (the certified counterexample + corrected `BoundedStrictlyLiveSupply`
docstring), `lrc14_upper_edge_refuted_kps_S127.py`/`.out`. Corrects
[[attacking-the-wall-the-witness-lives-at-the-bottom-kps-S127]]. See
[[the-residual-is-one-measure-floor-kps-S127]] and
[[describing-the-wall-the-scale-gap-is-the-separator-kps-S127]].*
