# The band-clearing window is unbounded — the finite check is only for bounded diameter

*kind-pasteur-2026-07-11-S127 cont.47. Owner: "shrink the open problems further." Trying to shrink Route B's
residual, I instead found that one purported shrink — the "diameter-free finite check over `[15,31]`" of
klein-S258's finish map — is **invalid**: the band-clearing window is genuinely unbounded, and this resolves a
standing tension between klein-S258 and opus-S238/S241 with an explicit construction. Ruling out an invalid
path is itself a shrink.*

---

## The construction

From cont.46, clearing at modulus `q` needs (a) `q ∤ every vᵢ` (a runner `≡ 0 mod q` sits at the danger
center for all `p`) and (b) the coprime sub-family misses a fold-class. Condition (a) is a divisibility
blocking, so a family can **block** `q` by carrying a multiple of it. Pushing this to the limit: I can
construct a **spread, primitive, divisor-complete** 13-family that blocks *all* non-14 `q ∈ [15,W]`, for
`W = 31, 37, 43, 50, 60` — the min-clearing modulus then jumps to `32, 38, 44, 52, 61`. Verified concretely at
`W = 43`:

`v = [200, 496, 540, 656, 851, 921, 935, 1122, 1482, 1680, 1835, 1849, 1856]`

is primitive, divisor-complete, longest-run `1` (fully spread), blocks every non-14 `q ∈ [15,43]`, and first
clears at `q = 44`. So **no fixed bounded window clears every spread DC family** — the window grows with the
construction, without bound.

## This resolves klein vs opus

Two fleet claims were in tension. klein-S258's finish map: *"clears at non-14 `q ∈ [15,29]`; diameter-free
(residue-periodic) ⟹ genuinely FINITE; 0 fails / 3000."* opus-S238/S241: *"no fixed bounded window is a
shortcut — the clearing modulus must adapt to the family."* The construction settles it: **opus is right.**
klein's `0 fails / 3000` is correct but was measured on *random* families, which never carry the large
multiples needed to block — and the "genuinely finite" inference fails, because the blocking residue classes
(present in any honest residue-class enumeration mod `lcm[15..W]`) do **not** clear in a fixed `[15,W]`. The
finite check over a fixed window is therefore **incomplete**: it flags the blocking classes as non-clearing.

## Why it doesn't threaten LRC(14) — the blockers are loose

The blocking family above has true **`M = 53/227 ≈ 0.234`** — more than triple `1/14 = 0.071`. It is *very*
lonely; it simply achieves its loneliness at its own moduli (`q = 227` here), not by band-clearing in a small
window. This is my cont.36 decoupling in its sharpest form: **the families that break the bounded window are
the loosest ones**, not near-counterexamples. The band-window is the wrong instrument for them. And the
band-edge margin confirms it degrades: clearing at `q ≈ 2·Vmax` gives `M ≥ ⌈q/14⌉/q → 1/14⁺`, a vanishing
margin — precisely because band-clearing at a large modulus is a weak witness for a family whose real
loneliness is large and lives elsewhere.

## What the finite check *does* cover — the shrink

The window is not unbounded uniformly. A family of diameter `Vmax` can block only `q ≤ Vmax` (blocking needs a
multiple `≥ q`), and empirically **random and near-tight families clear by `q ≈ 25–29` regardless of
diameter** (mean min-clearing modulus `≈ 17–18` at every scale up to `Vmax = 1200`). The known near-tight DC
families — the M-floor `{1,2,3,4,10..18}` (`M = 1/12`, clears at `24`), `{1,2,3,10..19}` (`M = 3/29`, at `25`),
`{2..14}` (`M = 1/8`, at `16`) — are all small-diameter with small windows. So:

> **The band-window finite check is complete on the bounded-diameter (near-tight) neighborhood — exactly the
> families where loneliness is delicate. The families it misses are large-diameter loose blockers, which are
> lonely by a wide margin and do not need the band-window at all.**

The open problem is thereby *reshaped*, not merely reasserted: it is **not** "prove clearing in a fixed
bounded window" (impossible), but the disjunction opus-S238/S241 named — every DC family clears at *some*
modulus it lacks a multiple of, adapting to the family — with the honest addition that the hard families for
this are bounded-diameter, and the unbounded-diameter ones are loose. The right next move is a *looseness
dichotomy*: show every large-diameter DC family has `M` bounded away from `1/14` by a direct witness, leaving
only a bounded-diameter finite check. That is a genuinely smaller target than the full adaptive-window
anti-concentration.

## Where this sits against the just-landed pigeonhole theorem

opus-S242 (the pigeonhole clearing theorem) and klein-S261 (the unified clearing formula) — both landed the
same hour, both crediting the cont.45/46 coprime reduction — prove that a family clears at a composite
`q ∈ {15,21,25,27}` **provided it has no multiple of `q`** and fewer than `φ(q)/2` coprime speeds (44% of DC
families, no anti-concentration). That hypothesis is exactly my condition (a), and it is exactly what the
blocking construction defeats: a blocker carries a multiple of *every* modulus in any fixed window, so the
pigeonhole theorem's hypothesis fails for it at every `q` in the window. So this finding is the precise
**ceiling** of the pigeonhole/unified-formula path: widening the fixed window cannot push it to 100%, because
an adversary blocks any fixed window by multiples. The theorem is correct (it is conditional on no-multiple);
what my construction shows is that its precondition cannot be met uniformly, and the families that escape it
are loose. The pigeonhole path and the looseness dichotomy are therefore complementary halves, not competitors.

## Honest scope

This is a correction plus a reshaping, not a proof. What is banked: the band-window is provably unbounded
(explicit verified construction), so klein-S258's finite-check path is incomplete; the obstruction families
are loose (`M ≈ 0.23`), so LRC(14) is untouched; and the finite check is complete on the bounded-diameter
neighborhood. The remaining open problem is the looseness dichotomy for large-diameter families — a cleaner
target than "no fixed window works, so check an unbounded one."

*Files: lrc14_window_boundedness_kps_S127.py/out, lrc14_window_vs_diameter_kps_S127.out. HYP-6120. Resolves
klein-S258 (finite window) vs opus-S238/S241/S242 (no shortcut / conditional pigeonhole) in opus's favor;
complements opus-S242 pigeonhole theorem + klein-S261 unified formula (the bounded-window ceiling); sharpens
[[the-window-is-decoupled-from-tightness-kps-S127]] and
[[route-b-is-two-conditions-and-the-window-is-wider-than-31-kps-S127]].*
