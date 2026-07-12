# The bounded window is decoupled from tightness — the hard core is loose

*kind-pasteur-2026-07-11-S127. Owner: "get an actual picture of the remaining hard part, and use that to
crack statements." I got the picture, and it inverts the expected one: the families that make the
bounded-window ruler search hard are not the near-counterexamples to LRC — they are the loosest, most
obviously-lonely families in the problem.*

---

## The expected picture, and why it's wrong

The bounded-window claim (HYP-6020/6025): every residual family has a `B5 > 0` ruler at a modulus `q ∈ [8, 43]`,
diameter-free, and the ~`3·10⁻⁴` "composite-core" (families failing every prime ruler) is where it's hardest.
opus-S231 and I both called this "LRC-adjacent" — the natural assumption being that the window-hard families
are the near-tight ones, the neighborhood of the conjecture's extremal locus.

Two measurements say the opposite.

**The composite-core is loose, not tight.** Sampling families that fail every prime ruler `11..43` and
measuring `M(v) = max_t min_i ‖v_i t‖`: mean `M = 0.184`, minimum `0.137` — *zero* near-tight, statistically
identical to random loose families (mean `0.187`). The tight value is `1/14 = 0.071`. These families are
lonely with nearly *double* the required margin. They are hard to *detect* with a small prime ruler — a
residue coincidence — not hard to be lonely.

**The tight families have the easiest ruler of all.** The genuine LRC-hard region — the tight AP `{1..13}`,
the Goddyn–Wong tight point `{1..11,13,24}`, and their neighbors — all have ruler `q = 14`. And by
**dilation-invariance** `B5(c·v, q) = B5(v, q)` for `gcd(c, q) = 1` (verified; the multiplier `p ↦ cp mod q`
just permutes the multipliers), *every dilate* inherits it: `3·AP`, `11·AP`, `101·AP`, `9999·AP` — all tight
(`M = 1/14`), all ruler `q = 14`, at any scale. The hardest families for the Lonely Runner are the *easiest*
for the window.

## Why the two facts fit together

The bounded window's status is governed by two rigorous invariances of `B5`:
- **Residue-periodicity** (`B5_congr_mod`, formalized cont.34): `B5(v, q)` depends only on `v mod q`.
- **Dilation-invariance** (this session): `B5(c·v, q) = B5(v, q)` for `gcd(c, q) = 1`.

Together they say the window-ruler status is a property of the family's *residue-and-scale class*, not its
loneliness margin. A near-tight family and a loose family in the same class share the same ruler. The tight
`AP`'s ruler `q = 14` is exactly its coarse-scale `1/14` witness, and dilation carries it everywhere. So
tightness never obstructs the window — it *is* the window's cleanest case.

What remains hard is the loose, residue-unlucky families: genuinely lonely, but whose residues mod every
small prime happen to defeat the prime `B5`-check, so a composite modulus is needed to *detect* the
loneliness. That is a covering/detection question about `B5` across residue classes — a bounded-modulus
combinatorial statement — with no near-counterexample anywhere near it.

## The statement this cracks

> **The bounded-window ruler existence is decoupled from LRC-tightness. It is a `B5`-detection-completeness
> statement over residue classes, not an LRC-margin statement. The LRC-hard tight families are handled
> trivially by `q = 14` (their coarse `1/14` witness, carried to all scales by dilation-invariance); the
> window-hard families are the loosest ones.**

This matters because it says the window claim is *not* blocked by the conjecture's difficulty. "LRC-adjacent"
was the wrong diagnosis. The obstruction — if there is one — is a residue-covering fact about a low-degree
alternating moment sum at moduli `≤ 43`, and the part of the problem that carries the Lonely Runner's actual
hardness sits entirely in the non-residual, isolated tight points (mac-mini THM-708/709), where `q = 14`
already answers it. The right next move is to attack the detection-completeness directly (a finite,
bounded-modulus covering argument, CRT-factored per HYP-6025), not to treat the window as a proxy for the
conjecture.

## The recurring shape

Once more the structure was hiding a decoupling. The wall was a decorrelation; the certificate was
`liveCount` minus a penalty; the extremal was a variance not a mean; and now the "LRC-hard core" of the
window is loose, while the hardest LRC families are the window's easiest. Each time, the thing that looked
like one hard object was two — a tractable part and a genuinely hard part pulling in opposite directions — and
the work is seeing which is which. Here the hard-looking part (the composite-core) is the easy one (loose,
detectable), and the easy-looking part (tightness) was never in the window at all.

*Files: this session's diagnostics; the two invariances (`B5_congr_mod` formalized cont.34; dilation-invariance
`B5(c·v,q)=B5(v,q)` verified, math-proof via multiplier reindex, formalizable). HYP-6030. Sharpens
HYP-6020/6025, corrects the "LRC-adjacent" framing (opus-S231). Extends
[[two-routes-one-ladder-kps-S127]].*
