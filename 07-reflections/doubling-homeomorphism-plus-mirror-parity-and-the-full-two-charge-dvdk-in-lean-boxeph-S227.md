# Doubling homeomorphism + mirror-parity; and the full two-charge DvdK in Lean

> **LRC CORRECTION (MISTAKE-230):** the claimed `chi(G_S)=0` descent is
> retracted. THM-2075 starts at the nonempty quotient-core safe set `G_C`, not
> at the empty full-row safe set `G_S`; the two original odd tails form a
> zero-safe-child outer step. Mirror parity survives, as do all independent
> two-charge DvdK Lean results in this note.

*boxeph-2026-07-21-S227. Owner: complete the GMC(2) formalization, and work the LRC math — doubling as a
continuous bijection, and the unique safe-child condition. Builds on codex THM-2073 (dyadic deletion tower /
unique safe-child law) + THM-2075 (doubling is a homeomorphism between consecutive safe sets, χ-invariant),
my S212/HYP-8845 (mirror-parity χ even), S224 (exact covering-min), and S226 (the two-charge DvdK Lean).
LRC verified in `04-computation/doubling_homeomorphism_meets_mirror_parity_boxeph_S227.py`; the Lean is
`GMC2DvdKTwoCharge.lean` (kernel-pure).*

## LRC: the doubling homeomorphism meets the mirror-parity

The owner's two prompts name exactly codex's just-committed frontier: **"making doubling a continuous
bijection" = THM-2075** ("doubling is a homeomorphism between consecutive safe sets in the dyadic tower"),
and **"the unique [safe-]child condition" = THM-2073's unique safe-child law**. My mirror-parity is the piece
that composes with them.

**The tower engine (verified).** With `δ=1/14`, `φ_Q(t)=min_q‖qt‖`, `G_Q={t:φ_Q≥δ}`, the doubling map
`D:t↦2t` satisfies `φ_{2C}(t)=φ_C(2t)` (since `‖(2c)t‖=‖c(2t)‖`), so `G_{2C}=D^{-1}(G_C)`. The safe set of
`2C` is the doubling-preimage of `C`'s — the engine of codex's dyadic tower, which lifts a safe set through
doubling and (THM-2073) keeps a **unique safe child** per binary lift.

**Doubling as a continuous bijection (verified).** `D` is 2-to-1 on `S¹` (preimages `t/2`, `t/2+½`); the
**binary address** `a∈{0,1}` (which half `[0,½)` or `[½,1)`) restricts `D` to a *homeomorphism* onto `S¹`
(THM-2075's addresses). The deck involution `τ:t↦t+½` swaps the two address sheets; on a single sheet `D` is
a genuine continuous bijection. (This `τ` is distinct from my mirror `ι:t↦1−t`, the reality/reversal `ℤ/2`.)

**χ-invariance + mirror-parity — the composition.** THM-2075 says component count, endpoint count, and
**Euler characteristic are invariant** under the doubling homeomorphism (per sheet). My S212/HYP-8845 says
for a covering set `G_δ` is `ι`-invariant, so `χ(G_δ)` is **even**. Verified: the dyadic-seam sets
`2{1,2,3}∪{1,3}` and `2{1..6}∪{1,7}` have `χ(G_{1/14})` = 8 and 24 — even. (Honest subtlety: the *full-circle*
`2C←C` doubling is 2-to-1, so raw `χ` *doubles* (`4→8`, `10→20`); THM-2075's *per-sheet* homeomorphism keeps
`χ` fixed. Either way a **disproof `χ=0` is preserved** — `0=2·0`.)

**The sharpened reduction (my contribution).** Combining THM-2073 (a dyadic-seam disproof descends to a
hereditarily-primitive **terminal core** in ≤8 levels) + THM-2075 (`χ` invariant down the tower) + my
mirror-parity (`χ` even):

> **A disproof of LRC(14) in the strict dyadic seam (`S=2C∪{x,y}`, THM-2061) has `χ(G_{1/14})=0`; `χ` is
> doubling-invariant, so the terminal core also has `χ=0`; and `χ` is even, so a disproof needs `χ=0`
> *exactly* (never 1). Hence Wall A (dyadic-seam case) ⟺ no hereditarily-primitive *terminal* core has
> `χ(G_{1/14})=0`.**

This composes two frontier theorems (codex's doubling tower + my mirror-parity) into a strictly smaller,
mirror-symmetric target: the terminal cores (one speed per 2-adic valuation, THM-2073's normal form), on
which the exact pair-sum covering-min (S224) and the `χ≥2` mirror criterion decide loneliness rigorously.
*"Making doubling a continuous bijection"* (THM-2075) is precisely what lets `χ` — and the disproof —
descend intact.

## GMC(2) formalization: the full two-charge DvdK, both orientations

Continuing S226, I extended `GMC2DvdKTwoCharge.lean` to **both sign orientations** of the two-charge
`DvdK1` — the elementary single-character case, now complete and kernel-pure:

- `exists_nonzero_ct_pair` — `f = c₀zᵖ + c₁z⁻ⁿ` (index 0 the `+p` charge): `CT(f^{p+n}) ≠ 0`;
- `exists_nonzero_ct_pair'` — the swapped orientation (index 0 the `−n` charge): same conclusion;
- `dvdk1_pair` — the `DvdK1`-shaped existential.

All three `#print axioms` report only `[propext, Classical.choice, Quot.sound]` — **sorry-free**. The proof
in each orientation is the unique-balanced-composition argument (`balanced_unique`/`balanced_unique'`): the
linear system `r₀+r₁=p+n`, `p·r₀=n·r₁` has a single integer solution, so `CT(f^{p+n})` is one uncancellable
multinomial term. This covers *any* injective straddling `Fin 2 → ℤ` up to relabeling — the whole
single-character face — DvdK-premise-free.

The general `DvdK1` (≥3 charges, where cancellation can occur) is codex's **THM-2067** (Galois orbit-product,
which superseded my retracted S222/S223 saddle route); its Lean formalization + the height package are the
remaining boundary for a fully DvdK-free GMC(2) in Lean.

## Scope

Two honest pieces of progress: (1) an LRC *reduction* composing codex's doubling homeomorphism (THM-2075)
and unique safe-child law (THM-2073) with my mirror-parity (S212) — the disproof descends the dyadic tower
with `χ=0` (even) to a hereditarily-primitive terminal core, verified in its parts; and (2) the *full*
two-charge `DvdK1` (both orientations) now kernel-pure in Lean, completing the single-character leaf of the
GMC(2) spine. Neither closes its conjecture: Wall A's terminal-core `χ>0` and the general `DvdK1` (THM-2067)
+ height package remain the open frontiers.

Links: HYP-8920, THM-2073 (codex), THM-2075 (codex), HYP-8845, THM-2067,
[[the-good-sets-reversal-symmetry-an-equivariant-mirror-parity-sharpening-of-the-chi-criterion-boxeph-S212]],
[[a-kernel-pure-lean-proof-of-the-two-charge-dvdk-seed-boxeph-S226]].
