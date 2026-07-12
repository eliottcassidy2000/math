# The three-gap disjunction collapses to "does a 13-block fit in the safe band"

*kind-pasteur-2026-07-11-S127 cont.44. Owner: "work the 3 gap disjunction." I worked it, and for the
AP corner the general three-gap machinery turned out to be overkill — the extreme case (a single
consecutive block, reached by two explicit multipliers) closes it, and the clearing condition is a
one-line modular interval.*

---

## Where the three-gap disjunction lives

The divisor-complete residual (`divisor-complete ⟹ M > 1/14`, = LRC(14) via THM-366) splits into two
halves the fleet has now separated cleanly:

- the **near-AP corner** (longest-AP ≥ `k−5` = 8) — closed by LEM-012 (klein-S196), the gap-splitting
  good-period bound;
- the **spread bulk** (longest-AP ≤ 7, ~99% of the class, opus-S237) — which opus-S238 showed needs a
  genuine, irreducible window disjunction (no fixed small-modulus shortcut).

opus-S236 isolated the AP's special structure: for an AP `{a + i·d}`, the residues `{(a+i·d)·p mod q}`
are *themselves an AP mod q* (common difference `d·p`), so clearing at `q` = **an AP on ℤ/q avoiding
the danger arc `{0, ±1}`** — a **three-gap / Steinhaus** statement. That is the "3-gap disjunction" for
the AP corner. This session makes it explicit, and finds it collapses.

## The collapse: the consecutive reduction

The three-gap theorem says the `13` points `{i·(dp) mod q}` cut the circle into arcs of at most three
lengths. The **extreme** case is `dp ≡ ±1`, where all `13` points are *consecutive* and there is a
single big gap. And that extreme case is reachable *by choice of multiplier*: since divisor-completeness
forces `gcd(d, q) = 1` for every `q` with prime factors `≤ 13` (see below), the multiplier `p ≡ d⁻¹
(mod q)` sends the AP to the `13` consecutive residues `{r, r+1, …, r+12}` with `r = a·d⁻¹ mod q`. Now
clearing is trivial to read off:

> **Consecutive-block lemma (exact, formalized kernel-pure).** For `16 ≤ q ≤ 28` the safe band is exactly
> `[2, q−2]`, and a length-13 consecutive block `{r, …, r+12}` avoids the danger arc `{0, ±1}`
> **iff `r ∈ [2, q−14]`** (the block fits without wrapping). The reflected multiplier `p ≡ −d⁻¹` clears
> iff `r ∈ [14, q−2]`.

So for the AP corner the whole three-gap apparatus reduces to: *does a 13-wide block fit inside the
`(q−3)`-wide safe band, at the offset `a·d⁻¹ mod q`?* Two multipliers, `p = ±d⁻¹`, and one modular
interval. Verified: this single disjunction over `p = ±d⁻¹`, `q ∈ [16, 27]`, clears **every** primitive
divisor-complete 13-term AP — 0 failures over 2161. The general three-gap (other `dp`) gives *more*
clearing multipliers, but the corner never needs them.

## Why it is provable, not just verified

Two structural facts turn the empirical 0/2161 into something with a proof-shaped spine:

1. **Divisor-completeness forces `d` coprime to every prime `≤ 13`.** If a prime `p ≤ 13` divided `d`,
   all `13` AP terms would share the residue `a mod p`; divisor-completeness demands a multiple of `p`
   among them, forcing `a ≡ 0`, hence *all* terms `≡ 0 (mod p)` — not primitive. So `q = 25 = 5²`,
   `26 = 2·13`, `27 = 3³` all have `gcd(d, q) = 1` and are **always usable**. (Only the primes `> 13` in
   the window — `17, 19, 23` — can share a factor with `d`.)

2. **`q = 27` is the sharp single modulus.** At `q = 27` the union `[2, 13] ∪ [14, 25]` is *everything
   except `{0, 1, 26}`*, so `q = 27` clears the AP **unless `a ≡ 0` or `±d (mod 27)`** — a
   codimension-one exception hitting only 11% of the corner (239 / 2161). Every one of those exceptions
   clears at a second modulus in `[16, 26]`.

So the pure-AP corner is closed by a proved elementary lemma plus a finite residue-covering disjunction
over `[16, 27]`, with `q = 27` alone carrying 89% and the coprimality guaranteeing the fallback moduli
`25, 26` are live. This is the sharp, explicit, Lean-ready form of opus-S236's "verified" AP result —
tighter window (`[16,27]` vs `[15,31]`), exact mechanism (`p = ±d⁻¹`), and a formalized core.

## What is formalized

`LRCThreeGapConsecutive.lean` (kernel-pure, `[propext, Classical.choice, Quot.sound]`):
`inBand_of_residue_mem_band` (a residue in `[2,q-2]` is safe for `16≤q≤28`) →
`bandCount_zero_of_consecutive_block` (a consecutive block in the band ⟹ `bandCount = 0`, the family
clears) → `liveCount_pos_of_consecutive_block` (⟹ a live multiplier, the hypothesis the band-edge margin
lemma consumes). The one hypothesis left abstract — that the residues *are* the block `r + i` — is exactly
"the consecutive reduction happened," dischargeable for an AP at `p ≡ d⁻¹` by `(a+i·d)·d⁻¹ ≡ a·d⁻¹ + i`.

## The recurring shape

Once more the hard-looking object had an easy spine. The "three-gap disjunction" sounds like it needs
the continued-fraction gap structure of an AP mod `q`; but the AP is so rigid that clearing collapses to
the single extreme gap — a consecutive block — and "does a block fit in a band" is a `≤`/`≥` pair that
`omega` closes. The structured corner was structured enough to be *elementary*. The genuine difficulty is
all in the other half: the spread bulk, where there is no AP to make consecutive, no `d⁻¹` to reach for,
and — by opus-S238 — no small set of moduli that suffices. That is the wall; this corner was a door.

*Files: lrc14_threegap_disjunction_kps_S127.py/out, lrc14_threegap_structural_kps_S127.out,
LRCThreeGapConsecutive.lean. HYP-6090. Sharpens opus-S236 (AP sub-case is three-gap), complements LEM-012
(near-AP gap-split, klein-S196) at the `L=13` endpoint, feeds opus-S235 (band-edge margin). Complements
mac-mini-S65cont44 (three-gap regularity is *why* the AP is the extremal pole) from the constructive side:
the same regularity is *how* the AP clears. Extends [[the-open-theorem-is-lem010-jstar-two-block-kps-S127]].*
