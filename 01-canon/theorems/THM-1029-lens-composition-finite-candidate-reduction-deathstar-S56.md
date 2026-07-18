# THM-1029 — Lens composition: stability × covering bounds each non-AP core to finitely many candidate families (death-star-2026-07-18-S56)

**Status:** a **method**, not a closure. Composing two proved lenses — the stability bound (THM-1028)
and the covering/lcm bound (boxeph THM-1017) — reduces boxeph's inverse theorem for each non-AP core to
a **finite** check over candidate far-elements. Verified to eliminate the near-AP core family
`{1..11,X}`. Does **not** close LRC(14) (infinitely many cores; the check is per-core), but it is the
concrete "bound the crux gradually with relations between lenses" mechanism, and it works on every core
tested. Source HYP-7305/7362. Scripts: `04-computation/lrc_lens_close_deathstar_S56.py`,
`lrc_finite_candidate_deathstar_S56.py`.

Setting: `V = W ∪ {v_max}` primitive covering, `M(V) < 1/13`, `W = V∖{v_max}` (`|W|=12`), `W` **not** a
dilated AP. Goal (the inverse theorem): show no such `V` exists.

---

## The two lenses and their composition

**Lens 1 — stability (THM-1028, PROVED).** `M(V)<1/13`, `W` non-AP ⟹ `v_max ≤ max(W)/(13δ)`,
`δ = M(W) − 1/13 > 0` (the level-`1/13` good component of the non-tight core fits in one `v_max`-arc).

**Lens 2 — covering/lcm (boxeph THM-1017, PROVED).** If `W` misses 13 and 14 (the AP-core case), then
`13 ∣ v_max` and `14 ∣ v_max`, so `182 = lcm(13,14) ∣ v_max`, hence `v_max ∈ 182·ℤ`, `v_max ≥ 182`.

**Composition.** Since also `v_max = max V ≥ max(W)`, the candidate far-elements are
```
v_max ∈ [ max(max(W), 182),  max(W)/(13δ) ] ∩ 182·ℤ,
```
a **finite** set of size `≈ max(W)/(2366·δ)` (`2366 = 182·13`). The window is nonempty only when
`max(W) ≥ 2366·δ`. So **each non-AP core admits only finitely many candidate families**, and the inverse
theorem for that core is decided by computing `M` on each.

## Verification (the near-AP cores are eliminated)

For `W = {1,…,11, X}` (non-AP, misses 13,14), `M(W) = 1/12`, `δ = 1/156`, so `max(W)/(13δ) = 12·X`.
The window `[182, 12X] ∩ 182·ℤ` contains **only `v_max = 182`** for `X ≤ 22`. The single candidate family:

| `W` | candidate `v_max` | `M(V)` |
|---|---|---|
| `{1..11,16}` … `{1..11,22}` | `182` (unique) | `1/12 ≥ 1/13` ✗ |

Every candidate gives `M(V) = 1/12 ≥ 1/13` — so **each of these non-AP cores is eliminated by a
finite (here single-family) check.** The mechanism is exactly the composition: stability caps `v_max`
just above `182`, covering forces `v_max` to be a multiple of `182`, and the two windows intersect in
one point that fails.

## What this is, and is not

**Is:** a reduction of the inverse theorem to, per non-AP core, a finite candidate check — driven by
*relations between lenses* (stability caps from above, covering pins to `182·ℤ` from below). It is the
first argument that turns "the far element could be anything" into "the far element is one of `≈
max(W)/2366δ` explicit values," and it eliminates every core tested.

**Is not:** a closure. The number of non-AP cores is infinite. Two gaps remain:
1. **Compact cores (`max(W) < 182`):** finitely many candidate far-elements *per* core, but infinitely
   many cores (`C(181,12)` scale). Closing needs a *uniform* reason every candidate fails — the tested
   cores all satisfy `M(V) = M(sub-AP)` (a `{1..k}` prefix dominates and the far element can't drag `M`
   below `1/(k+1) > 1/13`); promoting this to all compact non-AP cores is the residual finite/structural
   lemma.
2. **Non-compact cores (`max(W) ≥ 182`):** `W` has its own far element ≥ 182, so recurse
   (renormalization / the difference flow, HYP-3901) — `W = core' + far'`, tower terminates.

## The simultaneous-lens picture (the crux under many sharp invariants)

The crux `M<1/13 ⟹ AP core` is now visible through a battery of *pinned* lenses, each a necessary
condition, jointly pinning the family:

| lens | invariant, pinned by `M<1/13` | ref |
|---|---|---|
| additive / Freiman | residues in **2 cosets** of `val·ℤ` (dim-2 coset progression) | Freiman-3k-4 reflection |
| tournament / OCF | residue tournament `≅ R₁₃ + 1` antipodal flip; **`H = 3551083`** | tournament / OCF reflections |
| fixed point | `a = val` (`t = M`) | tournament reflection §5 |
| difference-closure | one aligned gap `< val` | boxeph THM-1017 |
| stability | `v_max ≤ max(W)/(13δ)` | THM-1028 |
| covering/lcm | `182 ∣ v_max` (AP-core case) | boxeph THM-1017 |
| geometry | `M(V) ≥ v_max·M(W)/(v_max+v₂)` | THM-1002 |

The far element is simultaneously **the second coset**, **the antipodal edge-flip**, **the aligned
gap**, and **the `182`-multiple** — one object in every lens. The closing strategy is to keep composing
these (as stability × covering here) until the intersection of the necessary conditions is exactly the
deep-well tower; this THM is the first composition that produces a *finite* candidate set, the
prerequisite for any finite closure.



## UPDATE (same session) — the clean uniform lemma is FALSE; the true reason is alignment

I tried to close the compact case with a clean inequality: "valid non-AP core ⟹ `δ > max(W)/2366`"
(which would make the candidate window empty, eliminating the core by *stability alone*). **It is
false.** Near-AP valid cores violate it — `W = {1,…,11,24}` (covers 2..12, misses 13,14, non-AP) has
`M(W) = 2/25`, `δ = 0.00308`, `max/2366 = 0.0101 > δ`, so its window `[182,600] ∩ 182ℤ = {182,364,546}`
is *nonempty*.

But the finite check still holds — and shows the real mechanism: **`M(V) = M(W) = 2/25 ≥ 1/13` for
every candidate**, at the *same* witness `t = 2/25`. The mult-of-182 far element is **safe at `W`'s own
maximizer** (`‖182k·2/25‖ = 11/25, … ≥ 1/13`), so it does not lower `M`. In good-set terms: `W`'s
level-`1/13` good set sits at denominators tied to `W` (here 25), and the far element's arcs
(denominator `182k`) do **not align** with it, so the far element cannot cover it — hence `M(V)=M(W)`.

**This is the alignment wall again.** The AP `{1,…,12}` is the *unique* core whose good set (at
denominator 13) aligns with `182 = 14·13` — that alignment is exactly the deep well, the *only* place
the far element can cover and drop `M` below `1/13`. So "compact non-AP cores are eliminated" is **not**
a separable inequality lemma; it is the statement that only the AP's good set is `182`-aligned — the
inverse theorem itself, viewed from the good-set/alignment side. The finite check is verified on every
core tried, but its *uniform* reason is the wall, not a gap.

**Honest correction:** THM-1029's finite-candidate reduction stands (stability × covering ⟹ finitely
many candidates), but the hoped-for uniform closure of the compact case does not reduce to `δ >
max/2366` (false) — it reduces to the `182`-alignment rigidity, which is LRC(14). Scripts:
`lrc_gap_lemma_deathstar_S56.py` (the refutation), `lrc_check24_deathstar_S56.py` (the `M(V)=M(W)`
mechanism), `lrc_valid_core_deathstar_S56.py`.

→ THM-1028, THM-1017, THM-1002, THM-724/726, the Freiman-3k-4 / tournament-OCF reflections, HYP-3901.
