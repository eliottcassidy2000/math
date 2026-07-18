# The soft Weyl bound closes the density route off-resonance — but its natural peel is maximally resonant, a Lonely Runner one level up

*boxeph-2026-07-18-S96. Aiming the soft Weyl bound at the density route (Route A), as directed. Outcome:
the soft bound genuinely closes the **off-resonance** part (my THM-886(II), proved), but the route's
natural far-element peel is **maximally resonant** (`Q_s = Θ(r²)`, no power-saving), so the density
route bottoms on the SAME resonance/cancellation as covering — confirming klein-S279's unification from
the density side, and correcting the finish-map's (pre-refutation) "any power-saving suffices." LRC(14)
not closed.*

## The target, and the tool

Density route (THM-729): `Q_s = Σ_ℓ |U_s(ℓw)|²/ℓ² = (2πw)²·[Riemann-discrepancy of the arc-set
autocorrelation at the w-grid]`; `Error = |S|/w`, `|S| = O(√Q_s)`. Closing needs `Q_s = o(r²)`
(any power-saving over crude `r²`) so `Error → 0`. The off-diagonal is the Weyl sum
`Σ_i w_i·e(−ℓw c_i)` over arc midpoints `c_i` (Farey fractions `a/v`, `v∈E`) — a soft equidistribution
statement, the *right* marriage for Weyl (S95).

## The soft Weyl bound closes the OFF-resonance part (proved)

My THM-886(II) is exactly the soft Weyl/geometric bound here. For the endpoint sum `S(m)` on `ℤ_P`,
each per-class active index set is a union of `R` runs, and a run's exponential sum is geometric:

> **`|S(m)| ≤ C_f + R/sin(π‖m/t‖)`** (proved, machine-refereed, 0 violations).

For `w` (i.e. `m = ℓw`) **not** near the mode lattice `tℤ`, `sin(π‖m/t‖)` is bounded below, so
`|S(m)| = O(R) = O(r)`, hence `Q_s = O(r)` and `Error → 0`. **The soft Weyl bound does its job — off
resonance the density row closes.** This is the concrete sense in which density is the "softer" side.

## But the natural peel is maximally resonant (`Q_s = Θ(r²)`)

The density route peels the far element `w = d` (its target families are `{≤6 bounded speeds} ∪ {d}`,
the unbounded-7th-speed case). But `d` **is** the mode-lattice generator: `ℓ·d ∈ dℤ` for **every** `ℓ`,
so the peel sits on the resonance at every harmonic. THM-886(IV/V):

- **The uniform bound is FALSE:** `max_w |S(w)|²/M` grows (14 → 226 across `t = 6 → 1200`), so HYP-6994
  is refuted (confirmed by klein-S314). `sup_w Q_s ≍ r²` — **sharp, no power-saving.**
- **At the resonant modes the teeth are full-strength (no cancellation):** THM-886(I),
  `S(ta) = S_frame(ta) + Σ_{c mod 7} e(ac/7)·N_c`, with `Σ_c |N_c| ≍ r` and **no `a`-decay** — the mode
  comb has full teeth at every `ta`, `a≢0 (mod 7)`. So the resonant contribution does not cancel; it is
  genuinely `Θ(r)` in `|S|`, `Θ(r²)` in `Q_s`.

So `Q_s = o(r²)` is **false uniformly**, and it fails precisely at the density route's own peel.

## The density route bottoms on the SAME resonance as covering — a Lonely Runner one level up

The classifier for the blow-up is **Diophantine, not arithmetic** (THM-886(V)): `w` is resonant ⟺ a
small multiple `ℓw` lands near the mode lattice `tℤ` — **`w` plays a Lonely Runner against the mode
lattice.** The far-element peel `w = d` is the extreme case (it generates the lattice), so it is
maximally resonant. Therefore:

> **The density route's remaining core is not a soft power-saving — it is the resonant far-element peel,
> a self-similar Lonely-Runner resonance one level up.** The soft Weyl bound removes everything except
> this resonance; the resonance is the *same* non-cancelling object the covering route bottoms on
> (klein-S279, "both routes bottom on the same multilinear cancellation"). This session confirms the
> unification **from the density side**, concretely.

**Correction to the finish-map (2026-07-13):** its "`Q_s = O(diam)` verified / any power-saving
suffices" is the *pre-refutation* optimistic state. klein-S314 and boxeph-S25 (both after 07-13)
refuted the uniform bound. Density is softer only *off* resonance; at the core it is the same crux.

## Net (honest)

- **Delivered:** the soft Weyl bound, aimed at the density route, closes it **off resonance**
  (`Q_s = O(r)` via THM-886(II), proved) — the right tool for the right (soft) part.
- **The obstruction:** the natural peel `w = d` is maximally resonant (`Q_s = Θ(r²)`, THM-886(V), sharp;
  resonant teeth full-strength, THM-886(I)); `o(r²)` is false uniformly. The soft Weyl bound cannot
  resolve the resonant peel.
- **The structure:** the density-route core is a Lonely Runner one level up (`w` vs the mode lattice
  `tℤ`) — the same resonance as the covering route. Both routes bottom on the resonance;
  self-similarity, not tractability, is what the density side reveals at the core.

LRC(14) is not closed. The soft Weyl bound is placed where it belongs and does its job (off-resonance);
the irreducible core is the resonant far-element peel, unified with the covering crux and self-similar.

Cross-links: [[weyl-is-the-wrong-tool-the-one-line-form-is-concentration-not-equidistribution-boxeph-S95]],
THM-886 (resonance law of Q_s), [[both-lrc14-routes-bottom-on-the-same-multilinear-cancellation-not-one-mollification-klein-S279]],
THM-729, [[density-needs-any-power-saving-not-the-sharp-bound-and-the-large-sieve-is-worse-klein-S281]].
