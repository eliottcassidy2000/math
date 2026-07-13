---
source: opus-2026-07-11-S255
status: CORRECTION + PROOF-AT-THE-TIGHT-POINT. Working S254's target M_core >= (182+s)/2379 shows the
  UNCONDITIONAL rigidity is FALSE (a covering counterexample: the dilated core 20*{1..11}u{143}), so S254's
  single-lemma reduction was over-stated. The escape closes instead by a beta case-split (beta = killer
  clearance at the core optimum): EASY (beta>=14/183) trivial via LRC(13); beta=0 the simple rigidity (exact);
  0<beta the full beta-balance beta*s+182*M_core >= 14(182+s)/183 -- all verified, 0 real counterexamples. And
  the tight beta=0 case (the deep well) is now RIGOROUS via S252: M_core=1/13 => interval => s=1 => equality.
tags:
  - lrc14
  - covering-min
  - correction
  - beta-case-split
  - deep-well-minimality
  - rigidity
  - proved-tight-point
---

# The single rigidity is false; the escape closes by a β-split; the tight point is proved

**opus-2026-07-11-S255.** Owner: work on `M_core ≥ (182+s)/2379`. Doing so both **corrects** S254 (the
unconditional rigidity is false) and **proves** the statement at the one point that matters — the deep well.

## Correction: the unconditional rigidity is FALSE

S254 reduced the single-killer covering-min to `M_core ≥ (182+s)/2379`. This is **false** for general covering
cores. Counterexample: the dilated core `20·{1..11} ∪ {143} = {20,40,…,220,143}` is **covering** with `182`,
has `M_core = 1/12` and binding speed `s = 20`, yet `(182+20)/2379 = 0.0849 > 1/12` — the rigidity **fails**.
But the family is **safe**: `M = 1/12 ≥ 14/183`, attained at `t = 1/240` where the **killer is non-resonant**
(`β = ‖182/240‖ = 29/120 ≫ 14/183`). So S254 wrongly demanded the balance-rigidity of an **easy-case** family
that clears via a different witness. The single-lemma reduction was over-stated.

## The correct closure: a β case-split (β = killer clearance at the core optimum)

Let `β = ‖182·t₀‖` at the core optimum `t₀ = a/q`. All three cases verified, **0 real counterexamples**:

- **EASY — `β ≥ 14/183`.** `M(family) ≥ min(M_core, β) ≥ 14/183`, since `M_core ≥ 1/13 > 14/183` (LRC(13)) and
  `β ≥ 14/183`. **Trivial.** This is where the dilated/spread cores live. *(603/603.)*
- **`β = 0`** (killer *exactly* resonant, `q ∣ 182`). The simple rigidity `M_core ≥ (182+s)/2379` is the
  **exact** balance requirement. *(140/140.)*
- **`0 < β < 14/183`.** Needs the **full β-balance** `β·s + 182·M_core ≥ 14(182+s)/183` — the head-start `β`
  makes the simple rigidity too strong here (it fails 14/157, the β-balance holds 157/157).

So the escape does **not** reduce to one core rigidity; it is a three-case witness argument, and the simple
rigidity is only the `β = 0` slice.

## Proof at the tight point: the deep well, via S252

The `β = 0` case is exactly where the deep well lives (`q = 13`; `182 = 14·13`, so `β = 0` automatically), and
there the rigidity is **provable**. For `q = 13`, `M_core = m/13`, and `m/13 ≥ (182+s)/2379 ⟺ s ≤ 183m − 182`:

- **`m = 1` (`M_core = 1/13`, the LRC(13) floor):** requires `s ≤ 1`. And `M_core = 1/13 ⟹` the core **is** the
  interval `{1..12}` (up to residue-preserving shift) because **`n = 13` is prime** (S252: prime ⟹ one tight
  pattern, no doubling) `⟹ s = 1 ⟹ req = 183/2379 = 1/13 = M_core`. **Equality.** So the deep well
  `{1..12,182}` is the **unique** family attaining the bound, exactly at `14/183`. **Proved.**
- **`m ≥ 2`:** `s ≤ 183m − 182 ≥ 184` — vast slack (verified, no violations).
- **Other `q ∣ 182` (2,7,14,26,91,182):** `M_core = m/q ≥ 1/13` (LRC(13)) forces `m ≥ q/13`, so `M_core` sits
  well above `1/13` unless it *equals* `1/13` — which forces the interval and hence `q = 13`. So these never
  bind.

## Net (honest)

Working the target flipped it: the **unconditional rigidity is false** (the `d=20` covering counterexample),
so S254's single-lemma reduction is **corrected** to a **β case-split** — easy (trivial via LRC(13)), `β = 0`
(simple rigidity), `0 < β` (full β-balance) — all verified with zero real counterexamples. And the **tight
point is now rigorous**: the deep well is the unique minimizer at `14/183`, via S252's prime-13 uniqueness
(`M_core = 1/13 ⟹` interval `⟹ s = 1 ⟹` equality). The remaining lemma is the **full β-balance for
`0 < β < 14/183`** (verified, unproved) plus multi-killer families — a smaller and correctly-scoped target
than S254's.

This is the useful shape of the result: the escape closes, the extremizer (deep well) is pinned by a real
proof, and the honest remaining work is a single β-dependent inequality — not the false unconditional rigidity.

→ opus-S254 (the rigidity, corrected here), opus-S253 (the balance), opus-S252 (prime-13 uniqueness — proves
the tight point), klein S267 (14/183 covering-min, verified), LRC(≤13) (`M_core ≥ 1/13`). Files:
`lrc14_rigidity_correction_beta_split_opus_S255.py` (+`.out`).
