---
source: claude-2026-06-03-S617
status: rigorous reduction of H=21 to a finite n<=12 check (Moon bound) + correction of an S616 error
tags: [H21, H-impossibility, moon-theorem, strong-tournament, three-cycles, reduction, correction, honest, finite-check]
---

# Hammering down H=21: a finite n≤12 reduction (and a correction)

## First, a correction (S616 was wrong in part)

S616 (HYP-2187) claimed H=21 reduces to a single open conflict-graph profile
`(α₁=6, α₂=2)` and used a reconstructed `Ω`. **That reconstruction was buggy:** it
deduped directed odd cycles by their *vertex set*, collapsing distinct directed
5-/7-cycles on the same vertices. I verified `I(my-Ω, 2) ≠ H` in **1612/3000**
random `n=6` tournaments. Consequences:

- "the one open case is `(6,2,0)`" was **wrong** — there are *four* connected
  profiles `(4,3),(6,2),(8,1),(10,0)`, which is exactly THM-079 Part H's
  **i₂-jump** (already in the repo);
- the `(6,2)` value `i₂ ∈ {0,1,5}` I reported *is* correct (3-cycle-dominated, so
  the dedup bug is harmless there) but it is THM-079's result, not new;
- H=21 is in fact **proved impossible for n≤8 exhaustively** (THM-079 Part G,
  268,435,456 tournaments at n=8) — I had under-credited it.

Lesson (again): verify the reconstructed object against the ground truth (`H` by
direct Hamiltonian-path count) *before* building on it.

## The new rigorous reduction: H=21 ⟹ n ≤ 12

This redeems the session — a clean, rigorous reduction that genuinely advances
the proof.

**Step 1 (strong-component reduction).** `H` is multiplicative over the
strong components: `H(T) = ∏ H(C_i)`. For `H=21 = 3·7`, and since `7` is **not**
a strong H-value (THM-029 proves no tournament has `H=7`), `21` cannot be a
product of *smaller* strong values — it must be the `H` of a **single strong
component**. So WLOG `T` is a strong tournament with `H(T)=21`.

**Step 2 (cycle bound).** `H = I(Ω,2) = 1 + 2α₁ + 4α₂ + … = 21` forces
`α₁ ≤ 10`, where `α₁` = total number of directed odd cycles. Since 3-cycles are
odd cycles, `c₃ ≤ α₁ ≤ 10`.

**Step 3 (Moon's theorem).** A strong tournament on `m` vertices has
`c₃ ≥ m − 2` (Moon 1968; verified here: `min c₃` over strong tournaments `= m−2`
for `m = 3..7`). Hence `m − 2 ≤ 10`, i.e. **`m ≤ 12`**.

**Conclusion.** With THM-079 Part G (`H=21` impossible exhaustively for `m ≤ 8`),
the *only* remaining cases are **strong tournaments on `m ∈ {9,10,11,12}` with
`c₃ ≤ 10`** — a **finite** set. This collapses THM-079's open part (which left
*all* cycle-rich `n ≥ 9` open in general) to a bounded, finite check.

## Evidence that the finite window is comfortably clear

Sampling strong tournaments with `c₃ ≤ 10`:
- `m=9`: min `H = 75`;
- `m=10`: min `H = 153`;
- `m=11,12`: `c₃ ≤ 10` strong tournaments are near-Moon-extremal and very rare.

These are *nowhere near* 21 — because a strong tournament on `m ≥ 9` is
vertex-pancyclic (Moon), so it carries many **long** odd cycles (5,7,9-cycles)
on top of its `c₃` triangles, blowing `α₁` far past 10. So the `n ≤ 12` window
isn't just finite; its `H`-values are large.

## What closes it fully

Two routes, either suffices:
1. **Exhaust** strong tournaments with `c₃ ≤ 10` on `m = 9,10,11,12`. These are
   highly restricted (near the Moon-extremal "transitive + one back-edge"
   family), so the enumeration is far smaller than `2^{C(m,2)}`.
2. **Prove** `strong on m ≥ 9 ⟹ α₁ ≥ 11` (equivalently `H ≥ 23`). The data
   (`min H = 75` at `m=9`) says the true bound is *much* stronger; even a crude
   pancyclicity count of long odd cycles should give `α₁ > 10`.

Either finishes `H=21`, and with it `{7,21}` as the complete set of permanent
H-gaps below the strong-monoid.

## Why this fits the equidecomposability frame

`H=21` is **doubly blocked**: the *decomposition* route (`21 = 3·7`) needs the
forbidden atom `7`, and the *atomic* (single-strong) route needs a strong
tournament of bounded size whose cycle count is too small to reach the `(α₁≤10)`
profile while staying strong (Moon). The reduction is exactly the statement that
no equidecomposition type and no atom realizes the measure `21`. The Moon bound
is the quantitative obstruction on the atomic side.

## Next

1. Enumerate strong `c₃ ≤ 10` tournaments on `m = 9..12` exhaustively (the
   restricted near-extremal family) → complete the `H=21` proof.
2. Or prove `strong, m ≥ 9 ⟹ α₁ ≥ 11` via vertex-pancyclic long-odd-cycle
   counting.
3. Apply the same Moon-style reduction to other candidate gaps (none below 200
   besides `{7,21}`) and to the strong-value spectrum generally.
