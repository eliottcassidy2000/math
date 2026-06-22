# The 3-set Venn Bonferroni-truncates the LRC wide bound at third order

*kind-pasteur-2026-06-22-S31t. Applying the owner's corrected 3-set-Venn understanding of the
recursion modes to the actual LRC(14) wide-bound proof. The Venn makes the multi-far cover a Newton
inclusion–exclusion whose 4-far-and-higher tail is negative — so the wide bound truncates at the
doublet + triple, both already controlled.*

## The setup the correction gives us
The owner's correction (S31s): the Legendre recursion is the **3-set inclusion–exclusion** over three
generators, and it **is** the LRC three-far coverage recursion (THM-548) — the three sub-tilings are the
three far runners, corners/edges/center are one/two/three-far Newton packets `T_r`, with
`p0(B∪far) = Σ_r T_r`. So the wide cover is literally a Newton/Möbius expansion in the far runners.

## What the expansion looks like (verified)
`lrc_venn_bonferroni_wide_kps.py`, `lrc_rfar_convergent_tail_kps.py` (base + up to 5 far):

| `r`-far | role | sign | magnitude (genuine-wide) |
|---|---|---|---|
| `T_1` | corner (single far) | `+` | **0** — a lone far runner cannot fill 6 sectors |
| `T_2` | **edge = A∩B (doublet)** | `+` | **DOMINANT / binding** (HYP-2797) |
| `T_3` | center (triple) | `+` | sub-dominant, `~0.4–0.5·T_2` (THM-557) |
| `T_4` | 4-far | **`−`** | `< T_3` |
| `T_5` | 5-far | `+` | `< |T_4|` |

Two structural facts, both verified across configs:
1. **`T_1 = 0` for genuine-wide** — coverage is irreducibly multi-far; the *binding* term is the
   **2-far doublet**, which is exactly the **Eisenstein `A∩B` edge** of the Venn (the even/degenerate
   mode). This is why HYP-2797's doublet is the maximizer: it is the largest Venn region.
2. **The tail alternates and decreases from `T_3`** (`T_3>0, T_4<0, T_5>0`, `|T_3|>|T_4|>|T_5|`), so the
   **`r≥4` tail is negative**.

## The consequence: a Bonferroni-3 upper bound
Because the `r≥4` tail is negative (alternating, decreasing), truncating the Newton series after third
order is an **upper bound**:

> **`p0(B∪far) ≤ T_1 + T_2 + T_3`** (the 1-far + doublet + triple).

VERIFIED on every genuine-wide config: `T_2+T_3 ≥ p0(full)` (`0.355 ≥ 0.293`; `0.372 ≥ 0.313`), and
`Σ_{r≤3} T_r < cap` (`0.424 < cap_9 = 0.494`). So the **infinite multi-far wide bound reduces to a
finite third-order truncation**: bound `T_2` (the doublet — CLOSED, THM-563/HYP-2797 periodic-Dedekind)
and `T_3` (the triple — sub-dominant, THM-557 / the R-tail HYP-2817 Mordell–Tornheim `12ζ(3)`), and the
4-far-and-higher corrections only *help* (they are negative).

## Why the Venn guarantees it
The center is the triple overlap `A∩B∩D ⊆ A∩B`, the edge — geometric containment forces
`|T_3| ≤ |T_2|`, and iterating (each higher packet sits inside the previous overlaps) forces the
magnitude decrease that makes the tail alternating-convergent. The **Möbius signs of the inclusion–
exclusion are exactly the Bonferroni signs**, so the partial sums bracket `p0` and the *odd*-order
truncation (`r ≤ 3`) is the upper bound. The owner's "even = degenerate (no lone `D`, no triple `G`)"
is the same statement seen from the doublet: the binding configuration lives on the even/edge mode, and
the odd/center mode is the first *correction*, not the leading term.

## Status and what it buys the proof
- **What's rigorous:** the Newton expansion `p0 = Σ T_r` (exact, inclusion–exclusion); `T_1 = 0` for
  genuine-wide (a single far runner leaves ≥1 sector uncovered, finite-checkable); the Bonferroni
  bracketing *given* the sign/decrease pattern.
- **What's verified-not-yet-proved:** the uniform sign+decrease `T_3 > |T_4| > |T_5| > …` (holds on all
  sampled genuine-wide configs; it is the coverage face of the repo's Dedekind ladder THM-557/HYP-2797,
  where the `r`-far packet is an `r`-fold Dedekind sum that decays). Proving that decay *uniformly* is
  the remaining analytic step — and it is now a **clean, finite, third-order target**: show
  `T_2 + T_3 ≤ cap` and `T_{≥4} ≤ 0`, instead of an unbounded multi-far sum.
- **Net:** the corrected Venn understanding turns the wide bound from "control an infinite tower of
  far-runner interactions" into "**the doublet binds, the triple corrects, and everything past third
  order subtracts**" — a Bonferroni-3 truncation that matches and structurally explains the repo's
  doublet-plus-convergent-tail route.

→ THM-548 (one/two/three-far), HYP-2797 (binding doublet), THM-557 (three-far sub-dominant), HYP-2817
(R-tail `12ζ(3)`), THM-563 (single-far periodic-Dedekind),
`the-three-modes-are-parity-stratified-and-lrc14-is-eisenstein-compose-legendre.md`, [[lrc14-thread]].
