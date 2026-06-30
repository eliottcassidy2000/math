# The descent's finite families: the WHOLE apex-7 floor is a 9-row table (cores of Z₇ up to the affine group AGL(1,7)), 8 rows clear the floor 4cos²(3π/7) by THM-590, and the entire LRC(14) wall is the ONE remaining row — that the descended odd-7 core is never the full Z₇ (the unique gap-0 boundary = the doublet's mirror, the disproof)

*opus-2026-06-30. Owner: work on the descent, come to know the finite families. The descent (mac-mini
HYP-3575) reduces the per-level floor `ρ_j` to the cyclotomic gap `g(O)=min_{k≠0}|Σ_{x∈O}ω^{kx}|²` of the
descended odd-7 core `O⊂Z₇`; THM-590 (PROVED) bounds `g`. Coming to know the cores collapses the whole
thing to ONE open row.*

## The finite families: 9 essential cores (the whole descent is a 9-row table)
The gap `g(O)` is invariant under the affine group `AGL(1,7)=Z₇⋊Z₇*` (order 42; translation and
multiplier leave `|Σω^{kx}|²` and its min unchanged). So the `2⁷=128` cores collapse to **9 orbits**:
| `|O|` | `g(O)` | orbit size | rep | role |
|---|---|---|---|---|
| **7** | **0** | 1 | `Z₇` | **BOUNDARY = the disproof (the only gap-0 core)** |
| 2 | `0.19806` | 21 | `{0,1}` | **THE FLOOR — the doublet (= my binding pair)** |
| 5 | `0.19806` | 21 | `{0,1,2,3,4}` | the quintuplet (doublet's complement) |
| 3 | `0.30798` | 21 | `{0,1,2}` | interval |
| 4 | `0.30798` | 21 | `{0,1,2,3}` | interval complement |
| 1 | `1` | 7 | `{0}` | singleton |
| 6 | `1` | 7 | `{0,…,5}` | singleton complement |
| 3 | `2` | 14 | `{0,1,3}` | **OPTIMAL — Fano/perfect-difference-set = octonion** |
| 4 | `2` | 14 | `{0,1,2,4}` | Fano complement |
- **5 cyclotomic values** `{0, 4cos²(3π/7), 0.30798, 1, 2}`, all in `Q(cos 2π/7)`; **complement-symmetric**
  `g(O)=g(Z₇\O)` (sizes pair `2↔5, 3↔4, 1↔6, 0↔7`).
- **THM-590 PROVES `g(O) ≥ 4cos²(3π/7)=0.198` for all 8 PROPER orbits** (finite, exact). The minimum is
  the **doublet** (and its mirror quintuplet); the maximum-gap (best-conditioned) is the **Fano/octonion**
  perfect-difference set (flat spectrum `λ_k=2`). `g=0` happens at **exactly one orbit: `O=Z₇`** (the full
  sum `1+ω+…+ω⁶=0`, by `Φ₇` irreducibility).

## The descent's binding core IS my binding pair
My razor's-edge binding pairs mod 14 descend to mod-7 doublets, all at the floor:
> `{1,13}→{1,6}`, `{5,9}→{2,5}`, `{3,11}→{3,4}` (mod 7), each `g=0.198`.
So **the descent's worst-case (binding) core = my binding pair = THM-578's doublet = the cyclotomic atom
= the genus-1 cusp form's support.** Five independent threads (binding pair, doublet, atom, cusp form,
3-cycle mirror) are the SAME size-2 core. The floor binds where the binding pair binds.

## The wall, reduced to ONE row
> **127 of the 128 cores (all 8 proper orbits) already clear the floor `0.198` — PROVED, finite (THM-590).
> The entire remaining LRC(14) wall is the single row `O=Z₇`:** the descended odd-7 core must never be the
> FULL `Z₇`. `O=Z₇` is the unique gap-0 orbit = the disproof boundary = "all residues active / fully
> resonant / no loneliness." The conjecture is exactly: **the descent never lands on `Z₇`.**
This is the cleanest the wall has ever been: not an inequality over an infinite family, not a measure
(which vanishes), not a sign (favorable, rank-0) — but a SINGLE FINITE EXCLUSION: the cusp (`m_R→0`) does
not degenerate the binding core to all of `Z₇`. The 2-adic descent (`14=2·7`, THM-580, the doubling-2
that `CV(H)~2/n` rehearses) peels to the odd-7 core; the wall is that this core stays proper.

## What "knowing the families" buys the descent
1. **The floor is a 9-row finite check, 8 rows done.** The proof obligation is one row (`O≠Z₇`), not an
   analytic bound — a finiteness the descent's transitivity (Z₇* spectrum-invariance) delivers.
2. **The binding core is the doublet (size 2), not `Z₇` (size 7).** The worst PROPER case is the doublet
   (`0.198`); `Z₇` is not a "binding" core, it is the degenerate all-active boundary. So the wall is a
   NON-DEGENERATION: at the cusp, the binding core has size 2, not 7.
3. **The Fano/octonion is the OPTIMAL row** (`g=2`) — the apex-7 octonion (HYP-3547) lives here as the
   flat-spectrum best case, the opposite end from the doublet floor.
4. **The exclusion is structural:** `g(O)=0` requires `Σ_{x∈O}ω^{kx}=0`, which by `Φ₇` irreducibility
   forces `O=Z₇` exactly — there is no "near-`Z₇`" core; the gap jumps from `0.198` (any proper) to `0`
   (only the full set). The floor is robust to within one discrete step of the boundary.

## Status
- **Computed/known (opus):** the 9 essential cores (AGL(1,7)-orbits), 5 values, complement symmetry, the
  floor=doublet / optimal=Fano / boundary=Z₇ structure; my binding pairs ≡ the binding doublets.
- **Proved (klein THM-590, finite):** 8 of 9 orbits clear `0.198`; `g=0` only at `Z₇` (by `Φ₇` irred.).
- **The wall, sharpened:** ONE row — the descended odd-7 core is never `Z₇` (cusp non-degeneration);
  conditional on mac-mini's `ρ_j=g(O)` descent reduction.
- **Open (the descent map):** characterize the descended core `O_j(S)` for covering `S` at the binding
  cusp and show `O_j ≠ Z₇` — equivalently the 2-adic-peeled odd-7 part of the binding structure is a
  proper subset (a doublet at the floor), never fully resonant.

Related: the f₁₄/14a-rank-0 reflection (the sign is favorable; the wall is the descent); the
roots-of-unity convergence; klein THM-590/HYP-3581/3585/3597 (the core landscape, finite families, the
infinite-family measure→0), mac-mini HYP-3575 (`ρ_j`=Z₇ Gram gap, the descent), THM-578 (doublet),
THM-580 (2-adic descent), my Dirac-comb + binding-pair reflections, OPEN-Q-108.
