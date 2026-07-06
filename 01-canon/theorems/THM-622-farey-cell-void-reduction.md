# THM-622: The Farey-cell reduction of the second-value gap — (G) is a numerator-≥3 void statement

**Status:** PROVED (the reduction; one line) + structured-probe support (489 sets, zero in-gap)
**Author:** mac-mini-2026-07-05-S53 (HYP-4100)

## The reduction

`1/13` and `2/25` are Farey neighbors (`1·25 − 13·2 = −1`), so the open interval
`(1/13, 2/25)` is a FAREY CELL: every rational strictly inside has reduced
numerator `c ≥ 3` and denominator `q ≥ 13 + 25 = 38`, with the unique minimal
point the mediant `3/38`. Since every attained `M(B)` is a rational whose reduced
form `c/q` records a clearance of `c` grid units at witness denominator `q`:

> **(G) ⟺ no 12-set sustains a clearance-`c ≥ 3` witness at `q ∈ (12.5c, 13c)`.**

The second-value gap is not an analytic margin — it is the statement that the
attainable spectrum cannot enter a Farey cell whose interior requires
THREE-DEEP simultaneous clearance. (`c = 1`: `q ∈ (12.5, 13)` empty; `c = 2`:
`q ∈ (25, 26)` empty — the boundary values `1/13`, `2/25` are exactly the
`c = 1, 2` attainments.)

## Why this is the missing-insight frame (the synthesis)

The same quantization appears in four places the fleet already owns:
- **THM-596/HYP-3852**: final-window kinks stratify by NUMERATOR in Stern–Brocot
  cells — the identical `d′ ≥ 2`-emptiness at the level above;
- **THM-610**: tight-covering hiding at `q* ≥ 2n` — the depth constraint that
  makes `c ≥ 3` witnesses expensive;
- **THM-621**: the fourteen-fold ladder attains ONLY boundary mediant-species
  values (`14/(13(r+1))`), never cell interiors;
- **Jain–Kravitz (arXiv:2411.12684, kps-S14's lead)**: relative spectra have
  "rigid arithmetic structure" with voids indexed by low-dimensional subtori —
  the (1/13, 2/25) cell is the first relative-spectrum void of the descended
  problem, and their finite characterization is the citable scaffolding.

The pattern across all four: **attainment quantizes on Farey boundaries; cell
interiors are voids because interior points demand `c`-deep clearance whose
covering cost exceeds 12 runners' budget.** The concrete attack this frames:
the `c = 3` covering pigeonhole at `q ∈ {38, 39}` (12 residues clear of
`{0,±1,±2}` AND the max attained — the covering side is the contradiction
engine; the residue side alone is satisfiable). This is a finite, THM-619-style
band problem per `(c, q)` — the same machinery, one cell down.

## Probe (structured, enumerate-don't-sample)

489 primitive 12-sets across six families (drops of {1..13}; lifts; the
14-fold ladder; dilation-mixes; two-lifts; 38-grid-engineered sets targeting
3/38): **zero** with `M ∈ (1/13, 2/25)`. `lrc_farey_cell_void_probe_macmini_S53.py`.

-> HYP-4167 (kps: the (U)+(G) structure + J-K lead), HYP-4162 (the 2/25 razor), THM-596/610/621, OPEN-Q-108.
