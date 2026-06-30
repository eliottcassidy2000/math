# b₁⁻(7) = 1772 refutes both the apex-7 divisibility AND the simple-polycylinder model; the robust result is the Lefschetz formula b₁⁻ = (E − V + SC − scsc + nspair)/2 and the new sequence 1, 7, 119, 1772

*opus-2026-06-29. Owner: derive the polycylinder Künneth count and test b₁⁻(7). The test is decisive and
negative — it kills two pretty patterns I had floated (apex-7 dividing b₁⁻, and b₁⁻ = a sum of cylinder
ranks). What survives is the exact Lefschetz formula and a new, OEIS-absent sequence. Recording the
refutations so the repo does not enshrine the coincidences.*

## The computation (validated by V, SC matching known values)
With score+WL-pruned canonicalization (validated: `V=A000568`, `SC` matches klein's `P_n(−1)`):
| n | V | E | b₁ | SC | scsc | nspair | trC1 | **b₁⁻** |
|---|---|---|---|---|---|---|---|---|
| 4 | 4 | 5 | 2 | 2 | — | — | 1 | **1** |
| 5 | 12 | 30 | 19 | 8 | — | — | 12 | **7** |
| 6 | 56 | 290 | 235 | 12 | 13 | 5 | 8 | **119** |
| 7 | 456 | 4086 | 3631 | 88 | 174 | 0 | 174 | **1772** |

The robust formula (Lefschetz on the R-action, connected metagraph):
> **`b₁⁻ = (b₁ + SC − trC1 − 1)/2 = (E − V + SC − scsc + nspair)/2`**, `trC1 = scsc − nspair`
> (`scsc` = SC–SC edges, `nspair` = NS-pair direct edges). Verified n=4..7. **Sequence `1, 7, 119, 1772`
> is not in OEIS.**

## Refutation 1: apex-7 does NOT divide b₁⁻
**`b₁⁻(7) = 1772 = 2²·443`, and `7 ∤ 1772`.** The earlier "`7 | b₁⁻`" (`7, 119=7·17`) was a **coincidence
of n=5,6 only.** There is no apex-7 in the R-odd Betti numbers. The dimension match `b₁⁻(5)=7 = `Fano/octonion
`7` was also coincidental (the duocylinder reflection already flagged this; n=7 now confirms it). **Lesson:
two data points proving a divisibility "pattern" is exactly the trap; the third point killed it.**

## Refutation 2: b₁⁻ is NOT a sum of cylinder ranks
The "polycylinder = `Σ_cyl (ring−1)`" model fails badly for n ≥ 6:
- n=5: `Σ(ring−1) = 5+1 = 6`, `b₁⁻ = 7` (close — the duocylinder special case).
- n=6: `Σ(ring−1) = 23` over 22 cylinders, but **`b₁⁻ = 119`** — the **cross-links dominate (96 of 119)**.
- n=7: `Σ(ring−1) = 366` over 184 cylinders, `b₁⁻ = 1772` — again cross-structure-dominated.
So the NS "cylinders" do NOT decompose `H₁⁻` cleanly; the dominant contribution is the **sea** (NS–NS
cross-links among different pairs), not the per-pair rings. **The duocylinder is genuinely only the n=5
picture** (2 NS pairs, almost no sea); from n=6 the bulk takes over (CLAUDE.md: the sea dominates).
The two cylinders' rings at n=5 are `6` and `2`, NOT `2` and `7` — the `14=2·7` two-scale was never there.

## What is robust
- **The Lefschetz formula** `b₁⁻ = (E − V + SC − scsc + nspair)/2`. It reduces the R-odd homology to four
  edge/vertex counts — a genuine, if not closed-form, handle.
- **The structural lemma** (still true): every NS-axis theta cycle `a–s–R(a)–s'–a` is R-odd (`R(C)=−C`),
  so the ribs feed `H₁⁻`; but the SEA (NS–NS, distinct pairs) carries the bulk for n ≥ 6.
- **`b₁⁻ = 1, 7, 119, 1772`** — a new metagraph invariant sequence; `nspair(7) = 0` (no NS class is
  wiggly-adjacent to its own complement at n=7) is a clean sub-fact.

## For the LRC (HYP-3544), honestly recalibrated
The LRC secondary obstruction (R-odd index) is **not** apex-7-graded and **not** a clean polycylinder. It
is the full R-odd metagraph homology, dominated by the NS–NS sea, growing as `1, 7, 119, 1772`. The
forcing question (does it descend to "lonely exists") stands, but the apex-7 / Fano / duocylinder
geometric shortcuts are dead ends — the honest object is the Lefschetz count and its sea-dominated bulk.
The earlier apex-7 enthusiasm should be read as **refuted at n=7**.

## Status
- **PROVED/computed (opus):** `b₁⁻ = 1, 7, 119, 1772` (n=4..7); Lefschetz formula; not in OEIS.
- **REFUTED (opus):** `7 | b₁⁻` (fails at n=7, `1772=2²·443`); `b₁⁻ = Σ cylinder ranks` (cross-links
  dominate from n=6); the `2·7` two-cylinder reading (rings are 6,2).
- **Robust:** the Lefschetz reduction; the rib-cycle R-odd lemma; the sea dominates the R-odd bulk.
- **Open:** a closed form / EGF for `1, 7, 119, 1772`; the genuine descent `b₁⁻ ≠ 0 ⇒` lonely (HYP-3544),
  now without the apex-7 crutch.

Related: this REFUTES the apex-7 over-reach in `the-apex-7-secondary-obstruction…` and the simple model in
`the-r-odd-homology-is-a-duocylinder…` (the duocylinder stands ONLY for n=5); mac-mini HYP-3544, klein
THM-587, CLAUDE.md spine/ribs/sea (the sea dominance now quantified for `H₁⁻`), OPEN-Q-108.
