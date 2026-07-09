# The lemniscate as a second-moment cue for the good-period existence bound

*kind-pasteur-2026-07-09-S95. The owner offered `(x²+y²)² = x²−y²` as "a strange source of
inspiration" for the a-priori `μ>1` (dissociated good-period existence) bound. This is an honest
account: the lemniscate is not a tool here — using its metric changes the problem — but it points,
correctly, at the functional that governs the bound. Three cues, one real, two thematic.*

---

## The setup it speaks to

The residual on the LRC(14) critical path (LEM-013): for a dissociated 13-cluster `E ⊂ [0,s]` and
ruler `Vmax ∈ [s+1, ⌊7s/6⌋]`, show a good period exists — some dilation `j` with
`maxgap{e_i·j mod Vmax} > Vmax/7`. kps-S95 reduced this to an **averaging** statement: if the
**mean** largest-gap over dilations exceeds `Vmax/7`, some dilation beats it (`max ≥ mean`).
Verified adversarially: `avg_j[maxgap]·7/Vmax ≥ 1.047 > 1` across all spreads.

The question the lemniscate answers: *what lower-bounds the average largest gap?*

## Cue 1 (real) — the lemniscate is the second-moment curve, and the gap bound is a second moment

Write the lemniscate in polar form: `(x²+y²)² = x²−y²` becomes

> **`r² = cos 2θ`.**

The *squared* radius equals a *doubled* angle's cosine. It is the canonical curve relating a second
moment (`r²`) to a frequency-doubled coordinate. That is exactly the shape of the elementary lower
bound behind the averaging route:

> **`maxgap ≥ (Σ_i gapᵢ²) / (Σ_i gapᵢ) = (Σ gap²)/Vmax`**   (contraharmonic mean ≤ max).

The largest gap is controlled from below by the **second moment of the gap distribution**. Dissociation
(few additive relations, low additive energy) makes the gaps *uneven* across dilations, which makes
`Σ gap²` large — a big gap opens. This is the same `Σ gap²` second moment the lemniscate foregrounds.
The cue is genuine: the good-period existence is a second-moment phenomenon, not a first-moment one
(the first moment `Σ gap = Vmax` is constant and says nothing; the mean gap `Vmax/13 < Vmax/7` is
below threshold — only the *second* moment sees the unevenness that opens a `>1/7` gap).

*Honest limit:* the pure second-moment bound is **necessary but not sufficient** — kps-S95 measured
`avg_j[Σgap²/Vmax]·7/Vmax ≈ 0.85 < 1`, so `(Σgap²)/Vmax` alone does not reach `Vmax/7`; the true
`maxgap` (which the average does clear, `1.047`) carries the rest. The lemniscate names the right
object; it does not by itself close the gap. Both inequalities are now in Lean
(`LRCArcCountExistence.maxgap_ge_contraharmonic`, `averaging_existence`).

## Cue 2 (thematic) — the doubling `2θ` and the dilation-multiplication maps

`r² = cos 2θ` has an intrinsic **angle doubling**, and the lemniscate carries complex multiplication
by `ℤ[i]` (`sl(iz) = i·sl(z)`). The LRC good-period machinery runs on exactly such
*multiplication maps*: the dilation `j ↦ 2j`, and klein's **`×7` collapse** (LEM-012) that folds a
near-`Vmax/7`-grid onto one point. The lemniscate is the geometry where "multiply the angle" is
built into the defining equation — the same operation the proof uses to manufacture a gap.

A curious arithmetic coincidence worth recording (not a claim): by Abel's lemniscate theorem, the
lemniscate's `n`-division (into `n` equal arc-length pieces) is constructible iff
`n = 2^a · (distinct Fermat primes)` — the exact Gauss condition. The odd primes that qualify are
`3, 5, 17, 257, 65537`. **`7` is the first odd prime that does not.** LRC(14) lives at the threshold
`1/7`; its hardness is the resistance of the `1/7`-covering to elementary closure. That the
lemniscate likewise "resists" `7`-division is very likely a coincidence — but it is the kind of
coincidence this project has learned to write down (cf. [[triangle_foundation]]: the constants
`π, e, γ, √2` all fell out of the staircase; the lemniscate constant `ϖ = 2∫₀¹dr/√(1−r⁴) ≈ 2.622`
is the "elliptic `π`", the natural next constant if the triangle foundation ever acquires an elliptic
refinement).

## Cue 3 (thematic) — the origin crossing is a `ℤ₂` quotient

The lemniscate's self-crossing at the origin identifies the two antipodal circle points `θ=π/2` and
`θ=3π/2` (the "crossing shortcut": Euclidean distance `0`). That is a `ℤ₂` involution quotient — the
same *complement symmetry* the project factors out everywhere (the merged metagraph `G_n/ℤ₂` is the
PRIMARY object; CLAUDE.md). For LRC the matching involution is the gap reflection `φ ↦ 1−φ` /
co-offset complement. The lemniscate offers a geometric realization of "antipodes merge," but it does
not (that I can see) simplify the gap count — the LRC quotient is already taken in the reduction.

## Verdict

The lemniscate is **inspiration, not tool**: its arc-length metric makes runner speeds
time-dependent (elliptic integrals), which changes the problem, and its Euclidean embedding breaks
the circular metric LRC depends on. But as a *cue* it earns its keep — it correctly identifies the
**second moment of the gap distribution** as the functional governing good-period existence (Cue 1,
now formalized), and it rhymes with the proof's **multiplication-map** and **`ℤ₂`-quotient** themes.
The actionable takeaway sits in Lean: existence follows from `mean maxgap > Vmax/7`, whose engine is
the second-moment lower bound `maxgap ≥ Σgap²/Vmax`. The lemniscate told us to look at `Σgap²`.

*Files: `04-computation/lrc14_avgmaxgap_apriori_kps_S95.py` (+ `.out`);
`lean/…/LRCArcCountExistence.lean` (`averaging_existence`, `maxgap_ge_contraharmonic`).
See LEM-013 (the existence margin), mac-mini-S62 / opus-S168 (arc-count, the "broken clock"),
[[triangle_foundation]].*
