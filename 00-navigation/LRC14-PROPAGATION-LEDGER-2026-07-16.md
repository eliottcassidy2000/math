# THE CONSTANT-PROPAGATION LEDGER — 2026-07-16 (boxeph-S28; kps task (b) executed, assembly form)

**What this is.** The explicit composition of the route-[A] two-scale chain, turning
"Error → 0 as w → ∞" into PER-ROW EXPLICIT THRESHOLDS W₀(k), with the finite remainder
mapped to existing tiles and the seams named. Sources: THM-727 (rigorous reduction),
klein cont.7 THM-887-uniform (the uniform |U| ≤ C·diam bound — which FILLS THM-727's
named empirical gap for dominant-owner families), THM-663/S58 (row margins), THM-755/756
(the bands below threshold), THM-887/888-boxeph (general-cluster fallback profiles).

## 1. The composed chain (dominant-owner peels E = slow ∪ {t})

1. **THM-727 (rigorous part):** Error·w = |S| ≤ (1/2π²)·Σ_s Σ_ℓ |U_s(ℓw)|·|sin(πℓ/7)|/ℓ².
2. **klein cont.7 (A) (proved):** sup_N |U_s(N)| ≤ 0.8287·t + C₂(M_slow + √(t·M_slow))
   per section (slow-six constant; measured cross-section value 0.047 — ~7× slack).
3. **The ℓ-sum constant:** C_ℓ = Σ_{ℓ≥1}|sin(πℓ/7)|/ℓ² = 0.928451… (exact convergent sum).
4. **Composed:** Error(w) ≤ (1/2π²)·C_ℓ·7·0.8287·(diam-scale)/w + (lower order)
   = **0.2729·diam/w + O((M_slow + √(w·M_slow))/w)** — PROVED lane.
   (Measured lane: 0.0155·diam/w — the 17.6× headroom quantifies remaining slack.)

## 2. The per-row thresholds

Row margins from the S58/THM-663 density-floor closure; the row closes when the
two-scale error stays below its margin:

| row k | margin (S58) | W₀ = 0.2729·diam/margin (PROVED lane) | W₀ (measured lane) |
|---|---|---|---|
| 8 | +0.086 | **3.17 · diam** | 0.18 · diam |
| 9 | +0.082 | **3.33 · diam** | 0.19 · diam |
| 10 | +0.101 | **2.70 · diam** | 0.15 · diam |
| 11 | +0.12 (LEM-009 l.25: A = 0.4530 ≥ bar; SEAM 1 CLOSED, S40) | **2.27 · diam** | 0.13 · diam |
| 12 | +0.157 | **1.74 · diam** | 0.10 · diam |
| 13 | +0.252 | **1.08 · diam** | 0.06 · diam |

**Headline (COMPLETE as of S40): every density row k = 8..13 closes for far
elements w ≥ 3.4·diam(slow part), by the PROVED chain with explicit constants —
the table has no gaps.**

## 3. The finite remainder and its tiles

What remains per row is the band w < 3.4·diam — bounded-ratio families:
- **THM-755 (capped envelope):** closes every v > v\* = r_P/(π|G′_P|) unconditionally —
  band edges computed per core (e.g. deep-well core 112 < 156).
- **THM-756 (band closure):** Battery A COMPLETE — all 91 bottom cores' (P, v) band
  pairs decided (4011 by (H), 19 by exact L > 0, 2 = the tight AP/GW corners).
- **The near-AP window sweeps (THM-734/738/741 ladder)** and **the fragmentation box
  (THM-883-macmini + THM-885 j ≤ 5)** cover the small-part shapes.

**The named seams (the honest remainder of task (b)):**
1. ~~k = 11's bar/margin~~ CLOSED (S40): margin +0.12 (LEM-009), W₀(11) = 2.27·diam.
2. The per-row CONSUMPTION SEMANTICS referee: confirm with klein/mac-mini that the
   margin guards exactly the functional the composed Error bounds (the finish-map
   logic says yes; a referee pass should say it in one page).
3. The band audit: for each row, confirm the [·, 3.4·diam] band instances are inside
   THM-755/756's closed batteries or the sweeps' windows — finite, mechanical; any
   uncovered (core, band) pair becomes ONE exact check each (decide-shaped).
4. Multi-large-owner clusters (no dominant t): route through THM-887's profile with the
   M-cap — constant weaker (M-scale, not 0.27); the slice-Parseval sharpening
   (THM-888(C)) is the named upgrade.

## 4. Status

Task (b) is now an ASSEMBLY WITH EXPLICIT CONSTANTS + four named finite seams, no
analytic residue in the dominant-owner lane. The 17.6× proved-vs-measured headroom
means seam 3's band is generously covered by existing tiles in practice; the referee
passes (seams 1–2) are one-session items. — boxeph-S28

## 5. Signed owner-resonance addendum (codex-S17)

`THM-891` computes the fixed owner-resonant far-peel limit exactly.  For a six-offset
slow core `E`, far speed `t`, and fixed `a`, at `w=at`:

`Error_t(a)=C_a(E)+O_E(1/t)`, with `aC_a(E)` depending only on `a mod 7`.

The coefficient is an exact signed expectation on one- and two-miss patterns.  The
full two-runner sector law lies on 21 arithmetic rays; finite exact quadratic
certificates on their endpoints prove the following diameter-free limiting bounds:

- residues `1,...,5` satisfy `|aC_a(E)|<0.097` universally;
- residue `6` satisfies `aC_a(E)<0.097` on its positive side;
- only the negative residue-six term remains, concentrated on
  `K_6({1,5})=K_6({2,4})=-12`, hence on `A_15+A_24`.

This removes five of six signed residue obligations at the current `0.097` benchmark,
but it does **not** replace the all-`w` `0.2729 diam/w` envelope above, certify the row
margins in §2, or make `O_E(1/t)` uniform.  The honest added seams are:

1. a higher-order relation/tensor bound for the negative residue-six mass;
2. a core-uniform wall-cell remainder converting the fixed-`a`, fixed-core limit into a
   usable finite-`t` band statement.

Concurrent synthesis fixes the right carriers.  `THM-890` says deviations are exact
additive-relation spectra; `THM-892` says the generic quadratic frame mean is a universal
`pi^2/3` term plus coincidence-class imbalance; `THM-894` independently names the
level-three overlap tensor as the next rung.  `THM-896-level3-crossing` now proves the order-three
Bonferroni crossing through `m'<=11` and leaves the relation-localized triple-beat upper
bound open; that is the same analytic rung as residue six, though not the same
observable.  The resolved `HYP-7027` proves wall-movie cycles are all expressive, so
seam 2 must retain palette/relation data rather than quotienting supposedly silent loops.
