---
source: opus-2026-07-09-S172
status: HONEST NEGATIVE (transport of the tournament OCF clean truncation to an a-priori resonant-sum
  bound FAILS) + the structural reason + an S171 correction + fleet convergence. Goal was to make the
  E_grid residual bound |R|<(6/7)^k a-priori by porting THM-076's clean Walsh-OCF truncation to the
  covering master law What=SUM prod g-hat (LEM-011). RESULT: it does NOT transfer to a crude absolute
  bound. (a) The per-COVERING sum SUM_n|What(n)| DIVERGES (shell ratios GROW: r=3/r=2 ~ 4-8) -- the
  OCF factorization is NOT absolutely summable. (b) The per-FREQUENCY sum SUM_m|What_1d(m)| (the real
  1-D Fourier coeff, = R_abs) IS < main (0.005-0.37) but is NOT provable by the rigorous BV bound
  |What_1d(m)|<=TV(W')/(2pi m)^2, because TV(W') ~ spread^2 (measured exponent 2.03, ~13*spread^2), so
  TV/(12 Vmax^2) ~ 1.1 (=8*main), USELESS -- the actual R_abs is 40-200x below by CANCELLATION. So
  R_abs<main is a genuine Mertens/L2 cancellation statement (converging kps-S98). THE STRUCTURAL REASON:
  the LRC covering over-covers (k/7=13/7>1) => ~26*spread order-swap breakpoints, each W' jump ~spread/2
  => TV~spread^2; the tournament OCF has NO over-covering, so IT truncates cleanly. The k/7>1 over-cover
  IS the transport obstruction (resolves the Agent-1 lead). NOT proof-critical: the fleet closes the good
  period by j=1 wraparound (non-strict, mac-mini-S64) + LEM-013 exhaustion + LEM-012, none needing |R|.
tags:
  - lrc14
  - good-period
  - egrid
  - walsh-ocf
  - mertens
  - honest-negative
  - over-covering
---

# The OCF transport fails because k/7 > 1

**opus-2026-07-09-S172.** Owner: work the transport of the OCF/decay truncation constant, then
formalize; pull often and use the fleet's work as signal.  The transport does not go through, for a
reason worth stating precisely.

## The goal and the two candidate objects

kps-S96's E_grid route closes the good period once `|R| < (6/7)^k`, `R = Σ_{Vmax|m} 𝒲̂₁(m)`,
`𝒲̂₁(m) = ∫₀¹ W e(−mx)dx`.  I hoped to make this a-priori by transporting the tournament Walsh-OCF
factorization (THM-076, which truncates cleanly) to LEM-011's covering master law
`𝒲̂(n) = (−1)^r(6/7)^{k−1−r}[∏b₀(n_i)](𝟙[σ=0]−c(σ))`.  Two absolute sums bound `|R|`:

- **per-covering** `Σ_n |𝒲̂(n)|` (the multi-dim OCF terms).  **DIVERGES.** Verified by shells (r = #
  nonzeros): `R_abs^{(r)}` *grows* with r (`r3/r2 ≈ 4–8`), so the per-coordinate geometric factor
  `(7/6)/π = 0.371` is beaten by the shell counts.  The OCF factorization is NOT absolutely summable
  (the L¹-divergence of opus-S154, seen term-by-term).  At small spread this exceeds `main` (7-str
  `1.7·main`, tight-AP `7.5·main`) — this is kps-S98's ">main".
- **per-frequency** `Σ_m |𝒲̂₁(m)| = R_abs` (the real 1-D coefficient, after collecting `n·e=m`).
  Converged, this is `< main` (`0.005–0.37`, `≈ |R_signed|`).  This is the right object; the collection
  `𝒲̂₁(m) = Σ_{n·e=m}𝒲̂(n)` is where the (free) cancellation lives.

## Why even the per-frequency bound is not a-priori: `TV(W') ~ spread²`

`R_abs < main` empirically, but to PROVE it a-priori the only rigorous handle is the BV bound
`|𝒲̂₁(m)| ≤ TV(W')/(2πm)²` (W continuous piecewise-linear), giving `R_abs ≤ TV(W')/(12 Vmax²)`.
The decisive measurement:

> **`TV(W') ~ spread^{2.03}` (`≈ 13·spread²`)**, so `TV/(12 Vmax²) ≈ 1.1` — a CONSTANT `≈ 8·main`,
> independent of spread.  The crude bound is useless; the actual `R_abs` (`0.005–0.03`) is `40–200×`
> below it.

My earlier hope — "only the `≤6` gaps exceeding `1/7` contribute to `W'`, so `TV = O(spread)`" — is
false: as `x` varies, the identity of those gaps swaps at `~26·spread` sorted-order breakpoints, each a
`W'` jump of `~spread/2`, so `TV ~ 13·spread²`.  Hence `R_abs < main` is **genuine phase cancellation**
among the `~spread²` jumps at the resonant frequencies — the Mertens/L2 wall (opus-S167, kps-S98), not
an absolute bound.  **This corrects my S171 claim** ("R_abs bounded absolutely, no cancellation"): the
per-frequency `R_abs` is `< main`, but *only by cancellation*; no crude absolute bound reaches it.

## The structural reason the transport fails: `k/7 > 1` (over-covering)

Why does the tournament OCF truncate cleanly while the covering one does not?  **The LRC covering
over-covers.**  Its `k = 13` arcs of length `θ = 1/7` have total length `k/7 = 13/7 > 1` — the circle
is covered `~1.86×` over.  Over-covering means `W = Σ(gap−1/7)_+` is usually `0` with rare deep
excursions, and the max-gap identity churns through `~26·spread` breakpoints ⟹ `TV(W') ~ spread²` ⟹
cancellation-dependence.  The tournament OCF (Ham-path covering) has no such over-covering, so its
collected Walsh coefficient has bounded variation and truncates cleanly (THM-076).  **`k/7 > 1` IS the
obstruction to transport** — exactly the "barely-covers wall" and the difference Agent-1's mining
flagged ("understand why their covering had no `k/7>1` obstruction").  Resolved: the over-covering is
the reason.

## Neither L1 nor L2 reaches it — the wall is arithmetic

The TV bound is an L¹ estimate.  The natural next move (my own S154 "L²-not-L¹") is Cauchy–Schwarz with
weight `m²`: `|R| ≤ 2√(Σ_{V|m}|m·𝒲̂₁(m)|²)·(Σ 1/(nV)²)^{1/2} ≤ 2√(E[(W')²]·ζ(2))/Vmax` (Parseval:
`Σ|m·𝒲̂₁(m)|² = ‖W'‖₂²`).  Tested:

> **`E[(W')²] ~ spread²`** (`√(E[(W')²])/spread ≈ 0.7–1.0`, constant), so **`B_L2 ≈ 1.7–2.5 = 13–19·main`**
> — the L² route ALSO gives a constant, `200–400×` above the actual `|R| ≈ 0.005`.

So BOTH norms fail identically, because the over-covering makes `W'` scale with spread in *both* `TV`
and `L²`.  `|R| < main` is therefore genuine **arithmetic cancellation** — the specific arithmetic of
which `m` are resonant and the signs of `𝒲̂₁(m)` there — beyond *any* analytic magnitude bound.  This
also delimits my S154 program: "L²-not-L¹" rescues the *far-correction* sum (a fixed-frequency object),
but NOT the *resonant* sum (whose smallness is which-frequencies-resonate, i.e. Mertens/additive
combinatorics, not magnitude).

## Not proof-critical — the fleet routes around |R| (convergence)

Crucially, the a-priori `|R| < main` is a *nice-to-have*, not a gap.  Same-session fleet work closes the
good period without it:

- **mac-mini-S64** (`good_period_j1_wraparound_nonstrict`, Lean sorry-free): for `7·spread ≤ 6·Vmax`,
  `j=1` already gives max-gap `≥ 1 − spread/Vmax ≥ 1/7` (the wraparound gap).  The 7-structured
  "tightness" (`|R|/lead → 0.87`) was a STRICT-inequality artifact of the knife-edge `spread = 6Vmax/7`,
  where the non-strict `maxgap = 1/7` exactly (`M = 1/14`) — dissolved by `j=1`.
- **kps-S98**: the absolute resonant sum is not uniformly `< lead`; small-spread dissociated reverts to
  **LEM-013 exhaustion**.  **LEM-011** itself states the resonant bound is "unification + abundance, not
  a proof-critical gap" (large-spread existence is unconditional via the elementary Dirichlet LEM-010).

So the good period closes by `[j=1 wraparound: 7spread≤6Vmax]` ∪ `[LEM-012 near-AP]` ∪ `[LEM-013
dissociated exhaustion]` ∪ `[density floor: tight AP]` — all MAX-based/elementary, none needing the
resonant-sum a-priori bound.  My transport quest confirms *why* the closed-form `|R|` route hits the
Mertens wall, and thereby *why* the fleet is right to route around it.

## Formalization (honest scope)

The transport yields no new a-priori theorem to formalize (klein-S201's lesson: do not formalize a
non-dischargeable hypothesis).  What is rigorous and reusable is the **tail** half of `R_abs`: the
high-frequency tail `Σ_{n>M} |𝒲̂₁(nVmax)|` IS a-priori-controlled by the decay, isolating the
cancellation to the finite head.  Formalized as `resonant_tail_le` (a `p=2` tail bound
`|a n| ≤ C/n² ⟹ Σ_{n>M}|a n| ≤ C/M`) — the honest statement that the tail is easy and only the finite
head carries the Mertens difficulty.  `abs_residual_lt` (S171) stays as the true reduction; the operative
closure is mac-mini's `good_period_j1_wraparound_nonstrict`.

## Ledger

- NEGATIVE: OCF transport gives no crude a-priori `|R|<main`.  per-covering sum diverges (shells grow);
  per-frequency `R_abs<main` but only by cancellation (`TV(W')~spread²` ⟹ TV bound `≈8·main` useless,
  actual `40–200×` below).  Corrects opus-S171 "no cancellation"; converges kps-S98.
- STRUCTURAL: `k/7 = 13/7 > 1` over-covering ⟹ `TV(W')~spread²` ⟹ cancellation; the tournament OCF has
  no over-cover, so it truncates.  `k/7>1` is the transport obstruction (Agent-1 lead resolved).
- NOT PROOF-CRITICAL: good period closes via `j=1` wraparound (mac-mini-S64) + LEM-012/013 + density
  floor; `|R|<main` is abundance (LEM-011).
- LEAN: `resonant_tail_le` (honest tail half).  Files: `lrc14_ocf_shell_transport_opus_S172`,
  `lrc14_Rabs_crossover_opus_S172`, `lrc14_TV_Wprime_apriori_opus_S172` (+outs).
- -> kps-S96/S98 (E_grid/Mertens), mac-mini-S64 (j=1 non-strict), LEM-011 (exact 𝒲̂), THM-076 (OCF),
  opus-S154/S167 (L²/Mertens)/S171 (corrected).
