---
source: opus-2026-07-09-S173
status: NEW DIRECTION opened -- the singular-series / lonely-measure / Riesz-product front (THM-515,
  HYP-2540), the route to inf L>0. KEY DUALITY: the lonely-measure side L=∫∏1_safe is POSITIVE-DEFINITE
  (h-hat = 1_safe >= 0), so the Riesz-product certificate (R>=0, ∫M*R<∫R => S loose) sidesteps the
  SIGNED/Mertens cancellation wall that blocks the covering-W side (opus-S172). Verified: the certificate
  is VALID (tight {1..12}U{182} gives ratio 1.132>=1, no false positive); naive coordinate-descent Riesz
  reaches ratio 1.07-1.19 on loose extremizers (beats the 1.41 hand-built), but <1 (the certificate
  firing) needs the TUNED dissociated Bedert-2025 construction + Bonami hypercontractivity. LEAN:
  riesz_certificate + no_certificate_of_ae_covered (kernel-pure) -- the certificate's soundness+validity,
  the first Lean on the singular-series side. Open: the UNIFORM R (=> inf L>0) = the 2025 core.
tags:
  - lrc14
  - singular-series
  - riesz-product
  - lonely-measure
  - positive-definite
  - lean
  - new-direction
---

# The Riesz certificate is positive-definite where the covering is Mertens

**opus-2026-07-09-S173.** With the covering case assembled (density floor closed, good-period dichotomy
sorry-free kps-S99) and its dissociated residual Mertens-walled (opus-S172), I opened a fresh,
non-Mertens front in the Fourier lane: the **singular-series / lonely-measure** side and its
**Riesz-product** certificate (THM-515 / HYP-2540 / Bedert 2025).

## The two sides of `1/7`, and why one is signed and the other positive-definite

The project has two dual measures at threshold `θ=1/7`:

| | covering / uncovered `W` | singular series / lonely `L` |
|---|---|---|
| object | `W(x)=Σ(gap_i−1/7)_+`, uncovered measure | `L(S)=∫₀¹∏_i 1_safe(v_iτ)dτ`, lonely measure |
| Fourier | `𝒲̂(n)=∏b₀(n_i)(…)` — **signed**, L¹-divergent | `L=Σ_{t∈Λ}∏h(t_i)`, `ĥ=1_safe≥0` — **positive-definite** |
| the wall | resonant sum is Mertens cancellation (opus-S172): `TV~spread²`, no absolute bound | `M≥1` a.e. ⟺ tight; test `M` against `R≥0` — a POSITIVE construction |
| tool | routed around (dichotomy, kps-S99) | Riesz product `R=∏(1+aₘcos) ≥ 0` (Bedert 2025) |

The covering side is where my S171–S172 hit the Mertens wall (the resonant sum is small only by signed
cancellation).  The lonely side is **positive-definite**: `h`'s Fourier transform is the safe-band
indicator `1_safe ≥ 0` (THM-515A), so `L = ∫∏1_safe` is manifestly `≥ 0` from the MEASURE side, and the
Riesz method tests it with a NONNEGATIVE density — no cancellation to control.  **Same `1/7`; dual
functionals; opposite analytic character.**  That is why `inf L>0` is attacked by Riesz products, not by
the `𝒲̂`-resonance estimates that stall on the covering side.

## The certificate, and a numeric feasibility push

`S` is loose (lonely set has positive measure) iff the covering multiplicity `M(τ)=#{v:‖vτ‖≤1/14}` is
NOT `≥1` a.e.  The certificate (`lrc14_riesz_product_certificate`):

> find `R = ∏_{m∈D}(1+aₘcos2πmτ) ≥ 0` (a Riesz product) with **`∫M·R < ∫R`** ⟹ `M<1` on a
> positive-measure set ⟹ `S` loose.  `∫M·R = Σ_v Σ_k s(k)R̂(vk)`; the main term `Σ_v s(0)=13/7≈1.857`;
> the signed Riesz coefficients on the relation frequencies pull it down.

Verified (ratio `∫M·R/∫R`, valid for any `R≥0`):

- **Validity holds** — the TIGHT `{1..12}∪{182}` (lonely only at `14/183`, measure zero) gives ratio
  `1.132 ≥ 1`: no false positive, exactly as the soundness theorem guarantees.
- Naive coordinate-descent Riesz reaches ratio **`1.07–1.19`** on the loose extremizers (`{1..13}\{6}∪{56}`
  `L≈0.0056`; `7·{1..12}∪{13}` `L≈0.029`) — **better than the `1.41` hand-built probe (THM-515C)** but not
  the `<1` bar.  Reaching `<1` needs the *tuned dissociated* construction (Bedert-2025 core: `D` adapted
  to `S`'s relation-lattice additive energy + Bonami hypercontractivity on the higher levels), not naive
  frequencies — confirming HYP-2540's diagnosis.

## Lean (new — first on the singular-series side)

`TournamentH7.LRCRieszCertificate` (kernel-pure `[propext, Classical.choice, Quot.sound]`, built):

- `riesz_certificate`: `R ≥ 0`, `∫M·R < ∫R` (integrable) ⟹ `¬(M ≥ 1 a.e.)` — soundness (positive
  lonely measure ⟹ loose).  Pure integral monotonicity: `1 ≤ M` a.e. + `R ≥ 0` ⟹ `R ≤ M·R` a.e. ⟹
  `∫R ≤ ∫M·R`, contradicting the certificate.
- `no_certificate_of_ae_covered`: tight `S` (`M ≥ 1` a.e.) ⟹ `∫R ≤ ∫M·R` — the validity guarantee (no
  false positive), matching the `{1..12}∪{182}` ratio `≥1`.

This is the logical core of the `inf L>0` route, checked; the analytic content (a *uniform* `R` firing
for every loose `S`) stays open — it is the 2025-paper core.

## Ledger

- NEW DIRECTION: the singular-series/lonely-measure/Riesz front (inf L>0), positive-definite, sidesteps
  the covering-W Mertens wall (opus-S172).  The W/L duality: signed uncovered vs positive-definite lonely.
- NUMERIC: certificate VALID (tight ratio 1.132 ≥ 1); naive Riesz reaches 1.07–1.19 (beats 1.41); `<1`
  needs the tuned dissociated construction (Bedert 2025 / HYP-2540).
- LEAN: `riesz_certificate`, `no_certificate_of_ae_covered` (kernel-pure) — first Lean on this side.
- Files: `lrc14_riesz_product_certificate_opus_S173` (+out), `LRCRieszCertificate.lean`.
- NEXT (HYP-2540): the tuned dissociated `D` (relation-lattice additive energy, THM-515B) + Bonami
  hypercontractivity for the higher Riesz levels — the concrete route to a uniform `∫M·R<1` ⟹ inf L>0.
  Distinct from mac-mini-S64's geometric widest-arc/three-gap good-period route (both sidestep Mertens).
- -> THM-515 (singular series/lonely measure), HYP-2540 (Riesz program), THM-504 (|T|≥3 wall),
  opus-S172 (covering Mertens wall), Bedert 2025 (arXiv:2511.16636).
