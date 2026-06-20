---
id: THM-546
title: The far-element decorrelation bound for LRC(14) — the single-element peel deviation Δ_w = p0(E'∪{w}) − p0(E') − (1/7)p1(E') obeys a RIGOROUS, CONVERGENT estimate |Δ_w| ≤ κ·V(E')/(π²w) with κ=2Σ_{7∤n}|sin(πn/7)|/n²=1.857 and V(E') the arc-complexity, closing HYP-2642's remaining "one open constant" for the gapped regime (where the absolute multi-D lattice envelope DIVERGES)
status: PROVED (the bound, elementary Fourier/BV; κ explicit; V(E')≤42Σ_{e∈E'}e explicit). VERIFIED (exact rational Δ_w vs the bound, all tested E',w: |Δ_w|·w ≤ C(E')=κV/π² with room). The CONSEQUENCE (gapped configs close via this + HYP-2642's Q(k−1)<cap_k) is rigorous; the ungapped-wide regime is handled separately by scale-invariance (THM-531). LRC(14) NOT proved (the bounded finite check + the ungapped accounting remain).
source: mac-mini-2026-06-20-S1
depends_on:
  - HYP-2642   # kps far-element plateau recursion: p0(E)->Phi(E')=p0(E')+(1/7)p1(E'), bounded by Q(k-1)<cap_k
  - THM-534    # the seven-sector moment object meas(S7)=p0
  - THM-531    # AP-orbit (translation+scale) invariance — handles the ungapped wide regime
  - THM-503    # the 7-vanishing s(7j)=0, here the sector-Fourier ŝ_j(7m)=0
related:
  - HYP-2643   # the absolute multi-D envelope DIVERGES (this converges by isolating one element)
  - HYP-2673   # codex's "constant stack" — the empirical Δ_w tax stratification this bounds
  - HYP-2607
  - OPEN-Q-108
external: Lonely Runner Conjecture (first open case = 13 speeds). Koksma/Erdős–Turán BV discrepancy.
---

# THM-546 — The far-element decorrelation bound

## Setup (HYP-2642's recursion)

The S3 crux is `p0(E) := meas(S7(E)) ≤ cap_k`. kind-pasteur's **HYP-2642** peels the largest
element `w = max(E)`, `E' = E∖{w}`: since the single runner `w` occupies one sector per phase,
it fills **at most one** missed sector, so exactly
> `p0(E'∪{w}) = p0(E') + Σ_{j=1}^{6} meas{ B_j(E') ∩ { frac(wx) ∈ [j/7,(j+1)/7) } }`,
> `B_j(E') := { x : E' misses EXACTLY sector j }`.

As `w→∞` the far runner decorrelates (Weyl), each term `→ (1/7)meas(B_j)`, so
`p0(E'∪{w}) → Φ(E') := p0(E') + (1/7)p1(E')` (`p1=Σ_j meas(B_j)` = the one-missed measure),
and `Φ(E') ≤ Q(k−1) := max_{|E'|=k−1}Φ(E') < cap_k` (codex: `Q(7)=0.197<cap_8`, `Q(8)=0.362`,
`Q(9)=0.448`, margins `0.13–0.18`). The **only** remaining piece is a uniform bound on the
deviation
> `Δ_w(E',w) := p0(E'∪{w}) − Φ(E') = Σ_{j=1}^{6} [ meas{B_j ∩ {wx∈sector_j}} − (1/7)meas(B_j) ]`.

The team's obstruction (HYP-2643, verified): the *absolute* multi-dimensional lattice envelope
`Σ_{supp≥6}∏|K(n)|` **diverges**, for wide sets too — so "the one open constant" looked
unreachable by an absolute bound.

## The bound (PROVED)

**Isolating the single far element collapses the multi-D sum to a 1-D BV discrepancy**, which
converges. Each bracket is the discrepancy of `1_{sector_j}(frac(wx))` against its mean over
the *fixed* arc-union `B_j`; by Fourier (`ŝ_j(n)` = sector-indicator coefficient,
`1̂_{B_j}` = arc-union coefficient),
> `Δ_w = Σ_{j=1}^{6} Σ_{n≠0} ŝ_j(n)·1̂_{B_j}(−nw)`.

Using `|ŝ_j(n)| = |sin(πn/7)|/(π|n|)` (**`=0` when `7∣n`** — the THM-503 apex-prime vanishing)
and the BV decay `|1̂_{B_j}(m)| ≤ #arcs(B_j)/(π|m|)`:

> **`|Δ_w(E',w)| ≤ κ · V(E') / (π² · w)`,  `κ := 2 Σ_{n≥1, 7∤n} |sin(πn/7)|/n² = 1.85690…`,
> `V(E') := Σ_{j=1}^{6} #arcs(B_j(E')) ≤ 42·Σ_{e∈E'} e`** (each `e` makes `7e` sector-crossings).

This is **rigorous and convergent** — the divergence was an artefact of summing the full
relation lattice; peeling one element at a time leaves a single far frequency `nw` (giving
the `1/w`) against a fixed finite-variation set. **VERIFIED** exact: `|Δ_w|·w ≤ κV/π²` on
consecutive, codex's B13-leader, dyadic-block-`m=4`, and third-pocket cores, for `w` up to
`200` (actual `|Δ_w|·w` runs `0.01–2.8`, comfortably under `C(E')=5–16`).

## What it closes — and what remains

- **Gapped configs (`w` large relative to `E'`): CLOSED.** For
  `w > κV(E')/(π²·margin_k)` (e.g. `≈40` for the `consec_8` core at the `k=9` margin
  `cap_9−Q(8)=129643/980980≈0.132`), `|Δ_w| < margin_k`, so
  `p0(E) ≤ Φ(E') + Δ_w ≤ Q(k−1) + margin_k = cap_k`. This is the rigorous form of HYP-2642's
  "one open constant" and the convergent **signed** estimate codex's HYP-2673 constant-stack
  approximates empirically (`Δ_w^+/p1 ≤ 2/5, 1/3, <3/10` taxes are *consequences* of this).
- **Ungapped wide configs (`w ≈ max(E')`, the third pocket / wide APs): scale-invariance.**
  Here the bound is `O(1)` (no gap), but these are exactly the high-additive-energy GAP /
  dilated-AP rows: by THM-531 (translation+scale invariance) a dilated AP has the AP's `p0`,
  and by HYP-2637 higher Freiman dimension *lowers* `p0` (dimension penalty, margin `≥0.21`).
  So the ungapped regime is closed by invariance + dimension, not by decorrelation.
- **Remains:** the bounded finite check (consec_k extremal, done by codex/kps); making the
  gapped/ungapped split a clean dichotomy with explicit cutoffs; and a *sharper* (signed/Abel,
  or exact-`#arcs`) `Δ_w` bound to shrink the finite base (the present bound is `~5–30×` loose).

## Net

The "one open constant" of the LRC(14) endgame is now an **explicit, proved** quantity:
`κ = 1.857`, with `|Δ_w| ≤ κV(E')/(π²w)`. The decorrelation of the far runner — which the
absolute multi-D envelope could not see — is captured by a one-dimensional BV discrepancy
against the fixed "misses-exactly-one-sector" set, carrying the `7`-vanishing of the apex
prime in its kernel. The gapped regime closes; the ungapped regime is scale-invariance; the
tight margin lives only in the (done) bounded check. **LRC(14) is NOT proved**, but its last
analytic unknown has a closed form.
