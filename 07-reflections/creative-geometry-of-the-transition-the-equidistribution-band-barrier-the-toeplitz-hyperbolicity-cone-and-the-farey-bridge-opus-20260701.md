# Creative geometry/topology of the covering-min transition (a(n)=n for n≥12): the EQUIDISTRIBUTION BAND-BARRIER (a character-sum mechanism — the construction is the interval with band-gap 3 for all n, disorder is penalized like (D/|S|)·ln|S|~1.5 ln n by Weil equidistribution, so beaters live in a thin low-discrepancy shell that collapses onto the pure interval), the covering as an integer program over the TOEPLITZ-PSD HYPERBOLICITY CONE (Gårding sense, HYP-2974) whose self-concordant barrier is my HYP-3769 rung ladder, and the FAREY-TESSELLATION BRIDGE that reconciles mac-mini's two "hyperbolics" (HYP-3780) — the CF/Stern-Brocot rung descent IS a geodesic in the Farey tessellation of the hyperbolic plane, which is BOTH the self-concordant-barrier hyperbolicity cone AND the apex-7 (2,3,7) hyperbolic geometry, glued by the Dedekind cocycle; plus two topological directions (Morse theory on S¹, Borsuk–Ulam ι-odd index)

*opus-2026-07-01. Owner: think character sums on the band + the barriers of hyperbolicity cones, and come up
with creative geometric/topological arguments. A brainstorm with one grounded mechanism (the band-barrier), one
precise cone (Toeplitz-PSD), one unifying bridge (Farey), and two topological leads — honestly marked.*

## 1. The equidistribution band-barrier (character sums on the band) — GROUNDED mechanism
klein-S60 (HYP-3778) verified the transition combinatorially (ILP); the open part is "no beater with speeds >4n."
A geometric/character-sum mechanism for it:
> A rung-k beater with `M=k/D` (D=k(n-1)+1) must, at every small modulus, keep its worst residue-gap small
> (a large gap = a deep hole = M too big). Tested: the **construction is the interval `{1,…,n-2}+killer`, with
> band-gap 3 for EVERY n** (it covers the radius-1 band `(n,2n-2]` uniformly — this is klein's HYP-3737); the
> beaters n=7..11 sit at **band-gap 4** (a thin "almost-interval" shell); generic spread sets scatter to
> band-gap 4–7. A **spread (disordered) set equidistributes mod a small modulus D'** (Weil / Erdős–Turán
> square-root cancellation of the exponential sums `Σ_v e(mv/D')`), giving max gap `~ (D'/|S|)·ln|S| ~ 1.5 ln n`,
> which **grows past the shell** as n increases. So the covering-min lives in a low-discrepancy shell around the
> interval; Weil equidistribution squeezes it; **at the transition the shell collapses onto the pure interval =
> the construction.**
Status: the *mechanism* is grounded (construction band-gap 3 ∀n; disorder penalized ~ln n), but the band-gap
metric alone is coarse (random sets also hit gap 4), so it is **not yet a proof** — the sharp version needs the
**full all-moduli discrepancy** (Erdős–Turán over q=2..2n simultaneously), which is the honest character-sum path
to klein's >4n residual. This is the elementary Ramanujan–Petersson (√V) tool (THM-546/HYP-2852), now aimed at
the band.

## 2. The covering as an integer program over the TOEPLITZ-PSD hyperbolicity cone (Gårding sense)
Precise convex-geometric home (complements mac-mini HYP-3780's abstract self-concordant framing):
> Covering `m_w(t)≥1` ⟺ (via Beurling–Selberg minorants) a **nonnegative trigonometric polynomial** condition ⟺
> the **Toeplitz moment matrix is PSD** (HYP-2974). The cone of nonneg trig polys of degree ≤N is a
> **spectrahedron** = a slice of the PSD cone `S^{N+1}_+`, and `det` is a **Gårding-hyperbolic polynomial**, so
> this is a genuine **hyperbolicity cone**. Its canonical self-concordant barrier `−log det(Toeplitz)`
> (parameter N+1) is the barrier theory; my HYP-3769 `1/M=(n-1)+1/rung` is the covering-side barrier residual.
So the **covering-min = an INTEGER point in this hyperbolicity cone** (the speeds are integers); the SDP/LP
relaxation is convex; **klein's ILP is exactly the integer program**, and the integrality gap = "which sets are
realizable" = the beaters. The construction is the cone's analytic-center-adjacent vertex (the interval, deepest
covering slack, band-gap 3); the transition = a **facial collapse** (for n≥12 no integer point near the floor
except the vertex). The modular value −1/12 is the relaxation's vertex value; the beaters are integer points that
undercut it only while the shell (§1) is open.

## 3. The FAREY-TESSELLATION BRIDGE — reconciling the two "hyperbolics" (refines HYP-3780)
mac-mini-S71 (HYP-3780) links "self-concordant barrier ⟺ hyperbolicity cone ⟺ apex-7 hyperbolic (2,3,7)." The
last ⟺ mixes two different notions of "hyperbolic": the **hyperbolicity CONE** (Gårding, convex optimization) and
**hyperbolic GEOMETRY** (negative curvature, the (2,3,7) Klein-quartic tessellation of ℍ). They are *not* the same
— but they are **bridged**, and the bridge is exactly the covering-min's own backbone:
> The **rung ladder / continued-fraction descent** `1/M=[0;n-1,rung]` (my HYP-3769; mac-mini's convex CF ladder)
> IS a **geodesic in the FAREY tessellation of the hyperbolic plane ℍ** — the ideal-triangle tessellation on
> which `SL₂(ℤ)` (the modular group) acts, and along which the **Dedekind sum is the transformation cocycle**
> (HYP-3773: the margin = the η-multiplier phase). The Farey/Stern–Brocot tree is simultaneously (a) the CF
> descent = the self-concordant-barrier reciprocity step (Güler; my HYP-3770), and (b) a tessellation of
> hyperbolic space = the apex-7 (2,3,7)/`X₀(14)` hyperbolic geometry.
So the two hyperbolics meet at the **Farey tessellation**: the hyperbolicity-CONE side is the convex/barrier
reading of the CF ladder; the hyperbolic-GEOMETRY side is its ℍ-tessellation reading; the **Dedekind/η cocycle
glues them** (the −1/12 anomaly is the modular cocycle on the tessellation). This is rigorous (the Farey
tessellation genuinely is both), and it upgrades mac-mini's ⟺ from a suggestive conflation to a real bridge.
(Aside: the general Gårding↔hyperbolic-geometry link is the **Lorentzian light cone** — the hyperboloid model of
ℍ sits inside the light cone, itself the hyperbolicity cone of the Lorentzian form; the Farey tessellation is the
`SL₂` arithmetic instance of exactly this.)

## 4. Two topological directions (leads, not developed)
- **Morse theory on the circle.** `g(t)=min_v ‖vt‖` is a piecewise-linear Morse function on `S¹`; `M=max g`; the
  holes are its local maxima. Three-distance ⟹ `g` has ≤3 slopes, so few critical points. The transition = a
  change in the critical-point structure (a hole merging/splitting as the set deforms). A discrete Morse/Euler
  count of holes vs n could give the beater count.
- **Borsuk–Ulam / equivariant index (OPEN-Q-108).** Covering is `ℤ/2`-equivariant under the antipode `ι:a↦−a`;
  the obstruction is the **ι-odd equivariant index** (klein-S56's genus cusp form f₁₄). This is the genuinely
  topological proof obstruction — the "residual" of HYP-3779 (the proof, not a(n)). The transition proof may be a
  vanishing-degree statement (odd index ⇒ a fixed hole survives), the topological form of the band-barrier §1.

## Synthesis
The transition has a convex-geometric skeleton (§2, hyperbolicity cone + integer program), a character-sum engine
(§1, equidistribution/Weil squeezing the low-discrepancy shell), a hyperbolic-arithmetic backbone (§3, the Farey
tessellation gluing the two hyperbolics via the Dedekind cocycle), and a topological obstruction (§4, the ι-odd
index). The **most actionable** is §1→§4: the band-barrier is the elementary (Weil/Erdős–Turán) route, and its
topological packaging is the Borsuk–Ulam ι-odd index — together the honest path to klein's open >4n residual.

## Status
- **Grounded (opus):** the band-barrier mechanism (construction=interval band-gap 3 ∀n; beaters band-gap 4;
  disorder ~1.5 ln n; the low-discrepancy shell) — a real mechanism, not yet a proof (band-gap coarse; needs
  all-moduli Erdős–Turán).
- **Precise (opus):** the covering = integer program over the Toeplitz-PSD hyperbolicity cone (Gårding);
  self-concordant barrier = HYP-3769; klein's ILP = the integer part; transition = facial collapse.
- **Bridge (opus, rigorous):** the Farey tessellation of ℍ is BOTH the CF/self-concordant ladder AND the apex-7
  hyperbolic geometry, glued by the Dedekind/η cocycle — refines mac-mini HYP-3780's two-hyperbolics ⟺.
- **Leads (opus):** Morse theory on S¹ (hole critical points); Borsuk–Ulam ι-odd index (the topological obstruction
  = OPEN-Q-108).
- **Open:** the all-moduli Erdős–Turán bound that closes klein's >4n residual; whether the ι-odd degree literally
  computes the beater count.

Related: HYP-3780/mac-mini-S71 (self-concordant⟺hyperbolicity cone; two-hyperbolics — refined here), HYP-3778/klein-S60
(transition ILP), HYP-3737 (radius-1 band), HYP-2974 (Toeplitz-PSD cone), HYP-3769/3770 (rung ladder=CF=self-concordant),
HYP-3773 (Dedekind=η cocycle), THM-546/HYP-2852 (√V=Weil), OPEN-Q-108 (ι-odd index). HYP-3781 (this).
Script: 04-computation/lrc_equidistribution_band_barrier_opus_20260701.py.
