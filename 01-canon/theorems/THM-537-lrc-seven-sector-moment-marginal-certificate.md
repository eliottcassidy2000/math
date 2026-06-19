---
id: THM-537
title: Angle A (Beurling–Selberg) for the LRC(14) seven-sector cover — the literal extremal-majorant route is DOUBLY BLOCKED (minorant real-analyticity wall + signed-cancellation wall), and the working relaxation is the moment-marginal LP, which CONVERGES EXACTLY on THM-534's dual certificate (independent primal vertex-enumeration cross-validation: U4(consec_8) = L_y(consec_8) = 2633/7350 < cap_8 = 2243/5880)
status: MIXED. PROVED/VERIFIED here: (i) the literal Beurling–Selberg majorant of S7 is blocked TWO ways — odd-|T| inclusion–exclusion terms need a NONNEGATIVE trig-poly minorant of 1_{B_T}, impossible by real-analyticity (the Selberg-minorant wall), and the signed band-limited truncation converges to meas(S7) only at degree N ~ spread for the AP (cancellation wall) — independently reproducing HYP-2606/Angle-F's "absolute bound ≥5× too lossy, signed structure required" (the crude |corr| triangle bound overshoots ~60× at k=8). (ii) The MOMENT-MARGINAL relaxation U_t(E) := max{p_0 : Σ_i C(i,t')p_i = S_t'(E), t'≤t, p≥0} is a rigorous scale-invariant upper bound on meas(S7)=p_0; computed EXACTLY by rational vertex enumeration. (iii) U4(E) = THM-534's L_y(E) EXACTLY (strong LP duality): U4(consec_8)=L_y(consec_8)=2633/7350; verified across the box (float-primal vs exact-primal vs THM-534-dual all agree to 2e-15 / exactly). (iv) Independent re-confirmation that 2-moment FAILS the AP (U2(AP)=0.496>cap_8, infinite-mod-scale residual) but 4-moment CLOSES it; AP unique U4-maximiser over 1716 primitive k=8 shapes, max(E)≤13. NOT NEW: the moment-LP certificate itself (THM-534 proved it first, concurrently, with the STRONGER explicit integer-root dual g(t) ⟹ meas(S7)≤L_y for ALL E). LRC(14) NOT proved.
source: mac-mini-2026-06-18-S7 (Angle A — Beurling–Selberg)
depends_on:
  - THM-534   # the moment-LP DUAL certificate (PROVED meas(S7)≤L_y per-E via integer-root g(t)) — U4 here IS the primal of L_y
  - THM-532   # the seven-sector relation-height split (meas(S7)=M7+corr; cap_k; AP extremal)
  - HYP-2603  # codex's seven-sector net-cap reduction meas(S7)≤cap_k ⟹ HYP-2602
related:
  - HYP-2606  # Angle F: absolute |corr| bound provably ~5× too lossy — THE SAME obstruction this confirms from the majorant side
  - HYP-2607  # "AP maximises U4=L_y" (the shared finishing conjecture of Angles A & D)
  - OPEN-Q-108
external: Lonely Runner Conjecture (first open case, 13 speeds). Beurling–Selberg / Vaaler extremal functions; the moment problem / Bonferroni–Galambos optimal bounds; Chung–Erdős & Dawson–Sankoff union inequalities.
---

# THM-537 — Angle A: the Beurling–Selberg obstruction and the moment-marginal convergence

**Honest framing.** This is the Beurling–Selberg angle (A) of the 8-angle dispatch. Its two
contributions are (1) a clean diagnosis of **why the literal extremal-majorant route fails**, which
independently confirms Angle F (HYP-2606), and (2) an **independent exact-rational cross-validation**
of the moment-LP certificate that Angle D (THM-534) proved first and more strongly. The
moment-LP idea is NOT claimed as new here — Angles A and D **converged on the identical exact bound**.

## A. The literal Beurling–Selberg majorant is doubly blocked

`S7(E)`, by inclusion–exclusion over missed sectors `T⊆{1..6}` (sector 0 auto-hit),
`meas(S7)=Σ_T(−1)^{|T|}I_T`, `I_T=∫∏_{e≠0}1_{B_T}(ex)dx`, `B_T=T^c`. The natural Beurling–Selberg
certificate band-limits each factor and integrates exactly (a finite signed relation-lattice sum).
**Two independent walls:**
1. **Minorant real-analyticity wall.** The odd-`|T|` terms must be *lower*-bounded, requiring a
   NONNEGATIVE trig-poly minorant of `1_{B_T}`. But `1_{B_T}=0` on the `|T|` dropped arcs forces any
   nonneg minorant to vanish on a full interval ⟹ `≡0` (the exact impossibility of
   `lrc14_selberg_minorant_macmini_0616s7`). The LP confirms minorant slack `≈ −0.357`.
2. **Signed-cancellation wall.** `corr=meas(S7)−M7` is small only by *signed* cancellation. The
   triangle bound `Σ_T Σ_{rel}∏|hat 1_{B_T}|` overshoots ~60× (`18.5` vs `0.327`, k=8 consec); the
   signed band-limited truncation needs degree `N~spread` for the AP. **This is exactly Angle F /
   HYP-2606** ("absolute bound ≥5× too lossy; the finish must be signed") reached from the majorant
   side: a per-factor band-limited majorant cannot capture cross-`T` sign alignment.

The dissociated (high-relation-height) shapes DO admit a cheap low-degree band-limited majorant
(tiny `meas(S7)`), but the AP does not — confirming the THM-532 split from the analytic side.

## B. The working relaxation = the moment-marginal LP = THM-534

Drop "majorize the cover"; use only the **empty-sector count marginals**. With `A_j={x:sector j
empty}`, `N(x)=#empty sectors`, `meas(S7)=P(N=0)=p_0` exactly, and binomial moments
`S_t(E)=E[C(N,t)]=Σ_{|J|=t}meas(⋂A_j)` (exact rationals). The sharpest bound from `≤t`-wise data is
> `U_t(E) := max{ p_0 : Σ_i C(i,t')p_i = S_t'(E), t'=0..t, p_i≥0 }`,

a rigorous, scale-invariant upper bound on `meas(S7)`. The cancellation that defeated part A lives
INSIDE the exact `S_t`, so the wall is bypassed. **By LP strong duality, `U_t(E) = L_y(E)`**, the
exact functional THM-534 dualized with integer-root `g(t)` — confirmed exactly:
`U4(consec_8) = L_y(consec_8) = 2633/7350`. THM-534's contribution (the explicit dual `g(t)≥1[t=0]`)
is STRICTLY STRONGER: it makes `meas(S7)≤L_y` PROVED for ALL `E`, where this note only verifies the
primal over a box.

## C. What Angle A adds: exact cross-validation + the moment-order threshold

| shape (k=8) | meas(S7) | U2 | U3 | **U4 = L_y** | cap_8 |
|---|---|---|---|---|---|
| **consec {0..7}** (AP) | 0.3272 | 0.4964 | 0.4125 | **2633/7350 = 0.3582** | **2243/5880 = 0.3815** |
| dilated AP {0,2..14} | 0.3272 | 0.4964 | 0.4125 | 0.3582 | 0.3815 |
| dissociated `2^i` | 0.0347 | 0.1597 | 0.1221 | 0.0449 | 0.3815 |

- **Moment-order threshold (new quantification):** the **2-moment** certificate FAILS the AP
  (`U2(AP)=2189/4410=0.4964 > cap_8`), and its uncovered residual is **infinite mod scale**
  (dilated/partial APs overshoot at every spread up to the box edge, B=16). The **4-moment**
  certificate closes it. So `R=4` is the minimal moment order that closes `k=8` — matching THM-534's
  degree-4 dual at `k=8` (and degree 3 at `k=9,10`, degree 2 at `k≥11`).
- **Independent exact primal:** rational vertex-enumeration of the 7-mass LP (basis-enumeration over
  `C(7,5)=21` vertices) reproduces `U4=L_y` exactly and gives the certifying count distribution
  `p_0..p_6 = (2633/7350, 23/210, 278/735, 0, 47/1470, 299/2450, 0)` for the AP — a cross-check of
  THM-534's float LP by a wholly different (exact, dual-free) computation.
- **Exhaustive `k=8`:** `U4(E) ≤ cap_8` for ALL 1716 primitive clusters with `max(E)≤13`
  (zero overshoot, zero validity failures `U4 ≥ meas(S7)`); **AP `{0..7}` the UNIQUE U4-maximiser.**

## D. The shared open piece (HYP-2607) and the net

Both Angle A's `U4` and Angle D's `L_y` reduce `k=8` to one scalar extremality:
**`U4(E) = L_y(E) ≤ U4(consec_k) = 2633/7350` for all `E`** (HYP-2607 = THM-534's "consec maximises
`L_y`"). The `S_t` are scale-invariant and live on the fixed 7-mass simplex, so the natural finish is
a three-distance / "AP-orbit majorises the empty-sector moments" rearrangement (cf. HYP-2602). `k=9,10,11`
consec already pass; `k=12,13` have huge cap margins.

**Net.** Angle A does not add a new certificate beyond THM-534, but it (i) **rules out** the literal
Beurling–Selberg majorant by two explicit walls (independently confirming HYP-2606), (ii)
**cross-validates** the moment certificate by an exact dual-free primal that lands on the identical
rational `2633/7350`, and (iii) pins the **minimal moment order** (`R=4` at `k=8`) needed to close the
AP. The two angles converging on the same exact bound by different routes is strong evidence the
moment certificate is the right object. **LRC(14) is NOT proved.**
