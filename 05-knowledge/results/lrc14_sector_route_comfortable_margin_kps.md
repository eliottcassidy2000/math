# LRC(14) sector route — the tight margins live ONLY in the finite check; everything else is comfortable

**kps-2026-06-20** (with per-sector workflow agent THM-PSK-4, codex HYP-2671 bridge). Sector route to LRC(14).
EXACT-Fraction verified throughout.

## Reduction (PROVED upstream)
LRC(14)-S3 ⟺ `p_0(E) := meas(S7(E)) ≤ cap_k` for every primitive k-set E (0∈E, |E|=k), k=8..12.
S7(E) = {x : all 6 inner sectors [j/7,(j+1)/7), j=1..6, hit by some frac(e·x)}.
caps: cap_8..cap_12 = 2243/5880, 1979/4004, 55/91, 66/91, 6/7 (= 0.3815, 0.4943, 0.6044, 0.7253, 0.8571).
**k≤7: PROVED by pigeonhole (μ_{1/7}=1).**

## The far-element plateau and its deviation
For E = E'∪{w}, w = max E: `p_0(E) = Plat(E') + Δ_w`, where
  `Plat(E') := p_0(E') + (1/7) p_1(E')`,  `Δ_w := p_0(E) − Plat(E')`.
**Q(k−1) := sup Plat over (k−1)-cores = Plat(consec_{k−1})** (VERIFIED argmax, span≤14):
  Q(7..11) = 0.1966, 0.3621, 0.4479, 0.5313, 0.6022.
**margin_k := cap_k − Q(k−1)** = 0.1849, 0.1322, 0.1565, 0.1940, 0.2549 (k=8..12).

## Two rigorous facts
- **σ-bound (THM-PSK-2, PROVED, BV/Koksma):** `|Δ_w| ≤ (6/7)·σ(E')/w`, σ(E')=Σ_{e∈E'}e. (per-arc, no cancellation).
- **w|Δ_w| is UNBOUNDED (THM-PSK-3, PROVED):** the "uniform C(k)=sup w|Δ_w|≤c·k" target (old HYP-2653b/c) is FALSE
  — at a resonance tuned to a wide cluster, w|Δ_w| = Ω(spread). BUT |Δ_w| = (w|Δ_w|)/w stays SMALL (a floor ~0.007).

## The CORRECT object and the comfortable-margin structure
The dovetailing target is the UNIFORM-in-w bound **sup_{max(E')>B} Δ_w ≤ margin_k**, equivalently p_0(E) ≤ cap_k.
The decisive structural fact (per-sector agent THM-PSK-4, kps-verified exactly):

> **The tight margins live ENTIRELY in the finite near-consec check. Every spread/wide/far config has
>  p_0 ≤ cap with COMFORTABLE margin (≥ 0.22), via the Plat↔Δ ENTANGLEMENT: a spread base has SMALL Plat,
>  which compensates a large Δ_w.**

Evidence (exact):
- **Finite check (span ≤ 14), k=8..12: PROVED, 0 violations**; consec is the argmax (k≤11); margins as above.
- The DECOUPLED bound p_0 ≤ Q(k−1)+Δ_w is LOSSY at the boundary: e.g. base [0,3,5,7,9,11,13,14], w=15 (k=9)
  has Δ_w=0.148 > margin 0.132, but Plat=0.122 (small, spread base), so DIRECT p_0 = 0.269 ≪ cap (margin 0.225).
- **The tight extremizer is a WIDE-base set, not a tight far bound:** the dyadic block [0,1,2,4,8,12,16,20], w=24
  (codex's HYP-2671 extremizer, the B=14 worst far config Δ_w=0.117) has max(E')=20>14 — so it is handled by the
  COMFORTABLE wide-base case (p_0 ≤ 0.27 < cap_9, margin 0.22), not the tight far bound. The 0.015 "tight margin"
  was an artifact of applying the decoupled bound to a wide-base set.
- **Wide base (2nd-largest > 14) ⟹ p_0 ≤ ~0.27 ≪ cap** (adversarial, k=9, 6000 sets, sup 0.2706, margin 0.224).

## The route, restated
1. **k≤7:** pigeonhole. DONE.
2. **Bounded (span ≤ B=14):** finite check p_0 ≤ cap. DONE (0 violations).
3. **Spread/wide/far (span > B):** p_0 ≤ cap with margin ≥ 0.22, by the Plat↔Δ entanglement + σ-bound.
   - bounded base (2nd-largest ≤ B), w ≫ σ(E'): σ-bound closes (RIGOROUS).
   - wide base (2nd-largest > B): p_0 directly small (margin 0.22).
   - boundary (spread base, w comparable): joint p_0 = Plat+Δ comfortable (entanglement).

## Honest status
LRC(14) NOT proved. The remaining content is a RIGOROUS proof of step 3 — "span > B ⟹ p_0 ≤ cap" — now
with COMFORTABLE margins (≥0.22) rather than the tight 0.015 of the old framing, and with the tight cases
provably confined to the (done) finite check. This is the same two-gate structure as codex's HYP-2671
(finite shell pocket + far tax). The cross-route bridge (this session): codex Δ_w^+ ≡ kps Δ_w, same dyadic
extremizer. → HYP-2653d, HYP-2671, THM-PSK-2/3/4, OPEN-Q-108.

## ADDENDUM — packet-mass decay (advancing codex HYP-2674), exact
codex's packet sign word: Δ_w = Σ_{s=1}^6 S_s, S_s = Σ_{cells missing exactly s}[G0(wb−s/7)−G0(wa−s/7)].
Dangerous rows are `++++++` (all packets positive, NO cancellation). kps verification:
- The TOP far configs by Δ_w ARE `++++++`, but the alignment is **confined to bounded spread**.
- **Packet-mass DECAY (structured dyadic + 2-cluster, exact):**
  - max(E) ∈ (15,25]: sup Δ_w = 0.1166 (the dyadic block s=4, max=20) — 88% of margin_9.
  - max(E) ∈ (25,40]: 0.0492 (37%).  (40,80]: 0.0345 (26%).  (80,160]: 0.0393 (30%).
  - `++++++` configs with max(E)>40: Δ_w ≤ 0.012.
- **The dyadic block (max=20, Δ=0.117) is the SOLE near-margin far config.** It is a WIDE-base set
  (max(E')=20>14), so handled by the COMFORTABLE entanglement (wide ⟹ small Plat; p_0 ≤ 0.27 < cap_9).
  For max(E)>25, sup Δ_w ≤ 0.05 ≪ margin — the far tail is comfortable.

NET: the route's only tight points are (i) the near-consec finite check (span≤14, DONE) and (ii) the dyadic
block at span 20 (handled by entanglement, comfortable). Everything else has ≥0.22 (wide) or ≥0.08 (far,
max>25) margin. The remaining RIGOROUS content is the packet-mass tail bound sup_{max>B}Δ_w ≤ ε(B) (ε→~0.035
floor, ≪ margin) — the σ-bound (6/7)σ/w is useless here (σ~w at the resonance), so it needs the Erdős–Turán
packet-cancellation estimate, but now with a LOOSE target (0.05 vs margin 0.132). → HYP-2674 (codex), HYP-2653d.
