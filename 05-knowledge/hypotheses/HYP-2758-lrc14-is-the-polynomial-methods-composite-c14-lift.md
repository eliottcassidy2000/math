---
id: HYP-2758
title: LRC(14) is OPEN because 14=2*7 is composite -- the repo's sector cover IS the analytic resolution of the polynomial method's intractable c=14 lift (arXiv:2604.23906)
status: OPEN; strategic reframing from the state-of-the-art paper -- a major connection
source: kind-pasteur-2026-06-21-S24
related:
  - OPEN-Q-108
  - HYP-2757   # the L7 curve atlas
  - HYP-2700   # mac-mini Z/7-coloring
  - THM-530    # the 1/7 pivot
---

# HYP-2758 - LRC(14) = the polynomial method's composite c=14 lift

## The external state of the art (arXiv:2604.23906, "Eleven, twelve, and thirteen lonely runners", Apr 2026)
Proves LRC for up to 13 runners (k=10,11,12). **k=13 (14 runners = repo LRC(14)) is EXPLICITLY LEFT OPEN.**
Their key tool (Prop 4.1, the POLYNOMIAL METHOD): for k+1 an ODD PRIME and p>k^2+k prime, any tuple with
`u_i ≡ i (mod p)` has the LR property. The proof builds `P(X)=prod_i(v_i i^{-1}+X)` as the indicator of a set
G in `Z_{k+1}^*`, then uses **Fermat's Little Theorem** (needs a PRIME field Z_{k+1}) to get a second indicator
`Q(X)=sum_{m in G}(1-(X-m)^k)`, and compares leading coefficients to force |G|=k.

## THE OBSTRUCTION (= why LRC(14) is the first open case)
**`k+1 = 14 = 2*7` is COMPOSITE.** Fermat's Little Theorem fails in the non-field ring `Z_14`, so Prop 4.1
**cannot apply**. Their paper states the bottleneck explicitly: without the polynomial shortcut they cannot
avoid the **`c=(k+1)=14` lift**, which "requires checking ~(13/2)^12 ~ 5e9 times more tuples" -- INTRACTABLE.

## THE CONNECTION (the repo IS the analytic resolution of the intractable c=14 lift)
- Their TWO-SCALE split `t = s/(k+1) + r/p` (Lemma 4.2): `s/14` is the p-independent SLOW scale, `r/p` the FAST
  scale. **This IS the repo's slow/fast reduction** (THM-525/527, the "slow time x").
- Their discretization `r_k(t) = (floor(14*frac(it)))_i` = the **14-SECTOR coloring**. The repo uses **7 sectors**
  -- exactly `14/2`, the COMPLEMENT-symmetry halving (a=-1 reverses sectors; THM-280/549). So the repo's
  Z/7-coloring (mac-mini HYP-2700) + complement = the CRT factorization **`Z_14 = Z_2 x Z_7`** of their Z_14.
- Their polynomial method handles the CANONICAL residue `(1,2,..,13)` = the repo's **CONSEC** case. Their
  intractable c=14 lift = the NON-canonical residues = the repo's **WIDE / multi-cluster** cases (OPEN-Q-108/L7).
> **So `meas(S7(E)) <= cap_k` (the repo's sector bound) is the ANALYTIC form of the paper's intractable
> computational `c=14` lift, performed through the apex prime 7 instead of by brute enumeration.**

## ACTIONABLE leads
1. **CRT the polynomial method onto Z_7.** Z_14=Z_2 x Z_7. The polynomial/Fermat argument WORKS in the prime
   field Z_7. Run their Prop-4.1 construction in Z_7 (and Z_2) and CRT-combine: does the prime-factor polynomial
   method certify the consec / canonical-residue case at the apex prime, bypassing the composite-14 failure?
2. **Their c=2 lift IS tractable; the gap is c=2 -> c=14.** The repo's complement halving already does the
   Z_2 (c=2) part. The remaining c=2 -> c=14 = the Z_7 lift = the 7-sector cover = OPEN-Q-108. Map their
   lifting sieve onto the 7-sector cover to see if the wide bound = a finite Z_7 lift after the polynomial
   reduction of the canonical residue.
3. **The canonical residue (consec) is the EASY half for them (polynomial method) but the repo treats it as the
   finite extremality (mac-mini W_a layers).** Reconcile: is consec-max provable via the Z_7 polynomial method?
4. The non-canonical residues that survive their sieve "are equivalents of (1,2,..,k)" -- so after sieving, ALL
   reduce to consec. The repo's domination (single-far = consec-base maximizer, kps-S23) is the analytic mirror.

## Status
Strategic reframing, not a proof. Confirms the repo's apex-prime-7 route is THE right tool for the composite-14
case the polynomial method cannot reach. NEXT: execute lead 1 (CRT the polynomial method to Z_7) -- if it
certifies the canonical residue at the apex prime, it could replace the intractable c=14 lift for consec, leaving
only OPEN-Q-108/L7 (the wide residue) for the sector cover. -> OPEN-Q-108, HYP-2757, HYP-2700, THM-530.
External: arXiv:2604.23906 (Rosenfeld-Trakulthongchai et al.), arXiv:2511.22427 (nine/ten), arXiv:2509.14111 (eight).
