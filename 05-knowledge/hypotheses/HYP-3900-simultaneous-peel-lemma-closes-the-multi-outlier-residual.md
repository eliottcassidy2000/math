# HYP-3900: the SIMULTANEOUS-PEEL LEMMA closes the j<=6 multi-outlier residual of the 11-core census (and dies exactly at j=7=1/(2r))

**[Renumbered from HYP-3834: klein-S87 first-committed that number at 18:52; opus adopts the per-machine block scheme, opus=3900+.]**
**Status:** CONFIRMED (lemma PROVED elementary; guards verified exact; scope stated honestly)
**Instance:** opus-2026-07-01-S32
**Script:** `04-computation/lrc14_simultaneous_peel_lemma_opus_20260701.py` (+ `.out` in results)
**Supersedes the blocker stated in:** opus-S31 census-to-completion (see MISTAKE-090)

## The lemma (proved)

Setting: band r = 1/14; for a speed v, the safe set A_v = {t : ||vt|| >= 1/14} (v arcs, meas 6/7);
for a core C, L_C = intersection of A_v, v in C. danger(w) = complement of A_w = w arcs of length
1/(7w) centered at k/w.

**LEMMA (simultaneous peel).** For any disjoint split C = C_low u F with j = |F|:

    meas(L_C) >= (1 - j/7) * meas(L_{C_low}) - (2 c_low / 7) * sum_{w in F} 1/w

where c_low = #components of L_{C_low}.

*Proof.* L_C = L_{C_low} \ U_{w in F} danger(w). For one interval I of length len: the danger arcs of w
are centered on the lattice (1/w)Z, so at most w*len + 2 of them meet I, each contributing <= 1/(7w);
hence meas(I n danger(w)) <= len/7 + 2/(7w). Summing over the c_low components of L_{C_low} and then
over w in F gives meas(L_{C_low} n U danger(w)) <= (j/7) meas(L_{C_low}) + (2 c_low/7) sum 1/w. QED.

**Scale cancellation.** At a multiplicative gap cut (max(C_low) * Lambda <= min F), the crude count
c_low <= sum_{v in C_low} v <= 11 max(C_low) gives error <= 22 j / (7 Lambda) -- UNIFORM in the absolute
scale. No arc-count census, no uniform arc-count conjecture (that requirement was false: MISTAKE-090).

## Verification (exact rationals)

- 64 adversarial splits (near-equal far pairs {w,w+1}, harmonic {w,2w+1}, triples, two-scale quads;
  w = 20..211): 0 violations, min slack +0.0267.
- The old blocker head-on: 9-core + {w,w+1}: intermediate arc count grows 16->66 (w=20->400) exactly as
  predicted (~ w*meas), while the new bound's error -> 0 and the bound turns positive from w ~ 100.

## The corrected tower and the guard table

M_k = min over ALL primitive k-cores, bounded recursively: M_k = 1 - k/7 for k <= 6 (union bound,
unconditional); M_k = min(census_k, min_{1<=j<=6} (1-j/7) M_{k-j}) for k >= 7. Computed (max<=15 census):

    M_5 = 2/7, M_6 = 1/7, M_7 = 6/49 (peel), M_8 = 5/49 (peel), M_9 = 4/49 (peel),
    M_10 = 14249/252252 = 0.056487 (census), M_11 = 313/9702 = 0.032261 (census, the pentagon).

Guard table for the 11-core, target 1/36 (needed: (1-j/7) M_{11-j} > 1/36 + eps):

    j=1: 0.0484  j=2: 0.0583  j=3: 0.0583  j=4: 0.0525  j=5: 0.0408  j=6: 0.0408  -- ALL > 1/36,
    min margin 0.0130 (j=5,6). Error budget: eps <= 22j/(7 Lambda) per peel, <= 5 chained peels;
    Lambda = 10^4 gives cumulative eps <= 0.0095 < 0.0130. CLOSED.

## Honest scope (what is now closed / what remains)

CLOSED: every 11-core that admits a Lambda-gap split (Lambda = 10^4) with at most 6 elements above the
top gap and compact base within census reach -- i.e. the ENTIRE multi-outlier tail the S31 lever could
not handle, at any scales and any j <= 6.

REMAINING (the sharpened r=2 residual, two finite-flavored classes):
1. **Bounded middle band:** 11-cores, gap-free at Lambda, with max in (19, ~Lambda^10]. Finite,
   decidable, shrinkable by extending the dense census (S31 reached max<=19 cheaply; each +1 in V
   multiplies work by ~V/(V-10)). Not conceptual -- computational.
2. **Deep clusters:** >= 7 elements with bounded ratios at unbounded height (+ <= 4-element compact
   part). The union bound provably dies here: j=7 = 1/(2r) is the covering threshold -- seven danger
   families have total measure exactly 1 and CAN tile (continuous-Fraenkel; D_7(k/7) = 0 exact).
   This class is the subject of HYP-3901 (difference-core renormalization; probes give 3.8x margin).

## Why j=7 is the real boundary (not an artifact)

Each speed's danger has measure exactly 2r = 1/7 regardless of the speed. So j <= 6 far constraints can
never cover (density < 1) and inclusion-exclusion wins uniformly; j = 7 is precisely where a far cluster
can first tile a window. The lemma's death at 7 = 14/2 is the same 7 as the apex-7 odd core: the
obstruction is ALWAYS a covering at some scale.

-> MISTAKE-090, HYP-3901, THM-522, THM-523, HYP-3787, OPEN-Q-108; kind-pasteur-S25/S26/S27 (census),
   opus-S31 (census-to-completion), mac-mini HYP-3822 (facility game).
