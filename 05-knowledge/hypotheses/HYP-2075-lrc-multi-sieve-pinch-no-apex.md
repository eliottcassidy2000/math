---
id: HYP-2075
status: METHOD/SYNTHESIS — multi-sieve implemented & measured; three assumptions challenged with data
source: opus-2026-06-02-S562
related:
  - HYP-2072
  - HYP-2059
  - HYP-2056
  - THM-369
---

# HYP-2075: the LRC multi-sieve — pair-sum moduli are the complete, apex-free primitive

Three sub-sieves (a config is caught if any finds a witness t, ‖v_i t‖≥1/n ∀i):
(D) division t=a/m, m∈{2..M}; (P) pinch t=a/(v_i+v_j) (pair-sum moduli, S557);
(F) recursive even-fold (n=14→n=7→…, doubling-preimage).

## Measured (`lrc_multi_sieve_recursive_s562.py`)
- composition over 329 random+loaded n=14 configs: div{2..14} 88.8%, fold 86.0%,
  div∪fold **95.1%**, **PINCH 100% (complete)**.
- **Assumption #1 (apex is a hard obstruction) — FALSE:** of 2366 configs with a
  multiple of 14, caught by m=14 alone 0%, by some m∈{2..13} **91%**. The apex is
  single-modulus; multi-sieve has no apex. (Refines S559/S561: the apex was an
  artifact of the single-corrector mechanism, not loneliness.)
- **Assumption #2 (small-integer moduli) — wrong scale:** division incomplete at
  any finite M (S551; residual min-modulus 23+ unbounded), but pinch (pair-sum
  moduli) complete with bounded COUNT O(k²) and no apex. Natural sieve moduli are
  v_i+v_j (the optimal witness, S557), not {1..M} nor the c=(k+1) lift.
- **Assumption #3 (flat sieve) — recurses:** the 2-adic even-fold catches 86%; its
  only miss is the odd-split residual (S554), covered by pinch; bottoms out at
  all-odd (t=½).

## Concrete proposal for the finite-check
Sieve at **pair-sum residues mod p** (t≡a/((v_i+v_j) mod p)) instead of / alongside
the c=(k+1) base lift: the optimal witness lives there (S557), count O(k²), no
zero-divisor apex. Recurse 2-adically (fold to proven LRC(13) on the even half).

## Honest scope
Exact witness-finders on samples, not an end-to-end check; "pinch complete" =
computes M(S) (LRC true). Value = de-mystification (no apex) + modulus reframe
(pair-sums), both implementable in the real pipeline. Not a proof of LRC(14).

## Independent scale-up + refinement (monad-compute-2026-06-02)
Stress-tested the pinch completeness claim at ~60× the original scale, integer
arithmetic (`t=a/q` safe iff `N·min((v·a)%q, q-(v·a)%q) ≥ q`), with an exhaustive
ground-truth witness search `a/q, q≤QMAX` for comparison.
- **Broad sweep (`lrc_pinch_completeness_stress_s_monad_compute.py`):** 7968 distinct
  n=14 configs (AP, V*, apex multiples-of-14, LCM-loaded, powers, geometric, random
  over speed ranges ≤20…≤1000), QMAX=400. **PINCH caught 7968/7968 (100%)**; every
  config also had a ground-truth witness with q≤400. 0 refutations.
- **Scarce-witness regime (`lrc_pinch_dense_adversarial_s_monad_compute.py`):** 12104
  distinct dense near-AP configs (13 speeds from {1..hi}, hi=14…24, + all single-
  coordinate AP perturbations), QMAX=600 — the regime where witnesses are rarest.
  **PINCH caught 12104/12104 (100%)**, 0 refutations, 0 witness-less configs.
- **Combined: 20072 distinct n=14 configs, 0 pinch misses.** Strongly corroborates
  the completeness (CAUGHT) claim well beyond the original 329.
- **REFINEMENT (corrects the "optimal witness is a pair-sum" gloss, S557):** the
  *minimal-denominator* witness is frequently NOT a pair-sum — it is often a small
  division modulus (observed q=4,6,7,8,10,12) strictly smaller than any pair-sum.
  So pinch finds **a** witness (complete) but **not** the smallest-q witness.
  Every dense config had a witness with q≤26 (max min-q over 12104 = 26). The
  precise statement is: pair-sum moduli form a *complete* (apex-free) witness
  family, not the *minimal-denominator* one — division moduli can be smaller.
- HONEST: still sample-based, not end-to-end; corroborates HYP-2075 completeness and
  sharpens the optimality sub-claim. No bearing on LRC(14) as a proof.

**See:** `07-reflections/lrc-multi-sieve-recursive-pinch-moduli-no-apex-s562.md`,
`04-computation/lrc_multi_sieve_recursive_s562.py` (+.out); HYP-2072 (two-tier CRT),
HYP-2059 (pinch r/s), HYP-2056 (even-fold), THM-369; S551, S557.
monad-compute stress scripts:
`04-computation/lrc_pinch_completeness_stress_s_monad_compute.py` (+.out),
`04-computation/lrc_pinch_dense_adversarial_s_monad_compute.py` (+.out).
