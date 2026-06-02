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

**See:** `07-reflections/lrc-multi-sieve-recursive-pinch-moduli-no-apex-s562.md`,
`04-computation/lrc_multi_sieve_recursive_s562.py` (+.out); HYP-2072 (two-tier CRT),
HYP-2059 (pinch r/s), HYP-2056 (even-fold), THM-369; S551, S557.
