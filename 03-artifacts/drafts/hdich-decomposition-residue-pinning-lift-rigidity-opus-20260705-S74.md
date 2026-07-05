# hdich decomposed: residue pinning (THM-593A) + the lift-rigidity trichotomy

**opus-2026-07-05-S74 (HYP-4097).** The n=12 rigidity leaf — `M(W) = 1/13 ⟹ W = c·{1..12}`,
the mathematical content of klein-S133's `hdich` — decomposes into corpus-existing machinery
plus finite checks. Verified exactly (`lrc14_hdich_lift_rigidity_opus_20260705_S74.py` + .out).

## Half 1 — residue pinning is free at the prime modulus

By primitivity + scale (klein-S131/S132) assume gcd(W) = 1; peel c so the target is: primitive
tight ⟹ W = {1..12} up to lifts. THM-593A (`exists_residue_one_of_tight`, in the corpus
mathlib-track since S34): a tight-from-above family with no multiple of q represents every unit
inverse-residue ±a⁻¹. **At q = 13 PRIME every nonzero residue is a unit**: 12 requirements
(the classes ±a⁻¹, a ∈ (ℤ/13)ˣ cover all 12 nonzero classes), 12 runners — pigeonhole makes the
residue multiset EXACTLY the full system {1,…,12}, each once. (No multiple of 13: a runner
≡ 0 would sit at distance 0 from the origin at every 13-witness, contradicting tightness — the
S73 chain, step (1).) So `hdich` reduces to: **no proper lift of the full-residue base stays
tight.**

## Half 2 — the lift-rigidity trichotomy (mechanism identified, verified 96/96)

A lift replaces r by r + 13k (residues unchanged, so every 13-witness value persists and
M ≥ 1/13 stays; the question is strictness). All 96 single lifts (r = 1..12, k = 1..8) are
STRICTLY loose; min M over all lifts = 1/12 (at r = 12, k = 1) — a uniform **rigidity gap
1/12 − 1/13 = 1/156** for single lifts. Three regimes, each matching existing machinery:

1. **Sieve regime** (most lifts): if r was the unique multiple of some modulus m ∈ {2..12} in
   the base and 13k ≢ 0 (mod m), the lifted family MISSES m — the sieve witness t = 1/m gives
   M ≥ 1/m ≥ 1/12 > 1/13. Verified signature: lifted r = 12, 11, 10 give M exactly 1/12, 1/11,
   1/10 at t = 1/12, 1/11, 1/10. Corpus lemma: `sieve_one_div`.
2. **Sieve-surviving lifts** (CRT classes of k, e.g. r = 7, k ≡ 0 mod 7 → 98 keeps the
   7-multiple): finitely many CRT classes per r; each is a parametric one-killer family over
   the 11-element base {1..12}∖{r}, which has M ≥ 1/12 by the LRC(12) CITATION — margin
   β = 1/12 over the 1/13 target. Large k closes UNIFORMLY by klein's one-window peel
   (`lonely_of_window_margin` at β = 1/12); small k is a finite check. Verified escapes: M =
   7/52, 2/19, 2/23, 14/117 — all ≥ 2/23 > 1/13, argmax denominators ∉ 13ℤ (the "better
   witness elsewhere" the window lemma constructs).
3. **Multi-lifts**: the same trichotomy one runner at a time — every lifted runner is either
   sieve-exposed, window-peelable against the remaining family's citation floor, or in the
   finite low-k census. This is the compressed-peel architecture one level down (the tower
   contracting again: n = 14 → 13 → 12 citations).

## The q = 5 contrast (why the mechanism is honest)

At q = 5 the lift 2 → 7 KEEPS tightness ({1,3,4,7}, THM-593-B3's perm-lift beater) — verified
here too. The difference: at q = 5 the base-minus-r family has M = 1/4 but the window/sieve
room is absorbed by the tiny unit group (4 runners, dense teeth); at q = 13 the 11-runner
remainder's citation floor 1/12 strictly beats 1/13 and the window lemma has room. Rigidity is
a LARGE-q phenomenon, and 13 is comfortably large — the proof above never needs exhaustion,
only the gap 1/12 > 1/13.

## What this hands the assembly (klein HYP-4096 / mac-mini HYP-4095)

`hdich` = Half 1 (one pigeonhole over THM-593A, Lean-ready — the lemma is already kernel-pure
in the corpus) + Half 2 regimes 1–2 (existing `sieve_one_div` + `lonely_of_window_margin` +
per-r finite checks) + regime 3 (induction wrapper). No new analytic input anywhere; the only
computations are the per-r CRT-class finite checks (≤ 12 × small). The rigidity gap 1/156 also
gives `hcorner` sharpening: near-tight 12-families sit at M ≥ 1/12 unless they ARE the AP —
the corner's β can be taken at 1/12, not 2/25, over the lifted classes.

-> THM-593A (opus S34), klein HYP-4093/4095/4096, mac-mini HYP-4089/4094, LRC(≤13) citation,
   MISTAKE-100/102 (this replaces the computational verification with a mechanism).
