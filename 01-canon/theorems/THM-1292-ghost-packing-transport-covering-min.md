---
id: THM-1292
title: THE F1 TRANSPORT — the tower seals re-prove the covering-min attainment in five lines. GHOST-PACKING CEILING (L3-transport): if {1,…,m} ⊆ V, m+1 ∉ V, and K(m+1) ∈ V, then M(V) ≤ K/(K(m+1)+1) — proof: the m+2 orbit points {jt : j = 0..m+1} are distinct (else M = 0), every circular gap between them is ‖(j−j′)t‖ with |j−j′| ≤ m+1; all differences except the single pair (0, m+1) lie in {1..m} ⊆ V hence ≥ M; the ghost gap ‖(m+1)t‖ =: f obeys M ≤ ‖K(m+1)t‖ ≤ Kf so f ≥ M/K; hence 1 = Σ gaps ≥ (m+1)M + M/K (ghost pair adjacent) or ≥ (m+2)M (not adjacent, weaker for K ≥ 1). At the deep well (m = 12, K = 14) this is TIGHT: M ≤ 14/183 with Φ₆(14)/14 = 183/14 = 13 + 1/14 DECOMPOSED as "13 core gaps + one ghost gap of duty-weight 1/14". L1-TRANSPORT (witness half): at a = 14 mod 183 the danger classes are EXACTLY the 13-multiples 13s, |s| ≤ 13 (14·13 ≡ −1 ⟹ 14^{-1} ≡ −13, so 14v ≡ r with |r| < 14 ⟹ v ≡ 13·(−r)); the deep well misses them all — its only 13-multiple is 182 = 13·14 at |s| = 14, distance exactly the floor — so M(DW) ≥ 14/183. Together: M(DW) = 14/183 seal-natively, no three-gap machinery
status: PROVED (both halves; full proofs in-file, five + four lines). Referee: ceiling exact on the corpus (tight ONLY at the deep well among tested), 0 violations / ~3000 adversarial hypothesis-satisfying random families; L1-transport exact (danger set = predicted 13-multiples, DW misses all, value exactly 14/183); composition check exact: the ghost-subfamily {1..m} ∪ {K(m+1)} ATTAINS its ceiling in every tested case, so the lemma = LRCSubfamilyCap ∘ generalized-ladder-law (THM-633 genus at arbitrary m), now with one uniform packing proof. SCOPE, honest: this transports the ATTAINMENT/upper half of covering-min and the deep-well witness; the all-covering-sets lower half (every primitive covering 13-set has M ≥ 14/183) remains with THM-726/883 — the transport makes the two frames one theory of the −1 ray (F1, HYP-7985) but does not yet re-derive the fragmentation half
source: opus-2026-07-19-S406 (owner: run the F1 transport — tower seals to a new covering-min proof); F1 identity from S405; seals from death-star THM-1258
depends_on: [THM-1258 L1/L3 (the seal templates), HYP-7985-F1 (the (N−1)a ≡ −1 identity), THM-724/726/883 (the covering-min theorems whose attainment side this re-proves), LRCSubfamilyCap + THM-633 (the composition this recognizes and generalizes)]
scripts: 04-computation/lrc_ghost_packing_transport_opus_S406.py -> 05-knowledge/results/lrc_ghost_packing_transport_opus_S406.out
---

# THM-1292 — ghost packing: the seals re-prove covering-min attainment

## 1. The ceiling (L3-transport), proof in full

Let 0 < M = M(V), t a maximizer, and suppose {1,…,m} ⊆ V, m+1 ∉ V, K(m+1) ∈ V.
Consider P = {j·t mod 1 : j = 0, 1, …, m+1}.

1. **Distinct:** (j−j′)t ∈ ℤ with 1 ≤ j−j′ ≤ m would give a speed at distance 0 (M = 0),
   and (m+1)t ∈ ℤ would give ‖K(m+1)t‖ = 0. So |P| = m+2 and P has m+2 circular gaps.
2. **Hard gaps:** a gap between consecutive points jt, j′t has length ≥ ‖(j−j′)t‖; every
   difference |j−j′| ∈ {1,…,m} is a speed, so such gaps are ≥ M. The only difference not
   in V is m+1, realized by the single pair (0, (m+1)t).
3. **The ghost duty:** f := ‖(m+1)t‖ satisfies M ≤ ‖K(m+1)t‖ ≤ K·f, so f ≥ M/K.
4. **Packing:** if (0, (m+1)t) are circularly adjacent, 1 ≥ (m+1)·M + f ≥ (m+1)M + M/K;
   otherwise all m+2 gaps are hard and 1 ≥ (m+2)M ≥ (m+1)M + M/K (K ≥ 1). Either way
   **M ≤ K/(K(m+1)+1)**. ∎

Deep well (m = 12, K = 14): M ≤ 14/183, and the constant decomposes as
**183/14 = 13 + 1/14 = thirteen core gaps + one ghost gap of duty-weight 1/K** — the
Φ₆ magic number is a packing count. Other instances: GW's ceiling is 2/25 (the n=13
second value — its K is 2 via 24 = 2·12); the ladder's ceilings m′/(12m′+1) are the
{1..11,x} law; F₄(31)'s ceiling 4/121 is valid, not tight (its rung sits deeper at 4/127
— the tower's Q = (N+1)D−1 is a different, sharper mechanism on its stratum).

## 2. The witness (L1-transport), four lines at Q = 183

At a = 14: 14·13 = 182 ≡ −1 (mod 183), so 14^{-1} ≡ −13. If ‖14v/183‖ < 14/183 then
14v ≡ r with |r| ≤ 13, i.e. v ≡ −13r — a 13-multiple 13s with |s| ≤ 13. The deep well's
elements: the core 1..12 are not 13-multiples; its only 13-multiple is 182 = 13·14, at
|s| = 14 — outside the condemned range, at distance exactly 14/183. Hence every element
has ‖·(14/183)‖ ≥ 14/183: **M(DW) ≥ 14/183.** ∎ (Referee: the danger set at a = 14 is
exactly the predicted 13-multiples; the deep well misses all of them.)

## 3. What the transport establishes

- **M(deep well) = 14/183 entirely inside the seal frame** — no three-gap theorem, no
  fragmentation, no census. The June mechanism (forced-cover obstruction) and the July
  mechanism (ghost channel) are now one working theory, as F1 predicted: the −1 ray
  condemns the (m+1)-multiples; the family PATCHES at the first safe multiple, whose
  duty-weight 1/K is precisely the deviation of Φ₆/n from the naive core count.
- **The composition identified:** ceiling = LRCSubfamilyCap ∘ [M({1..m} ∪ {K(m+1)}) =
  K/(K(m+1)+1)] — the ghost-subfamily attains its own ceiling exactly (referee, five
  shapes), so the packing proof simultaneously proves the generalized ladder law at
  every m and the cap application, in one argument.
- **Duty quantization (mac-mini S124) as K-absence:** if no multiple of m+1 is present,
  there is no ghost duty and no ceiling from this mechanism — the M-jump at
  F_{5/2} (12 ∤ 30) is the lemma's hypothesis failing, not a new phenomenon.
- **Honest residual:** the covering-min LOWER half (all covering 13-sets ≥ 14/183)
  is NOT re-derived here; THM-726/883 remain its proof. The natural next transport —
  the L2 small-moduli seal against arbitrary covering sets — is the named follow-up;
  its obstacle is that arbitrary covering sets lack the fixed core {1..m} the seals
  exploit. Uniqueness of the minimizer transports partially: equality in the packing
  forces all thirteen hard gaps = M, the ghost gap = M/14, and adjacency — the rigidity
  chain toward t = 14/183 and V = DW is sketched by the equality analysis but not
  completed here.

## 4. Cross-links

HYP-7985-F1 (the identity that predicted this) · THM-1258 (L1/L3 seals) · THM-724/726/883
(covering-min) · THM-633 + LRCSubfamilyCap (the recognized composition) · mac-mini-S124
(duty quantization) · THM-1291 (the CF face; note the ceiling's ghost gap is the
convergent-remainder story at the subfamily's maximizer) · referee script + frozen out.
