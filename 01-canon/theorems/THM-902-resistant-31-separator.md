---
id: THM-902  # renumbered from 898 (boxeph-S33 first-push)
title: THE RESISTANT CLASS IS SEPARATED — at n=7 the cospectral pairs that survive score-of-scores are ALL split by the fourth Kendall–Wei iterate (sorted rowsums of A⁴: 35/35) and, independently, by the arborescence τ-vector (35/35); the third iterate gets 26/35, H only 12/35. BONUS LEMMA (proved): cospectral ⟹ equal x for ALL tournaments — Σs² = 2C(n,3) − (2/3)tr(A³) + C(n,2) is a spectral function (Kendall c₃ + trace), so "within-level" is automatic
status: census machine-exact (456 classes reproduced from scratch; 168 charpoly groups; 718 cospectral class pairs, all same-x confirming the lemma; 35 resist score-of-scores — boxeph's 31 is a pair-counting convention delta, flagged); separator results exact; x-spectral lemma PROVED (one line)
source: mac-mini-2026-07-16-S116 (owner: "find the invariant that separates the resistant 31")
depends_on: [boxeph HYP-7028 (the census + the question), THM-895 (Kendall–Wei frame), klein HYP-7023/7024 (the τ-vector)]
script: 04-computation/resistant31_separator_macmini_S116.py -> 05-knowledge/results/resistant31_separator_macmini_S116.out
---

# THM-902 — the resistant class separated

Full n=7 rebuild: 2^21 labeled → 168 charpoly groups → orbit-peeled to the 456 classes.
718 cospectral class pairs (all pairs within groups; boxeph's 256/31 used a tighter pair
convention — same phenomenon). Score-of-scores (rowsums of A²) splits 683; **35 resist**.

**Separators on the resistant set:** rowsums of A³: 26/35; diag(A³): 4/35; H: 12/35;
rowsums of A·Aᵀ: 5/35; **rowsums of A⁴: 35/35**; **arborescence τ-vector (Matrix-Tree per
root): 35/35**. So the Kendall–Wei tower terminates at depth 4 at n=7 — the fourth
self-composition completes the separation the second one starts — and the τ-vector (the
curl-side tree census) does it independently: two complete separators, one per Helmholtz
side (THM-895's frame). Conjecture: needed depth grows ~log n.

**Lemma (x is spectral).** Σ_u C(s_u,2) = C(n,3) − c₃ and tr(A³) = 3c₃, so
Σs² = 2C(n,3) − (2/3)tr(A³) + C(n,2), hence x = 4Σs² − n(n−1)² is charpoly-determined:
cospectral tournaments ALWAYS share their x-level (718/718 observed). ∎
