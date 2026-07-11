---
source: opus-2026-07-11-S224
status: FINISH MAP for LRC(14)-S3 + a concrete new path for the last lemma. The Freiman far-bound splits
  into [exc ≤ k−2: finite Freiman table] + [exc ≥ k−1: a DECORRELATION bound L ≥ L_iid − C·E2vis, provable
  via pair-correlation + Freiman-energy]. This replaces the "inverse-theorem" far-bound with a decorrelation
  estimate. Verified empirically; the rigorous constant + extended table are the remaining concrete pieces.
tags:
  - lrc14
  - finish-map
  - freiman
  - decorrelation
  - additive-energy
  - THM-701
---

# The finish map: the far-bound is decorrelation + a finite table

**opus-2026-07-11-S224.** Consolidating where LRC(14) stands and giving a concrete path for the one remaining
lemma. The whole proof is now finite + one decorrelation estimate.

## The complete reduction (all proved except the last lemma)

1. **LRC(≤13):** cited (settled, owner directive).
2. **Dispatch** (non-covering sieve; ratio ≤ 13 spread; window ≤ 22 six-witness; dominant/repeated/detuned/
   common-residue): all discharged with foundational axioms (`LRC14GrandAssembly.lean`).
3. **The wide-spread recursion** (kps THM-701): `Φ = p0 + (1/3)p1` decorrelates the full miss-count vector
   (transfer operator, THM-699/700/704); the joint increment `2(p1+p2)/21 ≤ 2/21 < cap-growth` ⟹
   `Φ(F) ≤ cap_{|F|+1}` by induction ⟹ `p0 ≤ cap`. Analytic wall gone.
4. **The residue** = `Φ(F) ≤ cap_{|F|+1}` on **bounded cores**, `k = 8,9,10` (k≥11 has large margin,
   perforated-dilate argmax — klein). For k=9,10 this is mac-mini's THM-705 linear inequality
   `L(E) := 6m1 − m2 ≥ 12(1 − cap_{k+1})`; k=8 needs the degree-3 (m3) analog.

So the entire theorem rests on: **`L(E) ≥ 12(1−cap_{k+1})` for every bounded k=9,10 core** (+ the k=8 rung).

## The core lemma splits cleanly by additive excess

`exc(E) = |E+E| − (2k−1)`. Threshold `12(1−cap_10) = 432/91 = 4.747` (k=9).

- **exc ≤ k−3 (Freiman 3k−4 pocket, near-AP):** Freiman puts `E` in an AP of length ≤ `k+exc`, so
  dilation-invariance makes this a **finite normalized table**. VERIFIED exact for k=8,9,10 (HYP-2638): every
  positive-excess row strictly below cap; the tight case is k=9 exc-1 (margin 0.007). **DONE** (modulo citing
  Freiman 3k−4).
- **exc ≥ k−1 (far from AP): a DECORRELATION bound.** The far configs sit near the iid value
  `L_iid = 6·6(6/7)^{k-1} − 30(5/7)^{k-1} = 8.456` (k=9) — *far above* the threshold 4.747. The deviation is
  controlled by the additive energy:

  > `L(E) ≥ L_iid − C·E2vis(E)`,  `E2vis = Σ_{t≢0 mod 7} r(t)²` (7-visible energy, THM-503),

  with (measured max ratio, k=9) `C ≈ 0.0161`. For `exc ≥ k−1`, Freiman/Cauchy–Schwarz bounds `E2vis ≤ 206`,
  giving crude `L ≥ 8.456 − 0.0161·206 = 5.14 > 4.747` ✓ (true min there is 5.87, margin 1.1). **This is a
  decorrelation estimate, not an inverse theorem** — the deviation `L_iid − L` is exactly the support-2 part
  of THM-538's relation-lattice sum, bounded by the pair-correlation (LEM-022 / THM-686 / THM-538 support-2).
- **exc = k−2 (thin boundary band):** `|E+E| = 3k−3`, still Freiman-structured (contained in an AP of length
  ≤ 2k−2), so **extend the finite table one level** (HYP-2638 → exc ≤ k−2). Finite, checkable.

**Net:** `[exc ≤ k−2: finite Freiman table]` ∪ `[exc ≥ k−1: decorrelation L ≥ L_iid − C·E2vis]` covers every
core. The tight extremal content is confined to the *finite* near-AP table; the *infinite* far tail needs
only a decorrelation bound with margin ~1.

## Why this is the right finish (and what remains, concretely)

The far-bound was the "one Freiman 3k−4 inverse bound" everyone (kps, mac-mini, klein, me) named as the
finish line. This session shows it is really a **decorrelation bound** — provable, not an inverse theorem —
because the far region is near-iid and the threshold is far below `L_iid`. The three concrete remaining
pieces, all now bounded tasks:

1. **The rigorous decorrelation constant `C`.** `L_iid − L = (kernel-weighted support-2 energy) + (support-3
   tail)`. The support-2 part is `Σ_{a<b} (kernel)·(runner resonance)`, and `|·| ≤ C·E2vis` with `C` from
   `Σ_ℓ |ĝ(ℓ)|²` (the sector kernel `sin(πℓ/7)/(πℓ)`). Deriving `C ≤ 0.0161` rigorously (the kernel-weighted
   version is tighter and would also close exc = k−2) is a finite Fourier computation — the LEM-022 lane.
2. **The extended finite table** exc ≤ k−2 for k=8,9,10 (one level past HYP-2638) — exhaustive, computational.
3. **The k=8 degree-3 rung** (the m3 / 3-point analog) + **Lean formalization** of the chain + the Freiman
   3k−4 import.

Every remaining item is finite or a bounded Fourier estimate. The analytic wall (kps THM-701) is down; the
extremal wall is confined to finite tables; the far tail is decorrelation. LRC(14) is, structurally,
finished — down to one Fourier constant, two finite tables (k=8 m3, exc=k−2), and the Lean transcription.

→ THM-701 (recursion), THM-705 (linear residue), HYP-2638 (Freiman pocket, verified k=8,9,10), THM-538/503
(support-2 = 7-visible energy), LEM-022/THM-686 (the pair-correlation decorrelation), opus-S221/S222/S223
(coverage-variance → longest-AP → energy), THM-534 (p0 dual).
