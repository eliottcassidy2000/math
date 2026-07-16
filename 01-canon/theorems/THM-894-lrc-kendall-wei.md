---
id: THM-894
title: LRC'S KENDALL–WEI — the resonance classifier as a self-referential vector: for a speed set S, the RESONANCE MATRIX M[i][j] = ρ(x_i, x_j) (the exact sawtooth pair law; diagonal 2/13) over the cluster's OWN speeds has Perron pair (λ, w) with w = each speed's resonance centrality defined through all others simultaneously (the fixed point of the recursion into the second moment; residual 1e−16); (E) THE SPECTRAL EXCESS λ/⟨row⟩ is a new one-scalar classifier of feedback-structuredness: AP {1..13} MAXIMAL (1.0293) > mediant (1.0261) > floor (1.0216) > GAP (1.0181) > the U₀-rows (1.010) > 2AP−1 (1.0034) > near-AP (1.0006, essentially rank-1): the tight family amplifies its own resonance most, and the excess sees feedback structure the S321 coherence scalar (a first-moment sum) cannot; (S) the twins' AMBIENT rows separate already at level 2 (λ = 0.4668 vs 0.4565) — consistent with THM-893: only the folded shadows coincide; (N) honest negatives: the Perron profile is NOT uniform on the tight AP (0.052–0.102; the naive tightness=uniformity guess refuted) and the raw Perron mass largely shadows 1/speed (small speeds resonate broadly) — the spectral content is in the EXCESS and the profile DEVIATIONS
status: FORMALIZED + COMPUTED (8-family battery, exact ρ entries, machine-precision Perron; the excess ordering and the twins' separation referee-grade; negatives logged per BH)
source: opus-2026-07-16-S326 (owner: formalize LRC's Kendall–Wei, the resonance classifier as a self-referential vector; + arXiv:2607.14068)
depends_on:
  - THM-863 (the pair law ρ), THM-893 (the frames principle this instantiates)
related: [THM-891-codex (resonant peel — the other resonance-cancellation face), the Kikuchi weave reflection (this session), Kendall–Wei 1955 (the classical iterated-scores eigenvector)]
verification: 05-knowledge/results/lrc_kendall_wei_opus_S326.out
---

# THM-894 — LRC's Kendall–Wei

**The construction.** S = {x₁, …, x_k}; M[i][j] = ρ(x_i, x_j) (exact pairwise
overlap density, THM-863's corrected sawtooth law; M[i][i] = 2/13). M is
symmetric nonnegative; its Perron pair (λ, w) gives each speed a resonance
centrality DEFINED THROUGH the others' centralities — the owner's "the
resonance classifier is a vector over the cluster's own speeds," the LRC
avatar of Kendall–Wei iterated scoring, the fixed point of the recursion
into the second moment. Frames reading (THM-893): row sums are the
frame-dependent first moment; the Perron pair is the all-frames-at-once
invariant.

**The spectral excess.** λ/⟨row sum⟩ measures feedback beyond first moment:

| family | λ | λ/⟨row⟩ |
|---|---|---|
| AP {1..13} (tight) | 0.5460 | **1.0293** |
| mediant S₁₃ | 0.5227 | 1.0261 |
| floor ∖{6} | 0.4893 | 1.0216 |
| GAP | 0.5599 | 1.0181 |
| U₀-row (13,9) | 0.4668 | 1.0103 |
| U₀-row (17,13) | 0.4565 | 1.0094 |
| 2AP−1 | 0.4573 | 1.0034 |
| near-AP {21..32} | 0.4216 | 1.0006 |

The tight family maximizes its own resonance amplification; the E-twins
(AP vs 2AP−1, S321) are separated 1.0293 vs 1.0034; the GAP drops below the
coherent families (the dimension-vs-coherence tangle partially resolved: the
excess is a genuinely spectral, not first-moment, quantity).

**The twins.** The ambient 12-speed rows separate at level 2 (λ differs in
the third digit) — as THM-893 requires: the coincidence is a property of the
folded shadows only.

**Negatives (BH: 2 guesses, both refuted).** Tight ≠ uniform Perron profile
(the AP's resonance matrix is ratio-structured, not circulant); the raw
Perron mass largely shadows 1/speed. The information is in the excess and
the profile's deviations from the 1/speed shadow.

**Named next.** The excess under dilation/translation (it should inherit
saw's invariance signature — check); level-3 (the triple-overlap tensor's
top singular value: the Kikuchi ladder's next rung — see the weave
reflection); the excess along the well atlas.
