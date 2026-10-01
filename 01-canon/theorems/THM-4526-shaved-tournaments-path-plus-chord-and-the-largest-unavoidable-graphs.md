---
id: THM-4526
title: "Shaved tournaments. The Hamiltonian path plus the arc from its first to its last vertex lies in every n-tournament except exactly C3 and C3[1,C3,1] (equivalently: these are the only tournaments in which every Hamiltonian path closes into a Hamiltonian cycle), with an odd number of copies when n is even. The largest spanning oriented graphs contained in every n-tournament have u(n) = 1,2,4,6,8,9 arcs for n = 2..7 (unique for n <= 6: path + all span-3 arcs; 51 classes at n = 7, all inside the Paley heptagon), so kappa(7) = 12; and u(n) = Theta(n log n)."
status: >
  PROVED: Theorem A (all n; the odd case by a short constructive argument),
  Theorem C (bounds). PROVED by exhaustive computation: Theorem B (n <= 7).
  FINITE-EXACT: transversal multiplicities n <= 6; the n <= 6 embedding parity of
  path + span-3 arcs; T(n) not a power of 2 for 5 <= n <= 60.
  INDEPENDENT AUDIT: OWED.
  Notation: an n-shaving is a spanning oriented graph on n vertices contained in every
  n-tournament; u(n) = max arcs = C(n,2) - kappa(n) (kappa of HYP-3798); H_n = path
  v1 -> ... -> vn plus v1 -> vn; T(n) = A000568(n).
  (A) T has no copy of H_n iff T = C3 or C3[1,C3,1]; #copies = H(T) - n hc(T), odd for even n.
  (B) u = 1,2,4,6,8,9 (n = 2..7); the n <= 6 maximiser is unique (path + all i -> i+3);
      at n = 7 every shaving embeds in P7 (the only TT4-free class), no acyclic 10-arc
      subgraph of P7 is a shaving, 51 classes of 9-arc shavings; a minimum certificate is
      P7 plus 6 classes, all with |Aut| > 1.
  (C) (1/2 - o(1)) n log2 n <= u(n) <= C(n,2) - log2 T(n) <= log2 n!; the lazy-caterer
      formula 1 + C(n-2,2) fails both ways (kappa(7) = 12 > 11; kappa(n) < 1 + C(n-2,2)
      for 38 <= n <= 60 and all large n).
source: opus-2026-10-01-S15 (collatz-functional-uniqueness-20261001), owner's shaved 4-tournament prompt; the odd case of (A) was found by a proof-search subagent of the session and checked by the orchestrator and by an independent implementation (Check T3b)
depends_on:
  - Redei 1934 (H(T) odd); Camion 1959 (strong => Hamiltonian); Erdos-Moser (TT_(1+floor(log2 n)) in every n-tournament)
related:
  - 05-knowledge/hypotheses/HYP-3798-min-free-arcs-transversal-subcube-kappa.md (kappa; lazy-caterer formula, now exact only for n <= 6)
  - 05-knowledge/hypotheses/INDEX-HISTORICAL-THROUGH-2026-07-21.md (HYP-3805, opus S15 2026-07-01: kappa(7) = 12 'very likely' -> proved here)
  - 05-knowledge/hypotheses/HYP-3821-biquadratic-field-Q-sqrt-3-sqrt-7-Klein-four-and-sqrt21.md (excess law; kappa(7) = 12 confirmed)
  - 05-knowledge/hypotheses/HYP-3819-excess8-equals-4-proof-strategy-and-sqrt21-bridge.md (predicted excess(8) = 4)
  - 01-canon/theorems/THM-4524-selfie-tournaments-loop-gauge-arc-parity-and-the-odd-mallows-sloane-partner.md (a different 'shaving')
note: 05-knowledge/results/shaved_tournaments_unavoidable_cores_20261001.md
scripts: 04-computation/experiments/shaved_tournaments_20261001.py, 04-computation/experiments/shaved_tournaments_20261001_closing.c
output: 04-computation/experiments/shaved_tournaments_20261001.out (ALL CHECKS PASSED)
---

# THM-4526 — shaved tournaments

**Status: PROVED (A, C) + PROVED by exhaustive computation (B); independent audit OWED.**
Full statement, proofs and tables:
[`05-knowledge/results/shaved_tournaments_unavoidable_cores_20261001.md`](../../05-knowledge/results/shaved_tournaments_unavoidable_cores_20261001.md).

## Statement

The owner observed that every 4-tournament contains `A→B, B→C, C→D, A→D`. This is `H_4`: a Hamiltonian path plus
the arc from its first to its last vertex. Its 4 completions are the 4 classes, each exactly once.

**(A) The shape, for every `n`.** An `n`-tournament (`n ≥ 2`) contains `H_n` unless it is `C3` or `C3[1,C3,1]`.
The number of copies is `H(T) − n·hc(T)`, which is odd when `n` is even (Rédei). Equivalently, the tournaments in
which every Hamiltonian path closes into a Hamiltonian cycle are exactly `C3` and `C3[1,C3,1]`.

**(B) The largest shavings, `n ≤ 7`.** `u(n) = 1, 2, 4, 6, 8, 9` for `n = 2, …, 7`, so `κ(n) = 0, 1, 2, 4, 7, 12`.
- For `n ≤ 6` the maximiser is unique up to isomorphism: the Hamiltonian path plus all span-3 arcs (the owner's
  `H_4` at `n = 4`).
- At `n = 7` the Paley heptagon decides. It is the only class without a transitive 4-set, so every shaving lives
  inside it. No acyclic 10-arc subgraph of `P_7` is a shaving.
- There are 51 classes of 9-arc shavings at `n = 7`.

**(C) Growth.** `(1/2 − o(1)) n log₂ n ≤ u(n) ≤ C(n,2) − log₂ T(n) ≤ log₂ n!`.

## Proof sketch of (A), odd `n`, `T` strong

Take a Hamiltonian cycle `c_0 … c_{n−1}`. Call position `i` forward if `c_{i−1} → c_{i+1}`.
1. Two consecutive forward positions give the non-closing path `(c_i, c_{i+2}, …, c_{i−1}, c_{i+1})`.
2. Since `n` is odd, some two consecutive positions `i, i+1` are therefore backward.
3. The insertion lemma, applied to `c_i` in `T − c_i` and to `c_{i+1}` in `T − c_{i+1}`, either produces a
   non-closing path or forces `c_{i+1} ⇒ R ⇒ c_i → c_{i+1}`, where `R` is the rest.
4. In that case `(s, c_i, c_{i+1}, Y)` is non-closing for a suitable `s ∈ R` and Hamiltonian path `Y` of `R − s`,
   unless `R` is a single vertex or a 3-cycle.

Full proof: note §1.
