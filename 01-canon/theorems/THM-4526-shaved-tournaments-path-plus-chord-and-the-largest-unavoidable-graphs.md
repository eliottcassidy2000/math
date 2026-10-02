---
id: THM-4526
title: "Shaved tournaments. The Hamiltonian path plus the arc from its first to its last vertex (the bypass D(n,2)) lies in every n-tournament except exactly C3 and C3[1,C3,1] (Grunbaum 1971; independent proof here) (equivalently: these are the only tournaments in which every Hamiltonian path closes into a Hamiltonian cycle), with an odd number of copies when n is even. The largest spanning oriented graphs contained in every n-tournament have u(n) = 1,2,4,6,8,9,11 arcs for n = 2..8 (unique for n <= 6: path + all span-3 arcs; 51 classes at n = 7, all inside the Paley heptagon; 1617 classes at n = 8), so kappa(7) = 12 and kappa(8) = 17; and u(n) = Theta(n log n) (weaker than Linial-Saks-Sos 1983: u(n) ~ n log2 n)."
status: >
  PROVED: Theorem A (all n): a KNOWN theorem (Grunbaum 1971), proved independently here
  (Redei parity for even n; a short constructive argument for odd n). Theorem C (bounds):
  SUPERSEDED by Linial-Saks-Sos 1983. PROVED by exhaustive computation: Theorem B (n <= 8).
  FINITE-EXACT: transversal multiplicities n <= 6; the n <= 6 embedding parity of
  path + span-3 arcs; T(n) not a power of 2 for 5 <= n <= 60.
  INDEPENDENTLY AUDITED (2026-10-01, blind re-derivation, n <= 9 classes regenerated, all
  searches redone): SOUND; corrections applied (attributions, 5 minimum certificates at n = 7,
  wording; MISTAKE-555).
  Notation: an n-shaving is a spanning oriented graph on n vertices contained in every
  n-tournament; u(n) = max arcs = C(n,2) - kappa(n) (kappa of HYP-3798); H_n = path
  v1 -> ... -> vn plus v1 -> vn; T(n) = A000568(n).
  (A) T has no copy of H_n iff T = C3 or C3[1,C3,1]; #copies = H(T) - n hc(T), odd for even n.
  (B) u = 1,2,4,6,8,9,11 (n = 2..8); the n <= 6 maximiser is unique (path + all i -> i+3);
      at n = 7 every shaving embeds in P7 (the only TT4-free class), no acyclic 10-arc
      subgraph of P7 is a shaving, 51 classes of 9-arc shavings; a minimum certificate is
      P7 plus 6 classes, all with |Aut| > 1 (5 such certificates). n = 8: no forward 12-arc graph (C(28,12) = 30421755)
      is a shaving; 48571 forward 11-arc shavings = 1617 classes; kappa(8) = 17 as HYP-3819 predicted.
  (C) (1/2 - o(1)) n log2 n <= u(n) <= C(n,2) - log2 T(n) <= log2 n!; the lazy-caterer
      formula 1 + C(n-2,2) fails both ways (kappa(7) = 12 > 11; kappa(n) < 1 + C(n-2,2)
      for 38 <= n <= 60 and all large n).
  EXTENSION (thread session thread-2cjcob, 2026-10-01; INDEPENDENTLY AUDITED by a blind re-derivation with a
  different search): u(9) = 14, kappa(9) = 22
  (54 classes, all rigid; three hosts agree; a brute-force second proof of u(9) <= 14 and an
  exhaustive re-verification of the 54; the shave4 lane found the same 54 classes independently, THM-4533,
  whose audit covers u(9) >= 14, while the blind audit here re-derived u(9) <= 14 with a host-free search);
  Theorem A FINITE-EXACT for n = 10 (9733056 classes);
  the lazy-caterer formula is exact at n = 9 and too large for every n >= 38; the excess law of
  HYP-3817 holds for n = 3..9 and fails for every n >= 21 with n = 0, 2 (mod 3).
  Note: 05-knowledge/results/shaved_tournaments_proved_vs_conjectured_20261001.md.
source: opus-2026-10-01-S15 (collatz-functional-uniqueness-20261001), owner's shaved 4-tournament prompt; the odd case of (A) was found by a proof-search subagent of the session and checked by the orchestrator and by an independent implementation (Check T3b)
depends_on:
  - Redei 1934 (H(T) odd); Camion 1959 (strong => Hamiltonian); Erdos-Moser (TT_(1+floor(log2 n)) in every n-tournament)
  - Grunbaum, J. Combin. Theory Ser. B 11 (1971) 249-257 (Theorem A, first proof); Thomassen 1980 (>= n - 5 copies in strong T); Benhocine-Wojda, J. Graph Theory 7 (1983) 469-473 (all D(n,p))
  - Linial-Saks-Sos, Combinatorica 3 (1983) 101-104 (u(n) = n log2 n - O(n log log n); supersedes Theorem C)
related:
  - 05-knowledge/hypotheses/HYP-3798-min-free-arcs-transversal-subcube-kappa.md (kappa; lazy-caterer formula: exact for n = 3..6 and 9, fails at n = 7, 8 and for n >= 38; MISTAKE-557 corrected an earlier 'exact only for n <= 6')
  - 05-knowledge/hypotheses/INDEX-HISTORICAL-THROUGH-2026-07-21.md (HYP-3805, opus S15 2026-07-01: kappa(7) = 12 'very likely' -> proved here)
  - 05-knowledge/hypotheses/HYP-3821-biquadratic-field-Q-sqrt-3-sqrt-7-Klein-four-and-sqrt21.md (excess law; kappa(7) = 12 confirmed)
  - 05-knowledge/hypotheses/HYP-3819-excess8-equals-4-proof-strategy-and-sqrt21-bridge.md (predicted excess(8) = 4)
  - 01-canon/theorems/THM-4524-selfie-tournaments-loop-gauge-arc-parity-and-the-odd-mallows-sloane-partner.md (a different 'shaving')
note: 05-knowledge/results/shaved_tournaments_unavoidable_cores_20261001.md
scripts: 04-computation/experiments/shaved_tournaments_20261001.py, 04-computation/experiments/shaved_tournaments_20261001_closing.c, 04-computation/experiments/shaved_tournaments_20261001_u8.c, 04-computation/experiments/shaved_tournaments_n9_20261001.py, 04-computation/experiments/shaved_tournaments_n9_20261001.c (extension)
output: 04-computation/experiments/shaved_tournaments_20261001.out (ALL CHECKS PASSED), 04-computation/experiments/shaved_tournaments_n9_20261001.out (extension; --full, ALL CHECKS PASSED)
---

# THM-4526 — shaved tournaments

**Status: PROVED (A: a known theorem, Grünbaum 1971, with an independent proof; C: superseded by Linial–Saks–Sós
1983) + PROVED by exhaustive computation (B, n ≤ 8); INDEPENDENTLY AUDITED (2026-10-01, corrections applied).**
Full statement, proofs and tables:
[`05-knowledge/results/shaved_tournaments_unavoidable_cores_20261001.md`](../../05-knowledge/results/shaved_tournaments_unavoidable_cores_20261001.md).

## Statement

The owner observed that every 4-tournament contains `A→B, B→C, C→D, A→D`. This is `H_4`: a Hamiltonian path plus
the arc from its first to its last vertex. Its 4 completions are the 4 classes, each exactly once.

**(A) The shape, for every `n` (Grünbaum 1971).** An `n`-tournament (`n ≥ 2`) contains `H_n` (the bypass `D(n,2)`)
unless it is `C3` or `C3[1,C3,1]`.
The number of copies is `H(T) − n·hc(T)`, which is odd when `n` is even (Rédei). Equivalently, the tournaments in
which every Hamiltonian path closes into a Hamiltonian cycle are exactly `C3` and `C3[1,C3,1]`.

**(B) The largest shavings, `n ≤ 8`.** `u(n) = 1, 2, 4, 6, 8, 9, 11` for `n = 2, …, 8`, so
`κ(n) = 0, 1, 2, 4, 7, 12, 17`.
- For `n ≤ 6` the maximiser is unique up to isomorphism: the Hamiltonian path plus all span-3 arcs (the owner's
  `H_4` at `n = 4`).
- At `n = 7` the Paley heptagon is the principal obstruction. It is the only class without a transitive 4-set,
  so every shaving lives inside it. No acyclic 10-arc subgraph of `P_7` is a shaving; killing them needs `P_7`
  plus 6 more classes (5 minimum certificates, all made of symmetric classes).
- There are 51 classes of 9-arc shavings at `n = 7`.
- At `n = 8` no 12-arc graph is a shaving, and there are 1617 classes of 11-arc shavings. So `κ(8) = 17`, the value
  predicted by the excess law (HYP-3819, HYP-3821).

**(C) Growth.** `(1/2 − o(1)) n log₂ n ≤ u(n) ≤ C(n,2) − log₂ T(n) ≤ log₂ n!`. This is weaker than Linial–Saks–Sós
(1983): `u(n) = n log₂ n − O(n log log n)`.

**Extension, `n = 9` (2026-10-01, independently audited).** `u(9) = 14` and `κ(9) = 22`. The maximum 9-shavings are 54 classes,
all rigid, none with a Hamiltonian path; the shave4 lane found the same 54 classes independently
([THM-4533](THM-4533-redei-graphs-parity-of-shaved-tournaments-and-u9.md)). See
[`shaved_tournaments_proved_vs_conjectured_20261001.md`](../../05-knowledge/results/shaved_tournaments_proved_vs_conjectured_20261001.md),
which also proves that the excess law of HYP-3817 fails for every `n ≥ 21` with `n ≡ 0, 2 (mod 3)`.

## Proof sketch of (A), odd `n`, `T` strong

Take a Hamiltonian cycle `c_0 … c_{n−1}`. Call position `i` forward if `c_{i−1} → c_{i+1}`.
1. Two consecutive forward positions give the non-closing path `(c_i, c_{i+2}, …, c_{i−1}, c_{i+1})`.
2. Since `n` is odd, some two consecutive positions `i, i+1` are therefore backward.
3. The insertion lemma, applied to `c_i` in `T − c_i` and to `c_{i+1}` in `T − c_{i+1}`, either produces a
   non-closing path or forces `c_{i+1} ⇒ R ⇒ c_i → c_{i+1}`, where `R` is the rest.
4. In that case `(s, c_i, c_{i+1}, Y)` is non-closing for a suitable `s ∈ R` and Hamiltonian path `Y` of `R − s`,
   unless `R` is a single vertex or a 3-cycle.

Full proof: note §1.

**FORMALIZED 2026-10-01 (collatz-procgen-20260922 lean lane; Lean 4.30 core; package [`04-computation/lean/ProcgenSelfieEdim/`](../../04-computation/lean/ProcgenSelfieEdim/README.md); orchestrator-audited).**
- `every_four_contains_H4`.
- `avoids_H5_iff`: a 5-tournament avoids H_5 iff it is isomorphic to C3[1,C3,1].
- `copiesH_add`: #copies(H_n) + n hc = H.
- Theorem A for even n, through a from-scratch formal proof of Rédei's theorem (`redei`).
- The odd-n case for n > 5 and Theorems B and C are not formalized.

**LITERATURE + UPDATE 2026-10-01 (collatz-procgen-20260922 shave4 lane, THM-4533; orchestrator-audited).**
- **Theorem A is classical.** H_n is the oriented Hamiltonian cycle of block type (n-1, 1). In El Zein's proof of Rosenfeld's conjecture (arXiv:2204.11211, exactly 35 exceptions), the exceptions of this type are C3 and the 5-vertex class. El Zein credits them to Havet (JCTB 80 (2000)), and existence to Grunbaum (JCTB 11 (1971)). The insertion-lemma proof here is an independent proof, and the even-n parity statement is new in form only.
- **Theorem C's constant is settled.** Linial, Saks and Sos, "Largest digraphs contained in all n-tournaments", Combinatorica 3 (1983) 101-104, give u(n) = n log n - O(n log log n), so c = 1.
- **D70 resolved.** u(9) = 14 and kappa(9) = 22, matching the excess-law prediction; there are 54 classes of maximum 9-shavings.
- **D69 resolved.** Path + span-3 is Redei exactly for n = 4, 5, 6.
- **D72.** Rédei graphs are treated in THM-4533.
- **Independent audit.** Theorem A was verified independently for all classes with n <= 9, and its odd-n proof re-derived step by step with no gap found.
