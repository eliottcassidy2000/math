---
id: THM-4531
title: "The twisted Mallows-Sloane principle, and why OPEN-Q-060 has no natural bijection. For every S_n-invariant subspace U of the edge space, tournaments modulo reversal along U, up to isomorphism, are equinumerous with the even (untwisted) graphs in U-perp: U = 0 is Royle-Praeger-Glasby-Freedman-Devillers (2023), U = cuts is THM-4524 E1. A natural map is bijective on iso types iff a stabiliser-containment bipartite graph has a perfect matching. No natural map sends even Euler graphs to switching classes (n >= 3), and none sends switching classes bijectively onto even Euler graphs for 5 <= n <= 100. The same failure holds for tournaments versus even graphs (answering Royle et al.'s question in the equivariant sense for 5 <= n <= 15 and 16 larger n) and for Mallows-Sloane at every even n >= 4. For odd n, A049313(n) counts the tournaments with out-degree = in-degree (mod 4) everywhere. Higashitani-Ueyama's expectation s_{4,n} = t_{4,n} holds, and s_{l,n} = t_{l,n} for every l"
status: >
  PROVED + INDEPENDENTLY AUDITED:
  - Theorem G (twisted Mallows-Sloane principle; monomial-representation proof);
  - Lemma N1 (matching criterion);
  - Theorem O1 (no natural map from even Euler graphs to classes, n >= 3);
  - Theorem O2 at n = 5 (proof by hand);
  - Theorem O4 (Mallows-Sloane, all even n >= 4);
  - P1 (odd n: the canonical mod-4 Eulerian member; extends THM-1470);
  - P2 (bipartite Euler graphs are even);
  - P3 (s_{l,n} = t_{l,n} for all l, n; Brauer's lemma for finite abelian groups);
  - Theorem B (THM-479 branch integrality for n = 2 mod 4: N_lev(n) = #{classes with even-order Aut}/2);
  - S (odd n: self-converse switching classes on n correspond to self-converse tournaments on n - 1).
  PROVED by certificates (lane; the n = 5, 6 cases re-derived by the orchestrator): O2 for 5 <= n <= 100 (the block lemma;
  two forced classes per n); O3 for 5 <= n <= 15 and n in {22, 23, 26, 31, 32, 33, 38, 46, 47, 61, 62, 63, 71, 74, 83, 86}.
  FINITE-EXACT (lane): the exhaustive matchings for n <= 10 (CLASS) and n <= 9 (TOUR).
  CITED:
  - Royle, Praeger, Glasby, Freedman, Devillers, "Tournaments and even graphs are equinumerous", J. Algebraic Combin. 57
    (2023), arXiv:2204.01947, proving von Broemssen's conjecture;
  - Mallows-Sloane 1975;
  - Babai-Cameron 2000;
  - Higashitani-Ueyama arXiv:2409.10904v4 (Example 4.10 states s_{4,n} = t_{4,n} as an expectation without proof);
  - Seidel; THM-1470; THM-479.
  OPEN (HYP-9172): no natural bijection for every n >= 5 (CLASS and TOUR), which reduces to additive questions over
  twist-rigid block sizes; branch integrality for n = 0 (mod 4); a natural graph-type partner of A049313 for even n.
  Notation: "natural" means S_n-equivariant on labelled objects (a species morphism); even means sgn_X = 1 on Aut X, where
  sgn_X(g) = (-1)^#{edges of X reversed by g relative to a fixed orientation}. This is THM-4524's eps and Royle et al.'s
  sign.
source: collatz-procgen-20260922 session, tbij lane (2026-10-01; resumed after a reboot and a network drop), the bijective form of OPEN-Q-060 requested by the owner; audited and promoted by the session orchestrator 2026-10-01
depends_on:
  - 01-canon/theorems/THM-4524-selfie-tournaments-loop-gauge-arc-parity-and-the-odd-mallows-sloane-partner.md (E1)
  - 01-canon/theorems/THM-479-level-holonomy-and-branch-split.md (N_odd + N_lev)
  - 01-canon/theorems/THM-1470-even-tournaments-are-the-tournament-two-graph-theorem.md (score-parity law)
related:
  - 00-navigation/OPEN-QUESTIONS.md (OPEN-Q-060)
  - 05-knowledge/hypotheses/HYP-9172-no-natural-bijection-for-twisted-mallows-sloane-at-any-n-ge-5.md
note: 05-knowledge/results/procgen_tbij_20261001_natural_bijections.md
scripts: 04-computation/experiments/procgen_tbij_20261001_{run,lib,nat,twist,hu,blocks}.py
script_audit: 04-computation/experiments/procgen_tbij_20261001_orchestrator_check.py
output: 05-knowledge/results/procgen_tbij_20261001.out
output_sha256: 58e27d3ae9f28de41fffa59d9997af6d6f55511926979ec08cac8250ff3ff917
output_audit: 05-knowledge/results/procgen_tbij_20261001_orchestrator_check.out
output_audit_sha256: d450fc8faf49b8bc32e157cab67dd8c355b8758d0ac5326e1ed9e2ccdf0a0141
hash_basis: raw bytes
audit: >
  The orchestrator read and found sound:
  - Theorem G (the characters phi_X, X in U-perp, give a monomial basis of C[T_U]; Aut X acts on the line of phi_X by
    sgn_X; Frobenius reciprocity);
  - Lemma N1 (equivariance forces Stab(x) <= Stab(f(x)); a matching gives a well-defined map);
  - Theorem O1;
  - the n = 5 proof of O2;
  - P1 (via THM-1470's parity law);
  - P2 (orient each component across its bipartition);
  - P3 (the pairing <X_v, N> = row sum identifies modular Eulerian matrices with the dual of the switching classes;
    |A^g| = |A-hat^g| for a finite abelian group);
  - the structure of the Theorem B argument: Babai-Cameron Thm 5.1 makes 2-subgroups semiregular, hence of order <= 2
    when n = 2 mod 4, so rho = 1/2; tau has no fixed point by the oriented two-graph computation.
  Literature checked from the arXiv abstract pages: arXiv:2204.01947 is Royle-Praeger-Glasby-Freedman-Devillers
  (von Broemssen's conjecture); Higashitani-Ueyama arXiv:2409.10904v4 contains the sentence expecting
  s_{4,n} = t_{4,n} "though we do not have a proof", with t_4 = 1, 1, 3, 8, 62, 1760.
  Independent code (procgen_tbij_20261001_orchestrator_check.py; the lane's code was not read), 27 checks in 27 s:
  - the stabiliser-containment matching at n = 5 (2 class types, matching 1) and n = 6 (6 class types, matching 5);
  - no labelled class is S_n-fixed (n = 3..6);
  - the n = 5 hand-proof data (Euler graphs with an order-3 / order-5 automorphism: 4 / 3, only the empty one even);
  - bipartite Euler graphs are even (n <= 6);
  - #even graphs = #tournaments (n <= 5);
  - s = t by brute-force union-find for l = 3, 4 and n <= 4 (t_4 = 1, 1, 3, 8);
  - P1 on all switching classes for n = 3, 5, 7.
  THM-479's recorded N_lev(6) = 2, N_lev(8) = 15 and N_lev(10) = 280 agree with Theorem B and with the lane's n = 8
  analysis.
  The lane's runner was re-run in full (--full: 965 s, 39749 checks, ALL CHECKS PASSED); identical up to timing lines
  (peak RSS 634 MB, above the 500 MB lane guideline).
  Scope notes: the certificates for 7 <= n <= 100 (block lemma; nauty-checked non-isomorphism), the O3 certificates and
  the n <= 10 exhaustive matchings are lane-level.
---

# THM-4531 — the twisted Mallows–Sloane principle and the non-existence of natural bijections

**PROVED + INDEPENDENTLY AUDITED.** Full note: [procgen_tbij_20261001_natural_bijections](../../05-knowledge/results/procgen_tbij_20261001_natural_bijections.md).

## 1. One principle behind the even/odd counts

Fix a sign on graphs: `sgn_X(g)` is the parity of the edges of `X` that the automorphism `g` reverses. Call `X` *even* if this sign is trivial on `Aut X`.

**Theorem G.** For every relabelling-invariant subspace `U` of graphs, tournaments taken modulo reversal along members of `U` are equinumerous, up to isomorphism, with the even graphs orthogonal to `U`.
- **`U = 0`.** This is the 2023 theorem of Royle, Praeger, Glasby, Freedman and Devillers: tournaments and even graphs are equinumerous.
- **`U` = the cut space.** This is THM-4524's E1: switching classes of tournaments and even Euler graphs are equinumerous.

So E1 is the switching-class form of Royle et al.'s theorem and is proved by the same mechanism. THM-4524 is corrected accordingly.

## 2. Is there a one-to-one matching? (OPEN-Q-060)

Read "natural" as "given by a construction that respects relabelling". Then the answer is **no**:
- **From even Euler graphs to classes, never (`n >= 3`).** The empty graph is fixed by every permutation, and no switching class is.
- **From classes onto even Euler graphs, not for `5 <= n <= 100`.** At `n = 5` the two classes have symmetries of order 5 and 3. The only even Euler graph with such a symmetry is the empty graph, so both classes would have to map there.
- **Royle et al.'s own open question.** The same failure holds for tournaments versus even graphs, which answers their question negatively (in this sense) for `5 <= n <= 15` and 16 larger `n`.
- **Mallows and Sloane.** Their "unable to find such a correspondence" at even `n` is explained exactly: none exists, for every even `n >= 4`.

**For odd `n` there is a canonical object on the tournament side.** Every switching class contains exactly one tournament with out-degree ≡ in-degree (mod 4) at every vertex. So A049313 counts tournaments whose scores are all even (`n ≡ 1 mod 4`) or all odd (`n ≡ 3 mod 4`).

## 3. Two byproducts

- **Higashitani–Ueyama.** Their 2024–25 paper expects, without proof, that their switching classes and modular Eulerian matrices over `Z/4` are equinumerous. This holds, and it holds over every `Z/l`, by Brauer's lemma for finite abelian groups.
- **THM-479's branch split.** It gets a combinatorial meaning for `n ≡ 2 (mod 4)`: the level branch counts the pairs `{C, C + M_g}`.
