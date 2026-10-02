---
id: THM-4533
title: "Redei graphs: which spanning oriented graphs lie an odd number of times in every tournament. Deletion-reversal generates all parity relations among embedding counts. Path plus span-3 arcs is Redei exactly for n = 4, 5, 6 (D69), failing for all n >= 7 because of the smallest asymmetric tree. D_n = P_n + (0, n-2) + (1, n-1) is Redei for every odd n. Hereditarily symmetric trees have only parity-rigid orientations. The largest shaving on 9 vertices has 14 arcs: u(9) = 14, kappa(9) = 22"
status: >
  PROVED + INDEPENDENTLY AUDITED (statements and finite checks; the D_n proof and the trees theorem read at the level of
  their notes): B1 deletion-reversal and its completeness (B2); the short-certificate lemmas (B3); the necessary conditions
  (B4); heredity (B5); D1, the family D_n (odd n) and H_n (even n); G2 (orientations of C_n are all parity-rigid iff n is
  even); T (hereditarily symmetric trees); D3/D69 (path + span-3 is Redei exactly for n = 4, 5, 6); D4 (doubling).
  PROVED by exhaustive computation (lane; reproduced on re-run): u(9) = 14 and kappa(9) = 22, via an orderly search over
  the symmetric host C3[C3] (|Aut| = 81), confirmed with a second host; r(9) = 11. The orchestrator independently checked
  u(9) >= 14: the witness embeds in all 191536 classes.
  FINITE-EXACT: the Redei-class census 1, 1, 2, 5, 21, 43, 156, 220 for n = 1..8 (orchestrator re-derived n <= 5);
  r(n) = 0, 1, 2, 4, 6, 8, 9, 10, 11 for n = 1..9; 54 classes of maximum 9-shavings; the path + chords Redei lists for n <= 11.
  CITED (literature corrections for THM-4526, opus S15):
  - THM-4526 Theorem A is classical: the Hamiltonian cycle of block type (n-1, 1) in Rosenfeld's problem. El Zein,
    arXiv:2204.11211 (a proof of Rosenfeld's conjecture; exactly 35 exceptions), lists the exceptions of this type as C3
    and the 5-vertex class, credited to Havet 2000; existence goes back to Grunbaum 1971.
  - THM-4526 Theorem C's constant is settled: Linial-Saks-Sos, Combinatorica 3 (1983) 101-104, give
    u(n) = n log n - O(n log log n), so c = 1.
  OPEN (HYP-9173): whether r(n) is of order n log n or linear; "good iff hereditarily symmetric" for all trees; the path +
  chords lists for every n.
  Definitions: for a spanning oriented graph S on [n] and an n-tournament T, emb(S,T) is the number of bijections
  mapping S into T. S is Redei if emb(S,T) is odd for every T, and parity-rigid if emb(S, .) mod 2 is constant.
source: collatz-procgen-20260922 session, shave4 lane (2026-10-01; resumed after a reboot), the owner's "shaved down tournament as fundamental object ... as the size increases", building on opus S15's THM-4526 (directions D69-D72); audited and promoted by the session orchestrator 2026-10-01
depends_on:
  - 01-canon/theorems/THM-4526-shaved-tournaments-path-plus-chord-and-the-largest-unavoidable-graphs.md
  - 01-canon/theorems/THM-4524-selfie-tournaments-loop-gauge-arc-parity-and-the-odd-mallows-sloane-partner.md
related:
  - 05-knowledge/hypotheses/HYP-9173-redei-graphs-growth-and-trees.md
note: 05-knowledge/results/procgen_shave4_20261001_redei_graphs.md
scripts: 04-computation/experiments/procgen_shave4_20261001_run.py, procgen_shave4_20261001_{engine,closure,u9}.c
script_audit: 04-computation/experiments/procgen_shave4_20261001_orchestrator_check.{py,c}
output: 05-knowledge/results/procgen_shave4_20261001.out (+ _u9search.out, _r9search.out)
output_sha256: 9c0c842493ab24980a54a480cb1b0fb9923544c2a5f0d9055d2a40488ed52dac
output_audit: 05-knowledge/results/procgen_shave4_20261001_orchestrator_check.out
output_audit_sha256: a2f9dd837e58e4262c485cd78b7deef1f725a89edadd97832a3d922c5c4d9daf
hash_basis: raw bytes
audit: >
  Literature checked from primary pages:
  - arXiv:2204.11211 is El Zein's proof of Rosenfeld's conjecture ("with exactly 35 exceptions");
  - Linial-Saks-Sos, Combinatorica 3 (1983), states n log n - c1 n >= f(n) >= g(n) >= n log n - c2 n log log n.
  Independent code (procgen_shave4_20261001_orchestrator_check.{py,c}; the lane's code was not read), 14 checks:
  - path + span-3 is odd in every class for n = 4, 5, 6;
  - at n = 7 it is even in 174 of 456 classes and absent from exactly 2;
  - the near-transitive witness has exactly 2 copies;
  - D_5, D_7, D_9 are Redei on all classes (191536 at n = 9);
  - the 14-arc witness embeds in all 191536 9-classes;
  - the Redei census 1, 1, 2, 5, 21 and the maxima 0, 1, 2, 4, 6 for n <= 5, among 1, 2, 7, 42, 582 oriented-graph classes.
  The lane's runner was re-run (494 s, 149 checks, ALL CHECKS PASSED); identical up to timing lines.
  Scope: the non-existence of a 15-arc shaving on 9 vertices (u(9) <= 14) and r(9) = 11 rest on the lane's orderly searches,
  which reproduce S15's n = 7, 8 results.
---

# THM-4533 — Rédei graphs and the largest shaving on 9 vertices

**PROVED / FINITE-EXACT + INDEPENDENTLY AUDITED.** Full note: [procgen_shave4_20261001_redei_graphs](../../05-knowledge/results/procgen_shave4_20261001_redei_graphs.md).

## 1. The owner's 4-vertex object, as n grows

A *shaving* of order n is an oriented graph contained in every tournament on n vertices. The owner's object (a Hamiltonian path plus the arc from first to last) is the prototype.

**H_n itself (S15, THM-4526) is classical.** It is unavoidable except in C3 and one 5-vertex tournament. This is a case of Rosenfeld's problem on oriented Hamiltonian cycles (Grünbaum 1971; Havet 2000; El Zein 2022).

**The largest shavings.** Their size is  (Linial–Saks–Sós 1983), and on 9 vertices the maximum is 14 arcs.

## 2. Rédei graphs

Rédei's theorem says the Hamiltonian path lies an odd number of times in every tournament. Call any oriented graph with this property a *Rédei graph*. Every Rédei graph is a shaving, certified by parity alone.

- **One rule generates all parity relations.** Reversing an arc changes the count by the count of the graph with that arc deleted, mod 2.
- **D69 settled.** The path plus every span-3 arc is Rédei exactly for . It fails from 7 on, and the obstruction is the smallest asymmetric tree, the spider .
- **Infinite families.**
  -  = the path plus the two arcs  and , for every odd ;
  -  for every even ;
  - a doubling construction.
- **Parity certificates fall behind.** The largest Rédei graph has 11 arcs at , against 14 for the largest shaving.
