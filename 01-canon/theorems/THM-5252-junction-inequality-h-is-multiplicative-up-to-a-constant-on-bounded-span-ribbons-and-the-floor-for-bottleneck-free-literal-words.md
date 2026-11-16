---
id: THM-5252
title: "Junction inequality. (A) For a tournament T and a partition V = A + B with back arcs (from B to A) having heads A0 and tails B0: H(T[A]) H(T[B]) <= H(T) <= k+(A) k-(B) H(T[A]) H(T[B]), k+(A) = sum over R in A0 of |R|! prod 2^(outdegree in T[A]), k-(B) likewise with indegrees in T[B]. For span-w ribbon words: H(u) H(v) <= H(uv) <= K_w^2 H(u) H(v), K_3 = 26,425. (B) Every raw coordinate of the span-3 transfer vector of W is at most K_3 H(W). (C) A span-3 factor of m vertices whose minimal operator has rank >= 2 has a strong component with at least m - 12 vertices, so H >= 5^((m-13)/3); the least count g_m over such factors is supermultiplicative, rho_min = lim g_m^(1/m) = sup exists, and 5^(1/3) <= rho_min <= ((9+sqrt 113)/2)^(1/4). (D) rho_min > 5^(1/3) iff some g_m > 5^(m/3) iff the tree of slow factors is finite at some rate above 5^(1/3). (E) Exhaustive: among primitive periodic span-3 words of period 2-8 the only bottleneck-free cycles with growth below 1.7705 are 30 necklaces of period 4 (nonzero spectrum x^2 - 9x - 8) and 185 of period 8 (x^2 - 97x + 64, the same growth)"
status: >
  PROVED: (A), (B), (C), (D). (A) uses the abundancy inequality H(T) >= H(T[S]) H(T - S) (KNOWN: odd-cycle
  form of Grinberg-Stanley, THM-002; THM-5130 (H4)) and an iterated tail deletion; for one target vertex the tail
  deletion is the argument of chanpak_target_ribbon_audit_20261116/SHARP_GROWTH.md (another session, same day),
  which proves the sharp one-letter constant 39; the case of several targets, the two-sided statement for two
  blocks of arbitrary length, and (B) are added here. (C) uses Busch's lower bound for strong tournaments (CITED,
  as in THM-4810) and the elementary Lemma 2 of the note (an interior cut forces operator rank <= 1). (D) is an
  equivalence only: rho_min > 5^(1/3), which is HYP-9742 for literal span 3, is NOT proved.
  FINITE-EXACT: (E) (stable ranks exact with python-flint; growth values in floating point; pruning exact by
  abundancy); the observed junction constant (largest H(uv)/(H(u)H(v)) found: 3,965; exhaustive for |u| <= 6,
  |v| <= 4, beam search to |v| = 12); g_m = 1 for m <= 8, 3 for m = 9, 10, 11, 5 for m = 12, g_13 >= 9 (ranks
  modulo one prime); first six levels of the slow-factor tree at rate 1.72.
  NOT CLAIMED: Chan-Pak (OPEN); the literal span-3 lower-growth certificate (NOT built); density statements;
  sharpness of K_3^2.
  Not independently audited (owner directive 2026-11-13: ideas over audits).
session: mac-mini-2026-11-15-junction
date: 2026-11-15
source: 05-knowledge/results/chanpak_junction_inequality_20261115/README.md
scripts:
  - 05-knowledge/results/chanpak_junction_inequality_20261115/scripts/qm_ratio.py, qm_beam.py, cutcheck.py,
    gmin.py, badtree.py, cheap_periodic.py
related:
  - HYP-9802, HYP-9742, THM-5212, THM-5130, THM-4810, THM-4880, THM-002
  - 05-knowledge/results/chanpak_target_ribbon_audit_20261116/SHARP_GROWTH.md
  - 05-knowledge/results/chanpak_eight_braid_20261116/RETURN_GROWTH.md
---

# THM-5252 - the junction inequality

Source: [the note](../../05-knowledge/results/chanpak_junction_inequality_20261115/README.md).

`H(T)` is the number of Hamiltonian paths of a tournament.

**Lemma 1.** For distinct vertices `y_1, ..., y_r` of a tournament `D`, the number of tuples `(F_0, ..., F_r)` of
pairwise disjoint paths covering `D`, `F_0` possibly empty, `F_i` starting at `y_i`, is at most
`prod 2^(d+(y_i)) H(D)`.

**(A)** Let `V(T) = A + B`, `A0` the heads and `B0` the tails of the arcs from `B` to `A`. Then
`H(T[A]) H(T[B]) <= H(T) <= k+(A) k-(B) H(T[A]) H(T[B])`, with
`k+(A) = sum_{R in A0} |R|! prod_{y in R} 2^(d+_{T[A]}(y))` and `k-(B) = sum_{R in B0} |R|! prod_{z in R} 2^(d-_{T[B]}(z))`.
For words of span `w`: `H(u) H(v) <= H(uv) <= K_w^2 H(u) H(v)`, `K_w = sum_{R in {0..w-1}} |R|! prod_{i in R} 2^(w+i)`.

**(B)** In span 3 every raw coordinate of the 119-state transfer vector of `W` is at most `K_3 H(W)`, and
`H(W s) <= K_3^2 H(s) H(W)` for a continuation `s` of any length.

**(C)** If a factor `u` of `m` vertices has `rank A_u >= 2` in the minimal representation, then `T[u]` has a strong
component with at least `m - 12` vertices and `H(T[u]) >= 5^((m-13)/3)`. With `g_m` the least `H(T[u])` over such
factors: `g_(m+m') >= g_m g_m'`, `rho_min = lim g_m^(1/m) = sup_m g_m^(1/m)`, and
`5^(1/3) <= rho_min <= ((9 + sqrt 113)/2)^(1/4)`.

**(D)** Equivalent: `rho_min > 5^(1/3)`; `g_m > 5^(m/3)` for some `m`; for some `rho > 5^(1/3)` only finitely many
factors `u` of rank at least 2 have `H(T[u']) < rho^|u'|` for all their prefixes `u'`.

**(E)** See the table in section 4 of the note.

Proofs: sections 1 and 2 of the note.
