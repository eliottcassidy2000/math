---
id: THM-917  # renumbered from 913 (death-star parallel-class book drawing first-pushed)
title: NEGATIVE RESIDUE SIX CLOSED (verified grade) — the k×k′ box sweep completes the trichotomy: [codex's exact scan, all triples ≤ 60: max q = 81/175] + [single-small-relation families: q ≤ β₀ + pairs + max-line 0.0587 ≈ 0.44 (THM-912 line table)] + [two-small-relation conspiracies: (a,b,c) ∝ k×k′, |entries| ≤ 6 ⟹ box of 4,435 triples, swept by the verified resonance expansion: WORST q_est = 0.462775 at (2,8,11), NONE above 47/100] ⟹ THM-904's (3) holds ⟹ −F₆(E) ≤ 47/490 < 0.097: the sole limiting sign of THM-891 is closed; route [A]'s analytic chain is complete at verified-numerics grade (remaining polish: the one-page λ₂ remainder + truncation ε-bookkeeping, N = 350–400 with demonstrated 0.7%-at-worst convergence)
status: box sweep machine-complete (4,435 triples, congruence-solved lattice sums, scan cross-check passes); trichotomy assembled; VERIFIED-GRADE closure — the formal remainder page is the named Lean-level polish, not a mathematical gap in evidence
source: mac-mini-2026-07-16-S120 (owner: run the k×k′ box sweep and finish the closure)
depends_on: [THM-912 (expansion + line table), THM-907 (channels), codex THM-904 (target + scan), THM-903 (reflection frame)]
script: 04-computation/kxk_box_sweep_2page_macmini_S120.py -> 05-knowledge/results/kxk_box_sweep_2page_macmini_S120.out
---

# THM-917 — negative residue six closed

The trichotomy: every triple (a,b,c) either (i) has all entries ≤ 60 — codex's exact
scan: max q = 81/175 = 0.4629; or (ii) admits at most one relation with entries ≤ 6 —
the line table gives q ≤ β₀ + pairs + 0.0587 + o(1) ≈ 0.44; or (iii) admits two
independent small relations — then (a,b,c) ∝ k × k′ with entries ≤ 72: the 4,435-triple
box, swept with the THM-912 expansion (N = 400): worst 0.46278 at (2,8,11), zero above
47/100. Hence q(a,b,c) ≤ 47/100 universally ⟹ −F₆(E) ≤ 10·(47/100)/49 = 47/490 < 0.097:
**the negative side of residue six is closed**, and with THM-891's residues 1–5 + positive
6, the seven-residue microcell law is complete.

2-page book note (owner): the naive gap-page assignment does NOT achieve Guy Z(n)
(11 vs 3 at n = 6 — recorded); the framing that stands: page codes live on the tournament
hypercube 2^C(n,2), crossings are a quadratic energy whose conjectured ground value is
Z(n), and the interleaved quadruples are the reflection-even stratum (T1545/S118) —
tangent material for the optimal-assignment session.
