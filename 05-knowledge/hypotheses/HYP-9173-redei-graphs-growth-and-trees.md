---
id: HYP-9173
title: "Redei graphs (spanning oriented graphs with an odd embedding count in every n-tournament): (a) the largest has r(n) = o(u(n)) arcs, possibly linear in n; (b) REFUTED: a tree is good (every orientation parity-rigid) iff it is hereditarily symmetric (every subtree has a nontrivial automorphism); (c) for n >= 8 the Redei 'path + chords' graphs are exactly P_n and P_n + (0, n-1) (n even), resp. P_n and P_n + {(0, n-2), (1, n-1)} (n odd)"
status: >
  PARTIALLY REFUTED: (b) is false in its converse direction; (a) and (c) remain OPEN.
  (a) r(n) = 0, 1, 2, 4, 6, 8, 9, 10, 11 for n = 1..9 (FINITE-EXACT), against u(n) = 1, 2, 4, 6, 8, 9, 11, 14 (n = 2..9).
  The gap u - r = 0, ..., 0, 1, 3 opens at n = 8.
  (b) REFUTED in the converse direction by the nine-vertex twin-leaf tree in
  05-knowledge/results/twin_signature_parity_20261019/TWIN_SIGNATURE_PARITY.md. Hereditarily symmetric trees remain good, and
  asymmetric trees remain not good. The equivalence holds through eight vertices (FINITE-EXACT), where a tree is not
  good iff it contains S(1,2,3); nine is the minimal counterexample order. Every orientation of the counterexample has
  even embedding counts, so it is not Redei.
  (c) FINITE-EXACT for n <= 9; PROVED for n = 10, 11 by heredity plus explicit witnesses; reduces to finitely many explicit
  refutations per n.
source: collatz-procgen-20260922 session, shave4 lane (2026-10-01); promoted with THM-4533
related:
  - 01-canon/theorems/THM-4533-redei-graphs-parity-of-shaved-tournaments-and-u9.md
