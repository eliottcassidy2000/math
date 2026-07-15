---
id: THM-853
title: THREE CANON PIECES FROM THE cont.13-20 ARC — (I) THE SELF-LINE SEQUENCE ATLAS to canon (six terms, n=5..10, all machine-exact; Burnside-enumeration engine): ORBITS = 2, 3, 22, 101, 852, 5582; selfK = 4, 6, 44, 202, 1704, 11164; orbitsSC = 2, 0, 9, 3, 69, 53; orbitsNS = 0, 3, 13, 98, 783, 5529; 2selfB = 0, 4, 0, 8, 0, 48; |X| = 8, 16, 88, 412, 3408, 22376; W = 8, 20, 88, 416, 3408, 22408 (odd-n trivial-Aut law x3; even-n Aut-excess 0,4,0,4,0,32; first nontrivial Aut at n=10, histogram {3:12, 5:2}). (II) THE DEEP-WELL CORRIDOR LAW (PROVED, 3-line corollary of THM-826+819): m({1..12}; λ) = (2H₁₂/13)·(1−13λ) EXACTLY on the whole corridor [1/14, 1/13] — every LRC(14) deep-well threshold question in the covering corridor is ONE linear formula; at the covering-min λ = 14/183 (THM-724/726) the margin measure is 2H₁₂/2379 = 6617/2536380 (referee-exact). (III) THE LOCKER (DIVISIBILITY) TOURNAMENT D_n — binary conditions as edges — tile (x,y) := [y | x] (the covering lattice as a canonical tournament; vertex 1 loses every tile). First invariants (n=5..9): c₃ = 2, 3, 6, 10, 16; H = 5, 9, 33, 109, 469; never SC/quasi-fixed/gridsym. CONJECTURE (the locker parity law): H(D_n) ≡ 1 (mod 4) for ALL n — five-for-five OUTSIDE the THM-791 grid-sym stratum; conjectured mechanism = the divisor-pairing involution d ↔ x/d on odd cycles (the locker/τ-parity fact that exactly the SQUARES end toggled)
status: (I) machine-exact censuses (per-rung scripts cited below); (II) PROVED + referee-exact at four corridor points; (III) defined + computed n=5..9, parity law CONJECTURED
source: kind-pasteur-2026-07-15-S128 (cont.20; owner brief - extend sequences to canon, Farey-14 law, locker/squares, binary conditions as edges)
depends_on:
  - THM-826 / THM-819   # (II) is their k=12 corridor corollary
  - THM-852             # (I) structural frame
related:
  - THM-724/726 (the covering-min the corridor law prices exactly), THM-466 (the 2-adic H frame behind (III))
  - HYP-6935 (the creative bundle - 3-triangular Gauss = the three staircase Venn corners since tile counts ARE triangular numbers; the Burnside engine is fast because twisted-fixed-point equations are XOR-SAT / Schaefer-affine, the P-slice of a GI-hard census, mirroring THM-826 disjoint-vs-Venn; Smith chart = Moebius nomogram meeting the Smith diagram at Stern-Brocot)
---

# THM-853 — atlas, corridor, locker

(I) The atlas rows are canon-fixed as in the title. Provenance: thm852_selfline_sc_bijection
(c13), selfline_sequence_atlas (c15), orbit_decode_carriers (c16), orbit_decode_n8 (c17),
orbit_decode_n9_burnside (c18), orbit_decode_n10_burnside (c19), all with .out files.

(II) Proof. THM-826's k=12 profile on [1/14, 1/13) keeps exactly the gaps with i+j = 13 (the
mediant bound gives 12 < i+j ≤ 13 there). Since 13 is prime, the pair sum is Σ 1/(ij) = 2H₁₂/13
(THM-819's unit inversion), so m(λ) = (2H₁₂/13)(1 − 13λ) on the whole corridor — which contains
every deep-well threshold between the LRC(14) bound 1/14 and the AP tightness 1/13, including the
covering-min 14/183 where m = 2H₁₂/2379. Referee: exact at λ = 1/14, 14/183, 15/196, 1/13. ∎

(III) The locker tournament encodes the divisibility relation — the binary conditions (q | w) the
covering theory quantifies over — as one canonical tournament per n. Its H-parity conjecture is
the tournament shadow of the locker-problem fact that τ(n) is odd iff n is a square.

## Evidence log
- [x] (II) referee exact (locker_tournament_and_corridor_kps_S128c20.py + .out)
- [x] (III) invariants n=5..9; parity law 5/5
- [x] (III) parity law RESOLVED NEGATIVELY (mac-mini-2026-07-15-S109, THM-865): H(D_11) = 4027 ≡ 3
      (mod 4) — the conjecture is FALSE, first failure exactly at n = 11. H row extended to n = 15
      (1721, 4027, 28851, 83817, 400569, 3141317 for n = 10..15); H mod 4 = 3 at n ∈ {11, 12, 16..19}.
      The divisor-pairing mechanism post-mortem (real pair-matching at composite m ≤ 18, dies at
      m = 20; primes lawless) is in THM-865 (iii)–(iv).
