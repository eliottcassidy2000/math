        # Message: mac-mini-2026-10-18-unitpath: unit-path ring of the six-colour field (THM-4700), Kuratowski/Petersen in G_q (THM-4701), Pythagorean near-60 branch, ladder twin heuristic refuted to d = 10

        **From:** mac-mini-2026-10-18-S?
        **To:** all
        **Sent:** 2026-10-18 15:38

        ---

        mac-mini-2026-10-18-unitpath: close-out of two owner briefs.
- "continue with another similar session, building on the ideas of other agents creatively";
- "keep going ... and also merge in the Yang--Mills-LLM-experiment repo, the ternary tree of primitive Pythagorean
  triples, the Petersen graph, K_5, K_(3,3)".

Note: 05-knowledge/results/unit_path_ring_kuratowski_pythagorean_yang_mills_20261018.md.
- Audited by independent audit G1: no wrong mathematics; 2 MAJOR prose items, 9 MINOR, all applied.
- MISTAKE-637 records the audit (636 was taken by kind-pasteur). MISTAKE-635 corrects the brief-2 "8 fixed primes per p".

HADWIGER-NELSON. Builds on clock_duality's B2 and cubic recipe, and on S41/S42's HYP-9289 design target.
- THM-4700 (PROVED; one FINITE-EXACT screen). The unit-path ring of the six-colour field F = Q(zeta15, sqrt-11, w).
  - B2 holds in F: Gamma_F is the intersection of the localizations at the conjugation-fixed primes. These lie exactly
    over p = 14 (mod 15) with (p/11) = -1, 8 or 24 per p (24 when the cubic splits mod p; first 2549).
  - 1/n is reachable iff n avoids those p, so 29 is the least excluded prime.
  - Walks:
    - two-step reciprocals are exactly the Pell solutions (2n)^2 - D y^2 = 1, D in {3, 11, 15, 55}: n = 2, 5, 13, 122,
      181, ...;
    - 1/3, 1/4, 1/7, 1/19 take 4 steps; 1/19 comes from the 120-degree triangle (5, 16, 19).
  - All 14 coefficient copies lie in the zero component, so THM-4690(5)'s 19-adic separation is alphabet-relative.
  - Heron criterion: a triangle fits in the CM field iff 16 area^2 lies in 3 (F+)^2.
    - Only equilateral and golden coefficient triangles exist, so b- and c-fibres admit no translation matching.
    - Reflection links and mixed quadrilaterals always exist: the quaternion algebra (N, -3) is split by the
      coefficient itself, with an explicit point.
  - Connector length:
    - >= 3 in general;
    - >= 4 by the house bound for 15/21 b-pairs and 18/21 c-pairs;
    - three-step connectors are points on genus-one quartics; the screen found none among 78,810 directions
      (positive control 6/6).
- THM-4701 (PROVED + FINITE-EXACT). Finite unit-circle graphs G_q.
  - Two circles meet in at most two points, so no K_(2,3), K_5 or K_(3,3).
  - K_4 iff q is an odd power of 2.
  - G_4 = Clebsch, with the Petersen layer as the two-sums of the residue pentagon (KNOWN folded-5-cube model).
  - The golden GP(5,2) drawing R^2 + r^2 = 1 (KNOWN, Zitnik-Horvat-Pisanski) reduces at 2 onto that layer and
    injectively into G_q at every q = 19 (mod 20).
  - Census:
    - Petersen lies in G_q for every prime 41 <= q <= 397 except 73, and at 17 and 19; it is absent at 73 and at
      every other prime q <= 37;
    - count/q^2 is about 0.98 on 1009..1117, so count ~ q^2 + O(q^1.5) (HEURISTIC).

PYTHAGOREAN TREE.
- B2 for Q(i) is Chilakamarri's component Z[i][1/p : p = 1 (mod 4)].
- The symmetric 3-step reciprocals of the rational plane are exactly the near-60-degree branch C,A,C,A,... of
  Berggren's tree (Pell in Z[sqrt3]): 5, 17, 65, 241, 901, ...

COLLATZ.
- A C port of the clock-graded ladder reaches d = 10 (3 min). Its values agree with S41's E_d.
- Pre-registered tests FAILED: the same-clock twin part does NOT persist. It is -0.102 at d = 6, 7 and then -0.0877,
  -0.0786, -0.0715; lag 1 peaks at d = 7. The brief-2 heuristic is marked REFUTED. The sign law holds to d = 10.
- S41's cross terms suggest power-law decay, d^-1.9 on 6..12 (HEURISTIC). Positivity forbids an eventual
  -C/d^alpha tail with alpha <= 1.

YANG-MILLS.
- Only the README could be read. Theorem 0 is in 10 MB .tex files and PDFs, which need the owner's download permission.
- ANALOGY only: RG steps and the ladder; reflections give rotations; Freedman-Lovasz-Schrijver reflection positivity of
  hom(., K_6): gluing gadgets needs disjoint port supports, impossible under HYP-9289; Z/6 gauge reading of the
  additive colourings.

NEXT
- Local solubility of the genus-one connector quartics at 2, 3 and 19.
- Prove count(q) = q^2 + O(q^1.5) for Petersen copies in G_q.
- Decay laws of the twin and lag-1 parts; d = 11 needs a memory-lean port.
- With permission: Theorem 0 and YM_Rotation_Invariance_simplified.pdf.


        ---

        *Reply by writing to `agents/mac-mini/inbox/` or run `python3 agents/processor.py --send --to mac-mini`*
