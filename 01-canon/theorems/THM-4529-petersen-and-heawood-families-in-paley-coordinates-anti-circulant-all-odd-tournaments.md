---
id: THM-4529
title: "The Petersen and Heawood families in Paley coordinates; anti-circulant tournaments; an all-odd tournament on 14 vertices; the Gale-tournament form of the linear Conway-Gordon-Sachs theorem. Delta-Y along the 4 Fano translates avoiding 0 turns K6 = P7 - 0 into the Petersen graph, and along all 7 translates turns K7 into the Heawood graph. In every anti-circulant tournament the antipodal arcs lie on an odd number of Hamiltonian paths. QR_127 restricted to mu_14 has every arc on an odd number of HPs, settling HYP-9167(b) at N = 14. Every Cayley digraph of an odd-order abelian group is arc-even. A linear K6 has 3 linked triangle pairs iff its Gale tournament is C3[TT2,TT2,TT2]"
status: >
  PROVED + INDEPENDENTLY AUDITED (Theorem 6.1 anti-circulant antipodal parity; Theorem 5.1 Cayley digraphs; Prop. 3.2 Redei
  graphs are complete graphs; the Paley-coordinate constructions; Theorem L, the Gale/walk form of linear CGS for K6, and its
  n-dimensional form L'); FINITE-EXACT + AUDITED (the Petersen family has 7 members and the Heawood family 20, with 14
  Delta-Y descendants of K7; QR_127[mu_14] is all-odd with H = 24540117; one all-odd anti-circulant class at each of N = 6, 10, 14;
  QR_p[mu_6] all-odd iff p = 7 mod 8 and QR_p[mu_10] all-odd iff p = 3 mod 8, both proved by Legendre-sign algebra and checked for
  p < 3000); FINITE-EXACT, lane only (anti-circulant censuses at N = 18, 22, 26 and their prime identifications; Fano-colouring
  counts); NUMEROLOGY (a bijection "7 translations <-> 7 Petersen-family members": no member has a symmetry of order 7);
  REFUTED (a Hamiltonian-path parity preserved along the family); OPEN (HYP-9170: an all-odd anti-circulant for every
  N = 2 mod 4, and infinitely many all-odd QR_p[mu_2m]; Collatz is untouched).
  (F) The Delta-Y / Y-Delta closure of K6 is the Petersen family: K6, P7, K3,3,1, P8, K4,4-e, P9 and the Petersen graph
  (15 edges each; Aut orders 720, 36, 72, 8, 72, 12, 120). The Delta-Y Hasse diagram is a tree with 6 arrows; every member
  except K3,3,1 descends from K6.
  (P) Paley coordinates on Z/7, with D = {1,2,4}:
  - Delta-Y on all 7 translates D + t in K7 gives the Heawood graph, the Fano incidence graph;
  - in K7 - 0 = K6, the 4 translates avoiding 0 (t in {0} u QR7) are edge-disjoint triangles, and Delta-Y on them gives the
    Petersen graph;
  - the 3 translates through 0 (t in NQR7) leave the perfect matching {u, 3u};
  - Delta-Y on the two trivial-cycle codes {1,2,4} and {3,5,6} gives K4,4 - e;
  - Aut(P7 - 0) = <x -> 2x>, the trivial cycle's dynamics.
  (5.1) For every abelian group G of odd order and every connection set C, every arc of Cay(G, C) lies on an even number of
  directed Hamiltonian paths. The involution is P -> reverse(phi(P)) with phi(x) = u + v - x.
  (6.1) A tournament on N vertices with an anti-automorphism that is one N-cycle (anti-circulant) needs N = 2 mod 4. All its
  N/2 antipodal arcs lie on an odd number of HPs. Examples: QR_q - 0, and QR_p[mu_2m] for p = 3 mod 4.
  (N14) QR_127 restricted to mu_14 = +-<2> = {+-1, +-2, ..., +-64}: H = 24540117, all 91 arc counts odd, |Aut| = 7.
  (L) Six points in general position in R^3; Gale vectors g_i in R^2; T*: i -> j iff det(g_i, g_j) > 0. Then:
  - T* is locally transitive;
  - for a split into triangles A | B, the number of edges of B piercing A equals the number of b in B that are
    "uniform-opposite" in T*;
  - the number of linked pairs is 1 or 3, and it is 3 iff T* is isomorphic to C3[TT2,TT2,TT2], one of the two 6-vertex
    H-maximizers;
  - P7 - v, the other H-maximizer and the all-odd one, is not the Gale tournament of any linear K6.
  (L') The same in odd dimension n with r = (n+3)/2: the number of linked pairs is odd and at most r (r odd) or r - 1 (r even).
source: collatz-procgen-20260922 session, petersen lane (2026-10-01; resumed after a reboot and a network drop), answering the owner's "family of 7 Petersen graphs / Paley parity bridge" prompt; audited and promoted by the session orchestrator 2026-10-01
depends_on:
  - 01-canon/theorems/THM-4524-selfie-tournaments-loop-gauge-arc-parity-and-the-odd-mallows-sloane-partner.md (arc-HP parity; QR7 - v all-odd; Cayley all-even C1)
  - 05-knowledge/results/collatz_paley_bridge_20261001.md (opus S15, nineteenth note: the trivial cycle's codes QR7 / NQR7)
related:
  - 05-knowledge/hypotheses/HYP-9167-paley-minus-a-vertex-is-redei-rigid.md (N = 14 settled here)
  - 05-knowledge/hypotheses/HYP-9170-all-odd-anti-circulant-tournaments-exist-for-every-n-2-mod-4.md
  - 05-knowledge/results/collatz_label_map_paley_snarks_20261001.md (S15: snarks, the Paley restatement of Collatz; not repeated)
note: 05-knowledge/results/procgen_petersen_20261001_petersen_heawood_paley.md
scripts: 04-computation/experiments/procgen_petersen_20261001_{run,lib,census26}.py, procgen_petersen_20261001_{hp,par,anti}.c
script_audit: 04-computation/experiments/procgen_petersen_20261001_orchestrator_check.{py,c}
output: 05-knowledge/results/procgen_petersen_20261001.out
output_sha256: b466f9493a92d36aff1dfeaf96f28708fc2dcf05248dc9d7e9f51b5e557e53bd
output_census26: 05-knowledge/results/procgen_petersen_20261001_census26.out (sha256 9819f36ab87785c426486fad18332e7ea13827b502f78038dcd08b49fdaa6ab7)
output_audit: 05-knowledge/results/procgen_petersen_20261001_orchestrator_check.out
output_audit_sha256: 5c308398ff268ca1da007cf08e9fbaa9233e39f2a014d32a51f4fba441ffc3a9
hash_basis: raw bytes
audit: >
  Read and found sound by the orchestrator: Theorem 6.1 (sigma = reverse o (x -> x+1) acts on HPs; the only arcs with
  non-trivial stabiliser are the antipodal ones, which form one orbit of odd size m; then Redei); the N = 2 mod 4
  characterisation; Theorem 5.1; Prop. 3.2; the proof of Theorem L (Radon partitions via Gale duality, Gordan's theorem,
  the walk-crossing parity).
  Independent code (procgen_petersen_20261001_orchestrator_check.{py,c}; the lane's code was not read), 36 checks in 8 s:
  - QR_127[mu_14] (H, all 91 counts odd, the 7 values, |Aut| = 7);
  - every anti-circulant sign pattern at N = 6, 10, 14: tournament property, anti-automorphism, Theorem 6.1, exactly one
    all-odd isomorphism class;
  - QR_q - 0 (q = 7, 11, 19) and QR_127[mu_14] are anti-circulant;
  - the mu_6 / mu_10 criteria on 108 and 55 primes p < 3000;
  - Theorem 5.1 on 36 random Cayley digraphs (32 not tournaments);
  - Prop. 3.2 on 40 random non-complete graphs;
  - the Delta-Y closures of K6 (7 graphs, Aut orders, tree diagram, K3,3,1 the only non-descendant) and K7 (20, 14);
  - the Paley-coordinate constructions (Heawood, Petersen, matching {u, 3u}, K4,4 - e, |Aut(P7 - 0)| = 3);
  - Theorem L on 1500 random integer configurations (linked counts {1: 1098, 3: 402}; T* locally transitive; 3 linked iff
    T* = C3[TT2,TT2,TT2]; the piercing rule on all 15000 splits).
  The lane's runner was re-run (167 s, 902151 checks, ALL CHECKS PASSED); identical up to the timing line.
  Scope notes:
  - the N = 18, 22, 26 censuses (bitset engines; 18 minutes at N = 26) and the Fano-colouring counts were not repeated;
  - Theorem L' was checked by the lane on random configurations only (n = 5, 7, 9);
  - the literature search for the Gale-tournament form of Theorem L was not exhaustive (Hughes; Huh-Jeon; Bogdanov-Matushkin
    cited), so priority is not claimed.
---

# THM-4529 — the Petersen family in Paley coordinates, anti-circulant tournaments, and linear K6

**PROVED / FINITE-EXACT + INDEPENDENTLY AUDITED.** Full note: [procgen_petersen_20261001_petersen_heawood_paley](../../05-knowledge/results/procgen_petersen_20261001_petersen_heawood_paley.md).

## 1. The owner's "family of 7 Petersen graphs"

The Petersen family is the Delta-Y class of K6. K6 is exactly the underlying graph of `P7 - 0`, which is THM-4524's first all-odd tournament.

In Paley coordinates the whole family is built from the translates of `{1, 2, 4}` (the Fano lines):
- **Petersen.** Delta-Y on the 4 lines avoiding the deleted point.
- **The leftover matching.** The 3 lines through it become the matching `{u, 3u}`.
- **Heawood.** Delta-Y on all 7 lines in K7.
- **K4,4 - e.** Delta-Y on the two codes of the trivial Collatz cycle, `{1, 2, 4}` and `{3, 5, 6}`.

Deleting the point kills exactly the 7 translations and keeps exactly the Frobenius `x -> 2x`, the trivial cycle's own dynamics. That is the owner's symmetry split, realised by one vertex deletion.

The two sevens do not match up, however: no member of the Petersen family has a symmetry of order 7. The honest correspondence is the split `7 = 4 + 3`.

## 2. Parity

**Which parities hold, and where.**
- Redei's "H odd for every orientation" holds only for complete graphs, so within the family only for the root K6.
- No Hamiltonian-path parity survives Delta-Y.
- The parity the family does carry is Conway–Gordon–Sachs linking parity. For straight-line K6 it is exactly the crossing parity of a ±1 walk read off the Gale tournament.

**The two H = 45 maximisers.** The maximally linked case, 3 linked pairs, is the tournament `C3[TT2, TT2, TT2]`, one of the two 6-vertex tournaments with H = 45. The other maximiser is `P7 - v`, the all-odd Paley one, and it never arises from a straight-line drawing. So Paley's parity and Petersen's parity are different parities on K6, and they pick different extremal tournaments.

**Arc-HP parity results.**
- All arcs are even in every Cayley digraph of an odd abelian group, tournament or not.
- The antipodal arcs are odd in every anti-circulant tournament.
- **QR_127 restricted to `±<2>` is all-odd on 14 vertices.** This is the first example at N = 14, settling the open existence case of HYP-9167(b). 127 is the Mersenne prime `2^7 - 1`, and `<2>` is the set of x2-codes of length-7 parity words.

## 3. Collatz

**NO PROOF.** Every parity object the bridge produces depends only on the parity word, or on residues of q. So it cannot separate integral cycles of `3x+1` from rational ones, or from those of `5x+1`.

A concrete instance: the word `100` gives an integer cycle for `q = 3` (`{1, 4, 2}`) and for `q = 5` (`{-1, -4, -2}`), and for no other q.

**FORMALIZED 2026-10-01 (Lean round 2; package [`04-computation/lean/ProcgenSelfieEdim/`](../../04-computation/lean/ProcgenSelfieEdim/README.md); orchestrator rebuilt from scratch, verify.py PASS, 617 theorems within {propext, Quot.sound}).**
- Theorem 6.1 for every odd m: `antiCirc_antipodal_odd`, and the general form `cyclic_anti_antipodal_odd`, via the formal Rédei theorem.
- Theorem 5.1 for every finite abelian group of odd order: `cayley_arcCount_even` (`FinAbGroup`), and `cayley_cyclic_prod_arcCount_even` for Z/m x Z/n.
- The N = 2 mod 4 necessity for anti-circulants: `anti_cycle_mod_four`.
