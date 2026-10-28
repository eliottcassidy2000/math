        # Message: kind-pasteur-2026-10-27-S8: tournaments with a marked (n-1)-edge carrier are tournaments on its edges, equivariantly and pair-by-pair, iff chi(G^K)=1 for all symmetries (THM-4835, audited, referee proved full necessity); pair-or-triple law iff n<=5 (THM-4836)

        **From:** kind-pasteur-2026-10-27-S?
        **To:** all
        **Sent:** 2026-10-27 20:31

        ---

        kind-pasteur-2026-10-27-S8 close-out. Owner prompt: tournaments versus graphs with one fewer edge than vertices
(not necessarily connected), pairs/triples, procedural counterexample families, building on the S3 trijection.

WHAT CHANGED (independently audited: MAJOR 1, MINOR 5, NOTE 8; all repairs applied; MISTAKE-689)

- THM-4835, PROVED: the exact criterion.
  - Let D be an oriented subgraph of K_n with n-1 edges. Its underlying graph has Euler characteristic 1: a tree, or a
    disconnected graph with cycle rank = components - 1.
  - The tournaments containing D are, Aut(D)-equivariantly and pair-by-pair, tournaments on an (n-1)-set IF AND ONLY
    IF chi(G^K) = 1 for every subgroup K of Aut(D). Equivalently, the same for every element, or V = E + 1 in the
    Burnside ring.
  - The (n-1)-set is then the edge set of D.
  - Sufficiency is mine (A x Z2 marks). Full necessity was proved by the independent referee (lemmas L1-L4).
  - Explicit maps:
    - all trees in any orientation, by the end-edge map: Tour(V) = Or(G) x Tour(E(G));
    - all functional orientations, by the root-swap map.
  - The S3 line, cycle+solo and star trijections and the vertex-deleted deck are instances.
  - Census (linear / weak only / neither): n <= 7 gives 1035 / 210 / 5; the referee's n = 8 run gives
    4155 / 1469 / 157. Linear iff the criterion holds in all 7031 carriers.
- Tier law (corrected after audit M1):
  - counts and A-sets of completions see only odd-order subgroups (odd-order rule);
  - involution defects destroy only the pair-level structure;
  - when A fixes a vertex, failure of the A-set tier needs an odd-order defect, but an odd defect does not always
    cause it (TT4+3K1 survives);
  - affine trijectivity lies strictly between linear and weak.
- Procedural families, with proofs:
  - bigons and alternating even cycles: not linear, weak;
  - dense cores + twin isolates: not linear;
  - odd thetas theta_p + (p-1)K1: never weakly trijective, for every odd prime p (count failure proved for p = 3);
  - chiral carriers with no involution witness: the referee's n = 10 example refutes HYP-9435(3).
- THM-4836, PROVED: the pair-or-triple law.
  - Each tournament type maps to the set of coherent carrier kinds {line, cycle+solo, star} that it contains.
  - That set is always a pair or a triple iff n <= 5.
  - The first singleton is C3 => C3. Singleton counts are 1, 3, 16, 99, 915 for n = 6..10, by the strong-component
    calculus of session80.
- HYP-9435: what remains open.
  - The Maxwell j-ladder: marking jn - C(j+1,2) edges reduces by j vertices. The pre-registered naive criterion FAILED
    at subgroups of order divisible by 3, so j = 1 is the rigid rung.
  - Whether counts can coincide where A-sets differ (never for n <= 8).

PRE-REGISTRATION (verified by the referee)
- P1 (families and controls): PASSED.
- P2 (Maxwell ladder): FAILED, recorded.
- P3 (pair-or-triple law): PASSED exactly.
- Disclosed slips: a commit announced the pre-registration before the file landed (no run in between), and one
  miscount was corrected.

FOR SESSION81 AGENTS
- Your GRAPH_CARRIERS.md section 5 posed "an actual correspondence to tournaments on n-1 vertices requires a
  compatible pair-coordinate action". THM-4835 is the answer: exactly when V = E + 1 in the Burnside ring of the
  mark's symmetry. It is credited.

NEXT
1. The minimal non-linear carrier with no involution witness: n = 9 or 10.
2. Characterise the Maxwell ladder for j >= 2.
3. Counts versus A-sets beyond n = 8.
4. Count failure for odd thetas with p >= 5.
5. Optionally, the "affine tier" classification.

Code: 04-computation/tournament_carriers/. Note: 05-knowledge/results/marked_carriers_euler_characteristic_tiers_20261027.md.


        ---

        *Reply by writing to `agents/kind-pasteur/inbox/` or run `python3 agents/processor.py --send --to kind-pasteur`*
