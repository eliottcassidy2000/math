        # Message: mac-mini-2026-11-07-sevensphere: the pasted seven-point sphere holds with the anticomplementary vertex (THM-4910); signature simplex of a tournament, equicentric/skew-Hadamard, QR7 = Fano code (THM-4911); HYP-9500

        **From:** mac-mini-2026-11-07-S?
        **To:** all
        **Sent:** 2026-11-07 12:28

        ---

        mac-mini-2026-11-07-sevensphere: close-out. Owner directive (re-sent from 2026-10-21): tournament structure vs open questions; the Fano plane vs the pasted seven-point sphere; merge diverse ideas.

WHAT CHANGED
- THM-4910 (PROVED). The pasted seven-point sphere is TRUE when (C) is the anticomplementary vertex C* = V0 + 4(G - V0) = 4G - 3V0. The October lanes (session23b; S49 THM-4742) had read C = V0 and proved the claim false.
  - A, B, C*, D, E, F, H' and an eighth point X (the inverse of the facet centroid in the facet circumsphere) lie on THM-4742's S2, for every n-simplex (n >= 3, with 4 -> n+1) in the regime.
  - S2 is the sphere with power c(R^2 - |OP|^2) at the corners of the skew quadrilateral O G G' O', c = n+1, n, n-1, n. Vieta on its sides gives the pairs {D,E}, {F,B}, {C*,A}, {H',X}.
  - THM-4742's defect is the power of V0 along the median, V0A*V0C*.
  - Mod 2, C* is the PG(3,2) point {1,2,3}, the complement of V0's point. Audit J1 corrected my first reading: this congruence does not by itself select C*.
  - There is no flag extension for generic simplices (NEGATIVE, n = 4..6). There are real exceptions at n = 4 (audit J1).
- THM-4911 (PROVED). The signature simplex of a tournament, p_i = (I+M)e_i in {+-1}^n.
  - THM-1950's s(T) = 1/dist(0,hull)^2, r^2 = n - 1/s, and edges 4(1 + #separating vertices).
  - Equicentric (M^2u || u) iff every vertex is an endpoint of the same number of directed 2-paths (audit J1):
    - for odd n, iff regular (signed sub-Pfaffians of +-1 skew matrices are odd);
    - for non-strong T, iff R_a => R_b;
    - for almost-regular T, iff n = 2 mod 4;
    - lambda = n-1 mod 4.
  - Equilateral iff doubly regular (c = -1) or skew-Hadamard (c = 0). So the skew-Hadamard conjecture = orthogonal signature vectors at every n = 0 mod 4.
  - QR7's signature vectors are the 7 Fano lines. QR7, like every equicentric tournament, has no seven-point sphere.
  - The eight-point sphere transfers at vertices in the regime, with Pow(p_v) = 4(n - 1 + 2-paths at v). It does not transfer for some non-equicentric T, e.g. source => C3 => sink.
- HYP-9500.
  - (a) Strong equicentric tournaments exist at every even n >= 6: PROVED by audit J1 (substitution R[T2]; skew-Walsh towers at powers of 2).
  - (b) The lambda-values fill the progression {lambda = n-1 mod 4} between 1 or 3 and max ab: OPEN; FINITE-EXACT for n <= 10.
- AUDIT J1: PASS WITH CORRECTIONS, all applied (MISTAKE-732: generic-only negative, transfer scope, HYP-9500(a) mistyped, inverted complement).

DECISIVE EVIDENCE
- Exact symbolic proof (all n).
- Exact rational simplices n = 3..6 (all 8 labels distinct and on S2).
- Complete isomorphism-class censuses n <= 10: zero violations; equicentric counts 1, 2, 1, 5, 3, 14, 15, 165.
- QR7 = Fano code verified.
- Paley switching q = 7, 11, 19, 23, 31.

NEXT OBLIGATIONS
- HYP-9500(b): local search at n = 12, 16, 20; constructions for the middle lambda-values.
- Does the radius of S2(T, v) separate tournaments beyond the classical score data (Pow(p_v) is the Sonneborn-Berger-type count Delta_v)? Does it correlate with H(T)?
- S49's question 2 (the 2(n+1) spheres) is still open.
- One statement for both homotheties (the complementary nine-point circle and the anticomplementary sphere).

FILES
- 01-canon/theorems/THM-4910-...
- 01-canon/theorems/THM-4911-...
- 05-knowledge/hypotheses/HYP-9500-...
- 05-knowledge/results/sevensphere_anticomplementary_signature_simplex_20261107.md
- 05-knowledge/results/sevensphere_audit_J1_20261107.md
- 04-computation/experiments/sevensphere_anticomplementary_20261107/


        ---

        *Reply by writing to `agents/mac-mini/inbox/` or run `python3 agents/processor.py --send --to mac-mini`*
