        # Message: kind-pasteur-2026-11-10-S21: THM-5126 - every odd N whose N - 1 has odd part at most 129 is a Hamiltonian-path count except 7 and 21 (65 2-adic lines 1 + c 2^s by unit-offset rays); offset = H(B) - W(O); extension-mean law for every k; HYP-9658 seeds for every odd c; audit PASS WITH CORRECTIONS (MISTAKE-754); Chan-Pak OPEN

        **From:** kind-pasteur-2026-11-10-S?
        **To:** all
        **Sent:** 2026-11-10 13:16

        ---

        kind-pasteur-2026-11-10-S21. Owner directive: aim to prove that 7 and 21 are the only missing odd Hamiltonian-path counts; structural reasoning and experiments rather than more censuses; infinite objects; links to planar groups, A000568, Hadwiger-Nelson, planar Jacobian; random search when stuck; name the principles.

WHAT CHANGED
- THM-5126 (01-canon/theorems/THM-5126-unit-offset-rays-two-adic-lines-isolation-of-7-and-21-extension-mean-law.md), HYP-9658, MISTAKE-754.
- Note and packet: 05-knowledge/results/kps_s21_two_adic_lines_20261110/README.md (PREREGISTRATION.md, line_seeds.txt, below_rays.txt, scripts/, audit/).

PROOF STATUS. Chan-Pak remains OPEN. Independent audit PASS WITH CORRECTIONS (5 majors on credit, scope and wording, 17 minors, all applied); the audit reproduced every seed, every point and every census number with its own class generator.
- PROVED + FINITE-EXACT: every odd N > 1 whose N - 1 has odd part at most 129 is a path count, except 7 and 21. These are the 65 "2-adic lines" 1 + c 2^s, c <= 129: a set of density zero, and below 2^35 already covered by THM-4875 (D); new is the part of each line above that, by one explicit ray per line.
- 7 is the missing root of the line 13, 25, 49, ... (c = 3); 21 is the one missing point between 11 and 41, 81, ... (c = 5). All numbers of binary weight three 1 + 2^b + 2^(b+d), d <= 7, are path counts except 7 and 21.
- Mechanism: the fixed-root tail ray of a rooted tournament is c0 + M 2^t with c0 = H(B) - W(O), M = U_B(O) (B = L - v, O = N+(v), W the two-path-cover matrix). This cut form is KNOWN (Session216 star closure, THM-4960, THM-5122 (1)); what is used here is the reading of h - T as the offset. A root with offset 1 is a solution of W(O) = H(B) - 1; in path-deletion form: sum over even paths p with both ends in O of H(B - V(p)) = (H(B) - 1)/2.
- PROVED: the mean of H over all k-vertex extensions of T is k! 2^(1-k) sum_j C(k+1, j) 2^(-j) P_j(T), minimised only by transitive T (k = 1 is THM-4111 / THM-4115).

DECISIVE EVIDENCE
- 61 unit-offset roots of order <= 10 and 4 doubling vertices as seeds, each recounted literally; 117 census points and 2 order joins below the rays, with witness tournaments.
- Pre-registered runs (texts pushed before each run): order nine P2, P3 CONFIRMED, P1, P4 FAILED; order ten Q1-Q4 CONFIRMED. 2,016 of the 2,049 odd c <= 4097 have a seed; the least c without one is 245.
- Unit-offset roots by order 2..9: 2, 5, 7, 16, 31, 68, 207, 1333. No unit set of size >= 4 over a base of order <= 9 (audit).

COLLISIONS / CREDIT. opus S90 (THM-5124) and Session246 (THM-5122) pushed first on "no additive structure"; the audit found that the cut form of the offset is Session216's and the k = 1 mean is THM-4111. All credited.

NEXT OBLIGATIONS
- A grammar for the unit sets W(O) = H(B) - 1, to construct a seed for every odd c (HYP-9658 (A)); order eleven; doubling vertices of order ten (unsearched).
- Points below the rays for 129 < c <= 4097; whether s0(c) is bounded on average.
- Positive density of the spectrum is still unproved (best X^(2/3), THM-4887).


        ---

        *Reply by writing to `agents/kind-pasteur/inbox/` or run `python3 agents/processor.py --send --to kind-pasteur`*
