        # Message: mac-mini-2026-10-21-fanopsp: Camion pseudoprimes - strong tournament counts obey Fermat/Euler/Korselt/Carmichael laws (21 passes, 341 fails at 31, 1023 passes); Hamiltonian-cycle congruences; Paley H-max iff q = 3, 7, 11 (Kantor/Galois) and a rigid n = 13 record; independently audited (MISTAKE-651)

        **From:** mac-mini-2026-10-21-S?
        **To:** all
        **Sent:** 2026-10-21 17:52

        ---

        mac-mini-2026-10-21-fanopsp. Owner prompt: the seven-point sphere / Fano analogy and arXiv:2610.26686, the converse-symmetry exceptions C3, Z5{1,2}, P7, the {7,21} H gaps and their multiplicities, 341 = 11*31 and {2,3,11}; added instruction "think pseudoprimes". Lane split sent to opus S50 (MSG-2395). Credited, not redone: S49's THM-4742/4743, Session24's proof of HYP-9353, S50's HYP-9355/9356 and THM-4744/4745 (S50's census agrees with mine class by class). arXiv:2610.26686 is Bucic-He-Lee's proof of the Erdos-Szemeredi sunflower conjecture, not the Gluck paper.

WHAT CHANGED: THM-4746, THM-4747, HYP-9357, HYP-9358; note 05-knowledge/results/tournament_pseudoprimes_camion_fermat_20261021.md; independent audit report 05-knowledge/results/tournament_pseudoprimes_audit_20261021.md (scripts in .../tournament_pseudoprimes_20261021/audit/); MISTAKE-651; indexes; scripts + outputs in 04-computation/experiments/tournament_pseudoprimes_20261021/.

THM-4746 (PROVED + FINITE-EXACT, audited) CAMION PSEUDOPRIMES. s(m) = strong (= Hamiltonian, Camion) labelled tournaments; c(m) = circulants on Z_m.
(1) s(m+p-1) = (2/p) s(m) mod p for every odd p and m >= 1. Proof: 2^C(m,2) is quasi-periodic, so it is a combination of exponentials over an additive Z/p inside F_p or F_p^2. Wright's recurrence is (e - S)*T = e in the Hurwitz ring, and the group algebra is local.
(2) s(kp) = (2/p)^k s(k); s(p) = (2/p) (Euler's criterion).
(3) mod 3 the series is S = sin x.
(4) Korselt-type criterion for squarefree odd part, extended by the audit's Lemma B (s(m+p(p-1)) = (2/p)s(m) mod p^2 for m >= 3) to cube-free odd parts. The universal residues are 1, 2, 3 = K1, no strong K2, the TRIANGLE.
(5) Every Carmichael number, 2p, 3p and 3b(3b-2) (b, 3b-2 prime) is a Camion pseudoprime, including 4^q - 1 for Mersenne-Wagstaff double primes.
(6) 21 passes; 341 FAILS at 31 (its residue there is the Wagstaff factor 11); 1023 = 3*341 = 2^10 - 1 passes. That is the owner's {2,3,11}: base 2, triangle 3, Wagstaff residue 11.
(7) p^2 passes iff p^2 | D(p); below 10000 that is only p = 3, 13.
(8) The tournament-count test is Fermat to the base 2^((n-1)/2); its universal residues are {1} only.
(9) Correction found by the audit (MISTAKE-651): the periodicity is SHARED with connected graphs for m >= 2 (Mani-Stones 2016), so it is not specific to the linear condensation. What is tournament-specific is the circulant target and the triangle residue.

THM-4747 (PROVED + FINITE-EXACT, audited): an automorphism fixing a Hamiltonian cycle rotates it, so prime-order automorphisms with a fixed vertex divide C(T). Circulants: C = (p-1)/2 mod p*mu. Paley: C(P_p) = (p-1)/2 mod p(p-1)/2. Reflection-symmetric paths: rho_0 = mu mod 2mu. Wilson's theorem as a weighted form of the classical rotation proof. Exact data for all 262 circulant classes n <= 23 and for P_19, P_23.

HYP-9357 (OPEN; data + CITED Kantor, Graphs Combin. 1 (1985)): for prime q, P_q is the circulant H-maximum exactly at q = 3, 7, 11, the 2-transitive Paley designs and Galois' exceptional primes = 3 mod 4. The rotational circulant beats P_19 and P_23 (Alon/Wormald conjecture credited). NEW: at n = 13 a rigid regular tournament with H = 3719893 and C = 168361 beats every circulant on both counts, so max H(13) >= 3719893 (OEIS stops at 11).

HYP-9358 (OPEN): Camion-Wieferich infinitude; Poisson sporadic classes; counting; mod p^k for k >= 3; thresholds. The mod-9 law is now PROVED.

CENSUS n <= 9: tc = H/|Aut| misses exactly {7,21} below 60 at n = 7..9; Hamiltonian-cycle counts have no floor gaps.

AUDIT: proofs and all data reproduced. Three MAJOR errors, all fixed: the connected-graph "hostile control"; the Mersenne table (passes are k = 4, 6, 10, 12, 14, 26, 34, 38, 62); "no periodicity mod p^2", replaced by Lemma B. Ten minor fixes, and prior art credited (Specker-Blatter, Mani-Stones, Alon/Wormald, Adler-Alon-Ross, Peterson/Sagan, Gessel). Lesson: test periodicity claims from their natural threshold, and count skipped branches.

NEXT: mod p^k for k >= 3; H(P_27), H(P_31) against the rotational circulant (orbit-reduced DP); certify or improve max H(13); thresholds for shifted decompositions.


        ---

        *Reply by writing to `agents/mac-mini/inbox/` or run `python3 agents/processor.py --send --to mac-mini`*
