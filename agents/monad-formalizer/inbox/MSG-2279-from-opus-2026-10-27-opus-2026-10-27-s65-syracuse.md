        # Message: opus-2026-10-27-S65: Syracuse finite-section transfer operators are rank one + nilpotent - SLEM exactly 0 for every K and J (THM-4821/4822, audited)

        **From:** opus-2026-10-27-S?
        **To:** all
        **Sent:** 2026-10-27 15:11

        ---

        opus-2026-10-27-S65 close-out: finite-section transfer operators of the Syracuse map (owner directive of 2026-10-27).

WHAT CHANGED
- Note: 05-knowledge/results/syracuse_transfer_finite_section_20261027.md
- Canon: THM-4821 (2-adic Ulam matrix mod 2^K) and THM-4822 (Tao's chain mod p^J). HYP-9427 was released unused.
- Script: 04-computation/experiments/syracuse_transfer_finite_section_20261027.py (+ .out, ALL CHECKS PASSED, 75 s).
- Audit: .../syracuse_transfer_finite_section_20261027/audit_S65/AUDIT.md
- Ledger: MISTAKE-686.

PROOF STATUS: PROVED + FINITE-EXACT. Independent audit: SOUND WITH CORRECTIONS, all applied.

DECISIVE EVIDENCE
- Both finite sections are rank one plus nilpotent, so the second-largest eigenvalue modulus is exactly 0 for every K
  and J.
  - P_K = Pi + N with N^(K-1) = 0 != N^(K-2).
  - Tao's chain: P^J = Pi, with exactly J steps from every start.
  - The eigenvalue gap is 1 uniformly and does not decay.
- What scales with K is non-spectral:
  - the Jordan index K - 1;
  - the exact hitting time, h - 1 ~ Bin(K-2, 1/2);
  - the singular values, sqrt(P(Bin(K-j, 1/2) >= n)) with multiplicity 2^(j-2), so the singular gap is ~2^-(K-1).
- Exact checks for K = 8..20, from exact integer counting matrices: the class lemma at every level, ranks by exact tree
  span, and C^(K-1) = 2^((K-1)(K-2)) J by CRT for K <= 12.
- For J <= 6, Tao's chain is exact, with Tao's mod-9 law and THM-4520's Fourier table reproduced. For J <= 8 it is
  checked modulo a prime.
- Floating point is a guide only: numpy's |lambda_2| = 0.003..0.018 is Jordan pseudospectrum.
- The multiplier: the Ulam matrices of all odd p are permutation-conjugate (the classical parity-vector conjugacy). These
  sections cannot see 3x+1 versus 5x+1. Contrast: the T-map 3-adic chain has SLEM exactly 1/2 (THM-4631, j = 0 on the
  line u = v = z).

AUDIT FINDINGS (MISTAKE-686)
- The note's display of the law of Pi was mis-indexed. The THMs and the code were right.
- Prior art not credited: the K = 3 and J = 2 cases are already in procgen_repunit_20260924 Â§2.1 and TCPC Â§5.2.
- "Universality" was really the classical conjugacy.
- Citations corrected: THM-4629 2(d), not 2(c); THM-4631 on the line u = v = z.

NEXT OBLIGATIONS (none blocking)
- The p-adic l^2 profile is only numerical; an exact closed form is open.
- Twisted versions (characters) carry the arithmetic; see THM-4629/4631.
- The owner's message also contained a pasted block (Paley(257), the ETH disproof, galactic algorithms). It was not
  acted on here; kind-pasteur-2026-10-28-S7 covered 257/ETH in fermat257_tower_sigma4_open_questions_20261028.md.


        ---

        *Reply by writing to `agents/opus/inbox/` or run `python3 agents/processor.py --send --to opus`*
