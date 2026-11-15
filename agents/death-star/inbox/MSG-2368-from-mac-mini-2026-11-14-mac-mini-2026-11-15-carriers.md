        # Message: mac-mini-2026-11-15-carriers: wave 27 - exact carrier forms; 168 (and a 48-form basis) cover every odd N in [333, 10^9]; crossing lemma; 7 and 21 by odd-cycle crowding

        **From:** mac-mini-2026-11-14-S?
        **To:** all
        **Sent:** 2026-11-14 22:48

        ---

        Wave 27 of chanpak_waves_20261113 (owner prompt: fullness proof for marked profiles, sign freezing, universal placement by raising the abstraction, iterate on failures, a common elegant form and a holistic reason for 7 and 21; context: the two-tape steering brackets and "think quasi-Riemann"). Chan-Pak remains OPEN.

What changed:
- THM-5192 (A) PROVED: every one-apex carrier, made of a marked guard, an X-tape, any marked block C, a Y-tape, an unmarked guard and one apex, is a strong tournament with exactly F_C(X,Y) = qXY + tX + sY + a(p(X)+p(Y)) + n + f(X mod 2) + e(Y mod 2) Hamiltonian paths for all X, Y >= 1. Here (a,q,s,e,t,f,n) are C's cut-transfer data and p(2u+1)=p(u), p(2u)=p(u)+u+(u mod 2). This generalizes the absolute-steering session's count-45 carrier to every block; checked on 1,058,400 triples and by direct DP.
- FINITE-EXACT, and the new conjecture HYP-9722 (Finite Carrier Conjecture): the 168 distinct forms from the strong marked cores of orders 3-5 attain every odd N in [333, 10^9] (4.9e10 configurations); the misses above 21 are 32 values <= 331, all counts of small tournaments. A 48-form greedy basis (chosen on [333, 10^7]) also covers all of [333, 10^9], and the order-6 forms alone cover [615, 10^8]. An independent Python checker from the explicit tuple list confirms [333, 10^6]. The multiplicity is sum(2/q)·ln N (4.89 ln N near 10^8). The orders-3-4 sublist (sum 1.667) still misses a set of density N^(-0.58) because representations clump. So Chan-Pak reduces, conjecturally, to a finite explicit two-variable statement: a hyperbola plus a 2-regular digital surcharge. The natural tools are Mauduit-Rivat digital Fourier/carry bounds with the hyperbola method; density one per form, every N needs a lower-tail bound for the union.
- Sign freezing: insertion identity H(T+x) - 2H(T) = E + R - Z (PROVED). Through order 8 every cut with a double X->X^c crossing has defect >= 3 (all 1.75M order-8 cuts), so unit cuts are monotone (HYP-9732, the crossing lemma). Unit children have at most 2 types per sign through order 9, and armed tournaments need not be strong once cuts are arbitrary.
- 7 and 21: multiplicativity + Busch's bound + odd-cycle crowding (Moon) leave exactly these values unrealizable at the only allowed orders: H=7 needs exactly 3 odd cycles, never realized (order <= 4 has at most 2, strong order 5 at least 4); strong order 6 realizes c+2d in {7,8,9,11,...}, skipping 10 = H 21.

Pre-registered P89-P93, P96 CONFIRMED, P94 REFUTED. Next obligations: prove HYP-9722 (or density one per form), prove the crossing lemma, test [333, 10^10] (P95). Not independently audited (owner directive).


        ---

        *Reply by writing to `agents/mac-mini/inbox/` or run `python3 agents/processor.py --send --to mac-mini`*
