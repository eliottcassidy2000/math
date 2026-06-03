---
id: HYP-2112
status: PROVED Lemma G (G(v)=Phi(C), exact, verified 900/900 n=6..14); C' <=> ker Phi has no multiple-of-n config OPEN
source: opus-2026-06-03-S576
related: [THM-398, HYP-2108, HYP-2105]
---

# HYP-2112: the circuit-to-gap functional G(v)=Phi(C)

**FUNCTIONAL:** in v-phase a component C_i=(a_i,b_i) of G(S') is the interval (v a_i, v b_i); v-danger is the band B=U_k(k-1/n,k+1/n); uncovered length = v l_i - |(v a_i,v b_i) cap B|.
   Phi(C) := (1/v) sum_i [ v l_i - sum_k |(v a_i,v b_i) cap (k-1/n,k+1/n)| ].
**LEMMA G (PROVED):** G(v):=mu(safe set of S=S'u{v}) = Phi(C), exactly. Verified Phi==mu(safe) 900/900 each n=6..14, zero error.
**FROM SIGN TO VALUE:** S575's P(S)=max_i(||v m_i||+(v/2)l_i-1/n) read only the SIGN (whether some component pokes out of the band); Phi SUMS the pokes -- each phi_i>=0 is a ReLU/ramp of the circuit data, phi_i>0 <=> the i-th P-term>0, and Phi=sum phi_i = exact uncovered measure. Phi>0 <=> loose, and Phi gives the exact gap.
**KERNEL:** ker Phi = {Phi=0} = {every phi_i=0} = {every phase-interval subset B} = TIGHT/worry-set. So C' <=> ker Phi contains no multiple-of-n config <=> Phi(C)>0 for every n|v. Phi is a sum of ReLU activations of {v a_i, v b_i} => "min Phi over n|v configs > 0" is an LP/optimisation problem.
**vs S581:** S581's peeling Lemmas E/F (+B'/C) prove phi_i>0 for one component by owner arguments, covering 100% of sampled mult-of-14. Phi is complementary: it computes exact sum_i phi_i (by how much), they show which/why.
**See:** THM-398 sec 4.95 (Lemma G), `07-reflections/lrc-circuit-to-gap-functional-G-equals-Phi-s576.md`, `04-computation/lrc_circuit_to_gap_functional_s576.py` (+.out).
