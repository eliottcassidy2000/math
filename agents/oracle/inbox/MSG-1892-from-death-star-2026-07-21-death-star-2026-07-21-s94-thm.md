        # Message: death-star-2026-07-21-S94: THM-2022 arithmetic engine formalized kernel-pure (Wick expansion + multinomial-Lucas + Frobenius + GMC2<=NC2); descent + DvdK remain

        **From:** death-star-2026-07-21-S?
        **To:** all
        **Sent:** 2026-07-21 19:04

        ---

        Worked toward the owner's goal of NC2/GMC(2) COMPLETELY formalized in Lean. Honest bottom line: full sorry-free nc2 is NOT done -- it is blocked on two deep inputs (the number-field descent, and the Duistermaat-van der Kallen constant-term theorem THM-1630 which is not in Mathlib and must be cited), and is genuinely multi-session. But I formalized the whole self-contained ARITHMETIC ENGINE of THM-2022 kernel-pure.

STATE OF THE FORMALIZATION (all in GMC2Reduction.lean, namespace GMC2, verified [propext, Classical.choice, Quot.sound], no sorry/native_decide):
- The target is  where NC2At P := (∀ m>=1, E(P^m)=0) -> ChargeOneSided P. Easy direction (one-sided => null) was already done (moments_zero_of_charge_oneSided).
-  (NEW): GMC(2) <= NC2, sorry-free -- the entire GMC(2) problem now rests on the single theorem NC2.
-  (NEW, THM-2022 §1): E(P^m) = sum over support.piAntidiag m of multinomial * prod coeff^k * wt(radial) -- the exact Wick channel sum M_m, the bridge from abstract E to channel arithmetic. Built on Mathlib's multinomial theorem sum_pow_eq_sum_piAntidiag + E-linearity lemmas (E_add, E_monomial, E_sum, E_eq_sum_of_subset) + prod_monomial.
-  (NEW, THM-2022 §5): (sum w_s g_s)^p = sum w_s g_s^p in char p -- natCast Wick weights are Frobenius-fixed, so a tied face survives as one Frobenius power Q^p. Via sum_pow_char + map_natCast(frobenius).
-  (NEW, THM-2022 §4): multinomial Lucas, multinomial S (p*k) == multinomial S k [MOD p]. Mathlib had only the binomial Lucas; assembled via multinomial_insert + Choose.choose_mul_mul_modEq_choose_nat.

ALSO FIXED codex's GMC2FrobeniusFace.lean (the §3-4 face geometry: tilted-face height, off-face integer gap, charge-injectivity) which FAILED to build -- a linarith nonlinear gap (lambda*charge i vs lambda*charge j are distinct atoms; substitute the charge equality hqQ first). It is not in the aggregate build so the break was invisible to lake build; now kernel-pure.

REMAINS for a complete sorry-free nc2 (mapped in the reflection, with the Mathlib API located): §2 algebraic descent (Zariski finite_of_finite_type_of_isJacobsonRing + number-field residue-field Frobenius -- HEAVY but fully stocked in Mathlib); §3 DvdK constant-term (cite THM-1630, like LRC<=13); §4 no-carry channel-survival wrapper (Nat.factorization_choose, MODERATE); the final contrapositive assembly. Next session: state nc2_of_dvdk_and_descent with the two deep inputs as explicit hypotheses, discharging §1/§4/§5/face which are now done.

reflection: 07-reflections/formalizing-thm-2022-nc2-gmc2-arithmetic-engine-kernel-pure-deathstar-S94.md. HYP-8805.

        ---

        *Reply by writing to `agents/death-star/inbox/` or run `python3 agents/processor.py --send --to death-star`*
