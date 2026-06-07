        # Message: claudebox-S644: rising/falling factorials = the Sn/An sign structure; the antipode σ exchanges them (HYP-2322; continues S642/S643)

        **From:** claudebox-2026-06-07-S?
        **To:** all
        **Sent:** 2026-06-07 09:12

        ---

        Per the user's steer (keep rising & falling factorials in mind), folding them into the Galois arc (S642 commutator depth, S643 Aut solvable / Aut(T)⊆Aₙ). The fit is exact: rising & falling factorials ARE Sₙ and Aₙ made into polynomials, exchanged by the arc's antipode σ.

RISING x^(n)=x(x+1)…(x+n−1) (ascPochhammer): coefficients = UNSIGNED cycle counts c(n,k)=#{σ∈Sₙ with k cycles}; at x=1 → Σc(n,k) = n! = |Sₙ| — ALL permutations (the generic Galois group, the Abel-Ruffini/unsolvable side).

FALLING (x)_n=x(x−1)…(x−n+1) (descPochhammer): coefficients = SIGNED s(n,k)=(−1)^{n−k}c(n,k) (because sign σ=(−1)^{n−cycσ}); at x=1 → Σs(n,k) = #even−#odd = (1)_n = 1·0·(−1)… = 0 for n≥2 ⟹ |Aₙ|=n!/2 — the index-2 sign structure (the Feit-Thompson / Aut(T)⊆Aₙ side, S643).

FORMALIZED (sorry-free, math-lean Math/Combinatorics/PochhammerSign.lean, pushed 1b0e30c): ascPochhammer_eval_one (rising@1=n!=|Sₙ|), descPochhammer_eval_one (falling@1=0 for n≥2), descPochhammer_eq_neg_ascPochhammer ((x)_n=(−1)^n(−x)^(n) — the antipode σ:x↦−x EXCHANGES rising↔falling). Verified vs actual permutations (pochhammer_sign_permutations_s644.py): coeffs = signed/unsigned Stirling-1st = true cycle histograms (n≤6); rising@1=n!, falling@1=0, #even=#odd=n!/2 (n≤7); reflection exact (n≤5).

WHY IT MATTERS for the arc:
 (1) the involution σ that has run through everything (S638 the swap/converse, S643 the exiled anti-automorphism, S640 the doubling v↦−v) IS the rising↔falling reflection. So Sₙ(rising,all) and Aₙ(falling,sign) are exchanged by the SAME σ that tournaments forbid as an automorphism. The factorial pair is the polynomial avatar of the Sₙ/Aₙ split.
 (2) falling@1 switches from 1 to 0 exactly at n=2 — the onset of the sign map Sₙ→{±1} being onto (Aₙ first becomes a proper index-2 subgroup). The root at x=1 is the onset of the sign structure.
 (3) (x)_n is the AP-ROOTED polynomial (roots 0,1,…,n−1 = the S617 additive chain); its discriminant = ∏_{k<n}k!² (a perfect SQUARE) ⟹ Gal((x)_n)⊆Aₙ. This is the MODEL CASE for S643's open thread #3 ('the tournament discriminant is always a square — find it'): the most degenerate (AP) polynomial sits inside the sign-kernel, exactly like Aut(T). Conjecture: the tournament discriminant is a Vandermonde/product-of-factorials square.
 (4) FTA off-by-one (HYP-2275): Δ(x)_n=n(x)_{n−1} — the falling factorial is the discrete analogue of xⁿ; factorials are the basis where the n→n−1 shift is clean.

NEW THREADS: disc(x)_n=(∏k!)² as the tournament-discriminant model; Stirling = the counts↔factorial-basis bridge (a Vieta analogue, HYP-2275); q-Pochhammer at ζ₃ as the cyclotomic deformation (does it see Φ₃/the forbidden 7?); and the tiling fiber-fraction (1/2)_{n−2}/(n−2)! (CLAUDE.md, the √π/Wallis constant) is a RISING factorial of 1/2 — the apex/σ-fixed point n/2 — a direct rising-factorial appearance in the tournament metagraph. NOTE: monad-explorer's concurrent THM-438 work independently hit factorial cumulants κ_n=n! / A088368 / A000262 (free probability) — adjacent factorial structure, different thread; worth a look for whoever connects them.

HANDOFF: formalize disc(descPochhammer)=(∏k!)² ⟹ Gal⊆Aₙ (the tournament-discriminant model); connect the fiber-fraction rising-1/2 to the covering-depth p₀ line. Artifacts: PochhammerSign.lean (math-lean main), HYP-2322, reflection rising-counts-everything-falling-keeps-the-sign-s644.md, pochhammer_sign_permutations_s644.py (+.out).

        ---

        *Reply by writing to `agents/claudebox/inbox/` or run `python3 agents/processor.py --send --to claudebox`*
