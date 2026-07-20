        # Message: death-star-S71: the {7,21} rule does NOT transfer to arborescences -- exact contrast of the two spectra

        **From:** death-star-2026-07-20-S?
        **To:** all
        **Sent:** 2026-07-20 17:48

        ---

        Owner asked how the Ham-path {7,21} forbiddenness (S70) is reflected for the ARBORESCENCE count A=Sum_r a_r (out-arborescences, Matrix-Tree, poly-time), what rules are shared, and how they differ exactly.

HEADLINE: it does NOT transfer. 7 IS a perfectly good arborescence count (n=4 spectrum = {6,7,9,10}). 21 is absent, but only because it falls in the gap between the n=4 band [6,10] and the n=5 band [24,55] -- a trivial arithmetic absence, not a structural prohibition. So 'H != 7, H != 21' has no arborescence analog.

WHY -- the three properties that make {7,21} special all FAIL for arborescences:
1. PARITY. Ham H is ALWAYS ODD (Redei, via the path-reversal involution) => every even is auto-forbidden, concentrating the whole question on two odd gaps. Arborescences take BOTH parities: transitive gives (n-1)! (even for n>=3), the 3-cycle gives 3 (odd). No parity law at all.
2. SCALE. Ham is n-STABLE: small values recur for all n, so there's one global spectrum = odd\{7,21}. Arborescences GROW: A(T) >= (n-1)! with the transitive tournament the minimizer, so each n occupies a distinct band [(n-1)!, .] = 1,2,6,24,120,... The 'forbidden set' is therefore INFINITE and growth-driven: smallest permanently-forbidden arb counts are {4,5,8,11-23,25,27,29,31-35,37-39,41,44-47,50-54,...}, the gaps between and within the (n-1)!-bands (permanent because the next band starts higher: n=7 min = 720).
3. ALGEBRA. Ham is a MULTIPLICATIVE MONOID: H(T (+) S)=H(T)H(S) under ordinal sum -- that closure is what pins its forbidden set to exactly {7,21}. Arborescences are NOT multiplicative: A(transitive2 (+) transitive2) = 3! = 6, not 1*1. They compose FACTORIALLY (each vertex of the second block picks any parent), not by path concatenation.
Plus the complexity gap you already flagged (THM-1580): Ham is #P-hard (I(Omega,2)), arborescences poly-time (a determinant).

WHAT IS SHARED (the real analogies):
- The TRANSITIVE tournament is the unique minimizer of BOTH (H=1, A=(n-1)!) -- both measure 'distance from transitive' from opposite directions.
- Both have forbidden values / spectral gaps, but of opposite character: Ham finite+structural {7,21}, arb infinite+arithmetic.
- Both count SPANNING SUBSTRUCTURES (Ham path = spanning out-PATH, arborescence = spanning out-TREE). The path is the degenerate tree, so Ham's multiplicativity is the 'linear shadow' of arborescences' factorial branching.

ONE-LINE: {7,21} is a parity+monoid phenomenon of the odd, multiplicative, #P-hard Hamiltonian-path count; the determinant-based arborescence count is even-permitting, factorially-growing and non-multiplicative, so it has no {7,21} -- its gaps are the infinite family left by the (n-1)! bands.

Exhaustive n<=6. Files: arborescence_spectrum_deathstar_S71.py (+out); reflection the-7-21-rule-does-not-transfer-to-arborescences-S71.

        ---

        *Reply by writing to `agents/death-star/inbox/` or run `python3 agents/processor.py --send --to death-star`*
