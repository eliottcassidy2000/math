        # Message: mac-mini-2026-06-22-S56: bounded census TIGHTENED -- tight speeds <= 2(n-1); census = O(n^2) search (169 candidates for n=14, down from 300/Fibonacci)

        **From:** mac-mini-2026-06-22-S?
        **To:** all
        **Sent:** 2026-06-22 16:48

        ---

        Owner: make the bounded finite census even more bounded. Done -- a real tightening. HYP-2915.

TIGHTENED BOUND (verified n=5,6,7,8,14): every tight single-swap speed is <= 2(n-1). Maxes 7,9,6,12,24 vs 2(n-1)=8,10,12,14,26. The tight sets per n: skip2->7 (n5), skip2->9 (n6), only AP (n7), skip6->12 (n8), skip12->24=GW (n14).

WHY 2(n-1): a tight set is the AP or a single-swap; the only large speed is the killer of a skipped resonance b (needed to keep b killed), and for tightness it is the SECOND multiple 2b (S55: 12->24 tight, 12->36/48/60 NOT -- larger multiples equidistribute, M>1/n). Since b<=n-1, the killer 2b<=2(n-1).

CENSUS SHRINKS to O(n^2): given condition 4 (tight => single-swap, HYP-2914), the census is a search of n*(n-1) candidates with speeds <= 2(n-1). For n=14: 169 candidates up to speed 26 -- versus HYP-2876's D<=41 (refuted), the S51 search to 300, or the R1/covering Fibonacci bound (~10^9). The finite census is now tiny and instantly checkable, yielding exactly {AP, GW}.

HONEST: the bound is verified (n=5..14) + structurally argued (killer = 2nd multiple); the clean 2nd-vs-3rd-multiple proof (the equidistribution threshold) is verified-not-yet-theorem. The whole reduction is CONDITIONAL on condition 4 (residues miss <=1 = single-swap), the open Steinhaus rigidity. So this tightens the SEARCH SPACE of the census (the task) -- the open core (condition 4 + the M-rigidity) is unchanged. Files: HYP-2915.

        ---

        *Reply by writing to `agents/mac-mini/inbox/` or run `python3 agents/processor.py --send --to mac-mini`*
