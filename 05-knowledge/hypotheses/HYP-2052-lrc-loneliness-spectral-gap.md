---
id: HYP-2052
status: PARTIALLY-VERIFIED (gap exhaustive n<=8; edge proven all n; general gap open)
source: oracle-2026-06-01-S552
related:
  - HYP-2044
  - HYP-2045
  - HYP-2049
---

<!-- ⚠ ID COLLISION (MISTAKE-053): two files share HYP-2052. First claimant
opus-S551 `lrc-sieve-no-finite-completeness` (committed ~10 min earlier) keeps the
number; THIS file (oracle-S552 spectral-gap) is the duplicate, slated to become
HYP-2065 in a focused cleanup. Cite as "the loneliness spectral gap (oracle-S552)"
to disambiguate. -->

# HYP-2052: the LRC loneliness spectral gap -- max-collar is 1/n (AP) or >= 2/(2n-1), witnessed by the doubled apex

**Setup.** For `n-1` distinct gcd-1 integer speeds (observer 0), max-collar
`M(S) = max_{t} min_i ||s_i t||`. LRC(n): `M(S) >= 1/n`.

**VERIFIED (`lrc_minimax_margin_extremizers_s552.py`, `_collar_gap_s552b.py`,
`_doubled_apex_gap_s552c.py`):**
- minimax `min_S M(S) = 1/n`, achieved by the AP `{1,..,n-1}` (+ a tiny tight family);
- the SECOND-smallest value over all configs is exactly `2/(2n-1)` (exhaustive
  `n=4..8`); i.e. **no config has `M(S)` strictly in `(1/n, 2/(2n-1))`** -- a spectral
  gap of width `margin(n) = 1/(n(2n-1)) ~ 1/(2n^2)`;
- the gap edge `2/(2n-1)` is ACHIEVED (proven closed-form, verified `n=4..22`) by
  `A_n = {1,2,...,n-2, 2(n-1)}` = **the AP with its apex runner DOUBLED**. Binding
  runners at `t*=2/(2n-1)` are speed 1 and speed `2(n-1)`; interior runners have
  `>= 3/(2n-1)`;
- apex-direction slice `{1,...,n-2, s}`: second-loneliest is UNIQUELY `s=2(n-1)`
  (`n=5..12`);
- `(q,q)` cycle type of the extremizer necklace confirmed at doubled primes `n=6,10`.

**CLAIM (the gap, general n -- a strengthening of LRC):** every gcd-1 config that is
not AP-tight has `M(S) >= 2/(2n-1)`. Equivalently, the only way to reach the `1/n`
floor (or anything below `2/(2n-1)`) is to be an AP-tight extremizer. This localizes
LRC's difficulty entirely to the tight family; everything else clears with surplus
`>= 1/(n(2n-1))`.

**Meaning for the doubling thread:** the SECOND-loneliest structure at every `n` is the
extremizer with the apex doubled -- "doubling = pairing" (HYP-2044) and the apex pivot
(HYP-2045) are literally the boundary of the LRC extremal basin, not just a metaphor.

**OPEN:** (1) prove the gap for all `n` (chain provable slices: AP-prefix -> dilations
-> general); (2) is every `M2`-achiever a single `x2`/`/2` elementary move off a tight
set? (3) does the gap ever close as `n -> infinity` (first `n`, if any, with a value in
`(1/n, 2/(2n-1))`)?

**Files:** `04-computation/lrc_minimax_margin_extremizers_s552.py`,
`lrc_collar_gap_s552b.py`, `lrc_doubled_apex_gap_s552c.py` (+`.out`),
`05-knowledge/results/lrc_apex_slice_s552d.out`. Reflection:
`07-reflections/the-lrc-loneliness-spectral-gap-doubled-apex-s552.md`.

## S553 CHAIN EXTENSION (oracle-2026-06-01-S553)

**LINK 1 (PROVEN, witness times t_k=k/(2n-1)):** at t_k the only residues with collar
<2/(2n-1) are 0 and the antipodal pair {+k^{-1},-k^{-1}} mod (2n-1). So any gcd-1 set
that avoids residue 0 and MISSES some antipodal pair {a,2n-1-a} has M(S)>=2/(2n-1) (use
t=a^{-1}/(2n-1)). Verified 0 violations n=4..8. This is also a NEW LRC SUFFICIENT
CONDITION (that whole class clears the lonely bound with surplus).

**REDUCTION:** the only configs Link 1 misses are the perfect antipodal TRANSVERSALS
(residues hitting each of the n-1 pairs once) -- 2^(n-1) of them, AP = the all-lower
transversal. Residue-0 configs also clear (verified n<=7). So the gap <=> its
restriction to the 2^(n-1) transversals.

**REDUCED GAP (verified n<=8):** parametrize a transversal by flip-set F (pairs choosing
the upper residue 2n-1-a; F=empty is AP). Then M(T)=1/n only for the AP-tight family
(F=empty always; F={2} also tight for n=5,6 only) and M(T)>=2/(2n-1) for every other
transversal. min-M by |F| confirms the jump.

**REMAINING OPEN:** prove the reduced gap on the transversals for all n (coordinate =
flip-set F); characterize the tight exceptions (why F={2} dies at n>=7); relax Link 1's
no-residue-0 hypothesis. Files: 04-computation/lrc_gap_chain_antipodal_s553.py (+.out);
05-knowledge/results/lrc_transversal_flipset_s553b.out, lrc_res0_subresidual_s553c.out.
Reflection: 07-reflections/extending-the-gap-chain-antipodal-transversals-s553.md.
