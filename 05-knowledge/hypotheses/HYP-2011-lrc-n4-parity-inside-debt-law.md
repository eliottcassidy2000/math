---
id: HYP-2011
status: PARTIALLY-TRUE
source: oracle-2026-06-01-S531
related:
  - HYP-2004
  - THM-381
---

# HYP-2011: the n=4 parity law — LRC(4) closes on every odd-sum triple; n=4 is the last single-congruence case

**Mechanism.** n=4 covering character `g_k = -sin(πk/2)/(πk)` is supported ONLY on
odd `k` (`g_even = 0`). The inside-debt (order-3) term needs `k_a a+k_b b+k_c c=0`
with all `k_i` odd; mod 2 the LHS `≡ a+b+c`.

**RIGOROUS (verified 0 violations / 1336 triples, speeds≤22):**
> a+b+c ODD  ==>  no all-odd resonance  ==>  the n=4 inside-debt term is IDENTICALLY 0.

**Consequence — LRC(n=4) on odd-sum (verified, near-rigorous):** for a+b+c odd,
`|SAFE| = 1/8 + pairwise` (mean-field + 2-runner only). Bound: pairwise pair (x,y)
contributes only if x/gcd,y/gcd both odd, magnitude ≤ 1/(8xy); over odd-sum triples
`Σ 1/(xy) ≤ 0.875 < 1` so `|SAFE| > 0`. VERIFIED: all 752 odd-sum triples (≤22)
have `|SAFE| > 0`, min `1/18` at (1,3,9). => LRC(n=4) holds on ~half of all triples.

**Hard core isolated:** even-sum triples (2 odd + 1 even, primitive) carry the
active inside debt; AP `{1,2,3}` is the UNIQUE tight set (`|SAFE|=0`, speeds≤30) —
cleaner than n=5,6. Others bounded away (min positive `1/28`). Near-tight families
`(1,4k+2,4k+3)`, `(2,2j+1,2j+3)` collapse onto the AP.

**META (why n=4 is last clean case):** `g_k=0 iff (n/2)|k`; support residues mod
n/2: n=4 -> {1} (ONE class -> single parity law); n=6 -> {1,2}; n=8 -> {1,2,3}.
For n≥6 the support spans ≥2 classes so no single congruence on the speeds kills
the debt. The difficulty ramp is the width of the character support, not luck.

**Loneliness ladder:** L(n;s)=|SAFE| is an explicit arithmetic function per n;
L(3)=1/9+(2/9)χ_3(a)χ_3(b)/(ab); L(4,odd-sum)=1/8+pairwise(ψ_4)>0. LRC(n) <=>
L(n;s) >= 0 for all primitive s.

**OPEN:** (A) close even-sum n=4 (the mod-4 triple character sum; AP unique tight);
(B) multi-channel generalization of the parity law for n≥6.

**Files:** `04-computation/lrc_n4_parity_inside_debt_s531.py` (+.out). Reflection:
`07-reflections/lrc-n4-parity-law-and-the-loneliness-ladder-s531.md`. Parity law is
THM-worthy.
