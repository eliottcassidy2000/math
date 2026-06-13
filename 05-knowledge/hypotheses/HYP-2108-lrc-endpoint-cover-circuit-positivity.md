---
id: HYP-2108
status: PROVED Lemma D (exact criterion, 100% n=6..14) + P(S) functional + summed corollary; C' <=> P(S)>0 always OPEN
source: opus-2026-06-03-S575
related: [THM-398, HYP-2105, HYP-2104, HYP-2102]
---

# HYP-2108: endpoint-cover circuit positivity (Lemma D + P(S))

**LEMMA D (PROVED, verified 100% n=6..14):** component C_i=(a_i,b_i) of G(S') (midpoint m_i, length l_i), v=nw the multiple: C_i subset D_v <=> EXISTS j in Z with ||v a_i-j||<=1/n and ||v b_i-j||<=1/n (both endpoint v-phases near ONE integer) <=> ||v m_i|| <= 1/n - (v/2) l_i. S TIGHT <=> every component satisfies it.
**CIRCUIT:** the M component arc-indices j_i wind once (sum (j_{i+1}-j_i)=v); one integer v must SIMULTANEOUSLY put all M midpoints within 1/n phase of an arc centre j/v -- a closed circuit, not M separate problems.
**POSITIVITY:** P(S):=max_i(||v m_i||+(v/2)l_i-1/n); P>0 <=> loose, P<=0 <=> tight. Verified P>0 for every mult-of-n sample => C' <=> P(S)>0 on the whole multiple-of-n class (one scalar).
**PROVED COROLLARY (sum over circuit):** tight => sum_i ||v m_i|| <= M/n - (v/2)mu(G(S')) < M/n => avg midpoint v-phase < 1/n. Measured avg ~0.245 (generic 1/4) vs need 1/n (=0.071 @ n=14): positivity margin GROWS toward frontier; tightness needs v to resonate with ALL midpoints at once.
**UNIFIES:** B' = long component ((v/2)l_i>1/n => P>0); Lemma C = small owners can't share an arc centre. Residual = P>0 only via circuit simultaneous-resonance infeasibility (~1% large-owner short-component case).
**See:** THM-398 sec 4.75, `07-reflections/lrc-endpoint-cover-circuit-positivity-s575.md`, `04-computation/lrc_endpoint_cover_circuit_positivity_s575.py` (+.out).
