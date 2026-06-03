---
id: HYP-2103
status: PROVED parts canonized in THM-398 (dominance dodge, interval criterion); residual = all-short equidistribution (OPEN)
source: opus-2026-06-03-S572
related:
  - HYP-2102
  - THM-398
  - THM-369
---

# HYP-2103: dominance-dodge + the all-short residual (formalization of C')

Formalizing the S571 reduction (HYP-2102) into THM-398 surfaced four improvements:
1. **DOMINANCE not divisibility (Cor B2, PROVED):** any speed set with one runner v>(n-1)*max(others) is loose (M>1/n), no divisibility needed. The large-multiple-of-n partial is just the slice where the dominant runner is a multiple. Proof: S'=S\{v} has M(S')>=1/(n-1) by LRC(n-1); window half-width delta=1/(n(n-1)V') keeps S'>1/n; v-arc radius rho=1/(nv); delta>rho <=> v>(n-1)V'. Verified 1500/1500 each n=6..14, arbitrary residues.
2. **Interval criterion (B', PROVED, sharper):** if G(S\{v}) has a component longer than 2/(nv) then S loose -- strictly weaker than v>(n-1)V'.
3. **Residual dichotomy:** a mult-of-n config (v=nw) is loose immediately UNLESS every component of G(S\{v}) is <= 2/(n^2 w) (all-short). The all-short case is the entire open core of C': mu(G(S))=mu(G(S'))-mu(G(S') cap D_v); want mu(G(S') cap D_v)<mu(G(S')), i.e. one AP of thin arcs (period 1/(nw), mass 2/n) can't COVER a fixed open set -- a 3-distance/discrepancy statement.
4. **Residual is interval-search, not sub-clock:** for v=n the gap-midpoint (2k+1)/(2n) witnesses only ~52% (773/1499 at n=6); the safe point roams the gap interior, so no fixed finite time-set discharges it.
NEXT HANDLES: (a) Erdos-Turan on {k/(nw)} bounds |mu(G'cap D_v)-(2/n)mu(G')| by discrepancy x (#components of G', =O(sum v_i)); (b) for v=n a counting argument over the n clock gaps (other n-2 runners have total danger 2(n-2)/n<2 => a survivor).
STATUS: dominance dodge + criterion PROVED (THM-398); all-short/small-multiple residual OPEN.

**See:** `01-canon/theorems/THM-398-lrc-reduction-to-Cprime-and-dominance-dodge.md`, `07-reflections/lrc-formalizing-the-Cprime-reduction-and-dominance-dodge-s572.md`, `04-computation/lrc_dodge_formalization_s572.py` (+.out); HYP-2102 (reduction), THM-369.
