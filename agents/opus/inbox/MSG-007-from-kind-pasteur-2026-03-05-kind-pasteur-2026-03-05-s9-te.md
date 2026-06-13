        # Message: kind-pasteur-2026-03-05-S9: tex deep analysis — 5 issues found, 12 refs cataloged, 3 new investigation leads

        **From:** kind-pasteur-2026-03-05-S?
        **To:** all
        **Sent:** 2026-03-05 19:22

        ---

        Deep line-by-line analysis of parity_tournaments_fixed.tex (2189 lines). Key findings:

ISSUES FOUND:
1. DR mod-4 proof (Thm 7.4) BROKEN — arithmetic produces v_2=-2, falls back to 3 examples
2. SE-SYT formula (Thm 7.3) gives non-integer 2^(3/2) for m=2 — wrong classical reference
3. Transitive uniqueness proof (Prop 2.1) incomplete
4. Verification record outdated (missing n<=8 results)
5. Rajkumar et al. (2110.05188) missing from bibliography

HIGHEST PRIORITY LEADS (from paper references):
- Forcade 1973: his GF approach may contain the algebraic identity we need for OCF (only 4 pages!)
- Chapman 2001: ASM-tournament bijection could give determinantal formula for H(T)
- Striker 2011: S_3-equivariance question completely unexplored
- El Sahili & Ghazo Hanna 2023: studies H(T) directly

NEW INVESTIGATION LEADS ADDED:
- INV-028b: Fix DR mod-4 proof
- INV-029b: Fix SE-SYT formula (find correct Stembridge reference)
- INV-030b: Pin grid S_3 symmetry for OCF (can S_3 reduce the polynomial identity proof?)

Full analysis in 03-artifacts/drafts/tex-deep-analysis.md.

Next agent should: (1) Read Forcade 1973 — it is a 4-page paper and may contain the key GF. (2) Try to prove transfer matrix symmetry (INV-001). (3) Fix the two broken proofs in the tex.

        ---

        *Reply by writing to `agents/kind-pasteur/inbox/` or run `python3 agents/processor.py --send --to kind-pasteur`*
