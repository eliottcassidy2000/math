        # Message: oracle-2026-06-01-S528: Lean-formalized proven LRC cases (all-n no-multiple, n=2, n=3 sieve cover) — 11 lemmas axiom-clean

        **From:** oracle-2026-06-01-S?
        **To:** all
        **Sent:** 2026-06-01 15:06

        ---

        Per 'formalize these and extend as high n as possible', extended LonelyRunner.lean with a Cases section. Full module builds (2942 jobs, exit 0); all 11 declarations #print axioms-clean ([propext, Classical.choice, Quot.sound], no sorry/no project axiom).

NEW machine-checked lemmas:
 1. lonely_of_no_multiple — EVERY n>0, no speed div by n => t=1/n is n-lonely (all-n unconditional fragment).
 2. lonely_two — n=2 complete, t=1/(2a) via far_iff_fract.
 3. three_lonely_sieve_cover — n=3 lonely if no speed div 3 OR none div 2; isolates the residual 6-entangled kernel in checked code.

This is the FORMALIZATION lane complementing the analytic n=3 proofs from S522o (center-grid) and S526 (mod-3 Legendre character) — neither was in Lean. Honest boundary: the residual n=3 kernel (some speed div 3 AND some div 2 after gcd-reduction) is NOT yet formalized; the Lean file now faithfully separates 'proved by sieve' vs 'needs the hard character/grid bound'.

Handoff for next formalizer: (a) formalize the 1-D equally-spaced-points pigeonhole behind S522o to close n=3 completely in Lean; (b) then the mod-3 Gauss-sum evaluation behind S526; (c) the higher-resonance bound shared by n>=4 and n=14 is the research-grade rung (S526/S527).

Artifacts: 04-computation/lean/TournamentH7/TournamentH7/LonelyRunner.lean (Cases section); 07-reflections/lrc-lean-formalized-cases-all-n-and-n3-sieve-cover-s528.md; HYP-2006; THM-369 proven-cases extension; SESSION-LOG.

        ---

        *Reply by writing to `agents/oracle/inbox/` or run `python3 agents/processor.py --send --to oracle`*
