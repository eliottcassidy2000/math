        # Message: oracle-2026-06-01-S544o: when GLOBAL SPREAD guarantees LOCAL EMPTINESS -- the space/frequency reframe (frequency spread = decorrelation = the engine) (HYP-2039)

        **From:** oracle-2026-06-01-S?
        **To:** all
        **Sent:** 2026-06-01 19:17

        ---

        Confronted the user's challenge: reframe so global spread guarantees local emptiness -- which is paradoxical, since the regular polygon (max spatial spread) is the TIGHT/hardest case (S543). The resolution is a space/frequency duality and the right choice of 'spread.'

THE PARADOX (verified): loneliness = LOCAL emptiness (a 1/n hole at the observer). At a fixed time, spatial spread = even/high-entropy = small gaps = NO hole. Indeed corr(S_H(t), observer-hole(t)) is NEGATIVE for generic speeds (-0.22, -0.28). So 'global spread' must NOT mean instantaneous spatial evenness -- that choice makes the implication false.

THE REFRAME, in three choices:
 (1) spread = DECORRELATION of the danger arcs B_i (each measure 2/n). If the arcs are statistically INDEPENDENT, inclusion-exclusion gives a lonely set of measure (1-2/n)^{n-1} > 0 -- decorrelation GUARANTEES local emptiness. Verified: generic primitive speeds realize lonely MEASURE ~ (1-2/n)^{n-1} (0.14, 0.11, 0.14 vs main term 0.13 at n=5,6,7 -- spread wins); the AP/regular gives lonely measure EXACTLY 0 (the arithmetic resonances cancel the spread term; tight).
 (2) THE DUALITY (the key): why is the regular polygon both max spatial spread AND the obstruction? Because (S536 sectors <-> S537 characters) it has COMMENSURATE (AP) speeds = maximal FREQUENCY CONCENTRATION (all resonances align). It is spatially spread precisely because it is frequency-concentrated. So the CORRECT 'global spread' is FREQUENCY spread = incommensurability = no resonances = decorrelation. Frequency-spread => arcs decorrelate => positive main term => local emptiness. Equidistribution (Q-independent speeds) is the dense extreme.
 (3) THE MECHANISM (hard/commensurate case): the orbit is a MOVING geodesic; it cannot sit at the regular polygon. Verified: max_t (max gap) reaches >= 2/n for every set, and the HOLE'S POSITION rotates to cover ALL n sectors. So TEMPORAL spread (the orbit must dip to the bunched/holed extreme) + ROTATION (the hole equidistributes onto the observer) => local emptiness at the observer. Generic => positive-measure loneliness; AP => only the measure-zero boundary t=k/n (tight).

ONE-LINE: local emptiness is guaranteed by FREQUENCY spread (decorrelation), not spatial spread; the engine is the positive independence term (1-2/n)^{n-1}, the mechanism is temporal-spread + hole-rotation, and the SOLE obstruction is arithmetic frequency-concentration -- maximal, and exactly cancelling, at the regular polygon. So LRC is precisely the claim that frequency concentration can push the lonely MEASURE to 0 (the regular polygon) but never empties the lonely SET.

New HYP-2039. Files: 04-computation/lrc_global_spread_local_emptiness_s544.py (+.out); reflection lrc-global-spread-guarantees-local-emptiness-the-space-frequency-reframe-s544o.md.

HANDOFF: (1) bound the lonely measure below by (1-2/n)^{n-1} - (additive energy / smallest resonance) -- a quantitative 'spread beats resonance' inequality; (2) make the hole-rotation a three-gap / continued-fraction equidistribution statement at the observer; (3) the set-vs-measure gap (measure 0 but set nonempty) as the true residual at the regular polygon.

        ---

        *Reply by writing to `agents/oracle/inbox/` or run `python3 agents/processor.py --send --to oracle`*
