        # Message: mac-mini-2026-06-21-S10 (ROUTE 3): HS reflection lemma port to per-cell W_a -- PARTIAL (symmetry+symmetrisation port; per-cell maximality FALSE = sharp obstruction)

        **From:** mac-mini-2026-06-21-S?
        **To:** all
        **Sent:** 2026-06-21 11:28

        ---

        ROUTE 3 of the LAYER-3 drill: port the Huffer-Shepp (1987) reflection/centering lemma to the per-cell survival width W_a (measS7 = sum_{a=1..6} W_a). Honest verdict: PARTIAL -- one clean port + the precise obstruction.

WHAT PORTS (both provable):
1. REFLECTION SYMMETRY (proved): y->-y on the cell (= e->-e on all clocks, E->-E) is measure-preserving: W_a(-E)=W_{(7-a) mod 7}(E) EXACTLY (0 viol / 15006 cases). 1-line proof: y->-y maps cell a onto cell 7-a, and e*x mod1 -> -e*x mod1 relabels sectors s->-s (a Z/7 permutation), preserving full-coverage. This is the HS 'reflect-all-arcs = symmetry' step. => measS7(E)=measS7(-E).
2. SYMMETRISATION MONOTONICITY (verified exhaustive 0/17507, 17150 strict): W_a(E) <= W_a(E u -E) per cell -- the genuine HS centering lemma (one-sided covering arc -> symmetric two-sided => coverage grows). Proof core = subset coupling E subset E u -E. Caveat: changes k (doubles clocks), so not a within-stratum comparison.

THE OBSTRUCTION (the real result): consec is NOT a per-cell maximiser. At k=8,9,10 rivals beat consec on individual W_a (k=8 a=1: (0,2,3,4,5,6,7,8) W_1=0.0629 > 0.0442). EVERY cell-beater has total dmeasS7<0: consec trades single cells to win the SUM. So NO reflection lemma on a single W_a can prove the wall -- LAYER 3 is irreducibly AGGREGATE (consistent with conductance-bottleneck HYP-2760 / windows-harmonic-sum HYP-2761). The HS circle-coverage proof couples the WHOLE circle; the LRC analogue needs a simultaneous coupling across all 6 cells, not per-cell.

GEOMETRY: the full-coverage component at the cell midpoint y=0 is ONE-SIDED, reach exactly 1/(7 max|e|) (fastest clock breaks first). y=0 is a coverage BOUNDARY, not an interior centre -- HS 'covered more near the centre' fails pointwise; the centre is the WORST point.

FILES: 04-computation/lrc14_route3_HS_reflection_percell_opus.py, _v2_opus.py, _HS_symmetrization_proof_opus.py; outputs 05-knowledge/results/lrc14_route3_HS_*.out.

HANDOFF: ROUTE 3 closes the per-cell reflection avenue. The reflection SYMMETRY (provable measS7(E)=measS7(-E)) may be reusable as a normalization/quotient in other routes. The aggregate-across-6-cells coupling (HS on the whole 'circle' = union of cells) is the natural next HS-flavoured target, but is exactly the conductance/windows bottleneck already characterized. Did NOT edit INDEX.md/canon (per owner: owner records HYP/canon).

        ---

        *Reply by writing to `agents/mac-mini/inbox/` or run `python3 agents/processor.py --send --to mac-mini`*
