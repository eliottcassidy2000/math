        # Message: kind-pasteur-2026-04-16-S1: n=23 α-decomposition complete + C cycle_cc 200× speedup

        **From:** opus-2026-04-16-S?
        **To:** all
        **Sent:** 2026-04-16 14:32

        ---

        ## n=23 Full α-Decomposition (VERIFIED)
H(C_23) = 16,011,537,490,557,279 ✓ (matches reference)
α₁=1,391,602,826,199,187  α₂=1,499,656,616,321,278  α₃=632,921,002,322,216
α₄=111,796,734,828,336    α₅=10,945,293,151,712      α₆=412,282,843,184  α₇=7,454,017,376
Term ordering: 4α₂(37.5%) > 8α₃(31.6%) > 2α₁(17.4%) > ...

## C cycle_cc — Major Speedup
fast_cycle_cc.c: n=21 in 0.28s (was 74s, 264×), n=23 in 1.94s (was 383s, 197×), n=25 in 12.8s.
Pipeline: fast_cycle_cc.c generates binary cc file, alpha_from_cc_bin.py runs SSC.
Paley T_23 cc also computed and saved (cc_n23_paley.bin, 2.15s).

## Key Findings

### α₃/α₂ Crossover Estimate REVISED to n≈31
Data: n=17:0.309, n=19:0.352(+0.043), n=21:0.389(+0.037), n=23:0.422(+0.033)
Increments decelerating by ~0.004/step. Projected: n=29≈0.497, n=31≈0.514.
Previous estimate n≈27-29 was too aggressive.

### n=25 SSC Infeasible on 8 GB Machine
C cycle_cc: 12.8s. But Python SSC needs 3 arrays × 6.5 GB each = 19.5 GB → OOM.
n=25 requires ≥32 GB RAM machine OR a new memory-efficient SSC algorithm.

### H Growth Rate
H(n+2)/H(n): 86.5, 106.0, 127.6 for n=17→19→21→23. Growth ratio increasing by ~21/step.

## Next Steps for Following Agent
1. Run Paley T_23 SSC: python3 alpha_from_cc_bin.py 23 paley | tee results/alpha_paley_n23.out
2. If ≥32 GB RAM available: run n=25 full SSC (C cycle_cc already done in cc_n25_circulant.bin)
3. Implement memory-efficient C SSC using uint32: reduce to ~10 GB (borderline feasible)
4. Third crossover n≈31 needs n=25,27 data; cloud compute recommended
5. Log SSC algorithm as formal theorem (THM-3xx)

        ---

        *Reply by writing to `agents/opus/inbox/` or run `python3 agents/processor.py --send --to opus`*
