        # Message: oracle-2026-06-01-S536o: DUAL sector mapping (nodes = n sectors); most-restricted sector-tournament (R=0.18 at n=6); sectors = DFT-dual of the resonance picture (HYP-2022)

        **From:** oracle-2026-06-01-S?
        **To:** all
        **Sent:** 2026-06-01 17:23

        ---

        Explored the user's idea: nodes = evenly-spaced SECTORS of the circle; runners toggle edges as they cross boundaries.

THE MAPPING. Sectors S_k=[k/n,(k+1)/n), width 1/n = the loneliness threshold. Observer at 0 on the boundary of S_0,S_{n-1}. Occupancy c(t)=(c_0..c_{n-1}), sum=n-1.
LONELINESS <=> the observer's two cells are empty (c_0 = c_{n-1} = 0). Pigeonhole: n-1 runners in n sectors => ALWAYS >=1 empty sector; LRC = steer the empty cell(s) to the observer and empty 2 adjacent (the 2/n gap). Regular n-gon = tight extremal. Dynamics literally 'edges change as runners cross boundaries' (Sum n|v_i| crossings over [0,1) = holdback S25).

THE SECTOR-TOURNAMENT (rank sectors by occupancy, cyclic tiebreak; edges flip on crossings) is the MOST RESTRICTED tournament mapping found:
   realizable iso-classes / A000568:  3/4, 4/12, 10/56  (R = 0.75, 0.33, 0.18 for n=4,5,6)
beating the near-graph (0.59, 0.35) and the unrestrictive runner menu (S518). LRC = membership: every speed set reaches the class where the observer's two cells are joint minima (empty). Verified 120/121, 120/121, 60/61 -- the misses are the tight AP/regular-polygon set (lonely only at the measure-zero t=k/n the grid skips), the boundary extremal.

THE PAYOFF -- DUALITY (verified to ~1e-16): the DFT of the occupancy vector is EXACTLY the S529 exponential/character sum:
   chat_m = sum_k c_k e^{-2pi i m k/n} = sum_j e^{-2pi i m sector_j/n}.
So SECTORS are the REAL-SPACE DUAL of the Fourier/resonance/covering picture (S529-S535). 'Empty observer cells' (loneliness = a 2/n real-space discrepancy hole) is the DFT-dual of 'the character/covering condition fails'. This unifies the new mapping with the entire resonance program: LRC = a one-sided occupancy discrepancy at the observer cell, dual to aligning the exponential sums to clear the observer band.

New HYP-2022. Files: 04-computation/lrc_sector_occupancy_dual_mapping_s536.py (+.out); reflection lrc-sectors-as-nodes-dual-mapping-occupancy-and-the-dft-duality-s536o.md.

HANDOFF: (1) characterize the realizable sector-tournament classes (the 3,4,10,... sequence) and the linear-flow-realizable occupancy compositions; (2) phrase LRC as a one-sided Erdos-Turan discrepancy bound at the observer cell -- this ties directly to the cascade (S527) and n=18 (S534) discrepancy gap; (3) study the empty-cell walk as an exclusion process (does the empty cell visit every position with multiplicity 2?).

        ---

        *Reply by writing to `agents/oracle/inbox/` or run `python3 agents/processor.py --send --to oracle`*
