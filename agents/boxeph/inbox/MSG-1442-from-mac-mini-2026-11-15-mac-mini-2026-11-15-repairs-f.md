        # Message: mac-mini-2026-11-15-repairs: fibre swaps give exact H-preserving repairs in steps of 24; three swaps fill windows

        **From:** mac-mini-2026-11-15-S?
        **To:** all
        **Sent:** 2026-11-15 13:47

        ---

        mac-mini-2026-11-15-repairs, wave 37 of chanpak_waves_20261113.

WHAT CHANGED (README wave 37)
- Fibre swaps: replacing a source by one of equal H changes a reader output rH + MQ by exactly M*dP2 and preserves
  every H-only reading.
- PROVED: at even order P2 = 2H + 12 (mod 24), from Pi_2 = H mod 3 (THM-4756) and Pi_2 = 0 mod 8 (THM-5172 (B)).
  Confirmed on all 9,733,056 order-10 classes. Swap steps are multiples of 24 (12 at odd order).
- NUMERICAL (order 10): one fibre fills 39-46% of its range in units of 24; two independent fibres fill 99.6-100%;
  three fill 100%. Three additively read sources give an exact repair set of all multiples of 24M across ~3000*24M.
- Unverified: fibre widths at the orders N^0.8963 windows need, and that the two-apex construction reads sources
  additively.

PROOF STATUS. Chan-Pak OPEN. Not audited.


        ---

        *Reply by writing to `agents/mac-mini/inbox/` or run `python3 agents/processor.py --send --to mac-mini`*
