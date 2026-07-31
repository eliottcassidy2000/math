"""Extend lane D's exact greedy freeze table to the odd-ladder rates 1/3, 1/5
(corner-clocked per THM-2976) vs control 1/2. Uses lane D's Ledger verbatim."""
import importlib.util, sys
from fractions import Fraction as Fr

spec = importlib.util.spec_from_file_location(
    "laneD", __file__.replace("amm12592_oddladder_greedy_freeze_deathstar", "amm12592_carry_ledger_band_freeze_laneD_deathstar"))
laneD = importlib.util.module_from_spec(spec)
sys.modules["laneD"] = laneD
spec.loader.exec_module(laneD)

plist = [Fr(1,2), Fr(1,3), Fr(2,3), Fr(1285,2181)]
cps = [10, 20, 30, 40, 60, 80, 100, 120]

for gn, gd, D0, M in [(1,3,0,120), (1,5,0,120), (1,2,0,80)]:
    led, hist = laneD.run(gn, gd, D0, M, laneD.greedy, plist, cps,
                          f"oddladder gamma={gn}/{gd}")
    # envelope check at p=1/2: |D_M(1/2)| <= 2^-(M+1) required for feasibility
    for m in cps:
        if m > M: continue
        v = abs(hist[m][0])
        env = Fr(1, 2**m)  # (p^{M+1}+q^{M+1})/2 = 2^{-M-1}*2/2 = 2^-(M+1); use 2^-m as loose
        print(f"   gamma={gn}/{gd} m={m}: |D(1/2)|={float(v):.3e} vs 2^-(m+1)={float(Fr(1,2**(m+1))):.3e} "
              f"{'VIOLATED' if v > Fr(1,2**(m+1)) else 'ok'}")
