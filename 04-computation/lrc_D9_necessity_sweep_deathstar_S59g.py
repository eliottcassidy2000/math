#!/usr/bin/env python3
"""
death-star-2026-07-19-S59g -- HYP-7925 side: the D=9 NECESSITY sweep (lead xviii).

Direction under test: is N == 1 (mod L_9 = 30030) NECESSARY for the D=9 gate?
The branch mechanism says yes (each prime l in {3,5,7,11,13} kills a competitor
branch only when N == 1 mod l; parity needs N odd).  Since no N < 30031
satisfies N == 1 mod 30030, the prediction is: F_9(N) is NEVER in-window for
odd N in [26, 3000] with gcd(17, Q) = 1.  Any in-window hit REFUTES necessity
and reopens the gate law.

Ghost evaluator (gated in S59e/f); ~1450 families, each ms-scale.
"""
from fractions import Fraction as F
from math import gcd
import sys, time
sys.path.insert(0, '04-computation')
from lrc_D9_rung_N30031_deathstar_S59e import M_ghost

log = lambda s="": print(s, flush=True)

def main():
    D, p = 9, 17
    hits = []
    n = 0
    t0 = time.time()
    for N in range(27, 3001, 2):
        Q = (N+1)*D - 1
        if gcd(p, Q) != 1: continue
        n += 1
        Mg, qg, ag = M_ghost(N, D)
        if F(1, N+1) < Mg < F(2, 2*N+1):
            hits.append((N, Mg))
            log(f"   !! IN-WINDOW HIT F_9({N}): M = {Mg} -- NECESSITY REFUTED")
    log(f"D=9 necessity sweep: {n} odd N in [27,3000] with gcd(17,Q)=1 tested "
        f"({time.time()-t0:.0f}s)")
    log(f"in-window hits: {hits if hits else 'NONE -- necessity direction holds on the range'}")
    log("(no N in range is == 1 mod 30030, so zero hits = the L_9 congruence is "
        "necessary as far as tested; the branch mechanism predicts this exactly)")

if __name__ == "__main__":
    main()
