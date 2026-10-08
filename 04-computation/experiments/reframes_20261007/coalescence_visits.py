#!/usr/bin/env python3
"""Recurrence diagnostic for the debt walk: mean number of visits of the debt M to 1 (the identity of the multiplier
group) by time T, for chains that are still unabsorbed or not; absorbed chains stop. Prediction (Polya): rank 1 ~ sqrt T,
rank 2 ~ log T, rank >= 3 bounded. Also P(no merge by T).  Usage: NAME NSAMP TMAX SEED"""
import sys, random, math
from fractions import Fraction as Fr
from coalescence_phase import MAPS, check, modd
def run(name, nsamp, tmax, seed):
    d, br, desc = MAPS[name]; check(d, br)
    rnd = random.Random(seed)
    checkpoints = []
    T = 16
    while T <= tmax: checkpoints.append(T); T *= 4
    visits = {T: 0 for T in checkpoints}; alive = {T: 0 for T in checkpoints}
    for _ in range(nsamp):
        M, e = Fr(1), Fr(1); t = 0; nv = 0; absorbed = False; ci = 0
        while t < tmax:
            j = rnd.randrange(d)
            i = (modd(M, d) * j + modd(e, d)) % d
            mi, ri = br[i]; mj, rj = br[j]
            e = (mi * e + ri - Fr(mi, mj) * M * rj) / d
            M = M * Fr(mi, mj)
            t += 1
            if M == 1:
                nv += 1
                if e == 0: absorbed = True
            while ci < len(checkpoints) and t == checkpoints[ci]:
                visits[checkpoints[ci]] += nv; alive[checkpoints[ci]] += (not absorbed); ci += 1
            if absorbed:
                while ci < len(checkpoints):
                    visits[checkpoints[ci]] += nv; ci += 1
                break
    print(f"{name}: {desc}; {nsamp} samples")
    for T in checkpoints:
        print(f"  T={T:7d}  P(no merge) = {alive[T]/nsamp:.4f}   mean visits of debt to 1 = {visits[T]/nsamp:8.2f}   /sqrt(T) = {visits[T]/nsamp/math.sqrt(T):.3f}   /log(T) = {visits[T]/nsamp/math.log(T):.3f}", flush=True)
if __name__ == '__main__':
    run(sys.argv[1], int(sys.argv[2]), int(sys.argv[3]), int(sys.argv[4]))
