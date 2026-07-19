#!/usr/bin/env python3
"""
death-star-2026-07-19-S59f -- HYP-7915: the D=10 rung at N=510511 (pico-window)
+ the D=9 gate's second degree of freedom (N=60061, 90091).

Targets (all gate-predicted OPEN by THM-1257's law {N == 1 mod L_D, N != 1 mod 2D-1}):
  F_9(60061)   = {1..60059, 60061, 540540}     -> predicted  9/540557
  F_9(90091)   = {1..90089, 90091, 810810}     -> predicted  9/810827
  F_10(510511) = {1..510509, 510511, 5105100}  -> predicted 10/5105119
                  (window (1/510512, 2/1021023), width ~1.9e-12)

Machinery: THM-1271's e-channel reduction + THM-1270's ghost enumeration,
re-gated here (142/142 vs the full evaluator + the two big members).
The tower-closure/composite-skip lemma is verified at each target N
(gcd(2D'-1, Q_{D'}) > 1 for every lower rung D').
"""
from fractions import Fraction as F
from math import gcd
import sys, time
sys.path.insert(0, '04-computation')
from lrc_D9_rung_N30031_deathstar_S59e import M_ghost

log = lambda s="": print(s, flush=True)

def run_target(N, D):
    p = 2*D - 1
    Q = (N+1)*D - 1
    x = D*(N-1)
    lo, hi = F(1, N+1), F(2, 2*N+1)
    rung = F(D, Q)
    log(f"-- F_{D}({N}) = {{1..{N}}}\\{{{N-1}}} u {{{x}}}: p={p}, Q={Q} --")
    log(f"   window ({lo}, {hi}), width {float(hi-lo):.3e}; "
        f"N mod p = {N % p}; gcd(p,Q) = {gcd(p,Q)}")
    for D2 in range(3, D):
        p2, Q2 = 2*D2-1, (N+1)*D2-1
        g2 = gcd(p2, Q2)
        assert g2 > 1, f"lower rung D={D2} unexpectedly alive"
    log(f"   tower closure verified: every lower rung D'=3..{D-1} has "
        f"gcd(2D'-1, Q_D') > 1 (dead)")
    aw = (D * pow(p, -1, Q)) % Q
    t0 = time.time()
    dmin = min(min((v*aw) % Q, Q-(v*aw) % Q)
               for v in range(1, N+1) if v != N-1)
    dx = min((x*aw) % Q, Q-(x*aw) % Q)
    log(f"   L1 witness a = {aw}: min distance = {min(dmin, dx)} (want {D}) "
        f"[{time.time()-t0:.0f}s] -> M >= {D}/{Q}")
    t0 = time.time()
    Mg, qg, ag = M_ghost(N, D)
    dt = time.time() - t0
    ok = (Mg == rung) and (lo < Mg < hi)
    log(f"   EXACT (ghost, proof-backed): M = {Mg} at q={qg} [{dt:.0f}s]")
    log(f"   rung {rung}: {'ATTAINED' if Mg == rung else 'NOT ATTAINED'}; "
        f"in-window: {lo < Mg < hi}")
    log(f"   => {'PREDICTION CONFIRMED' if ok else 'PREDICTION FAILS -- investigate'}\n")
    return ok, Mg

def main():
    log("== HYP-7915: D=10 at N=510511 + D=9 at N=60061/90091 (death-star-S59f) ==\n")
    # re-gate (discipline): 142 rows + big members
    from lrc_singlefar_absorption_atlas_deathstar_S59 import M_exact_wit
    log("GATE re-validation:")
    mism = n = 0
    t0 = time.time()
    for D in (3, 4, 5, 6):
        p = 2*D - 1
        for N in range(3*D-1, 101):
            if N % 2 == 0: continue
            Q = (N+1)*D - 1
            if gcd(p, Q) != 1: continue
            fam = [v for v in range(1, N+1) if v != N-1] + [D*(N-1)]
            if M_exact_wit(fam)[0] != M_ghost(N, D)[0]:
                mism += 1
            n += 1
    big_ok = (M_ghost(211, 6)[0] == F(6, 1271)) and (M_ghost(2311, 7)[0] == F(7, 16183))
    log(f"   {n} table rows: {mism} mismatches; big members reproduced: {big_ok} "
        f"({time.time()-t0:.0f}s)\n")
    if mism or not big_ok:
        log("GATE FAILED -- aborting"); return

    results = []
    for (N, D) in ((60061, 9), (90091, 9), (510511, 10)):
        results.append((N, D, *run_target(N, D)))

    log("== SUMMARY ==")
    for N, D, ok, Mg in results:
        log(f"   F_{D}({N}): M = {Mg} -- {'CONFIRMED' if ok else 'FAILED'}")
    if all(ok for _, _, ok, _ in results):
        log("   All three predictions CONFIRMED -- out-of-sample hits #5, #6, #7;")
        log("   the D=9 gate is confirmed OFF the primorial diagonal (60061, 90091),")
        log("   and the tower reaches primorial 510510 in a pico-width window.")

if __name__ == "__main__":
    main()
