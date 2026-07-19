#!/usr/bin/env python3
"""
death-star-2026-07-19-S59g -- HYP-7925 main: the D=12 rung at N=9699691.

F_12(9699691) = {1..9699689, 9699691} u {116396280}; binder p = 23
(D=11's binder 21 = 3*7 is composite -- skipped, THM-1270 lemma);
Q = 116396303; N = 9699690 + 1 == 1 mod L_12 = 9699690 = 2*3*5*7*11*13*17*19,
N == 16 != 1 mod 23 -> gate predicted OPEN: M = 12/116396303, inside
W = (1/9699692, 2/19399383) of width ~5.3e-15 (FEMTO-scale), at ~10^7 speeds.

Machinery: THM-1271 (ex-1258) e-channel reduction + THM-1270 ghost enumeration
(M_ghost, gated 142/142 in S59e and re-gated in S59f -- re-gated once more here
on a spot subset for run-to-run integrity). Estimated 30-60 min.
"""
from fractions import Fraction as F
from math import gcd
import sys, time
sys.path.insert(0, '04-computation')
from lrc_D9_rung_N30031_deathstar_S59e import M_ghost
from lrc_singlefar_absorption_atlas_deathstar_S59 import M_exact_wit

log = lambda s="": print(s, flush=True)

def main():
    log("== HYP-7925 main: the D=12 rung at N=9699691 (death-star-S59g) ==\n")
    # spot re-gate (integrity): 12 sampled table rows + one big member
    log("Spot re-gate:")
    mism = 0
    for (N, D) in ((13,3), (31,4), (37,3), (61,4), (91,4), (43,3), (63,5),
                   (77,6), (85,3), (95,4), (99,6), (211,6)):
        Q = (N+1)*D - 1
        if gcd(2*D-1, Q) != 1 or N % 2 == 0: continue
        fam = [v for v in range(1, N+1) if v != N-1] + [D*(N-1)]
        if M_exact_wit(fam)[0] != M_ghost(N, D)[0]: mism += 1
    ok2311 = M_ghost(2311, 7)[0] == F(7, 16183)
    log(f"   sampled rows mismatches: {mism}; F_7(2311) reproduced: {ok2311}\n")
    if mism or not ok2311:
        log("   RE-GATE FAILED -- aborting"); return

    N, D = 9699691, 12
    p = 23
    Q = (N+1)*D - 1
    x = D*(N-1)
    lo, hi = F(1, N+1), F(2, 2*N+1)
    rung = F(D, Q)
    log(f"F_12({N}) = {{1..{N}}}\\{{{N-1}}} u {{{x}}}: p=23, Q={Q}")
    log(f"   window ({lo}, {hi}), width {float(hi-lo):.3e}; N mod 23 = {N % 23}; "
        f"gcd(23, Q) = {gcd(23, Q)}")
    for D2 in range(3, 12):
        assert gcd(2*D2-1, (N+1)*D2-1) > 1, f"lower rung D={D2} alive?!"
    log(f"   tower closure verified: all lower rungs D'=3..11 dead (gcd > 1)")
    aw = (D * pow(p, -1, Q)) % Q
    t0 = time.time()
    dmin = min(min((v*aw) % Q, Q-(v*aw) % Q)
               for v in range(1, N+1) if v != N-1)
    dx = min((x*aw) % Q, Q-(x*aw) % Q)
    log(f"   L1 witness a = {aw}: min distance = {min(dmin, dx)} (want 12) "
        f"[{time.time()-t0:.0f}s] -> M >= 12/{Q}")
    t0 = time.time()
    Mg, qg, ag = M_ghost(N, D)
    dt = time.time() - t0
    log(f"   EXACT (ghost, proof-backed): M(F_12({N})) = {Mg} at q={qg} [{dt:.0f}s]")
    inwin = lo < Mg < hi
    log(f"   rung {rung}: {'ATTAINED' if Mg == rung else 'NOT ATTAINED'}; in-window: {inwin}")
    log(f"   VERDICT: {'D=12 gate OPEN at N=9699691 -- EIGHTH out-of-sample tower confirmation, primorial 9699690, femto-window' if (Mg == rung and inwin) else 'PREDICTION FAILS -- investigate'}")

if __name__ == "__main__":
    main()
