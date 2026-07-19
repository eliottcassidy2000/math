#!/usr/bin/env python3
"""
death-star-2026-07-19-S59e -- HYP-7910: the D=9 rung at N=30031, via the
GHOST-ENUMERATION evaluator (THM-1258's e-channel reduction with an O(De)
clearance computation), plus the COMPOSITE-SKIP LEMMA verification.

TARGET: F_9(30031) = {1..30029, 30031} u {270270}; p = 17, Q = 270287;
gate (30031 == 1 mod 30030 = L_9 = 2*3*5*7*11*13; 30031 == 9 != 1 mod 17)
predicts OPEN -> M = 9/270287, inside a window of width ~5.5e-10.

GHOST ENUMERATION (replaces the O(N) clearance sweep): in the e-channel at
modulus S with multiplier a (gcd(a,S) = g'), the base elements at distance
exactly r from 0 are the representatives <= N of u == +-r*a^{-1} (mod S/g')
(none unless g' | r; if S/g' <= N the class contains S/g' -> distance 0 ->
dead).  Since N < S/g', each class has at most ONE representative <= N.  So
  c_eff := min( min{r <= De : some representative in base}, De )
computes the candidate's value min(c, De)/S in O(De) modular operations.
Exactness: value = min(c, De)/S needs only c ^ De, and any r < value-relevant
range is enumerated; u = N-1 (not in base) is skipped where it appears.

GATES (mandatory before trusting N=30031):
  G1  the ghost evaluator reproduces the FULL evaluator on all 142 odd-N
      table rows (D = 3..6, N in [3D-1, 100], gcd(2D-1, Q) = 1).
  G2  it reproduces the S59d values at (211, 6) and (2311, 7).

COMPOSITE-SKIP LEMMA (proved; verified numerically here): at tower-N
(N == 1 mod every odd prime < the level binder), any composite lower binder
p' = 2D'-1 has a prime factor l <= p'/3 in the primorial; N == 1 (mod l)
gives Q_{D'} = (N+1)D'-1 == 2D'-1 = p' == 0 (mod l), while l | D' would force
l | 1.  So the rung pair congruence p'*a == +-D' (mod Q_{D'}) is unsolvable:
composite-binder rungs are DEAD at every tower-N.  The same computation kills
every lower PRIME binder at tower-N (N == 1 mod p' => p' | Q_{D'}).
"""
from fractions import Fraction as F
from math import gcd
import sys, time
sys.path.insert(0, '04-computation')

log = lambda s="": print(s, flush=True)

def M_ghost(N, D):
    """exact M(F_D(N)) via THM-1258 reduction + ghost enumeration.
    Requires N odd, gcd(2D-1, Q)=1, N > 3D-2."""
    p = 2*D - 1
    Q = (N+1)*D - 1
    assert N % 2 == 1 and gcd(p, Q) == 1 and N > 3*D - 2
    x = D*(N-1)
    aw = (D * pow(p, -1, Q)) % Q          # Lemma-1 witness
    bn, bd, ba = D, Q, aw                 # rung floor
    Nm1 = N - 1
    for u in range(1, N+1):               # sums S = x+u, u in base; 2x is sealed
        if u == Nm1: continue
        S = x + u
        g = gcd(Nm1, S)
        E = S // N                        # e <= S/N (deleted-prefix packing)
        if g > E: continue
        M1 = S // g
        invb = pow(Nm1 // g, -1, M1)
        for e in range(g, E+1, g):
            De = D*e
            # prune: value <= De/S; skip if cannot beat current best
            if De * bd <= bn * S: continue
            a0 = (e // g) * invb % M1
            for s_a in (a0, (-a0) % M1):
                for k in range(g):
                    a = s_a + k*M1
                    if a == 0: continue
                    gp = gcd(a, S)
                    Mg = S // gp
                    if Mg <= N:           # class of 0 has representative <= N
                        continue          # distance 0 somewhere -> dead
                    ainv = pow(a // gp, -1, Mg) if gp > 1 else pow(a, -1, S)
                    # ghost enumeration: c = least r <= De with a base rep
                    c_eff = De
                    r = gp
                    while r < De:         # r must be multiple of gp
                        rr = (r // gp) * ainv % Mg
                        for cand in (rr, (Mg - rr) % Mg):
                            if 1 <= cand <= N and cand != Nm1:
                                c_eff = r
                                break
                        if c_eff == r: break
                        r += gp
                    if c_eff * bd > bn * S:
                        bn, bd, ba = c_eff, S, a
    return F(bn, bd), bd, ba

def main():
    log("== HYP-7910: the D=9 rung at N=30031 (death-star-S59e) ==\n")

    # --- G1: ghost evaluator vs full evaluator on the 142-row table ---
    from lrc_singlefar_absorption_atlas_deathstar_S59 import M_exact_wit
    log("G1: ghost evaluator vs FULL evaluator, all odd-N table rows:")
    mism = n = 0
    t0 = time.time()
    for D in (3, 4, 5, 6):
        p = 2*D - 1
        for N in range(3*D-1, 101):
            if N % 2 == 0: continue
            Q = (N+1)*D - 1
            if gcd(p, Q) != 1: continue
            x = D*(N-1)
            fam = [v for v in range(1, N+1) if v != N-1] + [x]
            Mfull, qf, af = M_exact_wit(fam)
            Mg, qg, ag = M_ghost(N, D)
            n += 1
            if Mfull != Mg:
                mism += 1
                log(f"   !! MISMATCH D={D} N={N}: full={Mfull}@{qf} ghost={Mg}@{qg}")
    log(f"   {n} rows, {mism} mismatches ({time.time()-t0:.0f}s)")
    if mism: log("   G1 FAILED -- abort"); return

    # --- G2: vs S59d pruned at (211,6) and (2311,7) ---
    log("G2: ghost evaluator at the two big members:")
    for (N, D, want) in ((211, 6, F(6,1271)), (2311, 7, F(7,16183))):
        t0 = time.time()
        Mg, qg, ag = M_ghost(N, D)
        tag = "OK" if Mg == want else "?!"
        log(f"   {tag} F_{D}({N}): M = {Mg} [{time.time()-t0:.1f}s] (want {want})")
    log("")

    # --- composite-skip + lower-rung-death verification at N=30031 ---
    N = 30031
    log("Tower closures at N=30031 (all lower rungs dead -- two-line lemma):")
    for D2 in (3, 4, 5, 6, 7, 8):
        p2 = 2*D2 - 1
        Q2 = (N+1)*D2 - 1
        g2 = gcd(p2, Q2)
        kind = "prime" if p2 in (5, 7, 11, 13) else "composite"
        log(f"   D={D2} (binder {p2}, {kind}): gcd(p, Q_D) = {g2} "
            f"{'-> rung pair p*a==+-D unsolvable, DEAD' if g2 > 1 else '-> ALIVE?!'}")
    log("")

    # --- THE RUN ---
    D = 9
    p, Q, x = 17, (N+1)*9 - 1, 9*(N-1)
    lo, hi = F(1, N+1), F(2, 2*N+1)
    log(f"== D=9 at N=30031: F_9 = {{1..30029, 30031, {x}}} ==")
    log(f"   p=17, Q={Q}, gcd(17,Q)={gcd(17,Q)}; window ({lo}, {hi}), "
        f"width {float(hi-lo):.3e}")
    aw = (9 * pow(17, -1, Q)) % Q
    fam_check = min(min((v*aw) % Q, Q-(v*aw) % Q)
                    for v in list(range(1, 30030)) + [30031, x])
    log(f"   L1 witness a = {aw}: min distance over all 30031 speeds = {fam_check} "
        f"(want 9) -> M >= 9/{Q}")
    t0 = time.time()
    Mg, qg, ag = M_ghost(N, D)
    dt = time.time() - t0
    rung = F(9, Q)
    log(f"   EXACT (ghost, proof-backed): M(F_9(30031)) = {Mg} at q={qg} [{dt:.0f}s]")
    log(f"   rung 9/{Q}: {'ATTAINED' if Mg == rung else 'NOT attained'}; "
        f"in-window: {lo < Mg < hi}")
    log(f"   VERDICT: {'D=9 gate OPEN at N=30031 -- FOURTH out-of-sample tower confirmation, primorial 30030' if (Mg == rung and lo < Mg < hi) else 'prediction FAILS -- see value'}")

if __name__ == "__main__":
    main()
