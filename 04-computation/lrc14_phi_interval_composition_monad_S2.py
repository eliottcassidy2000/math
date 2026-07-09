#!/usr/bin/env python3
r"""
lrc14_phi_interval_composition_monad_S2.py   (monad-explorer-2026-07-09-S2, HYP-5717 part 3)

THE PHI-INTERVAL COMPOSITION -- the repaired/final form of the low-anchor drift
composition, with the clean mass criterion.

For a gap (a, a+g) of the L-teeth at ruler index j (drift ratio r = S_L/Vmax), the
generalized drift-embed condition (klein-S205's proof with the midpoint replaced by a
free position phi, margin min(phi-a, a+g-phi) > 1/14 + r*phi) admits the phi-interval

    phi in ( (a + 1/14)/(1 - r),  (a + g - 1/14)/(1 + r) )   [valid iff lo < hi]

Each (j, phi) gives tau = (j + phi)/Vmax; the windows over j are DISJOINT tau-intervals.
Hence the MASS CRITERION (sufficient for a lonely instant):

    sum_j |validphi(j)| / Vmax  >  meas(G_P^c)  ( <= |P|/7 ),

and the exact checker intersects each valid window with G_P (rational intervals,
decidable) -- if ANY intersection is nonempty, Mreach(v) >= 1/14 is certified
(L-runners by the drift embed at (j, phi), P-runners by tau in G_P).

RUNS: the 2 Adversary-1 failures; the full 20-perturbation battery; a 60-set structured
sweep (dilated-AP patches, two-element perturbations); the mass-criterion ledger.
"""
import sys, random
from fractions import Fraction as F
from itertools import combinations
from math import gcd

_s = open('/home/bigo/math/04-computation/lrc14_clamp_port_composition_monad_S2.py').read()
exec(_s[:_s.rfind('if __name__')])

def is_primitive(v):
    g = 0
    for x in v:
        g = gcd(g, x)
    return g == 1

def bad_P_intervals(Pspeeds):
    """G_P^c as a merged list of rational intervals in [0,1)."""
    bads = []
    for vp in Pspeeds:
        for m in range(vp + 1):
            lo = (F(m) - F(1, 14)) / vp
            hi = (F(m) + F(1, 14)) / vp
            bads.append((max(lo, F(0)), min(hi, F(1))))
    bads.sort()
    merged = []
    for lo, hi in bads:
        if hi <= lo:
            continue
        if merged and lo <= merged[-1][1]:
            merged[-1] = (merged[-1][0], max(merged[-1][1], hi))
        else:
            merged.append((lo, hi))
    return merged

def subtract_intervals(win, bads):
    """win = (lo, hi); bads sorted merged; return list of remaining subintervals."""
    lo, hi = win
    out = []
    cur = lo
    for b0, b1 in bads:
        if b1 <= cur:
            continue
        if b0 >= hi:
            break
        if b0 > cur:
            out.append((cur, min(b0, hi)))
        cur = max(cur, b1)
        if cur >= hi:
            break
    if cur < hi:
        out.append((cur, hi))
    return [iv for iv in out if iv[1] > iv[0]]

def phi_composition(v, S_L):
    """Returns (fires, total_valid_mass, total_surviving_mass, witness_tau or None)."""
    v = sorted(v)
    Vmax = v[-1]
    E_L = [Vmax - vi for vi in v if Vmax - vi <= S_L]
    P = [vi for vi in v if Vmax - vi > S_L]
    if not E_L:
        return False, F(0), F(0), None
    r = F(S_L, Vmax)
    if r >= 1:
        return False, F(0), F(0), None
    bads = bad_P_intervals(P)
    tot_valid = F(0)
    tot_surv = F(0)
    witness = None
    for j in range(Vmax):
        teeth = sorted(set(F((e * j) % Vmax, Vmax) for e in E_L))
        nt = len(teeth)
        for t in range(nt):
            a = teeth[t]
            hi_t = teeth[t + 1] if t + 1 < nt else teeth[0] + 1
            g = hi_t - a
            if a + g > 1:
                continue
            plo = (a + F(1, 14)) / (1 - r)
            phi_hi = (a + g - F(1, 14)) / (1 + r)
            if phi_hi <= plo:
                continue
            # clip phi to [0, 1]
            plo2, phi2 = max(plo, F(0)), min(phi_hi, F(1))
            if phi2 <= plo2:
                continue
            tot_valid += (phi2 - plo2) / Vmax
            # tau window
            wlo, whi = (j + plo2) / Vmax, (j + phi2) / Vmax
            for slo, shi in subtract_intervals((wlo, whi), bads):
                tot_surv += shi - slo
                if witness is None:
                    witness = (slo + shi) / 2
    return witness is not None, tot_valid, tot_surv, witness

def check(v, name=""):
    v = sorted(v)
    Vmax = v[-1]
    pc = is_primitive(v) and is_covering(v)
    best = None
    for S_L in sorted(set([Vmax // 4, Vmax // 3, Vmax // 2, 2 * Vmax // 3,
                           int(0.8 * Vmax), Vmax - 1])):
        fires, tv, ts, wit = phi_composition(v, S_L)
        P = [vi for vi in v if Vmax - vi > S_L]
        crit = tv - F(len(P), 7)  # mass criterion slack (sufficient if > 0)
        if best is None or (fires and (best[0] is False or ts > best[2])):
            best = (fires, tv, ts, wit, S_L, crit)
        if fires and best[0]:
            pass
    fires, tv, ts, wit, S_L, crit = best
    status = "FIRES" if fires else "FAIL"
    print(f"  [{name}] prim-cov={pc} {status}: split S_L={S_L}, valid mass={float(tv):.4f}, "
          f"surviving={float(ts):.5f}, mass-crit slack={float(crit):+.3f}"
          + (f", tau={float(wit):.6f}" if wit else ""))
    return fires

if __name__ == "__main__":
    print("=" * 100)
    print("PART 1 -- the two Adversary-1 failures, repaired checker")
    print("=" * 100)
    f1 = [14, 28, 42, 56, 70, 83, 98, 112, 126, 140, 154, 168, 182]
    f2 = [14, 28, 42, 56, 70, 85, 98, 112, 126, 140, 154, 168, 182]
    r1 = check(f1, "fail-1 (83)")
    r2 = check(f2, "fail-2 (85)")
    sys.stdout.flush()

    print()
    print("=" * 100)
    print("PART 2 -- full 20-perturbation battery + two-element perturbations (phi-interval)")
    print("=" * 100)
    base = [14 * i for i in range(1, 14)]
    total = fails = 0
    failed_sets = []
    for idx in range(13):
        for d in (-2, -1, 1, 2):
            v = sorted(set(base[:idx] + [base[idx] + d] + base[idx + 1:]))
            if len(v) != 13 or not is_primitive(v) or not is_covering(v):
                continue
            total += 1
            Vmax = v[-1]
            got = False
            for S_L in sorted(set([Vmax // 3, Vmax // 2, 2 * Vmax // 3, int(0.8 * Vmax), Vmax - 1])):
                fires, _, _, _ = phi_composition(v, S_L)
                if fires:
                    got = True
                    break
            if not got:
                fails += 1
                failed_sets.append(v)
    print(f"  one-element perturbations: {total - fails}/{total} fire")
    rng = random.Random(20260709)
    total2 = fails2 = 0
    for trial in range(40):
        v = list(base)
        for _ in range(2):
            i = rng.randrange(13)
            v[i] = max(2, v[i] + rng.choice([-3, -2, -1, 1, 2, 3]))
        v = sorted(set(v))
        if len(v) != 13 or not is_primitive(v) or not is_covering(v):
            continue
        total2 += 1
        Vmax = v[-1]
        got = False
        for S_L in sorted(set([Vmax // 3, Vmax // 2, 2 * Vmax // 3, int(0.8 * Vmax), Vmax - 1])):
            fires, _, _, _ = phi_composition(v, S_L)
            if fires:
                got = True
                break
        if not got:
            fails2 += 1
            failed_sets.append(v)
    print(f"  two-element perturbations: {total2 - fails2}/{total2} fire")
    for v in failed_sets[:6]:
        print(f"    STILL FAILS: {v}")
    sys.stdout.flush()
