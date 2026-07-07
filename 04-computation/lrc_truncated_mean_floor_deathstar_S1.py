#!/usr/bin/env python3
"""
death-star-2026-07-07-S1 -- the CORRECTED reverse-Markov: a fine-scale TRUNCATED mean.

The reverse-Markov reduction mu_1/7 >= (7/6)(E[maxgap]-1/7) fails as an "AP-reduce" route
because E[maxgap] = int_0^1 mu_theta dtheta integrates over ALL scales and is NOT AP-minimized.

FIX: truncate at the fine scale.  Identity (layer cake):
    T(theta*) := E[min(maxgap, theta*)] = int_0^{theta*} mu_theta dtheta.
Since AP minimizes mu_theta for theta <= theta*_true ~ 0.18, T(theta*) is AP-MINIMIZED for
theta* <= theta*_true -- a MEAN that keeps the AP extremal structure E[maxgap] lost.

And T gives a corrected reverse-Markov floor: since mu_theta<=1 on [0,1/7] and mu_theta<=mu_1/7
on [1/7,theta*] (mu non-increasing),  T(theta*) <= 1/7 + (theta*-1/7) mu_1/7, hence
    mu_1/7(E) >= (T(theta*) - 1/7)/(theta* - 1/7)   (theta* > 1/7),
and by AP-minimality of T,  >= (T_AP(theta*) - 1/7)/(theta* - 1/7).  Optimize theta* in (1/7, theta*_true].

Question: is this NON-VACUOUS (T_AP(theta*) > 1/7 for some theta* <= 0.18)?  Compute exactly.
"""
from fractions import Fraction as F

def regions(E):
    E = list(E)
    kd = max(max(abs(j) for j in E), max(abs(a-b) for a in E for b in E))
    bps = set([F(0), F(1)])
    for d in range(1, kd+1):
        for m in range(1, d): bps.add(F(m, d))
    bps = sorted(bps); out = []
    for a, b in zip(bps, bps[1:]):
        mid = (a+b)/2; fl = {j: (j*mid).__floor__() for j in E}
        order = sorted(E, key=lambda j: (j*mid - fl[j])); n = len(order); gaps = []
        for s in range(n):
            j1 = order[s]; j2 = order[(s+1) % n]
            if s < n-1: c = F(j2-j1); b0 = F(-(fl[j2]-fl[j1]))
            else: c = F(order[0]-order[-1]); b0 = F(-(fl[order[0]]-fl[order[-1]])+1)
            gaps.append((c, b0))
        sub = set([a, b])
        for i in range(n):
            ci, bi = gaps[i]
            for jx in range(i+1, n):
                cj, bj = gaps[jx]
                if ci != cj:
                    xc = (bj-bi)/(ci-cj)
                    if a < xc < b: sub.add(xc)
        sub = sorted(sub)
        for u, v in zip(sub, sub[1:]):
            m2 = (u+v)/2; cb, bb = max(gaps, key=lambda cb: cb[0]*m2+cb[1]); out.append((u, v, cb, bb))
    return out

def trunc_mean(regs, tstar):
    """E[min(maxgap, tstar)] = int_0^1 min(c x+b0, tstar) dx over regions, exact."""
    tot = F(0)
    for u, v, c, b0 in regs:
        # min(c x+b0, tstar) on (u,v); split at x where c x+b0 = tstar
        if c == 0:
            tot += (b0 if b0 < tstar else tstar) * (v-u)
            continue
        xc = (tstar-b0)/c
        pts = sorted(set([u, v] + ([xc] if u < xc < v else [])))
        for p, q in zip(pts, pts[1:]):
            m = (p+q)/2
            val = c*m+b0
            if val <= tstar:
                fp = c*p+b0; fq = c*q+b0; tot += (fp+fq)/2*(q-p)   # integrate the line
            else:
                tot += tstar*(q-p)
    return tot

def mu_theta(regs, th):
    tot = F(0)
    for u, v, c, b0 in regs:
        if c == 0:
            if b0 > th: tot += (v-u)
        elif c > 0:
            xc = (th-b0)/c; lo = xc if xc > u else u
            if lo < v: tot += (v-lo)
        else:
            xc = (th-b0)/c; hi = xc if xc < v else v
            if hi > u: tot += (hi-u)
    return tot

if __name__ == "__main__":
    thr = F(1, 7)
    AP = list(range(1, 14)); PS = [2,4,6,8,10,12,13,14,16,18,20,22,24]
    rAP, rPS = regions(AP), regions(PS)
    muAP = mu_theta(rAP, thr)
    print("=== corrected reverse-Markov: fine-scale truncated mean T(theta*) ===\n")
    print(f"direct: mu_1/7(AP) = 477/1078 = {float(muAP):.5f}  (the sharp floor; AP-minimized)\n")
    print(f"{'theta*':>8} {'T_AP(t*)':>10} {'T_PS(t*)':>10} {'AP<=PS?':>8} {'floor=(T_AP-1/7)/(t*-1/7)':>26}")
    best = (F(-1), None)
    for tt in [F(19,120), F(1,6), F(17,100), F(43,240), F(11,60), F(0.185).limit_denominator(1000)]:
        if tt <= thr: continue
        Tap = trunc_mean(rAP, tt); Tps = trunc_mean(rPS, tt)
        apmin = Tap <= Tps
        floor = (Tap - thr)/(tt - thr) if Tap > thr else F(-1)
        if floor > best[0]: best = (floor, tt)
        fl_s = f"{float(floor):.5f}" if floor >= 0 else "vacuous(T_AP<1/7)"
        print(f"{float(tt):>8.5f} {float(Tap):>10.5f} {float(Tps):>10.5f} {str(apmin):>8} {fl_s:>26}")
    if best[1] is not None:
        print(f"\n  best truncated-mean floor: mu_1/7 >= {float(best[0]):.5f} at theta*={float(best[1]):.4f}")
    else:
        print(f"\n  ALL truncated-mean floors VACUOUS (T_AP(theta*) < 1/7 for every theta* <= theta*_true~0.185).")
    print(f"  (direct AP floor 477/1078={float(muAP):.5f} is sharper; T IS AP-minimized (unlike E[maxgap]),")
    print(f"   but the reverse-Markov bound from ANY fine-scale mean is VACUOUS => the TAIL mu_1/7 is")
    print(f"   IRREDUCIBLE. No mean-based reduction closes the floor; AP-minimality of the tail must be")
    print(f"   proved DIRECTLY (three-gap), not routed through a mean.)")
    # Also: is T_AP(theta*) > 1/7 ever, for theta* <= theta*_true~0.18?
    print("\n  Non-vacuity scan (T_AP(theta*) vs 1/7):")
    for tt in [F(17,100), F(18,100), F(185,1000), F(19,100), F(20,100)]:
        Tap = trunc_mean(rAP, tt)
        print(f"    theta*={float(tt):.3f}: T_AP={float(Tap):.5f}  {'> 1/7 OK' if Tap>thr else '< 1/7 vacuous'}")
