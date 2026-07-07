#!/usr/bin/env python3
"""
death-star-2026-07-07-S1 -- WHY does AP minimize mu_1/7 but not E[maxgap]?

Key identity: E[maxgap] = integral_0^1 mu_theta dtheta, mu_theta(E)=meas{x: maxgap>theta}.
AP minimizes mu_theta at theta=1/7 but E[maxgap] (the integral) is minimized by prim-sat.
=> AP-minimality of mu_theta must FAIL for some theta range (larger theta) whose contribution
   dominates the integral.  Find the crossover theta* where mu_theta(AP)=mu_theta(prim-sat).

If theta* > 1/7 with comfortable margin => mu_1/7 AP-min is robust in theta (reassuring).
This EXACTLY explains the E[maxgap] refutation and measures the floor's delicacy.
"""
from fractions import Fraction as F

def regions(E):
    E = list(E)
    kdenom = max(max(abs(j) for j in E), max(abs(a-b) for a in E for b in E))
    bps = set([F(0), F(1)])
    for d in range(1, kdenom+1):
        for m in range(1, d):
            bps.add(F(m, d))
    bps = sorted(bps)
    out = []
    for a, b in zip(bps, bps[1:]):
        mid = (a+b)/2
        fl = {j: (j*mid).__floor__() for j in E}
        order = sorted(E, key=lambda j: (j*mid - fl[j]))
        n = len(order); gaps = []
        for s in range(n):
            j1 = order[s]; j2 = order[(s+1) % n]
            if s < n-1:
                c = F(j2-j1); b0 = F(-(fl[j2]-fl[j1]))
            else:
                c = F(order[0]-order[-1]); b0 = F(-(fl[order[0]]-fl[order[-1]])+1)
            gaps.append((c, b0))
        subbp = set([a, b])
        for i in range(n):
            ci, bi = gaps[i]
            for jx in range(i+1, n):
                cj, bj = gaps[jx]
                if ci != cj:
                    xc = (bj-bi)/(ci-cj)
                    if a < xc < b: subbp.add(xc)
        subbp = sorted(subbp)
        for u, v in zip(subbp, subbp[1:]):
            m2 = (u+v)/2
            cbest, bbest = max(gaps, key=lambda cb: cb[0]*m2+cb[1])
            out.append((u, v, cbest, bbest))
    return out

def mu_theta(regs, theta):
    tot = F(0)
    for u, v, c, b0 in regs:
        if c == 0:
            if b0 > theta: tot += (v-u)
        elif c > 0:
            xc = (theta-b0)/c; lo = xc if xc > u else u
            if lo < v: tot += (v-lo)
        else:
            xc = (theta-b0)/c; hi = xc if xc < v else v
            if hi > u: tot += (hi-u)
    return tot

if __name__ == "__main__":
    fams = {
        "AP {1..13}": list(range(1,14)),
        "prim-sat 2*{1..12}+13": [2,4,6,8,10,12,13,14,16,18,20,22,24],
        "min-E[mg] {..29}": [1,3,5,6,7,8,9,10,11,13,15,20,29],
    }
    regs = {nm: regions(E) for nm, E in fams.items()}
    print("=== mu_theta(E) profile: AP vs prim-sat vs min-E[mg] ===\n")
    print(f"{'theta':>8} | " + " | ".join(f"{nm[:14]:>14}" for nm in fams) + " | AP is min?")
    thetas = [F(1,13), F(1,10), F(1,8), F(1,7), F(9,56), F(3,16), F(1,5), F(2,9),
              F(1,4), F(2,7), F(1,3), F(2,5), F(1,2)]
    for th in thetas:
        vals = {nm: mu_theta(regs[nm], th) for nm in fams}
        apmin = all(vals["AP {1..13}"] <= vals[nm] + F(0) for nm in fams)
        mark = "AP-min" if apmin else "NOT (AP not min)"
        print(f"{float(th):>8.4f} | " + " | ".join(f"{float(vals[nm]):>14.5f}" for nm in fams) + f" | {mark}")

    # locate crossover theta* where mu_theta(prim-sat) drops below mu_theta(AP)
    print("\n-- locating AP/prim-sat crossover (bisection on exact rationals) --")
    lo, hi = F(1,7), F(1,2)
    rAP, rPS = regs["AP {1..13}"], regs["prim-sat 2*{1..12}+13"]
    # ensure AP<=PS at lo, AP>PS at hi
    d_lo = mu_theta(rAP, lo) - mu_theta(rPS, lo)
    d_hi = mu_theta(rAP, hi) - mu_theta(rPS, hi)
    print(f"  at theta=1/7: mu_AP-mu_PS = {float(d_lo):+.5f} (AP {'<=' if d_lo<=0 else '>'} PS)")
    print(f"  at theta=1/2: mu_AP-mu_PS = {float(d_hi):+.5f} (AP {'<=' if d_hi<=0 else '>'} PS)")
    for _ in range(40):
        mid = (lo+hi)/2
        if (mu_theta(rAP, mid) - mu_theta(rPS, mid)) <= 0:
            lo = mid
        else:
            hi = mid
    print(f"  crossover theta* ~ {float(lo):.5f} = {float(lo)*7:.4f}*(1/7); "
          f"AP-min for mu_theta holds on theta in (0, ~{float(lo):.4f}), 1/7={1/7:.4f}")
    print(f"  margin: theta* - 1/7 ~ {float(lo)-1/7:+.5f}  (how far above 1/7 AP stays the minimizer)")
