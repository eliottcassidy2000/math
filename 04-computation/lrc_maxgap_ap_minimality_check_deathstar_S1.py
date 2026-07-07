#!/usr/bin/env python3
"""
death-star-2026-07-07-S1 -- Is the AP the true minimizer of E[maxgap]?

opus-S133 + kps-S58 both assert "AP minimizes E[maxgap]" and use
E[maxgap(AP_13)] = 93/440 ~ 0.2114 as the binding density floor.  BUT kps-S57's OWN
adversarial family [2,6,8,10,11,12,14,16,18,20,22,26,42] reportedly gives 0.203 < 0.2114.

opus-S133's Emaxgap_exact(E,kdenom) places breakpoints only at Farey m/d with d<=kdenom;
called with kdenom=k=13.  That is CORRECT for the AP (phase order changes at x=m/(e_i-e_j),
denominators <= k-1) but WRONG for a family with large speed gaps: [.,42] has order-changes
at denominators up to 40.  Fix: kdenom must be >= max|e_i - e_j| (= max(E) for positive sets).

Here: (1) corrected exact Emaxgap, validated vs high-res numeric; (2) settle AP-minimality.
"""
from fractions import Fraction as F

def Emaxgap_exact(E, kdenom=None):
    """EXACT E[maxgap of {frac(j x): j in E}], x~Unif[0,1]. kdenom must cover all
    order-change/wrap denominators: kdenom >= max(E) and >= max|e_i-e_j|."""
    E = list(E)
    if kdenom is None:
        kdenom = max(max(abs(j) for j in E),
                     max(abs(a-b) for a in E for b in E))
    bps = set([F(0), F(1)])
    for d in range(1, kdenom+1):
        for m in range(1, d):
            bps.add(F(m, d))
    bps = sorted(bps)
    total = F(0)
    for a, b in zip(bps, bps[1:]):
        mid = (a+b)/2
        fl = {j: (j*mid).__floor__() for j in E}
        order = sorted(E, key=lambda j: (j*mid - fl[j]))
        gaps = []
        for s in range(len(order)):
            j1 = order[s]; j2 = order[(s+1) % len(order)]
            if s < len(order)-1:
                c = F(j2-j1); b0 = F(-(fl[j2]-fl[j1]))
            else:
                c = F(order[0]-order[-1]); b0 = F(-(fl[order[0]]-fl[order[-1]])+1)
            gaps.append((c, b0))
        # maxgap(x)=max_s (c_s x+b0_s): sub-breaks at pairwise crossings inside (a,b)
        subbp = set([a, b])
        for i in range(len(gaps)):
            for jx in range(i+1, len(gaps)):
                ci, bi = gaps[i]; cj, bj = gaps[jx]
                if ci != cj:
                    xc = (bj-bi)/(ci-cj)
                    if a < xc < b:
                        subbp.add(xc)
        subbp = sorted(subbp)
        for u, v in zip(subbp, subbp[1:]):
            m2 = (u+v)/2
            cbest, bbest = max(gaps, key=lambda cb: cb[0]*m2+cb[1])
            fu = cbest*u+bbest; fv = cbest*v+bbest
            total += (fu+fv)/2*(v-u)
    return total

def Emaxgap_numeric(E, res=200_000):
    E = list(E); tot = 0.0
    for r in range(res):
        x = (r + 0.5) / res
        ph = sorted((e * x) % 1.0 for e in E)
        mg = ph[0] + 1.0 - ph[-1]
        for i in range(len(ph) - 1):
            g = ph[i+1] - ph[i]
            if g > mg:
                mg = g
        tot += mg
    return tot / res

if __name__ == "__main__":
    thr = F(1, 7)
    fams = {
        "AP {1..13}": list(range(1, 14)),
        "2*AP (M=1/14 tight)": [2*i for i in range(1, 14)],
        "prim-sat 2*{1..12}+13": [2,4,6,8,10,12,13,14,16,18,20,22,24],
        "S57-adversarial-.42": [2,6,8,10,11,12,14,16,18,20,22,26,42],
    }
    print("=== EXACT (corrected kdenom) vs NUMERIC E[maxgap];  1/7 = %.6f ===\n" % float(thr))
    print(f"{'family':>26} {'exact':>16} {'~exact':>9} {'~numeric':>9} {'vs AP 93/440':>13}")
    ap_val = F(93, 440)
    for nm, E in fams.items():
        ex = Emaxgap_exact(E)
        nu = Emaxgap_numeric(E)
        tag = "AP" if ex == ap_val else ("BELOW AP" if ex < ap_val else "above")
        print(f"{nm:>26} {str(ex):>16} {float(ex):>9.5f} {nu:>9.5f}  {tag:>13}")
    print(f"\n  AP exact = {ap_val} = {float(ap_val):.6f}")
    print(f"  Sanity: opus-S133's script called with kdenom=k=13 on the .42 family would MISS breaks.")
