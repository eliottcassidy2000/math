"""
Is inf_E E[gap@0] ABOVE or BELOW 1/7?  The decisive family is the inhomogeneous
AP {a + d j : j=0..12} (spread/detuned AP).  boxeph-2026-07-07-S1, HYP-4760.

klein-S153: inf E[gap@0] ~ 0.156 > 1/7  (single anchor closes it)
opus-S133 : spread cluster {6,11,..} gives E[gap@0]=0.134 < 1/7 (single anchor FAILS)
my random sample (205): inf ~ 0.1468 > 1/7 (missed the structured adversary)

Settle it EXACTLY on {a+dj}.  E[gap@0](E) = 2 E_x[min_i frac(e_i x)].
"""
from fractions import Fraction as F

def E_min_frac_exact(E):
    E = [e for e in E if e != 0]
    if not E: return F(0)
    bps = {F(0), F(1)}
    n = len(E)
    for i in range(n):
        ai = abs(E[i])
        for m in range(0, ai+1):
            bps.add(F(m, ai))
        for j in range(i+1, n):
            d = abs(E[i]-E[j])
            if d:
                for m in range(0, d+1):
                    bps.add(F(m, d))
    bps = sorted(b for b in bps if 0 <= b <= 1)
    tot = F(0)
    for a, b in zip(bps, bps[1:]):
        if a == b: continue
        mid = (a+b)/2
        best = None
        for e in E:
            fl = (e*mid).__floor__()
            val = e*mid - fl
            if best is None or val < best[0]:
                best = (val, e, fl)
        _, e, fl = best
        tot += e*(b*b-a*a)/2 - fl*(b-a)
    return tot

def E_gap0(E):
    return 2*E_min_frac_exact(E)

thr = F(1,7)
print(f"1/7 = {float(thr):.6f}")
print("="*66)
print("Inhomogeneous AP  {a + d*j : j=0..12}   (13 speeds)  -- E[gap@0] exact")
print("="*66)
best = (F(1), None)
rows = []
for d in range(1, 12):
    for a in range(0, 3*d+1):
        E = [a + d*j for j in range(13)]
        if 0 in E:  # a=0,j=0 -> 0 is a point; skip (co-offset convention differs)
            E = [e for e in E if e != 0]
            if len(E) < 13:
                # keep as 12-... no, keep only genuine 13 nonzero
                E = [a + d*j for j in range(1,14)]
        if len(set(E)) != 13: continue
        g = E_gap0(E)
        if g < best[0]:
            best = (g, (a, d, E))
        if float(g) < 0.148:
            rows.append((float(g), a, d))
rows.sort()
for g, a, d in rows[:14]:
    flag = "  <<< BELOW 1/7" if g < 1/7 else ""
    print(f"  a={a:3d} d={d:2d}  E[gap0]={g:.6f}{flag}   ({{a+dj}} = detuned AP)")
bg, (ba, bd, bE) = best
print(f"\nMIN E[gap0] over inhomog-AP grid = {str(bg)} = {float(bg):.6f}")
print(f"   at a={ba}, d={bd}, E={bE}")
print(f"   => inf over inhomog-AP is {'BELOW' if bg<thr else 'ABOVE'} 1/7   (margin {float(bg-thr):+.6f})")

print("\n" + "="*66)
print("Test opus-S133's specific claim: spread cluster near {6,11,...}")
print("="*66)
for E in [[6+5*j for j in range(13)], [6,11,16,21,26,31,36,41,46,51,56,61,66],
          [1+7*j for j in range(13)], [1+6*j for j in range(13)]]:
    print(f"  {E[:4]}...{E[-1]}  step {E[1]-E[0]}:  E[gap0]={float(E_gap0(E)):.6f}"
          f"{'  BELOW 1/7' if E_gap0(E)<thr else ''}")

print("\n" + "="*66)
print("CONCLUSION LOGIC:")
print("  If inf E[gap0] < 1/7  -> single-anchor floor FAILS; need multi-anchor /")
print("     spread-dichotomy (opus-S133).  klein-S153's single-anchor lead is WRONG.")
print("  If inf E[gap0] > 1/7  -> single anchor {0} already gives E[maxgap]>1/7,")
print("     nearly closing the crux (klein-S153 correct, opus's 0.134 an outlier).")
