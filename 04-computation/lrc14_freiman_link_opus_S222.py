"""
opus-2026-07-11-S222: the Freiman-stability link for the k=9,10 coverage extremality.

RETARGETED (klein S251: consec-Phi-extremality FALSE k>=11; TRUE k<=10). mac-mini THM-705: k=9,10 reduce to
ONE linear inequality each from the optimal deg-2 majorant q* = 1 - N/2 + N(N-1)/12:
   majorant = 1 - (1/2)m1 + (1/12)m2 <= cap_{k+1}   <=>   L(E) := 6 m1 - m2 >= 12(1 - cap_{k+1}),
and (identity) L = 6E[N] - E[N(N-1)] = 49/4 - [(E[N]-7/2)^2 + Var(N)], so the inequality is an UPPER bound
(E[N]-7/2)^2 + Var(N) <= 12 cap + 1/4.  consec should MINIMIZE L (k<=10).

THE FREIMAN QUESTION (the crux of the route): is L(E) controlled by the ADDITIVE STRUCTURE of the offsets,
so that "far from AP" (low energy / high doubling) => L large (bounded away from consec)?  Test the
correlation L vs additive energy / doubling / longest-AP, and whether consec is the argmin with a gap.
"""
from fractions import Fraction as F
from itertools import combinations

def sector(e,x): return int(((e*x)%1)*7)
def brk(E):
    Es=[abs(e) for e in E if e!=0]; bps={F(0),F(1)}
    for e in Es:
        for m in range(0,7*e+1): bps.add(F(m,7*e))
    return sorted(b for b in bps if 0<=b<1)
def m1m2(E):
    pts=brk(E)+[F(1)]; p=[F(0)]*7
    for i in range(len(pts)-1):
        a,b=pts[i],pts[i+1]
        if b<=a: continue
        occ=set(sector(e,(a+b)/2) for e in E); p[len(set(range(1,7))-occ)]+=(b-a)
    m1=sum(F(n)*p[n] for n in range(7)); m2=sum(F(n*(n-1))*p[n] for n in range(7))
    return m1,m2
def L(E):
    m1,m2=m1m2(E); return 6*m1-m2
def add_energy(E):
    # E(E) = #{(a,b,c,d) in E^4 : a+b=c+d} = sum_t r(t)^2, r(t)=#pairs summing to t
    from collections import Counter
    c=Counter(a+b for a in E for b in E)
    return sum(v*v for v in c.values())
def diff_doubling(E):
    return len(set(a-b for a in E for b in E))   # |E-E|
def longest_AP(E):
    s=set(E); best=1
    for a in E:
        for d in range(1,max(E)-min(E)+1):
            L_=1; x=a+d
            while x in s: L_+=1; x+=d
            if L_>best: best=L_
    return best

cap={9:F(1979,4004),10:F(55,91)}
for k,D in [(9,12),(10,12)]:
    print(f"=== k={k}, diam<=D={D}: L=6m1-m2 vs additive structure (consec should be argmin) ===")
    base=list(range(k)); Lb=L(base); thr=12*(1-cap[k+1])
    rows=[]
    for rest in combinations(range(1,D+1),k-1):
        E=[0]+list(rest); rows.append((float(L(E)), add_energy(E), diff_doubling(E), longest_AP(E), E))
    rows.sort(key=lambda r:r[0])
    print(f"  consec_{k}: L={float(Lb):.5f}  threshold 12(1-cap)={float(thr):.5f}  margin={float(Lb-thr):.5f}")
    print(f"  energy(consec)={add_energy(base)}  |E-E|(consec)={diff_doubling(base)}  longAP={longest_AP(base)}")
    print(f"  --- smallest L (consec should be #1): L | energy | |E-E| | longAP | set ---")
    for r in rows[:5]:
        tag=" <--consec" if r[4]==base else ""
        print(f"     L={r[0]:.5f}  en={r[1]:>4}  |E-E|={r[2]:>3}  lAP={r[3]:>2}  {r[4]}{tag}")
    print(f"  --- LARGEST L (should be low-energy / far-from-AP): ---")
    for r in rows[-3:]:
        print(f"     L={r[0]:.5f}  en={r[1]:>4}  |E-E|={r[2]:>3}  lAP={r[3]:>2}  {r[4]}")
    # correlation: does L decrease as energy increases? (Freiman link)
    import statistics
    Ls=[r[0] for r in rows]; ens=[r[1] for r in rows]
    n=len(rows); mL=statistics.mean(Ls); me=statistics.mean(ens)
    cov=sum((rows[i][0]-mL)*(rows[i][1]-me) for i in range(n))/n
    corr=cov/(statistics.pstdev(Ls)*statistics.pstdev(ens))
    print(f"  corr(L, additive energy) = {corr:.3f}  (strongly negative => Freiman link: high energy => low L)")
    print(f"  is consec argmin(L)? {rows[0][4]==base}   gap to 2nd = {rows[1][0]-float(Lb):.5f}\n")
