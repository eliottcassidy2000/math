#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""A FAMILY OF GRADERS on the Z_7 cores -- are they ONE axis or several? (klein-S22).

Each core O subset Z_p has an autocorrelation a(d)=#{x:x,x+d in O} and Fourier spectrum lam_k=|Ohat(k)|^2.
Build several spread/concentration graders (some standard, some invented for this problem) and rank-correlate
them to see whether the apex landscape has ONE concentration<->spread axis or several independent invariants.
  g      = gap = min_{k!=0} lam_k                         (THM-590; the WORST-mode concentration)
  F      = spectral FLATNESS = GM(lam_{1..p-1})/AM(...)   (global flatness, [0,1], 1=difference set)
  IPR    = (sum lam)^2 / (p * sum lam^2)                  (effective spread fraction, [0,1])
  Hs     = spectral ENTROPY (normalized, [0,1])           (Shannon of lam_k/sum)
  D      = difference-set DEFECT = Var_{d!=0}(a(d))       (NEW: distance from a perfect difference set; 0=diffset)
  godd   = Cayley ODD-GIRTH of Cay(Z_p, supp a)           (NEW, on-theme: shortest odd cycle; doublet->C_7=p)
"""
import math, cmath, itertools
import numpy as np
P=7; W=cmath.exp(2j*math.pi/P)
def autocorr(O):
    Os=set(O); return [sum(1 for x in Os if (x+d)%P in Os) for d in range(P)]
def spec(O): return [abs(sum(W**((k*x)%P) for x in set(O)))**2 for k in range(P)]
def gap(O): s=spec(O); return min(s[k] for k in range(1,P))
def flat(O):
    s=[spec(O)[k] for k in range(1,P)]
    if min(s)<=1e-12: return 0.0
    gm=math.exp(sum(math.log(x) for x in s)/len(s)); am=sum(s)/len(s); return gm/am
def ipr(O):
    s=spec(O); num=sum(s)**2; den=P*sum(x*x for x in s); return num/den
def sentropy(O):
    s=spec(O); tot=sum(s); ps=[x/tot for x in s if x>1e-12]
    H=-sum(p*math.log(p) for p in ps); return H/math.log(P)
def defect(O):
    a=autocorr(O)[1:]; m=sum(a)/len(a); return sum((x-m)**2 for x in a)/len(a)
def odd_girth(O):
    # Cayley graph on Z_7 with connection set = {d!=0 : a(d)>0} (symmetric). shortest ODD cycle.
    a=autocorr(O); S=[d for d in range(1,P) if a[d]>0]
    if not S: return math.inf
    # BFS bipartite-style: shortest odd closed walk = min over edges of (1+dist_even(u,v))... do brute small p
    adj={v:set((v+d)%P for d in S) for v in range(P)}
    best=math.inf
    for start in range(P):
        # BFS distances parity
        from collections import deque
        dist={}; dq=deque([(start,0)]); dist[(start,0)]=0
        while dq:
            u,par=dq.popleft()
            for w in adj[u]:
                np_=1-par
                if (w,np_) not in dist:
                    dist[(w,np_)]=dist[(u,par)]+1; dq.append((w,np_))
        if (start,1) in dist: best=min(best,dist[(start,1)])
    return best

shapes=[("doublet {0,1}",{0,1}),("pair-d2 {0,2}",{0,2}),("arc-3 {0,1,2}",{0,1,2}),
        ("diffset {1,2,4}",{1,2,4}),("singleton {0}",{0}),("arc-4 {0,1,2,3}",{0,1,2,3}),
        ("5-core {0,1,2,3,4}",{0,1,2,3,4}),("co-singleton",{0,1,2,3,4,5})]
print("="*100)
print(f"{'shape':<20}{'|O|':>4}{'gap':>8}{'flatF':>8}{'IPR':>7}{'entH':>7}{'defectD':>9}{'oddgirth':>9}")
print("="*100)
for nm,O in shapes:
    g=odd_girth(O)
    print(f"{nm:<20}{len(O):>4}{gap(O):>8.3f}{flat(O):>8.3f}{ipr(O):>7.3f}{sentropy(O):>7.3f}{defect(O):>9.3f}{(g if g<math.inf else 'inf'):>9}")

# rank-correlations over ALL 127 cores: do the graders agree (one axis) or diverge?
print("\n"+"="*100)
print(" RANK CORRELATIONS over all 127 cores (Spearman): which graders are the SAME axis vs distinct?")
print("="*100)
rows=[]
for r in range(1,P+1):
    for O in itertools.combinations(range(P),r):
        Os=set(O)
        gg=odd_girth(Os)
        rows.append([gap(Os),flat(Os),ipr(Os),sentropy(Os),defect(Os),(gg if gg<math.inf else 99),len(Os)])
A=np.array(rows,dtype=float)
names=["gap","flatF","IPR","entH","defectD","oddgirth","|O|"]
def spearman(x,y):
    rx=np.argsort(np.argsort(x)); ry=np.argsort(np.argsort(y))
    return np.corrcoef(rx,ry)[0,1]
print(f"\n{'':>10}"+"".join(f"{n:>9}" for n in names))
for i,n in enumerate(names):
    print(f"{n:>10}"+"".join(f"{spearman(A[:,i],A[:,j]):>9.2f}" for j in range(len(names))))
print("\n  (|corr|~1 => same axis; ~0 => independent grader. Look for graders that are NOT just g in disguise.)")
