#!/usr/bin/env python3
"""ts_tail_lemma_referee_kps_S128c38.py -- kind-pasteur S128 cont.38.
THE T_s TAIL LEMMA referee: for dissociated triples/quadruples (no relation with max coeff <= H),
the full-support sinc-product tail T_s = sum over relations of prod min(2lam, 1/(pi|h_i|))
obeys T_s <= C_s (1+ln(2+Vmax))^2 / H.  Compute exact truncated T_s for real dissociated
subsets at varying H; fit the envelope; verify the three-regime structure."""
import sys
from math import gcd, sin, pi, log
sys.stdout.reconfigure(line_buffering=True)
LAM=1/14
def c(h): return 2*LAM if h==0 else abs(sin(2*pi*h*LAM)/(pi*h))
def T3_exact(a,b,cc,H,BOX=4000):
    """full-support relation tail for triple (a,b,c): sum over h1 a + h2 b + h3 c = 0,
    all h_i != 0, max|h| > H, |h_i| <= BOX, of prod |c(h_i)|. Also return by-regime split."""
    g=gcd(a,b); a2,b2=a//g,b//g
    tot=0.0; reg={'smallt':0.0,'bigt':0.0}
    for t in range(-BOX,BOX+1):
        if t==0: continue
        if (t*cc)%g!=0: continue
        # solve h1 a + h2 b = -t c ; particular solution via extended gcd
        # ext gcd for a2, b2 (coprime)
        x0,y0,g2=None,None,None
        # extended euclid
        old_r,r=a2,b2; old_s,s=1,0; old_t2,t2=0,1
        while r:
            q=old_r//r
            old_r,r=r,old_r-q*r
            old_s,s=s,old_s-q*s
            old_t2,t2=t2,old_t2-q*t2
        inv_ok=old_r==1
        rhs=-(t*cc)//g
        x0=old_s*rhs; y0=old_t2*rhs
        # solutions: (x0 + m b2, y0 - m a2)
        # m range so |h1|,|h2| <= BOX
        mlo=int((-BOX-x0)/b2)-2; mhi=int((BOX-x0)/b2)+2
        for m in range(mlo,mhi+1):
            h1=x0+m*b2; h2=y0-m*a2
            if h1==0 or h2==0: continue
            if abs(h1)>BOX or abs(h2)>BOX: continue
            if max(abs(h1),abs(h2),abs(t))<=H: continue
            v=c(h1)*c(h2)*c(t)
            tot+=v
            if abs(t)<=H: reg['smallt']+=v
            else: reg['bigt']+=v
    return tot,reg
print("== dissociated triples: T3(H) vs C/H envelope ==")
triples=[(307,425,541),(671,944,1413),(541,800,2147),(1087,1943,3310),(313,741,1531)]
for (a,b,cc) in triples:
    # verify dissociation floor: smallest relation max-coeff
    best=10**9
    for h1 in range(-60,61):
        for h2 in range(-60,61):
            if h1==0 or h2==0: continue
            r=h1*a+h2*b
            if r!=0 and r%cc==0 and abs(r//cc)<=60 and r//cc!=0:
                best=min(best,max(abs(h1),abs(h2),abs(r//cc)))
    line="(%d,%d,%d) lam1>%s |"%(a,b,cc,best if best<10**9 else '60')
    for H in (5,10,20,40,80):
        t3,reg=T3_exact(a,b,cc,H,BOX=3000)
        line+="  H=%d: %.2e (t<=H: %.0f%%)"%(H,t3,100*reg['smallt']/t3 if t3>0 else 0)
    print(line)
print("== envelope check: T3(H)*H/(1+ln(2+Vmax))^2 should be ~bounded ==")
for (a,b,cc) in triples:
    line="(%d,%d,%d):"%(a,b,cc)
    for H in (5,10,20,40,80):
        t3,_=T3_exact(a,b,cc,H,BOX=3000)
        env=t3*H/(1+log(2+cc))**2
        line+="  H=%d: %.3e"%(H,env)
    print(line)
print("== budget check at the exhaustion horizon: sum over all C(13,3) triples of a real packet ==")
pk=[307,425,541,671,800,944,1087,1413,1943,2147,2570,3056,3310]
from itertools import combinations
for H in (10,30):
    tot=0.0
    for A in combinations(pk,3):
        t3,_=T3_exact(A[0],A[1],A[2],H,BOX=1200)
        tot+=t3
    print("  H=%d: sum T3 over 286 triples = %.4f ; x|E3|=24/49 -> %.4f  (budget 2052/16807 = 0.1221)"%(H,tot,tot*24/49))
print("DONE")
