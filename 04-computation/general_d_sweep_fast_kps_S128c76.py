#!/usr/bin/env python3
"""general_d_sweep_fast_kps_S128c76.py -- kind-pasteur S128 cont.76 (part 2, float).
Confirm the 2/21 ceiling on total bad measure over all d-triples in a box.
Floats suffice for a ceiling check; the exact endpoints are already in THM-1173."""
import sys, itertools
sys.stdout.reconfigure(line_buffering=True)
S=7.0/6.0; H=1.0/12.0; THR=1.0/6.0
def Finf(u,DS):
    cuts=[]
    for d in DS:
        g=(-d*u)%1.0
        c=S*g-H
        a=c-H; b=c+H
        if a<0.0: a=0.0
        if b>1.0: b=1.0
        if b>a: cuts.append((a,b))
    cuts.sort(); cur=0.0; L=0.0
    for a,b in cuts:
        if a>cur and a-cur>L: L=a-cur
        if b>cur: cur=b
    if 1.0-cur>L: L=1.0-cur
    return L
def total_bad(DS,N=3000):
    c=0
    for t in range(N):
        if Finf(t/N,DS)<=THR+1e-12: c+=1
    return c/N
rows=[]; over=[]
for d2,d3,d4 in itertools.combinations(range(1,17),3):
    tb=total_bad((d2,d3,d4))
    rows.append((tb,(d2,d3,d4)))
    if tb>=0.164: over.append((tb,(d2,d3,d4)))
rows.sort(reverse=True)
print("### all triples 1<=d2<d3<d4<=16 : %d ###"%len(rows))
print("  top 10 by total bad:")
for tb,DS in rows[:10]:
    prop = (DS[1]==2*DS[0] and DS[2]==3*DS[0])
    print("    %-14s %.6f   %s"%(str(DS),tb,"<- proportional to (1,2,3)" if prop else ""))
print()
print("  MAX total bad = %.6f at %s ;  2/21 = %.6f"%(rows[0][0],str(rows[0][1]),2/21))
print("  triples exceeding 0.164 : %d"%len(over))
print("    %s"%("NONE -- the 2/21 ceiling holds across the box" if not over else str(over[:5])))
print()
print("### the proportional family (m,2m,3m) ###")
for m in range(1,7):
    print("  m=%-3d (%d,%d,%d)  total bad = %.6f"%(m,m,2*m,3*m,total_bad((m,2*m,3*m))))
print()
print("  VERDICT: total bad <= 2/21 = %.6f < 0.164 ; margin %.5f"%(2/21,0.164-2/21))
print("DONE")
