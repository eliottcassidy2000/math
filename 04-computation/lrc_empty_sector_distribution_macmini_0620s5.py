#!/usr/bin/env python3
"""
lrc_empty_sector_distribution_macmini_0620s5.py  (mac-mini-2026-06-20-S5)

The coverage is p0=P(N=0), N=#empty nonzero sectors. The moment-LP (THM-534) bounds P(N=0) by
S_r=E[C(N,r)]. This script computes:
 (1) the exact empty-sector-count distribution p_t=P(N=t), t=0..6;
 (2) the REFLECTION structure: x->-x maps sector j -> 7-j, a MEASURE symmetry, so single-sector
     miss measures pair m_1=m_6, m_2=m_5, m_3=m_4. Verify, and check joint structure.
 (3) the per-sector miss measures m_j and whether QR{1,2,4} vs NQR{3,5,6} differ (cube-root probe).
This tells us what symmetry a refined moment-LP can exploit.
"""
import itertools, sys
from fractions import Fraction as F
sys.stdout.reconfigure(line_buffering=True)
def sector_of(p): return int((p%1)*7)
def bps(E):
    E=sorted(set(E)); b=set([F(0),F(1)])
    for e in E:
        if e==0: continue
        for m in range(0,7*e+1): b.add(F(m,7*e))
    return sorted(b)
def analyze(E):
    B=bps(E); pt=[F(0)]*7; mj=[F(0)]*7  # mj[j]=meas{sector j empty}
    for i in range(len(B)-1):
        x0,x1=B[i],B[i+1]
        if x1<=x0: continue
        secs=set(sector_of(e*((x0+x1)/2)) for e in E); w=x1-x0
        empty=[j for j in range(1,7) if j not in secs]   # nonzero empty sectors
        pt[len(empty)]+=w
        for j in empty: mj[j]+=w
    return pt, mj

rows={
 'consec_8':(0,1,2,3,4,5,6,7),
 'consec_9':(0,1,2,3,4,5,6,7,8),
 'true-wide leader':(0,4,6,8,10,12,14,15,16),
 'dyadic (0,1,2,4,8,12,16,20)':(0,1,2,4,8,12,16,20),
}
for name,E in rows.items():
    pt,mj=analyze(E)
    print(f"\n{name}  E={E}")
    print(f"  p_t (N=#empty, t=0..6): {[f'{float(x):.4f}' for x in pt]}   sum={float(sum(pt)):.4f}")
    print(f"  m_j (sector j empty), j=1..6: {[f'{float(mj[j]):.4f}' for j in range(1,7)]}")
    # reflection pairs {1,6},{2,5},{3,4}
    refl_ok = (mj[1]==mj[6] and mj[2]==mj[5] and mj[3]==mj[4])
    print(f"  reflection symmetry m1=m6,m2=m5,m3=m4: {refl_ok}  "
          f"(pairs: {float(mj[1]):.4f}={float(mj[6]):.4f}, {float(mj[2]):.4f}={float(mj[5]):.4f}, {float(mj[3]):.4f}={float(mj[4]):.4f})")
    qr=sum(mj[j] for j in [1,2,4]); nqr=sum(mj[j] for j in [3,5,6])
    print(f"  QR-sum m1+m2+m4={float(qr):.4f}  NQR-sum m3+m5+m6={float(nqr):.4f}  (cube-root probe: {'equal' if qr==nqr else 'DIFFER by '+str(float(qr-nqr))})")
