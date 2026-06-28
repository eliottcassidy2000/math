"""
lrc_joukowski_cyclotomic_dip_kps.py  (kind-pasteur-2026-06-27-S31ak)

Owner's connection: the IDEAL uniform PGF 1+z+...+z^6 (perfect 7-fold symmetry) has
roots at the 7th roots of unity; Joukowski w=z+1/z maps them to the de Moivre angles
{2cos(2pi j/7)} = {1.2470, -0.4450, -1.8019} (roots of the cyclotomic cubic x^3+x^2-2x-1).
So cap = the 7th-cyclotomic IDEAL; dip = deviation from 7-fold symmetry = Im(w) = the
resolvent's REAL-ROOTEDNESS DEFECT.

This LINEARIZES the Lee-Yang circle (HYP-3099): the coverage PGF Q(z)=sum q_t z^t has 6
zeros in 3 conjugate pairs near |z|=R; the R-normalized Joukowski w_j = z_j/R + R/z_j
sends them to 3 points near the real axis. If exactly on the circle => w real = 2cos(theta_j).
TESTS:
 (1) Re(w_j) vs the de Moivre angles {2cos(2pi j/7)} -- does the coverage approach the 7-fold ideal?
 (2) Im(w_j) = the off-circle/real-rootedness DEFECT -- does it track the dip and minimize at consec?
 (3) the zero-angles theta_j vs 2pi j/7 (RATIONAL APPROXIMATION of the 7-fold symmetry).
 (4) MERGE ferromagnetic: does the cyclotomic defect minimize at the ferromagnetic ground state?
"""
import sys, random
from fractions import Fraction as F
import numpy as np

DE_MOIVRE = sorted([2*np.cos(2*np.pi*j/7) for j in (1,2,3)])  # [-1.802,-0.445,1.247]

def sector_of(p): return int((p%1)*7)
def qdist(E):
    E=sorted(set(E)); b=set([F(0),F(1)])
    for e in E:
        if e==0: continue
        for mm in range(0,7*e+1): b.add(F(mm,7*e))
    b=sorted(b); q=[F(0)]*7
    for i in range(len(b)-1):
        x0,x1=b[i],b[i+1]
        if x1<=x0: continue
        cov=set(sector_of(e*((x0+x1)/2)) for e in E)
        t=7-len(cov)
        if 0<=t<=6: q[t]+=x1-x0
    return [float(x) for x in q]

def joukowski_analysis(E):
    q=qdist(E); coeffs=list(reversed(q))
    while len(coeffs)>1 and abs(coeffs[0])<1e-14: coeffs=coeffs[1:]
    if len(coeffs)<=1: return None
    roots=np.roots(coeffs); mags=np.abs(roots); R=np.exp(np.mean(np.log(mags)))
    # R-normalized Joukowski w = z/R + R/z
    W=roots/R + R/roots
    # group into 3 (conjugate pairs) by Re(w); take upper-half representatives
    order=np.argsort(W.real)
    Wr=W[order]
    # 6 values come as 3 conj pairs; dedupe by Re
    reps=[]; used=[False]*6
    for i in range(6):
        if used[i]: continue
        # find conjugate partner (same Re, opposite Im)
        reps.append(Wr[i]); used[i]=True
        for j in range(i+1,6):
            if not used[j] and abs(Wr[j].real-Wr[i].real)<0.05 and abs(Wr[j].imag+Wr[i].imag)<0.05:
                used[j]=True; break
    reps=sorted(reps,key=lambda w:w.real)[:3]
    Re=[w.real for w in reps]; Im=[abs(w.imag) for w in reps]
    defect=sum(Im)  # total real-rootedness defect = the dip signature
    # distance of Re parts to de Moivre angles
    dm_dist=sum(abs(Re[i]-DE_MOIVRE[i]) for i in range(min(3,len(Re))))
    return dict(R=R,Re=Re,Im=Im,defect=defect,dm_dist=dm_dist,q0=q[0])

if __name__=="__main__":
    sys.stdout.reconfigure(line_buffering=True)
    random.seed(7)
    print(f"de Moivre angles (ideal 7-fold, Joukowski of 7th roots of unity): {[round(a,4) for a in DE_MOIVRE]}")
    print(f"cyclotomic cubic x^3+x^2-2x-1\n")
    CAPS={8:2243/5880,9:1979/4004,10:55/91}
    DIPS={8:1081/76440,9:1/4004,10:0.0}
    for k in (8,9,10):
        E=tuple(range(k)); a=joukowski_analysis(E)
        print(f"=== consec_{k} (R={a['R']:.4f}) ===")
        print(f"  Joukowski Re(w) = {[round(x,4) for x in a['Re']]}  vs de Moivre {[round(x,4) for x in DE_MOIVRE]}  (dist {a['dm_dist']:.4f})")
        print(f"  Joukowski Im(w) = {[round(x,4) for x in a['Im']]}  total DEFECT={a['defect']:.5f}  dip_{k}={DIPS[k]:.5f}")
    # (2)+(4): does the cyclotomic defect track coverage / minimize at consec (ferromagnetic ground state)?
    print("\n=== defect vs coverage over k=8 clusters (does consec MINIMIZE the cyclotomic defect?) ===")
    consec=tuple(range(8)); ac=joukowski_analysis(consec)
    pool=[("consec",consec)]
    for _ in range(400):
        cfg=tuple(sorted([0]+random.sample(range(1,20),7)))
        if len(set(cfg))==8: pool.append(("rand",cfg))
    data=[]
    for nm,E in pool:
        a=joukowski_analysis(E)
        if a: data.append((a['defect'],a['q0'],E))
    data.sort()
    consec_rank=[i for i,d in enumerate(data) if d[2]==consec][0]
    print(f"  consec defect={ac['defect']:.5f} q0={ac['q0']:.4f}  rank={consec_rank}/{len(data)} (0=min defect)")
    print(f"  min-defect config: defect={data[0][0]:.5f} q0={data[0][1]:.4f} (is consec? {data[0][2]==consec})")
    defs=[d[0] for d in data]; q0s=[d[1] for d in data]
    md=np.mean(defs); mq=np.mean(q0s)
    cov=sum((d-md)*(q-mq) for d,q in zip(defs,q0s)); sd=(sum((d-md)**2 for d in defs))**.5; sq=(sum((q-mq)**2 for q in q0s))**.5
    print(f"  corr(cyclotomic defect, coverage q0) = {cov/(sd*sq) if sd*sq>0 else 0:+.3f}  (defect lowers coverage?)")
