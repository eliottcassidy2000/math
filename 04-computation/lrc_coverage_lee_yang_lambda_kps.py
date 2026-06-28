"""
lrc_coverage_lee_yang_lambda_kps.py  (kind-pasteur-2026-06-27-S31ah)

Deepen HYP-3099: the LRC coverage PGF Q(z)=sum q_t z^t has zeros on a near-circle;
the AP/consec maximizer is the MOST circular. Measure the phi^4 coupling
  lambda(E) = Var(|roots|/R)   (normalized off-circle variance, R=geomean radius)
and test:
  (1) does AP MINIMIZE lambda (most circular = phi^4 free/Gaussian limit)?
  (2) self-inversive/palindrome signature of zeros-on-circle: q_t * R^{2t} ~ q_{6-t} * R^{2(6-t)}?
      (a circle polynomial is self-inversive => the coefficient sequence is R-palindromic)
  (3) the circle radius R(k) and a closed-form guess;
  (4) does lambda track the coverage gap (cap - q_0)?
"""
import sys, itertools, random
from fractions import Fraction as F
import numpy as np

def sector_of(p): return int((p%1)*7)
def Ndist(E):
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

def lee_yang_signals(q):
    # Q(z)=sum_{t=0}^6 q_t z^t ; need q_6>0 for 6 finite roots
    coeffs=list(reversed(q))
    while len(coeffs)>1 and abs(coeffs[0])<1e-14: coeffs=coeffs[1:]
    if len(coeffs)<=1: return None
    roots=np.roots(coeffs)
    mags=np.abs(roots)
    R=float(np.exp(np.mean(np.log(mags))))     # geometric-mean radius
    lam=float(np.var(mags/R))                   # normalized off-circle variance (phi^4 coupling)
    # self-inversive palindrome residual: q_t R^{2t} vs q_{6-t} R^{2(6-t)} (scaled)
    deg=len(q)-1
    pal=0.0; npal=0
    for t in range(deg+1):
        lhs=q[t]*R**(2*t); rhs=q[deg-t]*R**(2*(deg-t))
        s=abs(lhs)+abs(rhs)
        if s>1e-12: pal+=abs(lhs-rhs)/s; npal+=1
    pal=pal/max(npal,1)
    return dict(R=R, lam=lam, pal=pal, q0=q[0], roots=roots, mags=mags)

if __name__=="__main__":
    sys.stdout.reconfigure(line_buffering=True)
    random.seed(3099)
    CAPS={8:2243/5880,9:1979/4004,10:55/91}
    for k in (8,9,10):
        consec=tuple(range(k))
        cap=CAPS[k]
        sig_c=lee_yang_signals(Ndist(consec))
        pool=[("consec",consec)]
        for _ in range(500):
            cfg=tuple(sorted([0]+random.sample(range(1,26),k-1)))
            if len(set(cfg))==k: pool.append(("rand",cfg))
        # measure lambda for all, find argmin
        data=[]
        for name,E in pool:
            s=lee_yang_signals(Ndist(E))
            if s is None: continue
            data.append((s['lam'],name,E,s))
        data.sort(key=lambda d: d[0])
        argmin=data[0]
        # is consec the argmin of lambda?
        consec_rank=[i for i,(lam,nm,E,s) in enumerate(data) if E==consec]
        print(f"\n=== k={k}  cap={cap:.4f} ===")
        print(f"  consec: q0={sig_c['q0']:.4f} (cap gap {cap-sig_c['q0']:+.4f})  "
              f"R={sig_c['R']:.4f}  lambda(phi4)={sig_c['lam']:.5f}  palindrome-resid={sig_c['pal']:.5f}")
        print(f"  consec lambda-RANK among {len(data)} configs: {consec_rank[0] if consec_rank else '?'} "
              f"(0 = MOST circular)")
        print(f"  global argmin-lambda: {argmin[1]} lam={argmin[0]:.5f} q0={argmin[3]['q0']:.4f} "
              f"(is consec? {argmin[2]==consec})")
        # correlation: lambda vs coverage gap (cap-q0)
        lams=[d[0] for d in data]; gaps=[cap-d[3]['q0'] for d in data]
        ml=np.mean(lams); mg=np.mean(gaps)
        cov=sum((l-ml)*(g-mg) for l,g in zip(lams,gaps))
        sl=(sum((l-ml)**2 for l in lams))**0.5; sg=(sum((g-mg)**2 for g in gaps))**0.5
        pear=cov/(sl*sg) if sl*sg>0 else 0
        print(f"  corr(lambda, cap-q0) = {pear:+.3f}  (lambda controls the coverage gap?)")
        # circle radius closed-form guesses
        print(f"  R(consec)={sig_c['R']:.5f}  R^6={sig_c['R']**6:.4f} = q0/q6 = {sig_c['q0']/Ndist(consec)[6]:.4f}")
        print(f"  consec root mags: {[f'{m:.4f}' for m in sorted(sig_c['mags'])]}")
