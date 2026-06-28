"""
lrc_k8_node_bimodality_newton_maclaurin_kps.py  (kind-pasteur-2026-06-27-S31ai)

THE LAST NODE (mac-mini HYP-3132): the entire LRC(14) proof reduces to the k=8
bounded-core dip. The dual L_yK8 = q0 + q6 + (1/10) q3 is the BIMODALITY functional
(mass at the two extremes N=0 [all covered] and N=6 [all empty]); consec must MAXIMIZE it.
The resolvent is the biquadratic phi^4 potential u^4 - 5u^2 + 4 (Galois-solvable).

Owner pointers merged here:
 - "Newton/Maclaurin quartic moment inequality extremal at the AP"
 - "k=8 dual is L_y = p0 + p6 + (1/10)p3, consec maximizes it"
 - "q0 = q6 R^6", Lee-Yang circle, phi^4 off-circle = the dip
 - information theory (entropy of N), everything as FUNCTIONS

For each k=8 cluster E, compute the empty-count distribution q_t (FUNCTION on N in {0..6})
and the FUNCTIONALS:
 - L_yK8 = q0 + q6 + (1/10)q3      (bimodality dual; consec maximizes?)
 - moments S1..S4, central moments, kurtosis kappa4/kappa2^2 (the phi^4 coupling)
 - Shannon entropy Hsh(N)          (information content)
 - Newton-Maclaurin: p_k = e_k/C(6,k) from the PGF; check p_k^2 >= p_{k-1} p_{k+1}
   and which inequality is TIGHT / extremal at consec
 - q6*R^6 vs q0 (Vieta check)
Find: which single functional does consec EXACTLY extremize? (the clean proof target)
"""
import sys, itertools, random
from fractions import Fraction as F
from math import comb, log2
import numpy as np

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
    return q

def functionals(q):
    qf=[float(x) for x in q]
    L_y = qf[0] + qf[6] + qf[3]/10.0                  # bimodality dual
    m1 = sum(t*qf[t] for t in range(7))               # E[N]
    m2 = sum(t*t*qf[t] for t in range(7))
    var = m2 - m1*m1
    mu3 = sum((t-m1)**3*qf[t] for t in range(7))
    mu4 = sum((t-m1)**4*qf[t] for t in range(7))
    kurt = mu4/var**2 if var>1e-12 else 0             # kappa4/kappa2^2 + 3
    exkurt = kurt - 3                                  # excess kurtosis (phi^4 sign)
    bimod = qf[0] + qf[6]                              # raw bimodality (the two wells)
    # Shannon entropy
    Hsh = -sum(p*log2(p) for p in qf if p>1e-15)
    # PGF roots, R, q6 R^6 vs q0
    coeffs=list(reversed(qf))
    while len(coeffs)>1 and abs(coeffs[0])<1e-14: coeffs=coeffs[1:]
    R=None; q6R6=None
    if len(coeffs)>1:
        r=np.roots(coeffs); mg=np.abs(r); R=float(np.exp(np.mean(np.log(mg))))
        q6R6=qf[6]*R**6
    # Newton-Maclaurin on the PGF Q(z)=sum q_t z^t (treat q_t as e-like): p_k=q_k/C(6,k)
    p=[qf[k]/comb(6,k) for k in range(7)]
    nm=[]  # p_k^2 - p_{k-1}p_{k+1} (>=0 for real-rooted, Newton)
    for k in range(1,6):
        nm.append(p[k]**2 - p[k-1]*p[k+1])
    return dict(L_y=L_y,m1=m1,var=var,exkurt=exkurt,bimod=bimod,Hsh=Hsh,R=R,q6R6=q6R6,q0=qf[0],
                nm=nm, mu4=mu4)

if __name__=="__main__":
    sys.stdout.reconfigure(line_buffering=True)
    random.seed(8)
    k=8
    consec=tuple(range(k))
    fc=functionals(qdist(consec))
    print(f"=== THE k=8 NODE: consec={consec} ===")
    print(f"  L_yK8(bimodality)={fc['L_y']:.5f}  bimod(q0+q6)={fc['bimod']:.5f}  exkurt={fc['exkurt']:.4f}")
    print(f"  E[N]={fc['m1']:.4f} var={fc['var']:.4f} entropy={fc['Hsh']:.4f}  R={fc['R']:.4f} q6R6={fc['q6R6']:.5f} q0={fc['q0']:.5f}")
    print(f"  Newton p_k^2-p_(k-1)p_(k+1) (k=1..5): {[f'{x:+.5f}' for x in fc['nm']]}")
    # sample competitors; find what consec extremizes
    pool=[("consec",consec)]
    pool.append(("even-AP",tuple(2*i for i in range(k))))
    for _ in range(3000):
        cfg=tuple(sorted([0]+random.sample(range(1,20),k-1)))
        if len(set(cfg))==k: pool.append(("rand",cfg))
    rows=[]
    for nm,E in pool:
        f=functionals(qdist(E)); rows.append((nm,E,f))
    def extremal(key, maximize=True):
        vals=[(f[key],nm,E) for nm,E,f in rows]
        vals.sort(reverse=maximize)
        top=vals[0]; consec_rank=[i for i,(v,n,E) in enumerate(vals) if E==consec][0]
        return top, consec_rank
    print(f"\n  Over {len(rows)} k=8 clusters -- does consec EXTREMIZE each functional?")
    for key,mx,desc in [('L_y',True,'bimodality dual (MAX)'),('bimod',True,'q0+q6 raw bimodality (MAX)'),
                        ('exkurt',True,'excess kurtosis kappa4 (MAX)'),('var',True,'variance (MAX)'),
                        ('Hsh',False,'entropy (MIN)'),('mu4',True,'4th central moment (MAX)')]:
        top,rank=extremal(key,mx)
        print(f"   {desc:32s}: consec rank={rank:>3} (0=extremal)  argext={top[1]}{'<--consec' if top[2]==consec else ''} val={top[0]:.5f}")
