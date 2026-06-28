"""
lrc_k8_two_new_angles_griffiths_caratheodory_kps.py  (kind-pasteur-2026-06-27-S31al)

Two NEW creative angles to attack the k=8 node (consec maximizes coverage/covariance):

ANGLE 1 -- GRIFFITHS/GKS FERROMAGNETIC MONOTONICITY. The empty-sector indicators are a
ferromagnetic spin system (HYP-3161, all 15 Cov>0 at consec). Griffiths-II: correlations
INCREASE with couplings. If there is a MONOTONE PATH from any config to consec along which
Sigma-kappa_2 (and coverage q0) increases, the ground-state extremality is a Griffiths
monotonicity. TEST: build paths config -> consec (move one speed at a time toward the AP),
check Sigma-kappa_2 is monotone non-decreasing.

ANGLE 2 -- CARATHEODORY-TOEPLITZ MOMENT EXTREMALITY. My Lee-Yang zeros-on-a-circle (HYP-3099)
IS the classical trigonometric moment problem: the coverage PGF's coefficients q_t are a
moment sequence; Caratheodory's theorem says the extremal positive measures with given moments
are ATOMIC on the circle (the atoms = the de Moivre angles). TEST: form the Toeplitz/Hankel
matrix of the coverage moments; is its smallest eigenvalue (the Caratheodory PSD margin)
extremal at consec? (consec = the extremal/rigid moment configuration.)
"""
import sys, itertools, random
import numpy as np
INNER=list(range(1,7))
def cells(E):
    E=sorted(set(E)); b=set([0.0,1.0])
    for e in E:
        if e==0: continue
        for mm in range(0,7*e+1): b.add(mm/(7*e))
    b=sorted(b); out=[]
    for i in range(len(b)-1):
        x0,x1=b[i],b[i+1]
        if x1-x0<1e-15: continue
        cov=set(int(((e*(x0+x1)/2)%1)*7) for e in E)
        out.append((cov,x1-x0))
    return out
def cov_sum_and_q(E):
    cl=cells(E); p={i:0.0 for i in INNER}; pij={}; q=[0.0]*7
    for cov,w in cl:
        empty=[i for i in INNER if i not in cov]
        for i in empty: p[i]+=w
        for a in range(len(empty)):
            for b2 in range(a+1,len(empty)):
                k=(empty[a],empty[b2]); pij[k]=pij.get(k,0.0)+w
        t=len(empty); q[t]+=w
    S2c=sum(pij.get((i,j),0.0)-p[i]*p[j] for i,j in itertools.combinations(INNER,2))
    return S2c, q

def caratheodory_margin(q):
    # moments c_m = sum_t q_t e^{-2pi i m t/7}? use the COVERAGE moments as a Toeplitz seq.
    # build Toeplitz from q_t treated as Fourier coeffs of a measure on the circle.
    q=np.array(q,dtype=float)
    # c_0..c_6 = q (real); Hermitian Toeplitz T[j,k]=c_{|j-k|}
    n=7; T=np.zeros((n,n))
    for j in range(n):
        for k in range(n):
            T[j,k]=q[abs(j-k)]
    ev=np.linalg.eigvalsh(T)
    return ev.min(), ev.max()

if __name__=="__main__":
    sys.stdout.reconfigure(line_buffering=True)
    random.seed(8)
    consec=tuple(range(8))
    Sc,qc=cov_sum_and_q(consec)
    print(f"consec_8: Sigma-kappa_2={Sc:.5f}  q0={qc[0]:.5f}")

    print("\n=== ANGLE 1: GRIFFITHS MONOTONICITY (path config -> consec, is Sigma-k2 monotone up?) ===")
    monotone_paths=0; total_paths=0; nonmono_examples=[]
    for trial in range(60):
        # random start config (with 0)
        start=tuple(sorted([0]+random.sample(range(1,18),7)))
        if len(set(start))!=8: continue
        # path: at each step move the speed that is 'farthest' toward making it consec={0..7}
        # greedy: sort current speeds, move the largest-deviating to its consec slot
        cur=list(start); path=[cov_sum_and_q(tuple(cur))[0]]
        target=list(range(8))
        for step in range(8):
            cur_sorted=sorted(cur)
            # find first index where cur_sorted != target, set it
            changed=False
            for idx in range(8):
                if cur_sorted[idx]!=target[idx]:
                    cur=cur_sorted[:idx]+[target[idx]]+cur_sorted[idx+1:]
                    # dedupe collisions by bumping
                    cur=sorted(set(cur))
                    while len(cur)<8:
                        x=max(cur)+1; cur.append(x)
                    cur=cur[:8]; changed=True; break
            if not changed: break
            path.append(cov_sum_and_q(tuple(sorted(cur)))[0])
        total_paths+=1
        diffs=[path[i+1]-path[i] for i in range(len(path)-1)]
        if all(d>=-1e-9 for d in diffs): monotone_paths+=1
        elif len(nonmono_examples)<3: nonmono_examples.append((start,[round(x,4) for x in path]))
    print(f"  monotone (non-decreasing) paths to consec: {monotone_paths}/{total_paths}")
    if nonmono_examples:
        print(f"  non-monotone examples (Griffiths path FAILS naively):")
        for s,p in nonmono_examples: print(f"    {s}: {p}")
    print(f"  => Griffiths monotonicity {'HOLDS on greedy paths (proof route viable)' if monotone_paths==total_paths else 'needs the RIGHT path/partial-order (naive greedy insufficient)'}")

    print("\n=== ANGLE 2: CARATHEODORY-TOEPLITZ margin (extremal at consec?) ===")
    cm_c=caratheodory_margin(qc)
    print(f"  consec Toeplitz(q) eigenvalues: min={cm_c[0]:.5f} max={cm_c[1]:.5f}")
    pool=[("consec",consec)]
    for _ in range(300):
        cfg=tuple(sorted([0]+random.sample(range(1,20),7)))
        if len(set(cfg))==8: pool.append(("r",cfg))
    data=[]
    for nm,E in pool:
        _,q=cov_sum_and_q(E); mn,mx=caratheodory_margin(q)
        data.append((mn,E))
    # is consec extremal (min or max) in the Caratheodory margin?
    smin=sorted(data); smax=sorted(data,reverse=True)
    rmin=[i for i,d in enumerate(smin) if d[1]==consec][0]
    rmax=[i for i,d in enumerate(smax) if d[1]==consec][0]
    print(f"  consec Caratheodory min-eigenvalue rank: {rmin} (from below) / {rmax} (from above) of {len(data)}")
    print(f"  => consec is {'EXTREMAL (rigid moment config)' if rmin==0 or rmax==0 else 'NOT extremal in the raw Toeplitz margin'}")
