"""
tournament_lrc_pgf_root_structure_kps.py  (kind-pasteur-2026-06-27-S31ah)

THE TWO MAPS (owner directive: "the whole PGF curve and its root structure",
Lee-Yang extremality, phi^4; "what maximizes LRC values relating to tournaments").

Stop measuring SINGLE values (H = I(Omega,2)). Measure the WHOLE partition
function and its ROOT STRUCTURE (Lee-Yang zeros), and locate the EXTREMIZER.

MAP 1 (tournament): for all tournaments n=5,6,7, the independence polynomial
I(Omega,x) = hard-core partition function (THM-485). Signals: roots, smallest
|root| (the Lee-Yang edge near 0), real-rooted?, and correlation with H.
BOLD HYP: H-maximizer (regular/Paley) has a root hugging 0 (Lee-Yang edge).

MAP 2 (LRC): for k-clusters, the coverage PGF Q(z)=sum q_t z^t (N=#empty inner
sectors). Signals: roots, and correlation with coverage q_0=P(N=0).
BOLD HYP: AP/consec maximizer sits at the Lee-Yang edge of the coverage PGF.
"""
import sys, itertools, random
from fractions import Fraction as F
import numpy as np

sys.path.insert(0, '04-computation')
from tournament_certificate_engine_kps import (all_tournaments, conflict_graph,
    independence_poly, H_value, random_tournament, scores)

# ---------- MAP 1: tournament partition-function root structure ----------
def alpha_and_roots(adj):
    m,E=conflict_graph(adj)
    a=independence_poly(m,E)
    while len(a)>1 and a[-1]==0: a.pop()
    # I(x) = sum a[k] x^k ; numpy.roots wants highest-degree first
    coeffs=list(reversed(a))
    if len(coeffs)<=1:
        roots=np.array([])
    else:
        roots=np.roots(coeffs)
    H=sum(a[k]*2**k for k in range(len(a)))
    return a, roots, H

def map1(n, sample=None, rng=None):
    print(f"\n===== MAP 1: tournament I(Omega,x) root structure, n={n} =====")
    print("  H | alpha-vector | #roots realrooted? | min|root| (Lee-Yang edge) | root nearest 0")
    rows=[]
    src = all_tournaments(n) if sample is None else (random_tournament(n,rng) for _ in range(sample))
    seen_alpha=set()
    for adj in src:
        a,roots,H=alpha_and_roots(adj)
        key=tuple(a)
        if key in seen_alpha: continue
        seen_alpha.add(key)
        if len(roots)==0:
            rows.append((H,a,0,True,float('inf'),None)); continue
        realrooted=np.all(np.abs(roots.imag)<1e-9)
        mags=np.abs(roots)
        minmag=float(mags.min())
        nearest=roots[np.argmin(mags)]
        rows.append((H,tuple(a),len(roots),bool(realrooted),minmag,complex(nearest)))
    rows.sort(key=lambda r:-r[0])
    Hmax=rows[0][0]
    for H,a,nr,rr,minmag,near in rows[:12]:
        nearstr=f"{near.real:+.4f}{near.imag:+.4f}j" if near is not None else "-"
        tag=" <== H-MAX" if H==Hmax else ""
        print(f"  H={H:>6} a={a} roots={nr} real={rr} min|r|={minmag:.4f} near0={nearstr}{tag}")
    # the key correlation: does H-max have the smallest min|root|?
    by_minmag=sorted([r for r in rows if r[4]<float('inf')], key=lambda r:r[4])
    print(f"  --> smallest min|root| config: H={by_minmag[0][0]} min|r|={by_minmag[0][4]:.5f} "
          f"(is it H-max? {by_minmag[0][0]==Hmax})")
    print(f"  --> H-max config min|root| = {[r[4] for r in rows if r[0]==Hmax][0]:.5f}")
    # real-rooted fraction
    rr_frac=sum(1 for r in rows if r[3])/len(rows)
    print(f"  --> distinct alpha-classes: {len(rows)}; real-rooted fraction: {rr_frac:.3f}")
    return rows

# ---------- MAP 2: LRC coverage PGF root structure ----------
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

def map2(k, rng, ntrials=60):
    print(f"\n===== MAP 2: LRC coverage PGF Q(z)=sum q_t z^t root structure, k={k} =====")
    consec=tuple(range(k))
    configs=[("consec",consec),("even-AP",tuple(2*i for i in range(k))),
             ("single-far",consec[:-1]+(21,))]
    for _ in range(ntrials):
        cfg=tuple(sorted([0]+random.sample(range(1,24),k-1)))
        if len(set(cfg))==k: configs.append(("rand",cfg))
    rows=[]
    for name,E in configs:
        q=Ndist(E)
        coeffs=list(reversed(q))  # highest power first
        # trim leading zeros
        while len(coeffs)>1 and abs(coeffs[0])<1e-12: coeffs=coeffs[1:]
        roots=np.roots(coeffs) if len(coeffs)>1 else np.array([])
        cov=q[0]
        mags=np.abs(roots) if len(roots) else np.array([np.inf])
        rows.append((cov,name,E,roots,float(mags.min()),float(mags.max())))
    rows.sort(key=lambda r:-r[0])
    covmax=rows[0][0]
    print("  coverage q0 | name | min|root| | max|root| | roots(mod) | all on circle?")
    for cov,name,E,roots,mn,mx in rows[:10]:
        modstr=",".join(f"{abs(r):.3f}" for r in roots[:6])
        oncirc = (mx-mn)<0.05*max(mx,1e-9) if len(roots) else False
        tag=" <== COV-MAX" if cov==covmax else ""
        print(f"  q0={cov:.4f} {name:10s} min|r|={mn:.3f} max|r|={mx:.3f} mods=[{modstr}] circle?={oncirc}{tag}")
    # does the coverage-max (consec) have extremal root magnitude?
    consec_row=[r for r in rows if r[1]=="consec"][0]
    print(f"  --> consec coverage={consec_row[0]:.4f} (max? {consec_row[0]==covmax}), "
          f"min|root|={consec_row[4]:.4f}, max|root|={consec_row[5]:.4f}")
    return rows

if __name__=="__main__":
    sys.stdout.reconfigure(line_buffering=True)
    rng=random.Random(485)
    for n in (5,6,7):
        map1(n, sample=(None if n<=6 else 4000), rng=rng)
    for k in (8,9,10):
        map2(k, rng)
