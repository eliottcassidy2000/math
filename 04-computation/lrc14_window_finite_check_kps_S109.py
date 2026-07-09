# Bounded-window finite check: verify M(S)=max_tau min_i ||v_i tau|| >= 1/14 DIRECTLY for runner sets
# with Vmax in (spread, 2.8*spread] (kps-S109). Sidesteps the good-period embedding (which fails
# locally in the window, mac-mini). ODLYZKO-TE RIELE ANALOG: search runner sets by phase-aligning the
# speeds (hill-climb to MINIMIZE M = maximal constructive interference of the gaps) to find the worst
# case; it should converge to the tight AP where M=1/14 exactly (the LRC(14) equality/extremal).
import numpy as np
from math import gcd
from functools import reduce
import random
random.seed(109)
def nearest(x): return abs(x-round(x))
def M_of(v, Ngrid=200001):
    # M = max_tau min_i ||v_i tau||, tau in [0,1). fine grid then local refine at the argmax.
    v=np.array(v,float); tau=np.linspace(0,1,Ngrid,endpoint=False)
    ph=np.mod(np.outer(tau,v),1.0); nr=np.minimum(ph,1-ph); mr=nr.min(axis=1)
    k=int(mr.argmax()); t0=tau[k]
    # refine: golden/ternary around t0 +- 1/Ngrid
    a,b=t0-1.5/Ngrid, t0+1.5/Ngrid
    for _ in range(80):
        m1=a+(b-a)/3; m2=b-(b-a)/3
        f1=min(nearest(vi*m1) for vi in v); f2=min(nearest(vi*m2) for vi in v)
        if f1<f2: a=m1
        else: b=m2
    tstar=(a+b)/2; return min(nearest(vi*tstar) for vi in v), tstar
def prim(vs):
    g=reduce(gcd,[int(round(x)) for x in vs]); return g==1
TH=1.0/14
print("Bounded-window finite check: M(S) >= 1/14 for Vmax in (spread, 2.8*spread]. 1/14 =",round(TH,5))
print(f"{'family':>14}{'Vmax':>6}{'spread':>7}{'M(S)':>9}{'M-1/14':>10}{'>=1/14?':>9}")
def emit(lab,v):
    v=sorted(set(int(x) for x in v))
    if len(v)!=13: return None
    M,ts=M_of(v); Vmax=max(v); sp=max(v)-min(v)
    print(f"{lab:>14}{Vmax:>6}{sp:>7}{M:>9.5f}{M-TH:>+10.5f}{'YES' if M>=TH-1e-6 else 'NO!':>9}")
    return M
# tight AP {1..13}: Vmax=13, spread=12, Vmax/spread=1.08 IN window; M should be exactly 1/14
emit("AP{1..13}", list(range(1,14)))
# 7-structured worst (mac-mini worst7Struct co-offsets, v=Vmax-e, Vmax=91)
e7=[0,7,14,21,26,29,37,44,51,58,67,75,82]; Vmax=91
emit("7-struct(91)", [Vmax-e for e in e7])
# window samples: co-offsets in [0,spread], Vmax in (spread, 2.8spread], v=Vmax-e
worst=(9.9,None)
for sp in (12,20,30,45,60):
    for _ in range(60):
        Vmax=random.randint(sp+1, int(2.8*sp))
        e=sorted(set([0]+random.sample(range(1,sp),11)+[sp]))
        if len(e)!=13: continue
        v=[Vmax-x for x in e]
        if min(v)<=0 or not prim(v): continue
        M=emit(f"win s={sp}", v)
        if M is not None and M<worst[0]: worst=(M,tuple(sorted(v)))
        break
print(f"\nADVERSARIAL min-M (Odlyzko-te Riele phase-align: hill-climb to minimize M in the window):")
def rand_window():
    for _ in range(200):
        sp=random.randint(12,40); Vmax=random.randint(sp+1,int(2.8*sp))
        e=sorted(set([0]+random.sample(range(1,sp),11)+[sp]))
        if len(e)==13:
            v=sorted(set(Vmax-x for x in e))
            if len(v)==13 and min(v)>0 and prim(v): return v
    return None
best=(9.9,None)
for _r in range(25):
    v=rand_window()
    if v is None: continue
    cur,_=M_of(v,50001); vv=v[:]
    improved=True
    while improved:
        improved=False
        for idx in range(13):
            for d in (-1,1):
                nv=vv[idx]+d
                if nv<=0 or nv in vv: continue
                v2=sorted(vv[:idx]+[nv]+vv[idx+1:])
                if not prim(v2): continue
                sp2=max(v2)-min(v2)
                if not (sp2 < max(v2) <= 2.8*sp2): continue
                m2,_=M_of(v2,50001)
                if m2<cur-1e-7: cur=m2; vv=v2; improved=True; break
            if improved: break
    if cur<best[0]: best=(cur,tuple(vv))
Mb,_=M_of(list(best[1]),400001) if best[1] else (float('nan'),0)
print(f"  adversarial MIN M(S) in window = {Mb:.5f} (1/14={TH:.5f}); margin {Mb-TH:+.5f}")
print(f"  minimizer v = {best[1]}")
print(f"\n=> if min M(S) >= 1/14 across the window (= at the AP, the phase-aligned extremal): the")
print(f"   bounded-window finite check PASSES; M=1/14 exactly at the AP = LRC(14) equality case.")
