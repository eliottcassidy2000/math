from fractions import Fraction as F
import functools, math
import numpy as np
from scipy.optimize import milp, LinearConstraint, Bounds
print=functools.partial(print,flush=True)
def norm(x):
    f=x-int(x); 
    if f<0: f+=1
    return min(f,1-f)
def breakpoints(S):
    S=sorted(set(S)); cand=set()
    for i in range(len(S)):
        for j in range(len(S)):
            for d in (S[i]-S[j],S[i]+S[j]):
                if d>0:
                    for k in range(1,d): cand.add(F(k,d))
    return cand
def lonely_witnesses(S, r):
    # all breakpoints where every speed is >= r from 0 (candidate fails to cover at radius <r)
    out=[]; 
    for t in breakpoints(S):
        if all(norm(v*t)>=r for v in S): out.append(t)
    return out
def beater_multicut(n,V,r,max_rounds=60):
    speeds=list(range(1,V+1)); P=len(speeds)
    rows=[[1.0]*P]; lb=[float(n-1)]; ub=[float(n-1)]
    for q in range(2,n+1):
        rows.append([1.0 if v%q==0 else 0.0 for v in speeds]); lb.append(1.0); ub.append(float(P))
    seen=set()
    for rd in range(max_rounds):
        A=np.array(rows)
        res=milp(c=np.zeros(P),constraints=LinearConstraint(A,lb,ub),integrality=np.ones(P),
                 bounds=Bounds(0,1),options={"time_limit":40})
        if not(res.success and res.x is not None): return None, rd, len(rows)
        S=[speeds[i] for i in range(P) if res.x[i]>0.5]
        w=lonely_witnesses(S,r)
        if not w:   # candidate has M<r everywhere -> beater
            return S, rd, len(rows)
        added=0
        for t in w:
            key=t
            if key in seen: continue
            seen.add(key)
            rows.append([1.0 if norm(v*t)<r else 0.0 for v in speeds]); lb.append(1.0); ub.append(float(P)); added+=1
        if added==0:  # all witnesses already cut but still lonely -> infeasible-ish; force re-solve won't help
            return None, rd, len(rows)
    return "TIMEOUT", max_rounds, len(rows)
for n in [13,14]:
    V=n*(n-1); thr=F(n,n*n-n+1)
    res,rounds,ncon=beater_multicut(n,V,thr)
    if res is None:
        print(f"n={n} V={V}: NO beater M<{thr} (infeasible, {rounds} rounds, {ncon} constraints) -> covmin={thr} speeds<=n(n-1) RIGOROUS")
    elif res=="TIMEOUT":
        print(f"n={n} V={V}: inconclusive ({rounds} rounds)")
    else:
        print(f"n={n} V={V}: BEATER {res} !!")
