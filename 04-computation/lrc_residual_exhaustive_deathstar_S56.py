# death-star-S56: EXHAUSTIVE enumeration of the boundary residual to close max<=34.
# Residual = valid non-AP core (covers 2..12, misses 13,14, max<=34) failing BOTH soft Weyl (C>464mu)
# and stability (delta<=max/2366). For max<=34, delta<=34/2366=0.0144 => M(W)<=0.0913 (near-tight).
from fractions import Fraction as F
from math import gcd
from itertools import combinations
import sys
POOL=[v for v in range(1,35) if v%13 and v%14]   # {1..34}\{13,14,26,28}, 30 elts
NEAR=F(1,13)+F(34,2366)  # M(W) threshold for the residual (0.0913)
def covers(fam): return all(any(v%q==0 for v in fam) for q in range(2,13))
def M_and_args(fam):  # returns (M, list of maximizer (a,q)) but bail if M>NEAR (return None)
    Q=2*max(fam)+2; best=F(0); args=[]
    for q in range(2,Q+1):
        for a in range(1,q):
            if gcd(a,q)!=1: continue
            m=min(min((v*a)%q,q-(v*a)%q) for v in fam)
            fr=F(m,q)
            if fr>NEAR: return None,None
            if fr>best: best=fr; args=[(a,q)]
            elif fr==best: args.append((a,q))
    return best,args
def C_mu(W,lvl=F(1,13)):
    ivs=[]
    for w in W:
        for k in range(0,w+1):
            lo=(F(k)-lvl)/w; hi=(F(k)+lvl)/w
            ivs.append((max(lo,F(0)),min(hi,F(1))))
    ivs=[(a,b) for a,b in ivs if b>a]; ivs.sort()
    mg=[]
    for a,b in ivs:
        if mg and a<=mg[-1][1]: mg[-1]=(mg[-1][0],max(mg[-1][1],b))
        else: mg.append((a,b))
    comps=[]; prev=F(0)
    for a,b in mg:
        if a>prev: comps.append((prev,a))
        prev=max(prev,b)
    if prev<F(1): comps.append((prev,F(1)))
    return len(comps), sum(float(b-a) for a,b in comps)
def is_AP(W): W=sorted(W); d=W[0]; return W==[d*i for i in range(1,13)]
def safe_at(vmax,args):  # is vmax safe (||vmax*a/q||>=1/13) at SOME maximizer of W?
    for (a,q) in args:
        m=(vmax*a)%q
        if 13*min(m,q-m)>=q: return True
    return False
residual=[]; counter=[]; near_cnt=0; done=0
for W in combinations(POOL,12):
    done+=1
    if done%5000000==0: print(f"  ...{done//1000000}M subsets, near-tight={near_cnt}, residual={len(residual)}",flush=True)
    if not covers(W): continue
    M,args=M_and_args(W)
    if M is None: continue          # M>0.0913 => stability handles, skip
    near_cnt+=1
    if is_AP(list(W)): continue
    C,mu=C_mu(list(W))
    delta=M-F(1,13)
    softweyl = (C<=464.4*mu); stability = (delta>F(max(W),2366))
    if softweyl or stability: continue   # covered by a rigorous lens
    # RESIDUAL: finite check via safe-at-maximizer for all candidates
    ub=int(max(W)/(13*delta)) if delta>0 else 0
    cands=[k for k in range(182,ub+1,182) if k>max(W)]
    allsafe=all(safe_at(vm,args) for vm in cands)
    residual.append((W,C,mu,float(delta)))
    if not allsafe: counter.append((W,cands))
print(f"DONE: {done} subsets, near-tight covering non-AP cores={near_cnt}",flush=True)
print(f"BOUNDARY RESIDUAL (fail both lenses): {len(residual)}",flush=True)
print(f"COUNTEREXAMPLES (a candidate not safe at maximizer): {len(counter)}",flush=True)
for w in residual[:20]: print("  residual:",list(w[0]),"C=%d mu=%.4f delta=%.5f"%(w[1],w[2],w[3]),flush=True)
if not counter: print("=> EVERY boundary residual core eliminated => max<=34 CLOSED (soft Weyl + stability + finite check).",flush=True)
else:
    for w in counter[:5]: print("  *** COUNTEREXAMPLE:",list(w[0]),w[1])
