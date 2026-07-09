import random
from math import gcd
def is_cov(S): return all(any(s%q==0 for s in S) for q in range(2,15))
def prim(S):
    g=0
    for s in S: g=gcd(g,s)
    return g==1
def margin(S):
    Vmax=max(S); E=sorted({Vmax-v for v in S}); best=0
    for j in range(1,Vmax):
        r=sorted({(e*j)%Vmax for e in E})
        g=max([r[i+1]-r[i] for i in range(len(r)-1)]+[Vmax-r[-1]+r[0]])
        m=7*g/Vmax
        if m>best: best=m
    return best
random.seed(2061)
print("ADVERSARIAL hill-climb: MINIMIZE the good-period margin over primitive covering 13-sets.")
print("(margin<=1 would mean a covering set with NO strict good period -- a real gap.)\n")
for CAP in [20, 30, 50, 100]:
    gbest=(9e9,None)
    for restart in range(10):
        # seed
        S=None
        for _ in range(20000):
            T=sorted(random.sample(range(1,CAP+1),13))
            if prim(T) and is_cov(T): S=T; break
        if S is None: continue
        m=margin(S); improved=True; it=0
        while improved and it<200:
            improved=False; it+=1
            for idx in range(13):
                for _ in range(20):
                    T=list(S); T[idx]=random.randint(1,CAP); T=sorted(set(T))
                    if len(T)!=13 or not prim(T) or not is_cov(T): continue
                    mt=margin(T)
                    if mt<m: S,m=T,mt; improved=True; break
        if m<gbest[0]: gbest=(m,list(S))
    m,S=gbest
    flag = "  <-- NO GOOD PERIOD (real gap!)" if m<=1 else ""
    print(f"  CAP={CAP:>3}: min margin found = {m:.4f}{flag}")
    print(f"           S = {S}")
print("\nmargin > 1 everywhere => the covering branch (the ONLY branch LRC(14) needs) has strict good periods.")
