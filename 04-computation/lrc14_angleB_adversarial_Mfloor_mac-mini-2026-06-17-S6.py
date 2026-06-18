from fractions import Fraction as F
import random
C=F(1,14)
def covering(S): return all(any(v%q==0 for v in S) for q in range(2,15))
def nrm(x):
    r=x-int(x); r=r+1 if r<0 else r; return r if r<=F(1,2) else 1-r
def g(S,t): return min(nrm(v*t) for v in S)
def cand(S):
    S=sorted(set(S)); Cc=set()
    for v in S:
        k=0
        while F(2*k+1,2*v)<=F(1,2): Cc.add(F(2*k+1,2*v)); k+=1
    for i in range(len(S)):
        for j in range(i+1,len(S)):
            for d in (S[i]+S[j],S[j]-S[i]):
                if d>0:
                    k=1
                    while F(k,d)<=F(1,2): Cc.add(F(k,d)); k+=1
    Cc.add(F(1,2)); return Cc
def Mexact(S): return max(g(S,t) for t in cand(S))
def gen_mixed(center, nsmall, jit):
    base=list(range(1,nsmall+1)); used=set(base); larges=[]
    needed=[q for q in range(2,15) if not any(b%q==0 for b in base)]
    for q in needed:
        k=round(center/q)+random.randint(-jit,jit); c=q*k
        while c in used or c in larges or c<=nsmall: k+=1; c=q*k
        larges.append(c); used.add(c)
    bl=sorted(set(base+larges)); hi=center+random.randint(0,30)
    S=list(bl)
    while len(S)<13:
        hi+=1
        if hi not in S: S.append(hi)
    S=sorted(set(S))
    return S[:13] if len(bl)>=13 else sorted(set(bl+[hi]))
random.seed(2024)
worstM=(F(99),None); tested=0; brk=0
for _ in range(6000):
    S=gen_mixed(random.randint(100,3000), random.choice([1,2,3,4,5,6,7]), random.choice([1,2,3,4]))
    S=sorted(set(S))
    if len(S)!=13 or not covering(S): continue
    tested+=1
    M=Mexact(S)
    if M<worstM[0]: worstM=(M,S)
    if M<C: brk+=1
print(f"ADVERSARIAL exact-M on {tested} clustered covering 13-sets")
print(f"  LRC breaks (M<1/14): {brk}")
print(f"  worst (smallest) M = {worstM[0]} = {float(worstM[0]):.6f}  at S={worstM[1]}")
print(f"  worst M >= 1/14? {worstM[0]>=C}")
