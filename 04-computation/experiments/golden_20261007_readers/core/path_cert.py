"""Follow one infinite parity word (the parity vector of a 2-adic integer y, here a negative-cycle point or -1) and
report the first depth s at which its class mod 2^s is certified by descent or by a backward-tree branch certificate
(ratio 3^{a-b} 2^{i-s} < 1, b <= a_s), using exact x_s mod 3^{a_s}."""
import sys
sys.setrecursionlimit(10000)
def T(x): return x//2 if x%2==0 else (3*x+1)//2
def word_of(y, L):
    w=[]; x=y
    for _ in range(L): w.append(x&1); x=T(x)
    return w
def lt1(e2,e3):  # 2^e2 3^e3 < 1
    # compare exactly using integers
    if e2<=0 and e3<=0: return not (e2==0 and e3==0)
    if e2>=0 and e3>=0: return False
    if e2<0: return 3**e3 < 2**(-e2)
    return 2**e2 < 3**(-e3)
def bsearch(z,p,e2,e3,budget):
    if lt1(e2,e3): return True
    if not lt1(e2+p,e3-p): return False
    if p==0 or budget[0]<=0: return False
    budget[0]-=1
    r=z%3
    if r==0: return False
    if r==2 and bsearch(((2*z-1)//3)%3**(p-1),p-1,e2+1,e3-1,budget): return True
    return bsearch((2*z)%3**p,p,e2+1,e3,budget)
def first_cert(w):
    a=0; c=0
    for s in range(1,len(w)+1):
        b=w[s-1]
        if b: c=3*c+2**(s-1); a+=1
        if 3**a < 2**s: return s,'descent'
        mod=3**a
        z=(c*pow(2,-s,mod))%mod if a>0 else 0
        budget=[2_000_000]
        if bsearch(z,a,-s,a,budget): return s,'branch'
        if budget[0]<=0: return s,'budget-exhausted'
    return None
for name,y in [('-1',-1),('-5',-5),('-7',-7),('-10',-10),('-17',-17),('-25',-25),('-37',-37),('-55',-55),('-41',-41),('-61',-61),('-91',-91),('-136',-136)]:
    w=word_of(y,150)
    print(name, ''.join(map(str,w[:33])), '...', first_cert(w), flush=True)
