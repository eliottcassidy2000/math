# Klapper-Goresky rational approximation algorithm: 2-adic complexity profile of a bit sequence.
import sys, random
def phi(h): return max(abs(h[0]),abs(h[1]))
def best_d(f,g):
    # odd d minimising phi(f + d g); candidates near the zeros of each coordinate
    cands=set()
    for (a,b) in ((f[0],g[0]),(f[1],g[1])):
        if b!=0:
            q=-a//b
            for dq in range(-2,3):
                cands.add(q+dq)
    cands={d if d%2 else d+1 for d in cands}|{d-2 for d in cands if d%2}
    best=None
    for d in cands:
        if d%2==0: continue
        h=(f[0]+d*g[0], f[1]+d*g[1]); v=phi(h)
        if best is None or v<best[0]: best=(v,d,h)
    return best[2]
def raa(bits, record_every=1):
    T=len(bits)
    k=next(i for i,b in enumerate(bits) if b)
    alpha=0
    for i in range(T):
        if bits[i]: alpha|=1<<i
    f=(0,2); g=(1<<k,1); prof=[]
    while k<T-1:
        k+=1
        mod=1<<(k+1)
        if (alpha*g[1]-g[0])%mod==0:
            f=(2*f[0],2*f[1])
        elif phi(g)<phi(f):
            h=best_d(f,g); f=(2*g[0],2*g[1]); g=h
        else:
            h=best_d(g,f); g=h; f=(2*f[0],2*f[1])
        prof.append((k+1, phi(g)))
    return g, prof
if __name__=="__main__":
    # validation: 2-adic expansion of -1/3 = 1/(1-4)*... ; of p/q rational -> complexity stops growing
    def expand(p,q,T):  # 2-adic digits of p/q (q odd)
        bits=[]; x=p
        for _ in range(T):
            b=x%2; bits.append(b); x=(x-b*q)//2
        return bits
    g,prof=raa(expand(-5,13,300)); print("rational -5/13 -> g =",g," log2 phi at end %.2f"%(prof[-1][1].bit_length()))
    random.seed(3); rb=[random.getrandbits(1) for _ in range(4000)]; rb[0]=1
    g,prof=raa(rb); import math
    print("random 4000 bits: log2 phi(g) =",prof[-1][1].bit_length()," (n/2 = 2000)")
