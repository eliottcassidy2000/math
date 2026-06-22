import random
def N(S,D):
    c=0
    for a in range(1,D):
        if all(min((s*a)%D,D-(s*a)%D)*14>=D for s in S): c+=1
    return c
# strong-component atoms: resonance graph G(S,D): s~s' if some small k,k' with k s + k' s' ~0 mod D (|k|,|k'|<=2)
def atoms(S,D):
    import sys
    n=len(S); par=list(range(n))
    def find(x):
        while par[x]!=x: par[x]=par[par[x]]; x=par[x]
        return x
    for i in range(n):
        for j in range(i+1,n):
            linked=False
            for k in range(-2,3):
                for kp in range(-2,3):
                    if k==0 and kp==0: continue
                    r=(k*S[i]+kp*S[j])%D
                    if min(r,D-r)*7<D:  # low-height relation (within 1/7)
                        linked=True; break
                if linked: break
            if linked:
                par[find(i)]=find(j)
    comps={}
    for i in range(n): comps.setdefault(find(i),[]).append(S[i])
    return list(comps.values())
primes=[83,89,97,101,103]
random.seed(3)
maxcover=0; worst=None
for _ in range(4000):
    S=random.sample(range(1,200),13)
    nz=[D for D in primes if N(S,D)==0]
    if len(nz)>maxcover: maxcover=len(nz); worst=(S,nz)
print("Multi-prime covering: max # of {83,89,97,101,103} with N=0 simultaneously, over 4000 13-sets:")
print("  max =",maxcover,"of",len(primes))
if worst: print("  worst S:",sorted(worst[0]),"covers (N=0) at",worst[1])
# atom structure of a bad case
S=list(range(1,12))+[13,84]
for D in [14,28,41,83]:
    print(f"  loosest at D={D}: N={N(S,D)}, #atoms(resonance)={len(atoms(S,D))}")
