from collections import Counter, defaultdict
from fractions import Fraction
from functools import lru_cache
from itertools import combinations, permutations


def decode(n, code):
    edge = [[False]*n for _ in range(n)]
    for k, (a, b) in enumerate(combinations(range(n), 2)):
        edge[a][b] = bool(code >> k & 1)
        edge[b][a] = not edge[a][b]
    return edge


def profile(A):
    n=len(A)
    full=(1<<n)-1
    end=[[0]*n for _ in range(full+1)]
    H=[0]*(full+1)
    H[0]=1
    for mask in range(1,full+1):
        for v in range(n):
            if mask >> v & 1:
                prev=mask ^ (1<<v)
                end[mask][v]=(sum(end[prev][u] for u in range(n) if A[u][v])
                              if prev else 1)
        H[mask]=sum(end[mask])
    return H[full], sum(H[mask]*H[full^mask] for mask in range(full+1))


def pack_polynomial(A):
    n=len(A)
    full=(1<<n)-1
    cycle=Counter()
    for d in range(3,n+1,2):
        for v in combinations(range(n),d):
            for sigma in permutations(v[1:]):
                t=(v[0],)+sigma
                if all(A[a][b] for a,b in zip(t,t[1:]+t[:1])):
                    cycle[sum(1<<j for j in t)]+=1

    @lru_cache(None)
    def term(mask):
        if not mask:
            return {0:1}
        low=mask & -mask
        out=Counter({k+1:v for k,v in term(mask ^ low).items()})
        for c, count in cycle.items():
            if c & low and c & mask == c:
                for k,v in term(mask ^ c).items():
                    out[k+1]+=2*count*v
        return out

    return term(full),cycle


def tick(condition, reason):
    if not condition:
        raise RuntimeError(reason)


def check(A):
    n=len(A)
    h,q=profile(A)
    f,c=pack_polynomial(A)
    tick(h==sum(f.values()) and q==sum(v*2**k for k,v in f.items()), "cycle packing")
    tau=sum(v for mask,v in c.items() if mask.bit_count()==3)
    value=Fraction(q,2**n)
    tick(1+Fraction(tau,2)<=value<=1+Fraction(tau,2)+Fraction(h-1-2*tau,16)
         <=Fraction(h+3,4), "envelope")
    tick(q >= (n+1)*h and q%3==((-1)**n*h)%3, "split or mod-3")
    if n>=3:
        c_top=sum(v for mask,v in c.items() if mask.bit_count()==n-(n%2==0))
        v=4 if n%2 else 8
        tick(q%(2*v)==v*(c_top%2), "cycle parity mod 2-power")
    return h,q


def strong(A):
    verts=set(range(len(A)))
    for a in (A, list(zip(*A))):
        reach={0}
        while new:= {v for u in reach for v in verts if a[u][v]}-reach:
            reach|=new
        if reach!=verts:
            return False
    return True


def main():
    profile_sets={}
    atoms={}
    total=0
    for n in range(1,6):
        prof=set()
        core=set()
        for c in range(1<<(n*(n-1)//2)):
            a=decode(n,c)
            h,q=check(a)
            prof.add((h,Fraction(q,2**n)))
            if n>1 and strong(a):
                core.add((h,Fraction(q,2**n)))
            total+=1
        profile_sets[n]=prof
        atoms[n]=core
    families={0:{(1,Fraction(1))}}
    for n in range(1,6):
        families[n]=set(families[n-1])
        for j in range(2,n+1):
            families[n].update((h*a,r*b) for h,r in atoms[j] for a,b in families[n-j])
        tick(families[n]==profile_sets[n], f"ordinal profile reconstruction n={n}")
    for n in (5,6):
        w=[]
        for c in (40,76):
            A=decode(5,c)
            if n==6:
                A=[x+[True] for x in A]+[[False]*6]
            w.append(check(A))
        tick(w == [(15,108*2**(n-5)),(15,120*2**(n-5))], f"sharp witnesses n={n}")
    A=decode(4,17)
    B=decode(5,548)
    AB=[row+[True]*5 for row in A]+[[False]*4+row for row in B]
    tick(profile(A)==(5,32) and profile(B)==(11,100) and profile(AB)==(55,3200), "cross-atom conflict")
    tick((3,Fraction(36,32)) not in profile_sets[5], "naive lattice survivor")
    D,E=decode(8,10760449),decode(8,18975246)
    tick(strong(D) and profile(D)==(223,4096), "223 dyadic toggle")
    tick(strong(E) and profile(E)==(599,8192), "599 dyadic toggle")
    F=decode(8,27280129)
    tick(strong(F) and profile(F)==(217,4096), "217 dyadic toggle")
    import math
    for r in range(3,7):
        for z in range(1,4):
            for v in range(1,4):
                mod=2**r*3**z*5**v
                jm=2**(r-2)*3**(z-1)*5**(v-1)
                period=2*3**(z-1)
                t=period*((r+period-1)//period)
                residues=set()
                for d in range(4):
                    for e in range(2):
                        for f in range(2):
                            W=6*2**(12*d+12*e+13*f+t)
                            U=217**d*223**e*599**f
                            for j in range((v+1)//2,(v+1)//2+jm):
                                val=(W*pow(100,j,mod)-U*pow(11,j,mod))%mod
                                residues.add(val)
                tick(residues=={a for a in range(mod) if math.gcd(a,mod)==1},
                     "joint residue decoder")
    for j in range(1,14):
        for d in range(4):
            for e in range(2):
                for f in range(2):
                    U=217**d*223**e*599**f
                    val=6*2**(12*d+12*e+13*f+9)*100**j-U*11**j
                    tick(math.gcd(val,2*3*5*11*U)==1, "selected-factor compatibility")
    print(f"PASS: {total} labelled profiles, 5 reconstructed layers, sharp gap pairs, three toggles, 36 joint residue systems.")

if __name__=="__main__":
    main()
