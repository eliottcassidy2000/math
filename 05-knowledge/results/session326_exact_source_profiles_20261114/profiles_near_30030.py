"""Complete H,Q source atlas at optional-slot budget <= floor(30047/99)."""
from itertools import permutations, product
from collections import defaultdict
from pathlib import Path
from json import loads
from inverse_carries import logical

def extend(A,cut):
    n=len(A)
    return tuple(r|(((cut>>i&1)^1)<<n) for i,r in enumerate(A))+(cut,)

def canon(A):
    groups=defaultdict(list)
    for i,r in enumerate(A): groups[r.bit_count()].append(i)
    best=None
    for runs in product(*(tuple(permutations(groups[d])) for d in sorted(groups))):
        p=tuple(x for run in runs for x in run)
        row=tuple(sum((A[i]>>j&1)<<k for k,j in enumerate(p)) for i in p)
        if best is None or row<best: best=row
    return best

def strong(A):
    for i in range(len(A)):
        rb=1<<i
        while True:
            b=rb
            for j,r in enumerate(A):
                if rb>>j&1: b|=r
            if b==rb: break
            rb=b
        if rb!=(1<<len(A))-1: return False
    return True

def hp_all(A):
    n=len(A); N=1<<n
    E=[[0]*n for _ in range(N)]; H=[1]+[0]*(N-1)
    for i in range(n): E[1<<i][i]=1
    for s in range(1,N):
        H[s]=sum(E[s])
        for i,e in enumerate(E[s]):
            if not e: continue
            W=A[i]&~s
            while W:
                b=W&-W; W-=b
                E[s|b][b.bit_length()-1]+=e
    return H

def profile(A):
    H=hp_all(A)
    return H[-1],sum(a*b for a,b in zip(H,H[::-1]))

def catalog(qmax):
    logical(qmax<=319) # orders 7 and greater are excluded by triangle lower bounds
    T={(0,)}; levels=[1]; prof={(1,2):(0,)}
    for n in range(2,7):
        T={canon(extend(A,c)) for A in T for c in range(1<<(n-1))}
        levels.append(len(T))
        for A in T:
            if strong(A):
                h,q=profile(A)
                if q<=qmax: prof.setdefault((h,q),A)
    logical(levels==[1,1,2,4,12,56])
    reached={(1,1):()}
    for b in range(1,qmax+1):
        for (h,q),L in tuple(reached.items()):
            if q!=b: continue
            for (a,c),X in prof.items():
                if q*c<=qmax: reached.setdefault((h*a,q*c),L+(X,))
    return reached

def local_parameters(A,S):
    n=len(A); M=(1<<n)-1
    def ends(R):
        E=[[0]*n for _ in range(M+1)]
        for i in range(n): E[1<<i][i]=1
        for q in range(1,M+1):
            for i,v in enumerate(E[q]):
                if v:
                    Z=R[i]&~q
                    while Z:
                        j=Z&-Z
                        Z-=j
                        E[q|j][j.bit_length()-1]+=v
        return E
    E=ends(A); P=ends(tuple(M^r^(1<<i) for i,r in enumerate(A)))
    H=[1]+[sum(e) for e in E[1:]]
    p=[sum(v for i,v in enumerate(e) if S>>i&1) for e in P]
    d=[sum(v for i,v in enumerate(e) if S>>i&1) for e in E]
    return H[M],sum(H[M^q]*p[q] for q in range(1,M+1)),sum(d[q]*p[M^q] for q in range(1,M))

def core_with_apex(I,m,cut):
    A=tuple(r|(((1<<m)-1)<<len(I)) for r in I)+tuple(
        ((1<<(len(I)+m))-1)^((1<<(len(I)+i+1))-1) for i in range(m))
    T=extend(A,cut)
    reach=[]
    for i in range(len(T)):
        r=1<<i
        while True:
            b=r
            for j,x in enumerate(T):
                if r>>j&1: b|=x
            if b==r: break
            r=b
        reach.append(r)
    C=tuple(i for i in range(len(T)) if reach[i]>>(len(T)-1)&1 and reach[-1]>>i&1)
    logical(len(I) in C)
    for i in set(range(len(T)))-set(C):
        logical(all(i==j or not ((reach[i]>>j&1) and (reach[j]>>i&1)) for j in range(len(T))))
    return tuple(sum((T[i]>>j&1)<<k for k,j in enumerate(C)) for i in C)

def main():
    bank=loads(Path(__file__).with_name('phase_table.json').read_text())
    P=catalog(30047//99)
    targets={30013,30017,30029,30031,30047}
    rows={}
    for t in targets:
        rows[t]=[(h,q,c) for h,q in P if q>1 and (t-99*q)%h==0 and (c:=(t-99*q)//h) in range(1,198,2)]
    print('atlas profiles (including empty):',len(P))
    print('exact relaxed target hits:',rows)
    logical(len(P)==35 and rows[30017]==[(23,296,31)] and bank['31']==0)
    logical(all(not rows[t] for t in (30013,30029,30031,30047)))
    logical({h for h,q in P if q==100}=={11} and {h for h,q in P if q==200}=={11})
    logical(profile((2,20,9,3,13))==(11,100))
    seeds=(
      (30013,(32,17,3,7,13,30,63),112,(37,470,104),6),
      (30017,(96,17,3,7,77,30,46),80,(115,468,50),6),
      (30029,(6,112,98,7,13,25,57),101,(151,947,426),5),
      (30031,(2,36,17,7,11,29,63),92,(29,472,206),6),
      (30047,(64,9,35,21,71,27,46),61,(113,953,562),5),
    )
    for t,A,cut,xyz,m in seeds:
        logical(local_parameters(A,cut)==xyz and xyz[1]*(1<<m)+xyz[0]-xyz[2]==t)
        T=core_with_apex(A,m,cut)
        logical(strong(T) and hp_all(T)[-1]==t)
    T=(462,480,378,114,419,385,176,12,200)
    logical(strong(T) and profile(T)==(1827,30000))
    print('five literal target repairs; strong source (1827,30000) verified')
if __name__=='__main__': main()
