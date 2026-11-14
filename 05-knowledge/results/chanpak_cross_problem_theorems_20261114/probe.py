from itertools import combinations, permutations
from math import gcd

def profile(a):
    n=len(a); dp=[{} for _ in range(1<<n)]
    for i in range(n): dp[1<<i][i]=1
    for m in range(1,1<<n):
        for i,j in dp[m].items():
            for k in range(n):
                if not (m>>k&1) and a[k][i]:
                    q=dp[m|1<<k];q[k]=q.get(k,0)+j
    h=[sum(q.values()) for q in dp];h[0]=1
    S=[sum(q.get(i,0)*h[((1<<n)-1)^m] for m,q in enumerate(dp))
       for i in range(n)]
    Q=sum(h[m]*h[((1<<n)-1)^m] for m in range(1<<n))
    return h[-1],Q,S

if __name__=='__main__':
    examples={}
    for n in range(3,6):
        es=list(combinations(range(n),2))
        for m in range(1<<len(es)):
            a=[[0]*n for _ in range(n)]
            for j,(u,v) in enumerate(es):
                a[u][v]=m>>j&1;a[v][u]=1-a[u][v]
            h,Q,S=profile(a)
            if h==1:continue
            ds={gcd(h,Q+sum(S[i] for i in range(n) if cut>>i&1)) for cut in range(1<<n)}
            if len(ds)>1 and h not in examples: examples[h]=(Q,S,tuple(sum(x<<i for i,x in enumerate(r)) for r in a))
    print(examples)

def seed(rows, mask):
    n=len(rows)
    spine=[[0]*(n+2) for _ in range(n+2)]
    for i in range(1,n+1):
        for j in range(1,n+1): spine[i][j]=rows[i-1]>>(j-1)&1
    for j in range(1,n+2): spine[0][j]=1
    for i in range(1,n+1): spine[i][n+1]=1
    cut=1|(mask<<1)
    g=[r+[int(not(cut>>i&1))] for i,r in enumerate(spine)]
    g.append([cut>>i&1 for i in range(n+2)]+[0])
    h,_,weights=profile(spine)
    N,*_=profile(g)
    return h,sum(s for i,s in enumerate(weights) if cut>>i&1),N

if __name__=='__main__':
    print("H=5 seeds:",[(u,seed((8,1,3,6),u)) for u in (0,10)])
    # Exhaust the finite orbit of two actual high append blocks, after
    # binary targeting. Enumerates states, not tournaments of each size.
    for rows, masks in [((8,1,3,6),(0,10)),((2,4,1),(0,1,7)),((16,1,3,7,14),(0,1,2,16))]:
        for mask in masks:
            h,S,N=seed(rows,mask)
            for m in range(1,140,2):
                start=(S%m,N%m); queue=[start]; seen={start}
                for p,y in queue:
                    for b in (1,(m+1)//2):
                        new=((2*b*(p+h)-2*h)%m,(y+b*(p+h))%m)
                        if new not in seen: seen.add(new);queue.append(new)
                d=gcd(h,S,m)
                assert {y for p,y in seen}=={y for y in range(m) if (y-N)%d==0}
    print("scalar congruence orbit replay passed")

def oddpack(a):
    n=len(a)
    cycles=[]
    for s in range(3,n+1,2):
        for sub in combinations(range(n),s):
            v,*rest=sub
            for perm in permutations(rest):
                path=(v,*perm,v)
                if all(a[i][j] for i,j in zip(path,path[1:])):
                    cycles.append(sum(1<<i for i in sub))
    coeff={(0,0):1}
    for c in cycles:
        for (mask,k), multiplicity in tuple(coeff.items()):
            if not (mask&c):
                key=mask|c,k+1
                coeff[key]=coeff.get(key,0)+multiplicity
    return coeff

def wrap_odds(a,m,r=1):
    n=len(a); nb=n*r+m
    b,d=nb,nb+1
    W=[[0]*(nb+2) for _ in range(nb+2)]
    for i,j in combinations(range(nb),2):
        W[i][j]=a[i%n][j%n] if j<n*r and i//n==j//n else 1
        W[j][i]=1-W[i][j]
    for i in range(nb): W[d][i]=W[i][b]=1
    W[b][d]=1
    f=oddpack(W)
    return [sum(v for (_,k),v in f.items() if k==i) for i in range(nb//3+2)]

def formula(a,m,r=1):
    from math import comb
    f=oddpack(a); n=len(a)
    A=[sum(v for (_,k),v in f.items() if k==i) for i in range(n//3+1)]
    K=[sum(v*2**(n-mask.bit_count()+k-i)*comb(k,i)
           for (mask,k),v in f.items() if k>=i) for i in range(n//3+1)]
    def pow_poly(P):
        coeff=[1]
        for _ in range(r):
            new=[0]*(len(coeff)+len(P)-1)
            for i,x in enumerate(coeff):
                for j,y in enumerate(P):new[i+j]+=x*y
            coeff=new
        return coeff
    x=pow_poly(A);y=pow_poly(K)
    return x+[0], [0]+[2**(m-1)*c for c in y]

if __name__=='__main__':
    controls=0
    def check_poly(a,m,r=1):
        lhs=wrap_odds(a,m,r)
        x,y=formula(a,m,r)
        rhs=[c+d for c,d in zip(x,y)]
        while len(lhs)>len(rhs) and lhs[-1]==0:lhs.pop()
        assert lhs==rhs,(a,m,r,lhs,rhs)
    for n in range(1,5):
        es=list(combinations(range(n),2))
        for m in range(1<<len(es)):
            a=[[0]*n for _ in range(n)]
            for j,(u,v) in enumerate(es):
                a[u][v]=m>>j&1;a[v][u]=1-a[u][v]
            for tail in (1,2):
                check_poly(a,tail)
                controls+=1
    check_poly([[0,1,0],[0,0,1],[1,0,0]],1,2)
    print(f'{controls+1} literal odd-cycle pack controls passed')

if __name__=='__main__':
    for m in range(1,7):
        n=m+2;b,a=m,m+1
        W=[[0]*n for _ in range(n)]
        W[b][a]=1
        for i in range(m):
            W[a][i]=W[i][b]=1
            for j in range(i+1,m):W[i][j]=1
        h,q,_=profile(W)
        assert (h,q)==(1+2**m,2**(m+2)+2*(3**m-1))
    print("six direct Fermat-rank profile controls passed")
