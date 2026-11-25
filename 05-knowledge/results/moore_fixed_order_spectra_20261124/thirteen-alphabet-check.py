from itertools import combinations
from functools import lru_cache

def tournament(n,m):
    A=[0]*n
    for b,(i,j) in enumerate(combinations(range(n),2)):
        if m>>b&1:A[i]|=1<<j
        else:A[j]|=1<<i
    return A

def hp(A):
    n=len(A);d=[[0]*n for _ in range(1<<n)]
    for i in range(n):d[1<<i][i]=1
    for s in range(1,1<<n):
        for i,v in enumerate(d[s]):
            if v:
                for j in range(n):
                    if (A[i]&~s)>>j&1:d[s|1<<j][j]+=v
    return sum(d[-1])

def pack(A):
    n=len(A);allm=(1<<n)-1;c=[0]*(1<<n)
    def walk(r,i,m):
        if m.bit_count()>=3 and m.bit_count()%2 and A[i]>>r&1:c[m]+=1
        todo=A[i]&~m&~((1<<r)-1)
        while todo:
            b=todo&-todo; todo^=b
            walk(r,b.bit_length()-1,m|b)
    for r in range(n):walk(r,r,1<<r)
    @lru_cache(None)
    def f(m):
        if not m:return (1,)+(0,)*(n//3)
        v=m&-m
        z=list(f(m^v))
        for t,d in enumerate(c):
            if d and t&v and t&m==t:
                g=f(m^t)
                for k in range(1,len(z)):z[k]+=d*g[k-1]
        return tuple(z)
    return f(allm)

def inflate(A,x):
    n=len(A);B=A[:];B.append(B[x])
    for i in range(n):
        if B[i]>>x&1:B[i]|=1<<n
    B[x]|=1<<n
    return B

def seam(n,m,x,i,j,want):
    A=tournament(n,m)
    assert min(A[x].bit_count(),n-1-A[x].bit_count())==1
    B=A[:];B[i]^=1<<j;B[j]^=1<<i
    assert [hp(A),hp(B),hp(inflate(A,x)),hp(inflate(B,x))]==want

def invariant_control():
    # two moving 5-orbits and one fixed site, with different cross voltages.
    A=[0]*11
    for c in (0,5):
        for i in range(5):
            for d in (1,2):A[c+i]|=1<<(c+(i+d)%5)
    for i in range(5):
        for j in range(5):
            if (j-i)%5 in (0,1):A[i]|=1<<(5+j)
            else:A[5+j]|=1<<i
        A[i]|=1<<10
        A[10]|=1<<(5+i)
    r=pack(A)
    for a in range(5):
        assert sum(v*a**j for j,v in enumerate(r))%5==pow(1+2*a,2,5)

def cyclic_models():
    from itertools import permutations
    for p in (7,11,13):
        D={i*i%p for i in range(1,p)} if p%4==3 else set(range(1,7))
        A=[sum(1<<((i+j)%p) for j in D) for i in range(p)]
        assert hp(A)=={7:189,11:95095,13:3711175}[p]
        if p==7:
            good=[v for v in permutations(range(p)) if
                  all(A[i]>>j&1 for i,j in zip(v,v[1:]))]
            assert sum(not(A[v[-1]]>>v[0]&1) for v in good)==21
    assert 95095//55==1729 and 3711175//13==285475
    assert 1**3+12**3==9**3+10**3==7*13*19==1729

if __name__=='__main__':
    p=13
    assert len({pow(2,j,p) for j in range(p-1)})==12 and pow(2,12,p*p)==40
    for k in (1,2,3):
        q=p**k
        r={b for j in range(12*13**(k-1)) for b in
           ((12*pow(2,j,q)-1)%q,(54*pow(2,j,q)-3)%q)}
        assert len(r)==q
    seam(6,12,1,2,4,[11,13,23,25])
    seam(8,75588141,1,2,5,[53,51,107,105])
    vals=[]
    for m in (815318702,815320750):
        A=tournament(9,m)
        r=pack(A)
        v=hp(A)
        assert v==sum(a*2**i for i,a in enumerate(r))
        vals.append((m,v,r))
    assert [x[1] for x in vals]==[1729,1731]
    assert 815318702^815320750==1<<11
    invariant_control()
    cyclic_models()
    print('exact taxicab seam and cycle-packing coefficients:',vals)
