from math import isqrt

def hs():
    A=[[0]*50 for _ in range(50)]
    def edge(i,j): A[i][j]=A[j][i]=1
    for i in range(5):
        for x in range(5):
            edge(5*i+x,5*i+(x+1)%5)
            edge(25+5*i+x,25+5*i+(x+2)%5)
            for j in range(5): edge(5*i+x,25+5*j+(x+i*j)%5)
    assert all(sum(a)==7 for a in A)
    assert all(sum(A[i][k]*A[k][j] for k in range(50))+A[i][j]==(7 if i==j else 1)
               for i in range(50) for j in range(50))
    return A

def rank_mod(M,p=101):
    A=[[x%p for x in row] for row in M]
    N=len(A[0]);j=0
    for c in range(N):
        q=next((i for i in range(j,len(A)) if A[i][c]),None)
        if q is None: continue
        A[j],A[q]=A[q],A[j]
        inv=pow(A[j][c],-1,p)
        A[j]=[(x*inv)%p for x in A[j]]
        for i in range(j+1,len(A)):
            z=A[i][c]
            if z: A[i]=[(a-z*b)%p for a,b in zip(A[i],A[j])]
        j+=1
        if j==len(A): break
    return j

def fixture():
    A=hs(); root=0; hub=1
    reps=[i for i,v in enumerate(A[root]) if v and i!=hub]
    labs=[i for i,v in enumerate(A[hub]) if v and i!=root]
    grid=[[next(v for v in range(50) if A[i][v] and A[x][v]) for x in labs] for i in reps]
    assert len({v for row in grid for v in row})==36
    cells=[3,4,6,11,12,13,14,15,17,18,20,21,23,24,26,27,28,29,30,33,34,35]
    F=[divmod(k,6) for k in cells]
    assert len({i for i,j in F})==len({j for i,j in F})==6
    seen={0}
    while True:
        ext={i for i,j in F if j+6 in seen}|{j+6 for i,j in F if i in seen}
        if ext<=seen: break
        seen |= ext
    assert len(seen)==12
    O=sorted(set(range(50))-{grid[i][j] for i,j in F})
    GP=[[5*(A[i][j]+3*(i==j))-1 for j in O] for i in O]
    GM=[[10*(-A[i][j]+2*(i==j))+1 for j in O] for i in O]
    # Integer multiples of the two projectors; modular ranks lower-bound rational ranks.
    return rank_mod(GP),rank_mod(GM),len(cells)-12+1

def h(v):
    if not v:return 1
    u=v//2
    return h(u) if v&1 else h(u)+u+(u&1)

def q(v):
    if not v:return 0
    u=v//2
    return q(u)+(v&1)*((u+1)//2)

def rev(v,L):
    o=0
    for i in range(L):
        o=2*o+(v&1)
        v//=2
    return o

def prime(v):
    return v>=2 and (v==2 or (v%2==1 and all(v%d for d in range(3,isqrt(v)+1,2))))

def counts(B):
    v=B//2
    f=1662*v-58*q(v)+427
    return f,f+2

def test():
    r=fixture()
    assert r==(26,21,11),r
    print('HS frozen ranks and residual cycle rank:',r)
    for B,n in [(2684355418,2203103445089),(2684355568,2225910104699)]:
        a,b=counts(B)
        assert (a,b)==(n,n+2) and prime(a) and prime(b)
        assert B.bit_length()+9==41
        print('prime seam:',B,a,b)
    for v in range(1,2**16):
        L=v.bit_length()
        assert 8*abs((q(v)<<(L+1))-v*rev(v,L))<=3*L*(1<<(L+1))
        assert h(v)==v-2*q(v)+1-(v&1)
    for d in range(1,513):
        for a in range(64):
            for c in (2*((a+1)//2+d)-1,2*((a+1)//2+d)):
                low=h(a)+h(c)-1+2*a*c+2*((a+1)//2)
                assert low>=h(d-1)
                assert (low==h(d-1))==(a==0 and c==2*d-1)
    print('digit error and one-apex optimum checks passed')

if __name__=='__main__':test()
