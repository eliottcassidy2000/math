from itertools import permutations,product

def q4(size,Q):
    """Normalized by 2 h_A h_B."""
    return sum(prod(sum((t<i,j<t<k,t>l)) for t in range(size) if t not in Q)
               for i,j,k,l in permutations(Q) if j<k)

def prod(it):
    ans=1
    for a in it:ans*=a
    return ans

def bundled(size, swaps):
    return sum((-1)**sum(b)*q4(size,tuple(sorted(swaps[j][i] for j,i in enumerate(b))))
              for b in product((0,1),repeat=4))

if __name__=='__main__':
    for swaps in [((0,1),(2,3),(4,5),(6,7)),((0,1),(3,4),(6,7),(9,10)),
                  ((0,7),(1,6),(2,5),(3,4)),((0,2),(1,3),(4,6),(5,7))]:
        print(swaps,bundled(max(map(max,swaps))+1,swaps))
    for x in [(0,1,2,3),(0,2,4,6),(1,2,3,4)]:
        print(x,q4(8,x))

def kernel(a,b,c,d):
    p,q,r=b-a-1,c-b-1,d-c-1
    return 2**(p+r)*(4*3**q+3)+2*(p==0)*2**r+2*(r==0)*2**p+int(p==r==0)

def matchings(a):
    if not a:yield ()
    else:
        for j in range(1,len(a)):
            for s in matchings(a[1:j]+a[j+1:]):yield ((a[0],a[j]),)+s

def oracle():
    from itertools import combinations
    from random import Random
    n=10
    for q in combinations(range(n),4):
        assert q4(n,q)==kernel(*q)
    for u in combinations(range(n-1),4):
        if any(y-x<2 for x,y in zip(u,u[1:])):continue
        assert bundled(n,tuple((x,x+1) for x in u))==kernel(*(x-i for i,x in enumerate(u)))
    b=[bundled(8,s) for s in matchings(tuple(range(8)))]
    assert len(b)==105 and min(b)==8 and max(b)==184
    print('105 interleaved perfect-match controls: nonzero; min/max',min(b),max(b))
    rng=Random(432)
    for _ in range(40):
        n=10;a_len=2;b_len=2;m=n-a_len-b_len-2
        z=m+4
        # a,b arbitrary small transitive tiles; only separator contacts vary.
        A=list(range(a_len));V=list(range(a_len,a_len+m));B=list(range(a_len+m,a_len+m+b_len))
        base=[0]*n
        for i in range(z):
            for j in range(i+1,z):base[i]|=1<<j
        for u in A+B:
            for v in (z,z+1):
                if rng.randrange(2):base[u]|=1<<v
                else:base[v]|=1<<u
        def count(v,flip):
            a=list(base)
            a[z+flip]|=1<<(z+1-flip)
            for i,j in enumerate(V):
                for port in (z,z+1):
                    if v>>i&1:a[j]|=1<<port
                    else:a[port]|=1<<j
            N=1<<n
            d=[[0]*n for _ in range(N)]
            for j in range(n):d[1<<j][j]=1
            for s in range(1,N):
                for j,h in enumerate(d[s]):
                    if not h:continue
                    r=a[j]&~s
                    while r:
                        x=r&-r;r-=x;d[s|x][x.bit_length()-1]+=h
            return sum(d[-1])
        for flip in (0,1):
            diff=sum((-1)**x.bit_count()*count(x,flip) for x in range(1<<m))
            assert diff==2*kernel(0,1,2,3),(diff,flip)
    print('80 complete Boolean cubes checked against literal path DP')

if __name__=='__main__':
    oracle()
