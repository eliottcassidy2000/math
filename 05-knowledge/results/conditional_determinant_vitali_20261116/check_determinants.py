from itertools import product

def rev(t,m):
    return int(f'{t:0{m}b}'[::-1], 2) if m else 0

def r(t,m):
    return (1<<m)-1-rev(t,m)

def tape(t,m):
    return tuple(3*int(c) for c in f'{t:0{m}b}')

def addresses(w):
    a,b=zip(*(divmod(x,2) for x in w))
    L=lambda bits: sum(e<<i for i,e in enumerate(bits[::-1]))
    x,y=L(a),L(b)
    return x,y,r(x,len(w)),r(y,len(w))

def gap(w):
    a,b=zip(*(divmod(z,2) for z in w))
    def F(j,i):return (1-a[i])*b[j]-(1-b[i])*a[j]
    return (b[0]-a[0]+b[-1]-a[-1]
            +sum((1<<(i-j-1))*F(j,i) for i in range(len(w)) for j in range(i))
            +sum(F(i+1,i) for i in range(len(w)-1)))

def fast(w):
    def step(old,a,b):
        v=[0]*12
        for f,l,m in product(range(2),range(2),range(3)):
            x=old[6*f+3*l+m]
            v[6*a+3*l+m]+=x
            if a or m:v[6*f+3*l+1+b]+=x
            if b or l:v[6*f+3+m]+=x
        return v
    x=y=[1]+[0]*11
    for z in w:
        a,b=divmod(z,2)
        x=step(x,a,b); y=step(y,b,a)
    return x[0]+x[3]+x[1]+x[4]+y[1]+y[4], y[0]+y[3]+y[1]+y[4]+x[1]+x[4]

p=(0,0,0,1,1,2,0)
assert addresses(p)==(2,12,95,103)
assert gap(p)==8
checked=0
for j in range(1,6):
    m=4*j+7
    seen=set()
    for eps in product(range(4),repeat=j-1):
        t=25*16**j+sum(e*16**i for i,e in enumerate(eps,1))
        T,S=5*t+2,4*t+1
        w=tape(T,m)+p+tape(r(S,m),m)
        assert w[:2]==(3,3) and w[-2:]==(0,0)
        a,b,c,e=addresses(w); D=b*c-a*e
        Z=sum((0,2,1,3)[v]*16**(j-1-i) for i,v in enumerate(eps,1))
        assert D==(339*4**m)//2-800*2**m-42*2**(m+7)*Z
        assert D>158*4**m and D not in seen
        seen.add(D)
        h,i=fast(w)
        assert 3**(len(w)-4)<=i<h<=2*3**len(w)
        assert h-i==2 and gap(w)==2
        checked+=1
    assert len(seen)==4**(j-1)
print(checked,'literal three-run replay pairs; formula and distinctness passed.')
for j in range(5,27):
    m=4*j+7
    inds=list(range(j-1,0,-4))
    prev=None
    for eps in product(range(2),repeat=len(inds)):
        t=25*16**j+sum(e*16**i for i,e in zip(inds,eps))
        T,S=5*t+2,4*t+1
        w=tape(T,m)+p+tape(r(S,m),m)
        hi,lo=fast(w)
        assert hi-lo==2
        if prev is not None:assert lo>prev
        prev=hi
print('Separated subfamily strictly ordered, with disjoint absolute intervals, at j=5,...,26.')
# Small-count application: abridge neither the paired tape search nor the literal witnesses.
assert 5**8>69**3
found=0
for n in range(7):
    for w in product(range(4),repeat=n):
        assert set(fast(w))!={67,69}
        found+=1
assert found==5461

def literal_count(rows):
    n=len(rows)
    dp=[[0]*n for _ in range(1<<n)]
    for i in range(n):dp[1<<i][i]=1
    for mask in range(1,1<<n):
        for j in range(n):
            for i in range(n):
                if not mask>>i&1 and rows[j]>>i&1:
                    dp[mask|1<<i][i]+=dp[mask][j]
    return sum(dp[-1])

def induced(rows,S):
    return [sum((rows[i]>>j&1)<<k for k,j in enumerate(S)) for i in S]

R=[86,124,24,113,96,69,4]
assert literal_count(R)==67
assert literal_count(induced(R,[0,1,3,4,6]))==3
for (a,b),h in [((2,5),69),((5,6),71)]:
    R[a]^=1<<b
    R[b]^=1<<a
    assert literal_count(R)==h
print('All 5461 small tapes exclude {67,69}; raw subset DP gives 67→69→71 off the mutual-arc model.')
