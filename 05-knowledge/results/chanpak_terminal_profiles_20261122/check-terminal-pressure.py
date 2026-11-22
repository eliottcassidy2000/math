from fractions import Fraction

def Z(v):
    d,u,w=v
    return d+2*u,w,d+u

# Prove the rational upper bounds on the step ratios, without floats.
assert 3**24 < 2**23 * 5**24
assert 3**24 * 3**23 < 2**24 * 5**23
for i in range(3):
    e=tuple(int(i==k) for k in range(3))
    assert all(a>=2*b for a,b in zip(Z(Z(Z(e))),e))
base=10**6
rhs=base**24 * 3**23
def ceil_entry(c):
    low,high=0,base
    while high-low>1:
        k=(low+high)//2
        if k**24*c**23>=rhs: high=k
        else: low=k
    assert high**24*c**23>=rhs
    return high

a,b=(3,0,0),(2,0,1)
rows=[]
for j in range(1,21):
    a,b=Z(a),Z(b)
    vals=[]
    for m in range(1,21):
        c=min(sum(t)+(1<<m)*(sum(t)+t[1]) for t in (a,b))
        vals.append(ceil_entry(c))
    rows.append(sum(vals)+2*vals[-1])
bound=Fraction(2*sum(rows)+3*sum(rows[-3:]),2*base)
assert bound<Fraction(999,1000)
print("exact entire pressure upper bound:",bound,"< 999/1000")

def profile(A):
    n=len(A);N=1<<n
    p=[[0]*n for _ in range(N)]
    for j in range(n):p[1<<j][j]=1
    for s in range(1,N):
        for i,y in enumerate(p[s]):
            R=A[i]&~s
            while R and y:
                d=R&-R;R-=d;p[s|d][d.bit_length()-1]+=y
    H=[1]+list(map(sum,p[1:]))
    return H[-1],sum(H[j]*H[N-1-j] for j in range(N))
lanes=[(1,2)]+[profile(A) for A in ((2,8,3,5),(18,4,9,19,6),(2,16,11,19,37,15))]
assert lanes==[(1,2),(5,32),(11,100),(19,256)]
den,num=200,199
assert 3**den < 2**num*5**den
assert 3**den*3**num < 2**den*5**num
def round_lane(c):
    lo,hi=0,base
    goal=base**den*3**num
    while hi-lo>1:
        d=(lo+hi)//2
        if d**den*c**num>=goal: hi=d
        else: lo=d
    assert hi**den*c**num>=goal
    return hi
ra,rb=(3,0,0),(2,0,1)
rr=[]
for j in range(1,21):
    ra,rb=Z(ra),Z(rb)
    B=0
    for h,Q in lanes:
        vals=[round_lane(min(h*sum(t)+(Q<<k)*(sum(t)+t[1]) for t in (ra,rb))) for k in range(20)]
        B+=sum(vals)+2*vals[-1]
    rr.append(B)
bound=Fraction(2*sum(rr)+3*sum(rr[-3:]),2*base)
assert bound<1
print('four padded-lane pressure upper bound:',bound,'< 1')

def strong(A):
    n=len(A);full=(1<<n)-1
    for i in range(n):
        S=1<<i
        while True:
            new=S
            for j in range(n):
                if S>>j&1:new|=A[j]
            if new==S:break
            S=new
        if S!=full:return False
    return True
def rows(n,m):
    v=[0]*n
    for i in range(n):
        for j in range(i+1,n):
            if m&1:v[j]|=1<<i
            else:v[i]|=1<<j
            m>>=1
    return v
profiles=set()
for n in range(3,6):
    for m in range(1<<(n*(n-1)//2)):
        A=rows(n,m)
        if strong(A):profiles.add(profile(A))
assert profiles=={(3,12),(5,32),(9,84),(11,100),(13,104),(15,108),(15,120)}
assert 20**den<23**num and 60**den<70**num
ra,rb=(3,0,0),(2,0,1)
by_lane=[[],[],[]]
for j in range(1,21):
    ra,rb=Z(ra),Z(rb)
    for rec,(h,Q) in zip(by_lane,lanes[:3]):
        vals=[round_lane(min(h*sum(t)+(Q<<k)*(sum(t)+t[1]) for t in (ra,rb))) for k in range(20)]
        rec.append(sum(vals)+2*vals[-1])
B=[Fraction(2*sum(rr)+3*sum(rr[-3:]),2*base) for rr in by_lane]
bound=B[0]+Fraction(20,19)*B[1]+Fraction(60,59)*B[2]
assert bound<1
print('all ordinal gap parts through order five, prime-height pressure:',bound,'< 1')
