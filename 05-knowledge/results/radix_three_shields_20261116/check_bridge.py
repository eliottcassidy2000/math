"""Run python3 check_bridge.py. Verify physical fixtures independently of search."""
from itertools import permutations

LEFT = (
  (0,(94,60,56,112,96,65,6)),
  (3,(96,65,3,39,15,22,60)),
  (4,(72,97,3,6,15,29,60)),
  (2,(30,60,120,80,96,73,3)),
)
RIGHT_START = ((0,(2,4,1)), (3,(8,1,3,6)))

def require(ok):
    if not ok: raise RuntimeError("invalid certificate")

def test(r):
    n=len(r)
    require(all(r[i] < 1<<n and not (r[i] >> i & 1) for i in range(n)))
    require(all((r[i]>>j&1) ^ (r[j]>>i&1) for i in range(n) for j in range(i)))
    for i in range(n):
        seen={i}
        for _ in range(n):
            seen |= {k for j in seen for k in range(n) if r[j]>>k&1}
        require(len(seen)==n)

def profile(r,v):
    test(r)
    n=len(r)
    a=[0,0]
    h=0
    # Each permutation and each cut encodes the ordered root path followed by
    # its ordered complement path, with no condition on the two seam vertices.
    for p in permutations(range(n)):
        missing=[i+1 for i in range(n-1) if not r[p[i]]>>p[i+1]&1]
        if not missing: h+=1
        if p[0]==v and len(missing)<=1:
            for k in missing or range(1,n+1): a[k%2]+=1
    return h,*a

def count(r):
    test(r)
    n=len(r); f=[[0]*n for _ in range(1<<n)]
    for s in range(1,1<<n):
        for j in range(n):
            if s>>j&1:
                f[s][j]=(1 if s==1<<j else sum(f[s^(1<<j)][i] 
                       for i in range(n) if r[i]>>j&1))
    return sum(f[-1])

def join(a,x,b,y):
    n,m=len(a),len(b)
    c=[row | (((1<<m)-1)<<n) for row in a]
    c.extend(sum((b[j]>>i&1)<<(n+j) for j in range(m)) for i in range(m))
    c[x] ^= 1<<(n+y)
    c[n+y] |= 1<<x
    return c

for j,(x,a) in enumerate(LEFT):
    require(profile(a,x)==(109,131,141+j))
    for k,(y,b) in enumerate(RIGHT_START):
        require(profile(b,y)==((3,1,2),(5,3,4))[k])
        value=count(join(a,x,b,y))
        require(value==(1133+2*j,2439+6*j)[k])
        print(j,k,value)

SAME_CONTINUATION = (
    (7,(28,37,56,194,106,201,135,23),(657,540,617)),
    (6,(26,44,97,148,102,201,139,23),(657,541,616)),
    (7,(28,41,82,164,106,197,139,23),(657,542,615)),
)
for j,(x,a,p) in enumerate(SAME_CONTINUATION):
    require(profile(a,x)==p)
    u=[(a[i]&((1<<x)-1))|((a[i]>>(x+1))<<x) for i in range(len(a)) if i!=x]
    # Strongness is incidental here; these particular complements are strong.
    require(count(u)==153)
    y,b=RIGHT_START[0]
    r=join(a,x,b,y)
    require(count(r)==5365+2*j)
    r=join(r,x,b,y)
    require(count(r)==56217+6*j)
    print("same continuation",j,5365+2*j,56217+6*j)
