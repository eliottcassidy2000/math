from collections import Counter
from fractions import Fraction
from math import factorial, gcd, prod
from itertools import combinations

def partitions(n, lo=1):
    if not n:
        yield ()
    for a in range(lo,n+1,2):
        for t in partitions(n-a,a):
            yield (a,)+t

def count(n, power=1, center=0):
    def term(l):
        r = sum((a-1)//2 for a in l) + sum(gcd(a,b) for a,b in combinations(l,2))
        d = prod(a**k*factorial(k) for a,k in Counter(l).items())
        return Fraction(2**(power*r+center*len(l)),d)
    return sum(map(term,partitions(n)))

for n in range(4,8):
    t, s = count(n), count(n//2,2,n%2)
    print("converse",n,t,s,(t-s)/2)
    assert (t,s) == ((4,2),(12,8),(56,12),(456,88))[n-4]

br = {3:(144224,144226),4:(141420,141423),5:(137972,137974),6:(134800,134802)}
for n,(l,h) in br.items():
    assert l**n < n*100000**n < h**n
assert br[3][1]-br[4][0] < 2810
assert br[5][0]-br[6][1] == 3170
assert br[4][0]-br[5][1] > br[5][1]-br[6][0]
print("MRB root comparisons verified")

# F25 = F5[w]/(w^2+w+1), index a+5b.
def add(x,y):
    return (x%5+y%5)%5 + 5*((x//5+y//5)%5)
def mul(x,y):
    a,b=x%5,x//5
    c,d=y%5,y//5
    return (a*c-b*d)%5 + 5*((a*d+b*c-b*d)%5)
def pw(x,n):
    r=1
    for _ in range(n):
        r=mul(r,x)
    return r
def pts(q,n):
    t=[(1,x,y) for x in range(q) for y in range(q)] + [(0,1,y) for y in range(q)] + [(0,0,1)]
    return {v for v in t if add(add(pw(v[0],n),pw(v[1],n)),pw(v[2],n))==0}
fibers=Counter(tuple(pw(x,2) for x in v) for v in pts(25,6))
assert set(fibers)==pts(25,3) and len(pts(5,3))==6
assert Counter(fibers.values())=={2:9,4:27}
assert [(-h**3)%5 for h in (32,96,960,5280,640320)]==[2,4,0,0,0]
print("F5 elliptic points 6 | F25 elliptic points",len(fibers),"| Hermitian points",sum(fibers.values()),"| fibers",Counter(fibers.values()))

# (side,horizontal,vertical)
coords=[(s,h,x) for s in range(2) for h in range(5) for x in range(5)]
index={v:j for j,v in enumerate(coords)}
def adj(a,b):
    s,h,x=a
    t,i,y=b
    if s==t:
        return h==i and (x-y)%5 in ((1,4) if s==0 else (2,3))
    if s: return adj(b,a)
    return (y-x-h*i)%5==0
G=[{j for j,v in enumerate(coords) if adj(u,v)} for u in coords]
assert all(len(ns)==7 for ns in G)
assert all(len(G[i]&G[j])==(0 if i in G[j] else 1) for i,j in combinations(range(50),2))
R=[index[(s,(-h if s==0 else h)%5,-x%5)] for s,h,x in coords]
fixed=[i for i,r in enumerate(R) if i==r]
assert len(fixed)==6
assert sorted(len(set(fixed)&G[i]) for i in fixed)==[1]*5+[5]
base=list(range(50))
for i in fixed:
    R[i]=len(R)
    R.append(i)
    base.append(i)
assert all(R[R[i]]==i and R[i]!=i for i in range(56))
Gt=[{j for j in range(56) if base[j] in G[base[i]]} for i in range(56)]
assert all({R[j] for j in ns}==Gt[R[i]] for i,ns in enumerate(Gt))
assert all((Gt[i]==Gt[j])==(base[i]==base[j]) for i,j in combinations(range(56),2))
reps={i for i in range(56) if i<R[i]}
def makepath(path,left):
    if not left:
        return path
    for i in sorted(left-Gt[path[-1]], key=lambda j:(len(left-Gt[j]),j)):
        r=makepath(path+[i],left-{i})
        if r:
            return r
    return None
first=makepath([0],reps-{0})[::-1]
spine=first+[R[i] for i in first[::-1]]
assert len(set(spine))==56 and all(j not in Gt[i] for i,j in zip(spine,spine[1:]))
pos={i:k for k,i in enumerate(spine)}
def arc(i,j):
    return (pos[i]<pos[j]) != (j in Gt[i])
assert all(arc(i,j)!=arc(R[i],R[j]) for i,j in combinations(range(56),2))
print("HS fixed star:",fixed)
print("56 site spine:",spine)
print("edges in blown-up graph:",sum(map(len,Gt))//2)
