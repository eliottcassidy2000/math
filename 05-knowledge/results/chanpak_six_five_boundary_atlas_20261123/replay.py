from itertools import combinations
from fractions import Fraction as F
V=range(6)
boundary={tuple(sorted((i,(i+1)%6))) for i in V}
D=set(combinations(V,2))-boundary
def crossing(a,b):
    x,y=a;u,v=b
    return len({x,y,u,v})==4 and ((x<u<y) != (x<v<y))
T=[S for S in combinations(sorted(D),3) if not any(crossing(a,b) for a,b in combinations(S,2))]
assert len(T)==14
long={(0,3),(1,4),(2,5)}
assert all(sum(d in S for S in T)==4 for d in long)
assert sum(not (set(S)&long) for S in T)==2
def mirror(d):return tuple(sorted((1-x)%6 for x in d))
assert [d for d in sorted(long) if mirror(d)==d]==[(2,5)]
assert all({mirror(d) for d in S}!=set(S) for S in T if (2,5) in S)
W=('0','10','11','12','2')
H=(19,61,63,65,23)
p=tuple(F(1,3**len(w)) for w in W)
for w,h in zip(W,H):
    mid=sum(F(int(c),3**(i+1)) for i,c in enumerate(w))+F(1,2*3**len(w))
    assert F(h,3**len(w))==6+2*mid
assert sum(t*F(h,3**len(w)) for w,h,t in zip(W,H,p))==7
R=[]
for a in (0,2):
  for b in (1,3):
    c=[F(1,4),F(1,16),F(1,16),F(1,16),F(1,4)]
    c[2*a]=F(1,2);c[b]=F(1,8)
    assert sum(c)==1 and sum(c[1:4])==F(1,4)
    assert sum(abs(x-y) for x,y in zip(c,p))/2==F(13,72)
    R.append(c)
q=tuple(sum(row[i] for row in R)/4 for i in range(5))
assert q==(F(3,8),F(3,32),F(1,16),F(3,32),F(3,8))
assert sum(abs(x-y) for x,y in zip(q,p))/2==F(1,12)
print("hexagon, frontier centers, reflection, and binary TV checks pass")

from collections import Counter
from itertools import permutations

def pairs(n):
    return list(combinations(range(n),2))
P=set(pairs(5))-{(0,1),(1,2),(2,3),(3,4),(0,4)}
mark={(0,2),(0,3),(1,4)}
pent=[A for A in combinations(sorted(P),2) if not crossing(*A)]
assert len(pent)==5
assert Counter(len(set(A)&mark) for A in pent)=={1:4,2:1}
assert next(A for A in pent if len(set(A)&mark)==2)==((0,2),(0,3))
seeds=(1,19,61,63,65,23)
X=[None]+[F(3,2)*(F(h,3**len(w))-7) for h,w in zip(seeds[1:],W)]
assert X==[None,-1,F(-1,3),0,F(1,3),1]
R5=[None]+[F(x).numerator*pow(F(x).denominator,-1,5)%5 for x in X[1:]]
assert R5==[None,4,3,0,2,1]
def inv5(x):
    if x is None: return 0
    if x==0: return None
    return 3*pow(x,-1,5)%5
assert [inv5(x) for x in R5]==R5[3:]+R5[:3]
pe=[(a,b) for a,b in combinations(pairs(5),2) if not set(a)&set(b)]
assert len(pe)==15
colors=Counter(next(iter(set(range(5))-set(a)-set(b))) for a,b in pe)
assert colors==dict.fromkeys(range(5),3)
print('six marked pentagon fillings, real branch labels, projective half turn, and 15 Petersen edges pass')

def counts(p):
    t=next(a for a in range(1,p) if pow(a,(p-1)//2,p)==p-1)
    def add(u,v): return ((u[0]+v[0])%p,(u[1]+v[1])%p)
    def mul(u,v):
        return ((u[0]*v[0]+t*u[1]*v[1])%p,(u[0]*v[1]+u[1]*v[0])%p)
    def count(es):
        square=Counter(mul(x,x) for x in es)
        n=1
        for x in es:
            x2=mul(x,x)
            n+=square[mul(x,mul(add(x2,(-1,0)),add(mul((9,0),x2),(-1,0))))]
        return n
    n1=count([(x,0) for x in range(p)])
    n2=count([(x,y) for x in range(p) for y in range(p)])
    a=p+1-n1
    b=(n2-p*p-1+a*a)//2
    return n1,n2,1-a+b-p*a+p*p
assert counts(5)==(6,6,16)
assert counts(7)==(8,46,48)
print('curve counts at 5,25 and 7,49 and Jacobian orders 16,48 pass')

def H(W):
    n=len(W); z=n
    def arc(i,j):
        if i==z: return W[j]=='O'
        if j==z: return W[i]=='I'
        return i<j
    return sum(all(arc(i,j) for i,j in zip(p,p[1:]))
               for p in permutations(range(n+1)))
ctx=['W','WI','WII','OW','OOW','OWI','WOI']
def row(w,columns):
    return [H(ctx[i].replace('W',w)) for i in columns]
def det(A):
    if not A: return 1
    return sum((-1)**i*u*det([w[:i]+w[i+1:] for w in A[1:]])
               for i,u in enumerate(A[0]))
assert det([row(w,range(7)) for w in ('O','I','OO','OI','IO','II','OOI')])==128
assert det([row(w,(0,1,2,3,5,6)) for w in ('O','OO','OI','OOI','OIO','OII')])==64
assert det([row(w,(0,1,3,5,6)) for w in ('OI','OOI','OII','OOII','OIOI')])==64
print('literal apex probe determinants: 128, 64, 64')
