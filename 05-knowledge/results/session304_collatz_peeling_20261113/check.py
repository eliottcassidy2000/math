"""Replay bounded controls for Session304 (proofs are in Markdown)."""
from itertools import permutations

def D(n,s):
    y=3*n+s
    return y//(y&-y)
def odd(n):
    return n//(n&-n)
def plan(n):
    s=1 if n%4==1 else -1
    y=D(n,s)
    if y!=7:return [y]
    y=D(n,-s)
    if n==9:return [y,D(y,1)]
    z=D(y,-s)
    if n==19:return [y,z]
    if n==37:return [y,z,D(z,1)]
    return [y,z,D(z,-s)]

def H(w):
    a=b=0
    h=1
    for x in map(int,w):
        if not x:h+=a+b
        a,b=2*a+x,x
    return h

def literal(w):
    n=len(w)
    def edge(i,j):
        if i==n:return w[j]=='1'
        if j==n:return w[i]=='0'
        return i<j
    return sum(all(edge(i,j) for i,j in zip(p,p[1:])) for p in permutations(range(n+1)))

def check(b):
    if not b:raise ValueError('failed')

for n in range(3,100002,2):
    if n in (7,21):continue
    p=plan(n)
    check(299*p[-1]<=253*n and all(v not in (7,21) for v in p))
    check(all(y in (D(x,1),D(x,-1)) for x,y in zip([n]+p,p)))
f={299}; minimum=299
for _ in range(3):
    f={D(x,s) for x in f for s in [-1,1]}-{7,21}
    minimum=min(minimum,min(f))
check(minimum==253)
for k in range(1,17):
    x=2**k
    cases=[('100'+'1'*k+'00','1000'+'1'*k+'010',15*x+3,45*x+7),
           ('101'+'0'*k+'10','1001'+'0'*k+'110',15*x+1,45*x+1)]
    for a,b,u,v in cases:
        it=iter(b)
        check(all(e in it for e in a))
        check(H(a)==u and H(b)==v and v==3*u-2)
        if k<=1:check(literal(a)==u and literal(b)==v)
print('PASS arithmetic, sharpness, and physical charts')
