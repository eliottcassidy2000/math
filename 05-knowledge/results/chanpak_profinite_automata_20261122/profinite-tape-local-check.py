from collections import deque
from math import gcd
from functools import lru_cache

def run(state,word,q):
    y,f,b=state
    for ch in word:
        c=int(ch)
        f=(f+(1-c)*(y+b))%q
        y=(2*y+c)%q
        b=c
    return y,f,b

def check(q):
    m=q
    while not m%2:m//=2
    dq=q//m
    t=max(1,dq.bit_length()-1)
    while pow(2,t,m)!=1%m:t+=1
    loops=[b+'0'*t for b in ('0','10','110','1010')]
    start=run((0,0,0),'10011'+'0'*t,q)
    Q=deque([start]); seen={start}
    while Q:
        state=Q.popleft()
        for w in loops:
            v=run(state,w,q)
            if v not in seen:seen.add(v);Q.append(v)
    assert len(seen)==m*m*max(1,dq//2),(q,len(seen))
    assert all(y%dq==0 and b==0 for y,f,b in seen)
    for suffix in ('','11001','0001','10000'):
        outputs=set()
        for state in seen:
            y,f,b=run(state,suffix,q)
            outputs.add((43+81*y+25*f+25*b)%q)
        expected={x%q for x in range(1,2*q,2)}
        assert outputs==expected,(q,suffix)

def graph(y):
    core=(4,5,24,19,35,15)
    right=bin(y)[2:][::-1].translate(str.maketrans('01','10')) if y else ''
    tape='1000001'+right
    n=len(tape)
    rows=[((1<<n)-1)^((1<<(i+1))-1) for i in range(n)]
    for i,r in enumerate(core): rows[i]=(rows[i]&~63)|r
    new=0
    for i,b in enumerate(tape):
        if b=='1':new|=1<<i
        else:rows[i]|=1<<n
    return rows+[new]

def h(rows):
    @lru_cache(None)
    def d(A,i):
        B=A^(1<<i)
        if not B:return 1
        return sum(d(B,j) for j in range(len(rows)) if B>>j&1 and rows[j]>>i&1)
    return sum(d((1<<len(rows))-1,i) for i in range(len(rows)))

def strong(rows):
    for a in range(len(rows)):
        seen=1<<a
        while True:
            new=seen
            for i in range(len(rows)):
                if seen>>i&1:new|=rows[i]
            if new==seen:break
            seen=new
        if seen!=(1<<len(rows))-1:return False
    return True

for y in range(32):
    low=graph(y)
    high=low[:];high[1]^=1<<3;high[3]^=1<<1
    _,f,b=run((0,0,0),bin(y)[2:],10**9)
    N=43+81*y+25*f+25*b
    assert not y or (strong(low) and strong(high)),(y,low,high)
    assert (h(low),h(high))==(N,N+2)
for q in list(range(1,65))+[72,81,96,125,128,189,225,256,324]:
    check(q)
print('31 literal strong seams plus the zero-address sink case; joint state and four terminal cylinders at 73 moduli')
