from itertools import product
from math import isqrt
from replay import rows,dp,packets

WORDS=[(8850,54),(8881,54),(8881,53)]
QUAD=[(528449,m) for m in (99,66,93,94)]

def paint(x):
    return [int(c) for c in bin(x)[2:]]
def p(x):
    a=0
    while x:
        k,b=divmod(x,2)
        if not b:a+=k+(k&1)
        x=k
    return a

def forward(x,y,j):
    return 432*x*y+129*x+258*y+15*(x&1)+30*(y&1)+45*(p(x)+p(y))+135+2*j

def build(n,c,m,x,y):
    L,R=paint(x),[1-b for b in paint(y)[::-1]]
    C=[(m>>i)&1 for i in range(n)]
    size=len(L)+n+len(R)
    a=[((1<<size)-1)^((1<<(i+1))-1) for i in range(size)]
    block=((1<<n)-1)<<len(L)
    for i,v in enumerate(rows(n,c)):
        j=i+len(L)
        a[j]=(a[j]&~block)|(v<<len(L))
    marks=L+C+R
    for i,b in enumerate(marks):
        if not b:a[i]|=1<<size
    a.append(sum(b<<i for i,b in enumerate(marks)))
    return a

def row(A,b,h,N):
    # Invert A*y+h*p(y)+b*(y&1)=N over y>=0, least significant bit first.
    if N==0:
        yield 0
    elif N>=A:
        for z in row(2*A+h,h,h,N):
            if z:yield 2*z
        for z in row(2*A,0,h,N-A-b):
            yield 2*z+1

def invert(N):
    if N<999 or N%2==0:return
    for x in range(1,isqrt(N//432)+1):
        for j in range(3):
            C=135+2*j+45*p(x)
            for y in row(432*x+258,30,45,N-C-129*x-15*(x&1)):
                if y>0:
                    assert forward(x,y,j)==N
                    yield x,y,j
            for y in row(432*x+129,15,45,N-C-258*x-30*(x&1)):
                if y>x:
                    assert forward(y,x,j)==N
                    yield y,x,j

if __name__=='__main__':
    for n,B,expected,profile in [(6,WORDS,[135,137,139],(45,432,0,0)),(7,QUAD,[151,153,155,157],(27,564,0,0))]:
        # The complete (s,t,e,f) collision for WORDS is checked by full_trial.py.
        for (c,m),v in zip(B,expected):
            q,a,_=list(packets(rows(n,c)))[m]
            if n==7: assert q==profile
            assert a==v
    for x,y in [(1,1),(2,3),(3,2),(3,5),(5,3),(9,13)]:
        checks=[dp(build(6,c,m,x,y))[0][-1] for c,m in WORDS]
        assert checks==[forward(x,y,j) for j in range(3)]
        for N in checks: assert set(invert(N))
        print(x,y,checks)
    for x in [1,2,3,5,11]:
        v=[dp(build(7,c,m,x,x))[0][-1] for c,m in QUAD]
        base=564*x*x+537*x+27*(x&1)+54*p(x)+151
        assert v==list(range(base,base+8,2))
        print('mirror',x,v)
    cap=20000
    want=set()
    for x in range(1,cap//432+1):
        for y in range(1,cap//(432*x)+1):
            for j in range(3):
                N=forward(x,y,j)
                if N<=cap:want.add(N)
    have={N for N in range(999,cap+1,2) if any(invert(N))}
    assert want==have
    print('exact inverse coverage through',cap,len(have))
