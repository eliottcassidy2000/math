from itertools import combinations
from math import factorial, comb
from fractions import Fraction as F

D=(0,1,2,5)
assert {(i-j)%11 for i in D for j in D}==set(range(11))
assert 3*2<10

def lower(m,d,k,q,c,bits):
    r=m-k
    L=m-2*q
    B=q-2*d-2-k
    b=r-4*(d-1)-2
    E=(d-1)*r+1
    assert B>=2*k and L>=k and b>0
    mu=k+F((d-1)*k*(k-1),L-1)
    share=1-mu/(c*k+1)
    assert share>0
    num=(d-1)**((d-1)*m)*factorial(r)*b**E*comb(L,k)*(factorial(B)//factorial(B-2*k))*share.numerator
    den=d**((d-2)*m)*(b+1)**E*share.denominator
    assert num**c>(den**c*factorial(c)**k)<<(bits*c)

lower(560,10,46,209,4,5640)

e0=F(27182818,10**7)
assert e0<sum(F(1,factorial(i)) for i in range(14))
w=tuple(F(z,1000) for z in (1000,1000,1000,1000,1209,1471,1790,2178,2650,3225,3924,4775,5811,7071,8605,10472,12744,15510,18875,22972))
wa=F(25585,1000)
def ceiling(mate):
    scale=10**6
    def ceil(z):
        z=z*scale/e0
        return -(-z.numerator//z.denominator)
    top=ceil(56*sum(wa/x for x in w)+5*wa)
    for i,x in enumerate(w):
        early=[x/y for y in w[:i]]
        val=sum(early)+55*sum(x/y for y in w[i+1:])+x/wa
        exc=sum(max(y-1,0) for y in [*early,x/wa])
        for a in range(56):
            r=val+55-a
            ex=exc
            if mate and a==0:
                r-=1
                if i:
                    y=x/w[i-1];r+=y;ex+=y-1
            if mate and a==55:
                r+=1
                if i+1<len(w):r-=x/w[i+1]
            top*=ceil(r+F(1413,500)+F(1087,500)*x+F(ex,50))
    den=scale**1121
    bit=top.bit_length()-den.bit_length()
    while top>=den<<bit:bit+=1
    print('upper',mate,bit)
    assert bit==7566
ceiling(False)
ceiling(True)

V=[(s,i,x) for s in (0,1) for i in range(5) for x in range(5)]
def edge(a,b):
    s,i,x=a
    t,j,y=b
    if s==t:return i==j and (y-x)%5 in (1+s,4-s)
    if s:return edge(b,a)
    return (x+i*j-y)%5==0
def shift(a):
    s,i,x=a
    return (0,i,(x-i)%5) if s==0 else (1,(i+1)%5,x)
for a,b in combinations(V,2):assert edge(a,b)==edge(shift(a),shift(b))
A=lambda i,x:(0,i,x)
B=lambda i,x:(1,i,x)
rows=[A(0,4)]+[B(i,0) for i in range(5)]
cols=[A(0,2)]+[B(i,1) for i in range(5)]
S=[next(a for a in V if edge(a,r) and edge(a,c)) for r in rows for c in cols]
assert len(set(S))==36
g=[S.index(shift(a)) for a in S]
es={tuple(sorted((i,j))) for i,j in combinations(range(36),2) if edge(S[i],S[j])}
assert len(es)==90

def check(order,mate,endpoint):
    n=len(order)
    out=[0]*n
    for i,j in combinations(range(n),2):
        back=edge(S[order[i]],S[order[j]])
        if j==i+1:
            assert not back
            back=mate
        if back:out[j]|=1<<i
        else:out[i]|=1<<j
    dp=[0]*(n*(1<<n))
    for i in range(n):dp[(1<<i)*n+i]=1
    allbits=(1<<n)-1
    for a in range(1,allbits):
        rem=a
        while rem:
            bit=rem&-rem; rem-=bit
            i=bit.bit_length()-1
            v=dp[a*n+i]
            if not v:continue
            z=out[i] & (allbits^a)
            while z:
                b=z&-z;z-=b
                dp[(a|b)*n+b.bit_length()-1]+=v
    h=sum(dp[-n:])
    # Permanent of augmented graph via fixed row expansion, one value per column subset.
    M=(1<<(n+1))-1
    apx=1<<n
    row=[v|apx for v in out]+[allbits]
    f=[0]*(M+1);f[0]=1
    for a in range(M):
        if not f[a]:continue
        z=row[a.bit_count()]&(M^a)
        while z:
            b=z&-z;z-=b
            f[a|b]+=f[a]
    return dp[-n+endpoint],h,f[-1]

spines=(( [0,1,4,10,9,13,15,20,22,24,25,26,35,30],10,
           ((519088,3480779,12304822),(752938,5796107,20737356)) ),
        ( [6,7,8,9,10,11,12,13,14,15,16,17,19,18,20,21,22,23],15,
           ((596935418,5571834411,24670179110),(813914947,7744560549,34510679111)) ))
O=set(spines[0][0])
seen={e for e in es if all(i in O for i in e)}
F0=set()
for _ in range(5):
    F0|=seen
    seen={tuple(sorted((g[i],g[j]))) for i,j in seen}
assert F0==es
for order,i,truth in spines:
    res=tuple(check(order,mate,i) for mate in (False,True))
    print(len(order),res)
    assert res==truth

retained=set(spines[1][0])
covered=set()
for e in es:
    a,b=e
    for _ in range(5):
        if a in retained and b in retained:
            covered.add(e)
            break
        a,b=g[a],g[b]
assert covered=={e for e in es if min(e)>=6}
