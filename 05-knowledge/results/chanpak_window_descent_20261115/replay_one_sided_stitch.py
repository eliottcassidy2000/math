"""Literal subset counts for opposing sewn palettes and one full seam."""
from itertools import product

def ordinal(x,y):
    n=len(x)
    return tuple(v|(((1<<len(y))-1)<<n) for v in x)+tuple(v<<n for v in y)

def build(r,am,p,q):
    def tape(w):
        n=len(w)
        return tuple(((1<<n)-1)^((1<<(i+1))-1) for i in range(n))
    x=ordinal(ordinal(tape(p),r),tape(q))
    m=int(p[::-1] or '0',2)|(am<<len(p))|(int(q[::-1] or '0',2)<<(len(p)+len(r)))
    return x,m

def apex(r,m):
    n=len(r)
    return tuple(x|(0 if (m>>i)&1 else 1<<n) for i,x in enumerate(r))+(m,)

def tables(r):
    n=len(r); f=(1<<n)-1
    dp=[[0]*n for _ in range(f+1)]; H=[1]+[0]*f
    for K in range(1,f+1):
        s=K
        while s:
            b=s&-s; s^=b; i=b.bit_length()-1; L=K^b
            w=r[i]&L; v=0
            while w:
                z=w&-w; w^=z; v+=dp[L][z.bit_length()-1]
            dp[K][i]=v if L else 1
        H[K]=sum(dp[K])
    return dp,H

def state(r,m):
    dp,H=tables(r); f=len(H)-1
    S=d=0
    for i in range(1,f+1):
        v=sum(x for j,x in enumerate(dp[i]) if (m>>j)&1)*H[f^i]
        S+=v; d+=v if i.bit_count()%2 else -v
    return H[-1],tables(apex(r,m))[1][-1],S,d

A=((17,'10',6,'11'),(23,'',6,'1011'),(27,'1010',7,''),(29,'1',4,'011'),(43,'101',4,'1'))
B=((233,'100',9,'0'),(255,'',0,'111000'),(259,'',8,'11000'),(261,'',6,'1000'),(267,'10',2,'00'))
D=[]
for r,bank,profile in [((2,4,1),A,(3,129,3)),((6,8,10,1),B,(5,280,0))]:
    P={}
    for N,p,m,q in bank:
        x,mark=build(r,m,p,q)
        if state(x,mark)!=(profile[0],N,*profile[1:]):
            raise AssertionError((r,N,p,m,q,state(x,mark)))
        P[N]=(x,mark)
    D.append(P)
U={5*a+3*b:(a,b) for a,b in product(D[0],D[1])}
if not all(x in U for x in range(910,919,2)) or 920 in U or len({x%30 for x in U})!=15:
    raise AssertionError(U)
a,b=U[914]; x,m=D[0][a]; y,l=D[1][b]
f=(1<<len(y))-1
y=tuple(f^(1<<j)^r for j,r in enumerate(y)); l^=f
z=ordinal(ordinal(x,(0,)),y); mz=m|(l<<(len(x)+1))
H=tables(apex(z,mz))[1][-1]
if H!=36765*2+269:
    raise AssertionError(H)
print('Ten literal one-sided palette members verified')
print('25 outputs, 15 odd residues mod 30; offsets 910,912,914,916,918 present; 920 absent')
print(f'Full literal order-{len(z)+1} tournament: H={H} at height one')
