from functools import lru_cache
from collections import defaultdict

def ending(A):
    n=len(A); F=(1<<n)-1
    out=[[0]*n for _ in range(F+1)]
    for i in range(n):out[1<<i][i]=1
    for S in range(1,F+1):
        for i,x in enumerate(out[S]):
            if x:
                B=A[i]&~S
                while B:
                    b=B&-B;B^=b
                    out[S|b][b.bit_length()-1]+=x
    return out
def count(A):
    return sum(ending(A)[-1]) if A else 1
def rev(x,n):
    return int(f'{x:0{n}b}'[::-1],2)
def g(x):
    h=1
    while x:
        x,b=divmod(x,2)
        if not b:h+=x+(x&1)
    return h

def profile(A,U):
    n=len(A);F=(1<<n)-1
    ee=ending(A)
    ss=ending([r^F^(1<<i) for i,r in enumerate(A)])
    h=[1]+[sum(row) for row in ee[1:]]
    f=sum(h[s]*h[F^s] for s in range(F+1))
    P=sum(ss[s][i]*h[F^s] for s in range(F+1) for i in range(n) if U>>i&1)
    Q=sum(ee[s][i]*h[F^s] for s in range(F+1) for i in range(n) if not U>>i&1)
    S=sum(ee[F][i] for i in range(n) if U>>i&1)
    T=sum(ss[F][i] for i in range(n) if not U>>i&1)
    K=count(tuple(r|((((U>>i)&1)^1)<<n) for i,r in enumerate(A))+(U,))
    E=sum(ee[s][i]*ss[F^s][j] for s in range(F+1) for i in range(n) for j in range(n) if U>>i&1 and U>>j&1)
    assert K-P+S==h[F]-E
    return h[F],f,P,Q,S,T,K

def physical(A,U,x,y):
    left=f'{x:b}'
    right=''.join('1' if b=='0' else '0' for b in f'{y:b}'[::-1])
    l,n,r=len(left),len(A),len(right); F=(1<<(l+n+r))-1
    C=[]
    for i,b in enumerate(left):
        C.append((F^((2<<i)-1))|(((b=='0')*1)<<(l+n+r)))
    for i,R in enumerate(A):
        C.append((R<<l)|(((1<<r)-1)<<(l+n))|((((U>>i)&1)^1)<<(l+n+r)))
    for j,b in enumerate(right,l+n):
        C.append((F^((2<<j)-1))|(((b=='0')*1)<<(l+n+r)))
    return C+[sum((b=='1')<<i for i,b in enumerate(left))|(U<<l)|sum((b=='1')<<(l+n+i) for i,b in enumerate(right))]

from itertools import product

@lru_cache(None)
def info(A,U):
    h,f,P,Q,S,T,K=profile(A,U)
    E=h+P-S-K
    return (h,f,P,Q,S,T,E)

def combine(A,B):
    ha,fa,pa,qa,sa,ta,ea=A
    hb,fb,pb,qb,sb,tb,eb=B
    return (ha*hb,fa*fb,pa*fb+ha*pb,qb*fa+hb*qa,ha*sb,ta*hb,
            hb*ea+ha*eb+pa*(fb-hb-qb)+sa*(hb-tb))

def fold(blocks):
    from functools import reduce
    return reduce(combine, [info(*a) for a in blocks])

def readout(C,x,y):
    h,f,P,Q,S,T,E=C
    return f*x*y+Q*x+P*y+h*(g(x)+g(y)-2)+h+P-S-E+T*(x&1)+S*(y&1)

@lru_cache(None)
def literal(C,Ubits):
    E=ending(C)
    S=ending([r^((1<<len(C))-1)^(1<<i) for i,r in enumerate(C)])
    I=(1<<len(C))-1
    def c(mask,r,inside):
        return sum(r[mask][i] for i in range(len(C)) if bool(Ubits>>i&1)==inside)
    H=[[0]*3 for _ in range(3)]
    for m in range(I+1):
        w=I^m
        L=[int(m==0), c(m,E,False),c(m,E,True)]
        R=[int(w==0), c(w,S,True),c(w,S,False)]
        for i,j in product(range(3),repeat=2):
            H[i][j]+=L[i]*R[j]
    return H

def direct(blocks,x,y):
    # A pair of old paths crossing forward through the ordinal sum can be
    # joined at the apex iff its final first endpoint and initial second
    # endpoint are each status 0 (absent) or 1 (valid).
    def wing(w):
        n=len(w)
        return tuple(((1<<n)-1)^((2<<i)-1) for i in range(n)),sum((b=='1')<<i for i,b in enumerate(w))
    wL=wing(f'{x:b}')
    wR=wing(''.join(str(int(b)^1) for b in f'{y:b}'[::-1]))
    pairs={(0,0):1}
    for B in [wL,*blocks,wR]:
        M=literal(*B)
        Q=defaultdict(int)
        for (i,j),k in pairs.items():
            for a,b in product(range(3),repeat=2):
                Q[a or i,j or b]+=k*M[a][b]
        pairs=Q
    return sum(pairs[i,j] for i,j in product(range(2),repeat=2))

Egy=[
    (((2,20,9,35,77,23,47),17),((2,20,33,7,45,75,31),33)),
    (((130,172,169,209,135,217,151,0),130),((0,9,147,101,171,199,23,75),9)),
    (((2,12,17,37,139,151,63,79),130),((2,12,17,37,139,151,63,79),36))]
for pair,h in zip(Egy,[65,109,161]):
    a,b=(info(*o) for o in pair)
    assert a[0]==h and a[:-1]==b[:-1] and b[-1]==a[-1]+2
v1=[Egy[0][1],Egy[1][0],Egy[2][0]]
v2=[Egy[0][0],Egy[1][1],Egy[2][1]]
for marked_filler in [[], [((2,0),1)], [((2,0),1),((0,),0)]]:
    b1=v1[:1]+marked_filler+v1[1:]; b2=v2[:1]+marked_filler+v2[1:]
    state=[fold(b) for b in [b1,b2]]
    assert state[0][:-1]==state[1][:-1] and state[1][-1]==state[0][-1]+2
    for x,y in [(1,1),(2,7),(30,11),(231,314)]:
        a,b=[direct(blocks,x,y) for blocks in [b1,b2]]
        assert a==readout(state[0],x,y) and b==readout(state[1],x,y) and a-b==2
        if (x,y)==(1,1): print('Egyptian literal counts',a,b,'fillers',len(marked_filler))

a9=((0,133,9,147,135,159,191,5),139)
a45=((10,20,33,22,37,11),12)
def v3(x):
    assert x
    v=0
    while x%3==0: x//=3;v+=1
    return v
def lower(word):
    c=[0]+[2 if q==9 else 1 for q in word]+[0]
    d=[0]+[1 if q==9 else -1 for q in word]+[0]
    L=max([2]+[c[i]+c[j]+sum(d[i+1:j]) for i in range(len(c)) for j in range(i+1,len(c))])
    return max(0,2*len(word)-L)
for k in range(2,9):
    for word in product((9,45),repeat=k):
        B=[a9 if q==9 else a45 for q in word]
        C=fold(B)
        X=readout(C,2,3)
        assert X==direct(B,2,3) and v3(X)>=lower(word)
        if k>=3 and all(q==9 for q in word):assert v3(X)==k-2
        # arbitrary one-way filler, transitive h=f-unit and C3 plural.
        for fill in [((2,0),1),((2,4,1),3)]:
            F=B[:1]+[fill]+B[1:]
            Y=readout(fold(F),3,2)
            assert Y==direct(F,3,2) and v3(Y)>=lower(word)
            if k>=3 and all(q==9 for q in word):
                assert (v3(Y)==k-2)==(info(*fill)[1]%3!=0)
print('Mixed valuation and pure-h9 exactness replayed with two legal filler types.')
