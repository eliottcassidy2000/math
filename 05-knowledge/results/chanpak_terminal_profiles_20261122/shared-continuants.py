from itertools import combinations,product
from collections import deque
from math import prod

def trans(n):
    return tuple(((1<<n)-1)^((1<<(i+1))-1) for i in range(n))

def rows(n,m):
    a=list(trans(n))
    for i,(x,y) in enumerate(combinations(range(n),2)):
        if m>>i&1:a[x]^=1<<y;a[y]^=1<<x
    return tuple(a)

def path_read(A):
    n=len(A);N=1<<n
    p=[[0]*n for _ in range(N)]
    for i in range(n):p[1<<i][i]=1
    for m in range(1,N):
        for i,v in enumerate(p[m]):
            if not v:continue
            e=A[i]&~m
            while e:
                b=e&-e;e-=b;p[m|b][b.bit_length()-1]+=v
    return [1]+list(map(sum,p[1:]))

def Z(A,qmax):
    h=path_read(A);N=len(h);z=[1]+[0]*(N-1);ret=[z[-1]]
    for q in range(qmax):
        w=[0]*N
        for m in range(N):
            sub=m
            while True:
                w[m]+=h[sub]*z[m^sub]
                if not sub:break
                sub=(sub-1)&m
        z=w;ret.append(z[-1])
    return tuple(ret)

def inflate(gaps):
    n=1+sum(1+len(a) for a in gaps)
    rows=list(trans(n));j=0
    for a in gaps:
        m=len(a);k=j+m+1
        rows[j]^=1<<k;rows[k]^=1<<j
        for u,r in enumerate(a):
            for v in range(u+1,m):
                if not r>>v&1:
                    rows[j+1+u]^=1<<(j+1+v)
                    rows[j+1+v]^=1<<(j+1+u)
        j=k
    return tuple(rows)

def advance(state,p,q):
    d,u,v=state
    a,b,c=p[q+1],p[q-1],p[q]
    return a*d+2*a*u+(a-b)*v,b*v,q*c*d+q*c*u+(q*c-b)*v

def transfer(gaps,q):
    X=(0,0,q)
    for A in gaps:X=advance(X,Z(A,q+1),q)
    return X

def strongly(A):
    n=len(A);full=(1<<n)-1
    for i in range(n):
        reach=1<<i
        while True:
            S=reach
            for j in range(n):
                if reach>>j&1:S|=A[j]
            if S==reach:break
            reach=S
        if reach!=full:return False
    return True

def old(bound):
    # First zero erases at the level of states (repo BACKCHAIN); omit it.
    seen={(0,0,1)};q=deque(seen);heights={1}
    while q:
        d,u,v=q.popleft();F=d+u+v
        W=[]
        if F>1:W.append((d+2*u,v,d+u))
        k=1
        while F+(1<<k)*(d+2*u+v)<=bound:
            W.append(((1<<k)*(d+2*u+v),0,F));k+=1
        for a in W:
            S=sum(a)
            if S<=bound and a not in seen:heights.add(S);seen.add(a);q.append(a)
    return heights

def tests():
    base=((),(0,),(2,0),(2,4,1),(2,8,3,5))
    for k in range(4):
        for P in product(base,repeat=k):
            if sum(map(len,P))+k+1 > 11:continue
            A=inflate(P);Zs=Z(A,3)
            for q in (1,2,3):
                x=transfer(P,q)
                assert sum(x)==Zs[q] and min(x)>=0,(P,q,x,Zs)
    print('literal generalized H,Z(2),Z(3) transfer checks pass',flush=True)

def find():
    O=old(100001);print('old',len(O),flush=True)
    assert len(O)==13263
    first={};first2={};misses=[];bank={}
    for n in range(1,6):
        for m in range(1<<(n*(n-1)//2)):
            A=rows(n,m);v=Z(A,3);h,c,a=v[1:]
            bank.setdefault(v,(n,A))
            if h+c not in O and strongly(A):misses.append((h+c,n,A,v))
            k=h+c
            if k in first and first[k][1]!=h:
                if not first[k][3]: print('same H differing G',k,first[k][:3],(n,h,A),flush=True);first[k]=(n,h,A,True)
            elif k not in first:first[k]=(n,h,A,False)
            k=h,c
            if k in first2 and first2[k][1]!=a:
                if not first2[k][3]: print('same H,G differing Z2',k,first2[k][:3],(n,a,A),flush=True);first2[k]=(n,a,A,True)
            elif k not in first2:first2[k]=(n,a,A,False)
    print('module-only first missing',sorted(set(misses))[:25])
    print('bank',len(bank))
    ext={}
    for p,(n,A) in bank.items():
      for q,(m,B) in bank.items():
        for z in (0,1):
            sp=(A,B) if z==0 else (A,(),B)
            H=sum(transfer(sp,1))
            if H not in O and strongly(inflate(sp)):
                ext.setdefault(H,(n+m+3+z,sp))
    print('two-module missing',len(ext),sorted(ext.items())[:14])
def verify_examples():
    K=(2,8,3,5)
    assert Z(K,3)==(0,5,32,117)
    assert strongly(inflate((K,))) and Z(inflate((K,)),1)[1]==37
    O=old(100001)
    assert len(O)==13263 and 37 not in O
    for m in (1,2,3,4,5,8,9,10):
        P=(trans(m),(),K);d,u,v=transfer(P,1);F=d+u+v
        assert F==74*(1<<m)+69 and v==10*(1<<m)+5 and F not in O
        if m<=2:assert strongly(inflate(P)) and Z(inflate(P),1)[1]==F
    X=(2,4,65,23,39,15,59);Y=(2,8,11,97,15,23,55)
    assert all(strongly(a) for a in (X,Y))
    assert Z(X,2)[1:]==(45,732) and Z(Y,2)[1:]==(37,740)
    assert [sum(transfer((a,),1)) for a in (X,Y)]==[777,777]
    assert [transfer((a,),1)[2] for a in (X,Y)]==[45,37]
    assert sum(transfer((Y,(),(0,)),1))-sum(transfer((X,(),(0,)),1))==8
    X=(80,13,25,97,74,23,166,63);Y=(128,129,27,51,99,71,15,124)
    assert all(strongly(a) for a in (X,Y))
    assert Z(X,3)[1:]==(345,5640,45225) and Z(Y,3)[1:]==(345,5640,44505)
    assert Z(inflate((X,)),2)[2]-Z(inflate((Y,)),2)[2]==1440
    for P in ((X,),(Y,),((0,),X,(),K),((0,),Y,(),K)):
        assert all(x>=0 for x in transfer(P,2))
    assert transfer(((0,),X,(),K),1)==transfer(((0,),Y,(),K),1)
    assert sum(transfer(((0,),X,(),K),2))>sum(transfer(((0,),Y,(),K),2))
    print("new absolute old-gap exclusions, equal-height boundary separation, and second-slot checks pass")

if __name__=='__main__':
    tests();verify_examples()
