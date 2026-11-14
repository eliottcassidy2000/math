"""Two-coordinate literal digit bank and nonunit-residue tail compiler."""
from dataclasses import dataclass
from itertools import product
from random import Random
from replay_f1 import Core,audit,path_arrays,core29,make
from compile_absolute import binary_h as g
from replay_positional import A,puncture,rev
FIRST=(None,3806,3222,3488,2734)
LATER=(6430,4746,6932,5848,5534)
F0=4096; F=8192
def profile_menu():
    H,S=path_arrays(A)
    _,E=path_arrays(tuple(15^(1<<i)^R for i,R in enumerate(A)))
    result=[]
    for U in range(16):
        L=[i for i in range(4) if U>>i&1]; V=[i for i in range(4) if not (U>>i&1)]
        P=sum(S[s][i]*H[15^s] for s in range(1,16) for i in L)
        Q=sum(E[s][i]*H[15^s] for s in range(1,16) for i in V)
        Ss=sum(E[-1][i] for i in L);T=sum(S[-1][i] for i in V)
        R=[r|((0 if i in L else 1)<<4) for i,r in enumerate(A)]+[U]
        result.append((P,Q,Ss,T,path_arrays(R)[0][-1]))
    return result

def extract():
    Pout=profile_menu()
    first={p:{} for p in FIRST[1:]}
    later={p:{} for p in LATER}
    for D,G,useK in [(7,first,True),(8,later,False)]:
        for a in range(D+1):
            b=D-a
            for x,z,u in product(range(1<<a),range(1<<b),range(16)):
                p,q,s,t,K=Pout[u]
                P=(1<<b)*(32*x+p)+5*z
                if P not in G:continue
                S=5*(z%2) if b else s
                T=5*(1-(x>>(a-1))) if a else t
                if S or T:continue
                y=(1<<b)-1-rev(z,b)
                Q=(1<<a)*(32*y+q)+5*((1<<a)-1-rev(x,a))
                K=K+5*(g(x)+g(y)-2)+q*x+p*y+32*x*y+t*(x%2)+s*(y%2)
                rows,mask=puncture(a,b,x,z,u)
                digit=Core(rows,mask,5,32<<D,P,Q,S,T,K)
                G[P].setdefault((K if useK else Q)%5,digit)
        assert all(len(v)==5 for v in G.values())
    assert tuple(p%5 for p in FIRST[1:])==(1,2,3,4)
    assert tuple(p%5 for p in LATER)==(0,1,2,3,4)
    return first,later

def class_inverse(w,k):
    # Separate all k actual P labels from their prescribed mod-5^k weighted sum.
    classes=[];z=w-F0*pow(F,k-1)*10
    for j in range(k):
        r=z*pow(F,-(k-1-j),5)%5
        P=(FIRST if j==0 else LATER)[r]
        z=(z-pow(F,k-1-j)*P)//5
        classes.append(P)
    return classes

def prefix(classes):
    out=[10]
    for j,P in enumerate(classes):
        out.append((F if j else F0)*out[-1]+pow(5,j)*P)
    return out

def read(sources):
    p=prefix([d.P for d in sources])
    H=9
    for j,d in enumerate(sources):
        H=5*H+5**j*(d.K-5)+p[j]*d.Q
    return p[-1], H

# Tail blocks with h coprime to 5 and f admitting valuations 1 and 2.
tail1=core29[0]
tail2=Core((24,1,3,6,14),0,11,100,0,0,0,0,0) # only H and F2 used
def tail(k,t):
    q,r=divmod(t,2)
    # 5*floor(t/2)+6*(t mod 2) <= 5*floor(k/2)+6.
    L=5*(k//2)+7
    s=L-5*q-6*r   # at least one singleton so the final reader is nonempty
    h=11**q*29**r; f=100**q*320**r*2**s
    return h,f,q,r,s

def solve(Hres,Rres,k,groups,actual=False):
    """Physical witness for any residue pair; exact pair membership if actual=True."""
    assert k>=1
    M=5**k
    RR=Rres%M;t=0
    if RR==0:t=k
    else:
        while RR%5==0: t+=1; RR//=5
    h,f,q,r,s=tail(k,t)
    if actual:
        P,z=divmod(Rres-M*h,f)
        if z or P<=0 or P%5==0:return None
    elif t==k:
        P=1
    else:
        W=5**(k-t)
        P=RR*pow(f//5**t,-1,W)%W
    c=class_inverse(P,k)
    p=prefix(c)
    if actual and P!=p[-1]:return None
    z=Hres-(f-h)*p[-1]
    if actual:
        z,rem=divmod(z,h)
        if rem:return None
    else:
        z=z*pow(h,-1,M)%M
    ds=[]
    for j,Pj in reversed(list(enumerate(c))):
        bank=(groups[1] if j else groups[0])[Pj]
        if j:
            d=bank[z*pow(p[j],-1,5)%5]
        else:
            d=bank[(z+5)%5] # P0=10 kills Q locally
        z=(z-5**j*(d.K-5)-p[j]*d.Q)//5
        ds.append(d)
    if actual and z!=9:return None
    return ds[::-1],t

def output(k,t,ds):
    P,H=read(ds);h,f,*_=tail(k,t)
    return h*H+(f-h)*P,5**k*h+f*P

def literal(k,t,ds):
    h,f,q,r,s=tail(k,t)
    A2=Core(tail2.masks,0,0,0,0,0,0,0,0)
    A1=Core(tail1.masks,0,0,0,0,0,0,0,0)
    return make(ds+[A2]*q+[A1]*r,10,(1<<s)-1)

def minimum_census():
    from collections import defaultdict
    menu=profile_menu();out=[]
    for D in range(9):
        G=defaultdict(lambda:[set(),set()])
        for a in range(D+1):
            b=D-a
            for x,z,u in product(range(1<<a),range(1<<b),range(16)):
                p,q,s,t,K=menu[u]
                S=5*(z%2) if b else s
                T=5*(1-(x>>(a-1))) if a else t
                if S or T:continue
                P=(1<<b)*(32*x+p)+5*z
                y=(1<<b)-1-rev(z,b)
                Q=(1<<a)*(32*y+q)+5*((1<<a)-1-rev(x,a))
                K=K+5*(g(x)+g(y)-2)+q*x+p*y+32*x*y+t*(x%2)+s*(y%2)
                G[P][0].add(K%5);G[P][1].add(Q%5)
        out.append(tuple({p%5 for p,S in G.items() if len(S[i])==5} for i in range(2)))
    return out

if __name__=="__main__":
    groups=extract()
    candidates={d for g in groups for tbl in g.values() for d in tbl.values()}
    for d in candidates:audit(d)
    for D in (tail1,tail2):
        H,_=path_arrays(D.masks)
        assert (H[-1],sum(x*y for x,y in zip(H,H[::-1])))==(D.h,D.f)
    rng=Random(4032026)
    for k in range(1,7):
        M=5**k
        controls=product(range(M),repeat=2) if k<=3 else ((rng.randrange(M),rng.randrange(M)) for _ in range(600))
        for H,R in controls:
            D,t=solve(H,R,k,groups)
            e,p=output(k,t,D)
            assert (e%M,p%M)==(H,R)
            assert solve(e,p,k,groups,actual=True)==(D,t)
            assert len(literal(k,t,D))==12*k+5*(k//2)+11
    M=minimum_census()
    assert M[:7]==[(set(),set())]*7
    assert M[7]==({1,2,3,4},{1,2,3})
    assert M[8]==(set(range(5)),set(range(5)))
    print("PASS:",len(candidates),"literal digit profiles; complete pairs through 5^3; exact joint inverse")
