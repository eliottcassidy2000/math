from math import prod,comb
from itertools import permutations,combinations
from collections import Counter
from importlib import import_module
c=import_module('e61-cyclic93')
M=c.M
D=(8,5,17,6,11)
B=(
(1321,775,5,(256,289,67,23,295,141,59,95,236)),
(1655,3641,3,(258,12,393,305,295,71,31,123,226)),
(8723,1321,7,(768,269,73,561,551,391,59,351,604,230)),
(1359,1599,4,(258,288,283,403,99,141,303,87,176)),
(155,427,3,(2,40,35,21,39,201,31,95)),
(287,151,5,(10,20,25,82,161,15,55,111)),
)
U=tuple(x for *_,x in B)+((0,),)*33

def inflate(D,L,z=0):
    A=[];v=[]
    for i in range(len(D)):
        v.extend((i,j) for j in range(len(L) if i==z else 1))
    for i,j in v:
        A.append(sum(1<<a for a,(b,t) in enumerate(v) if
          (b!=i and D[i]>>b&1) or (b==i==z and L[j]>>t&1)))
    return tuple(A)

def recipes(D,new=None):
    q=len(D);V=tuple(range(1,q));out=[]
    for outer in permutations(V):
      for tag in range(1,1<<q):
        P=[]
        for j in range(q):
            if tag>>j&1:P.append(0)
            if j+1<q:P.append(outer[j])
        require=None; forbid=None;good=True
        for u,v in zip(P,P[1:]):
            goodedge=bool(D[u]>>v&1)
            if new is not None and (u,v) in ((0,new),(new,0)):
                span=('E' if u==0 else 'S')
                if goodedge:forbid=span
                elif require is None:require=span
                else:good=False;break
            elif not goodedge:
                good=False;break
        if good:out.append((tag.bit_count(),require if require else forbid,require is not None))
    return Counter(out)

def cover(Z,r):
    return sum((-1)**(r-j)*comb(r,j)*Z[j] for j in range(1,r+1))

def rooted(P,i,u,r,mode):
    _,s,e=P[i]
    if mode=='S':
        def Q(j):return s[j-1][u]*prod(x[0][j-1] for x in P[:i])*prod(x[0][j] for x in P[i+1:])
    else:
        def Q(j):return e[j-1][u]*prod(x[0][j] for x in P[:i])*prod(x[0][j-1] for x in P[i+1:])
    return sum((-1)**(r-j)*comb(r-1,j-1)*Q(j) for j in range(1,r+1))

def evaluate(rec,P,Z,i=None,u=None):
    a=0
    for (r,mode,req),t in rec.items():
        v=cover(Z,r)
        if mode is not None:
            s=rooted(P,i,u,r,mode)
            v=s if req else v-s
        a+=t*v
    return a

def controls():
  for L in (((2,4,1),),((2,8,3,5),(0,)),((0,),(2,4,1),(0,))):
    P=tuple(c.read(X,len(D)) for X in L)
    Z=tuple(prod(X[0][r] for X in P) for r in range(len(D)+1))
    A=inflate(D,c.ordsum(L));n=sum(map(len,L))
    assert c.read(A,1,False)[1]==evaluate(recipes(D),P,Z)
    for a,b in combinations(range(1,len(D)),2):
        assert c.read(c.flip(A,n+a-1,n+b-1),1,False)[1]==evaluate(recipes(c.flip(D,a,b)),P,Z)
    for a in range(1,len(D)):
        rec=recipes(D,a)
        d=0
        for i,X in enumerate(L):
            for u in range(len(X)):
                got=evaluate(rec,P,Z,i,u)
                assert c.read(c.flip(A,d+u,n+a-1),1,False)[1]==got,(L,a,i,u)
            d+=len(X)

def scan(L):
    P=tuple(c.read(X,len(D)) for X in L);n=sum(map(len,L))
    Z=tuple(prod(X[0][r] for X in P) for r in range(len(D)+1))
    H,Q=Z[1:3];target=H+6*Q
    assert H==M and 6*Q==(M<<61) and target==c.TARGET
    assert recipes(D)==Counter({(1,None,False):13,(2,None,False):6})
    def all_changes():
      for a,b in combinations(range(1,len(D)),2):
        yield ('outer',a,b),evaluate(recipes(c.flip(D,a,b)),P,Z)-target
      for a in range(1,len(D)):
        rec=recipes(D,a)
        for i,X in enumerate(L):
          for u in range(len(X)):
            yield ('contact',a,i,u),evaluate(rec,P,Z,i,u)-target
      for i,X in enumerate(L):
        F,*_=P[i]
        for u in range(len(X)):
          for v in range(u+1,len(X)):
            R=c.read(c.flip(X,u,v),2,False)
            yield ('internal',i,u,v),H//F[1]*(R[1]-F[1])+6*Q//F[2]*(R[2]-F[2])
          for j,Y in enumerate(L[i+1:],i+1):
              for v in range(len(Y)):
                yield ('cross',i,j,u,v),c.weighted_cross(P,i,j,u,v,1)+6*c.weighted_cross(P,i,j,u,v,2)
    count=Counter();nearest={};seams=[]
    for tag,y in all_changes():
      name=tag[0];count[name]+=1
      if (d:=nearest.get(name)) is None or (abs(y-2),y,tag)<d:nearest[name]=(abs(y-2),y,tag)
      if y==2:seams.append(tag)
    assert sum(count.values())==comb(n+len(D)-1,2)
    return H,Q,count,nearest,seams

if __name__=='__main__':
    controls()
    for h,odd,v,A in B:
        assert c.strong(A)
        assert c.read(A,2,False)[1:]==(h,odd<<v)
    h,q,counts,nearest,seams=scan(U)
    A=inflate(D,c.ordsum(U))
    assert len(A)==90 and c.strong(A)
    assert all(any((D[i]&m) not in (0,m) for i in range(5) if not m>>i&1)
               for m in range(1,31) if m.bit_count()>1)
    assert not seams and tuple(counts[k] for k in ("outer","contact","internal","cross"))==(6,344,209,3446)
    print('strong on',len(A),'vertices: H=',h+6*q)
    print('arc counts',dict(counts))
    print('nearest to +2',nearest)
    print('seams:',seams)

def class_contacts():
    cores=tuple(x for *_,x in B)
    p=tuple(c.read(A,3) for A in cores+((0,),))
    L=len(cores);mx=1<<L
    Z=tuple(prod(X[0][r] for X in p[:L])*r**33 for r in range(4))
    ratios=tuple(tuple(prod(p[i][0][r] for i in range(L) if S>>i&1) for r in range(4)) for S in range(mx))
    powers=tuple(tuple(q**j for j in range(34)) for q in range(4))
    Cs=[cover(Z,r) for r in range(4)]
    recipes90=[]
    for a in range(1,5):
        d=-c.TARGET;s=[0]*4;e=[0]*4
        for (r,mode,req),t in recipes(D,a).items():
            if not req:d+=t*Cs[r]
            if mode:
                add=s if mode=='S' else e
                for q in range(1,r+1):add[q]+=t*(1 if req else -1)*(-1)**(r-q)*comb(r-1,q-1)
        recipes90.append((d,s,e))
    nn=[0]*4;near=[None]*4;hits=[]
    for i in range(L+1):
        avail=(mx-1)^(1<<i) if i<L else mx-1
        n_pad=33-(i==L)
        Z0,s0,e0=p[i]
        before=avail
        while True:
            after=avail^before
            for ell in range(n_pad+1):
                rr=n_pad-ell
                za=tuple(ratios[before][q]*powers[q][ell] for q in range(4))
                zb=tuple(ratios[after][q]*powers[q][rr] for q in range(4))
                for u in range(len(s0[0])):
                    S=(0,)+tuple(s0[q-1][u]*za[q-1]*zb[q] for q in range(1,4))
                    E=(0,)+tuple(e0[q-1][u]*za[q]*zb[q-1] for q in range(1,4))
                    for a,(d,s,e) in enumerate(recipes90):
                        y=d+sum(s[q]*S[q]+e[q]*E[q] for q in range(1,4))
                        nn[a]+=1;t=(abs(y-2),y,i,before,ell,u)
                        if near[a] is None or t<near[a]:near[a]=t
                        if y==2:hits.append((a,i,before,ell,u))
            if before==0:break
            before=(before-1)&avail
    return tuple(nn),near,hits

if __name__=='__main__':
    counts,minimum,hits=class_contacts()
    assert counts==(59776,)*4 and not hits
    print('class-wide mixed exact cases per outer vertex:',counts)
    print('class-wide nearest contacts to +2:',minimum)
    print('class-wide +2 contacts:',hits)
