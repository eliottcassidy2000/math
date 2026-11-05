#!/usr/bin/env python3
"""Exact independent physical controls. Standard library only; -O-safe gates."""
from itertools import combinations, permutations, product
from fractions import Fraction
from math import comb
import json

def require(test, why):
    if not test: raise RuntimeError(why)

def trim(p):
    p=list(p)
    while len(p)>1 and not p[-1]: p.pop()
    return p

def add(*ps):
    q=[0]*max(map(len,ps))
    for p in ps:
        for j,c in enumerate(p):q[j]+=c
    return trim(q)

def scale(p,c):return trim([c*v for v in p])
def mul(p,q):
    r=[0]*(len(p)+len(q)-1)
    for i,a in enumerate(p):
        for j,b in enumerate(q):r[i+j]+=a*b
    return trim(r)
def mon(k):return [0]*k+[1]
def shift(p,a):
    q=[0]*len(p)
    for n,c in enumerate(p):
        for j in range(n+1):q[j]+=c*comb(n,j)*a**(n-j)
    return trim(q)
def value(p,x):
    a=0
    for c in reversed(p):a=a*x+c
    return a

def D(p):return add(mul(mon(2),p),mul(mon(1),add(shift(p,1),scale(shift(p,-1),-1))))
def E(p):return mul(mon(1),add(mul([-2,0,1],p),mul([1,-2],shift(p,-1)),mul([1,2],shift(p,1))))
def basis(k,parity):
    p=mon(1 if parity else 2)
    for j in range(1,k+1):p=mul(p,[-j*j,0,1])
    return p

def coordinates(p,parity):
    p=trim(p);c=[0]*((len(p)-1-(1 if parity else 2))//2+1)
    for k in range(len(c)-1,-1,-1):
        b=basis(k,parity);c[k]=p[len(b)-1] if len(p)>=len(b) else 0;p=add(p,scale(b,-c[k]))
    require(p==[0],('bad parity',p,parity))
    return trim(c)

def tournament(n,mask):
    a=[0]*n
    for b,(i,j) in enumerate(combinations(range(n),2)):
        u,v=(i,j) if mask>>b&1 else (j,i);a[u]|=1<<v
    return tuple(a)
def transitive(n):return tuple(sum(1<<j for j in range(i+1,n)) for i in range(n))
def validate(a):
    n=len(a)
    require(all(not(a[i]>>i&1) and a[i]<(1<<n) for i in range(n)), 'diagonal')
    require(all(((a[i]>>j)^(a[j]>>i))&1 for i,j in combinations(range(n),2)), 'tournament')
def join(a,b):
    n=len(a);m=len(b);bm=((1<<m)-1)<<n
    return tuple(x|bm for x in a)+tuple(x<<n for x in b)
def cyclic(a,b,c):
    na,nb,nc=map(len,(a,b,c));bm=((1<<nb)-1)<<na;cm=((1<<nc)-1)<<(na+nb);am=(1<<na)-1
    out=tuple(x|bm for x in a)+tuple((x<<na)|cm for x in b)+tuple((x<<(na+nb))|am for x in c)
    validate(out);return out
def wrapper(a):return cyclic(a,(0,),(0,))
def edge_wrapper(a):return cyclic(a,transitive(2),(0,))

def h_subsets(a):
    n=len(a);N=1<<n;dp=[[0]*n for _ in range(N)];h=[0]*N;h[0]=1
    for i in range(n):dp[1<<i][i]=1
    for mask in range(1,N):
        h[mask]=sum(dp[mask])
        for i,v in enumerate(dp[mask]):
            if not v:continue
            todo=a[i]&~mask
            while todo:
                bit=todo&-todo;todo-=bit;j=bit.bit_length()-1;dp[mask|bit][j]+=v
    return h

def cover_polynomial(a):
    n=len(a);h=h_subsets(a);c=[[0]*(n+1) for _ in h];c[0][0]=1
    for mask in range(1,len(h)):
        bit=mask&-mask;s=mask
        while s:
            if s&bit:
                for k in range(1,n+1):c[mask][k]+=h[s]*c[mask^s][k-1]
            s=(s-1)&mask
    return trim(c[-1])
def bin_polynomial(a):
    c=cover_polynomial(a);fall=[1];out=[0]
    for k,v in enumerate(c):
        out=add(out,scale(fall,v));fall=mul(fall,[-k,1])
    return out

def literal_cover(a):
    n=len(a);c=[0]*(n+1)
    for word in permutations(range(n)):
        for cut in range(1<<(n-1)):
            pieces=[];s=0;valid=True
            for j in range(n-1):
                if cut>>j&1:pieces.append(word[s:j+1]);s=j+1
                elif not a[word[j]]>>word[j+1]&1:valid=False;break
            if not valid:continue
            pieces.append(word[s:])
            if pieces==sorted(pieces):c[len(pieces)]+=1
    return trim(c)

def det(m):
    m=[[Fraction(v) for v in r] for r in m];d=Fraction(1)
    for j in range(len(m)):
        p=next((i for i in range(j,len(m)) if m[i][j]),None)
        if p is None:return 0
        if p!=j:m[j],m[p]=m[p],m[j];d=-d
        z=m[j][j];d*=z
        for i in range(j+1,len(m)):
            z1=m[i][j]/z
            for k in range(j,len(m)):m[i][k]-=z1*m[j][k]
    require(d.denominator==1,'nonintegral determinant');return int(d)

def discriminant(p):
    n=len(p)-1;dp=[j*p[j] for j in range(1,len(p))];f=list(reversed(p));g=list(reversed(dp));m=[]
    for j in range(n-1):m.append([0]*j+f+[0]*(n-2-j))
    for j in range(n):m.append([0]*j+g+[0]*(n-1-j))
    return (-1)**(n*(n-1)//2)*det(m)//p[-1]

def odd_cycle_subsets(a):
    n=len(a);N=1<<n;cy=[0]*N
    for root in range(n):
        d=[[0]*n for _ in range(N)];d[1<<root][root]=1
        for mask in range(1<<root,N):
            if mask&((1<<root)-1) or not mask>>root&1:continue
            if mask.bit_count()>=3 and mask.bit_count()%2:cy[mask]=sum(d[mask][v] for v in range(n) if a[v]>>root&1)
            for v,c in enumerate(d[mask]):
                if not c:continue
                todo=a[v]&~mask&~((1<<(root+1))-1)
                while todo:
                    bit=todo&-todo;todo-=bit;d[mask|bit][bit.bit_length()-1]+=c
    return cy

def packing_polynomial(a):
    n=len(a);cy=odd_cycle_subsets(a);p=[[0]*(n//3+1) for _ in cy];p[0][0]=1
    for mask in range(1,len(cy)):
        bit=mask&-mask;p[mask]=p[mask^bit][:];s=mask
        while s:
            if s&bit and cy[s]:
                for k in range(1,len(p[mask])):p[mask][k]+=cy[s]*p[mask^s][k-1]
            s=(s-1)&mask
    return trim(p[-1])

def power(p,n):
    q=[1]
    for _ in range(n):q=mul(q,p)
    return q

def run():
    report={};census=0;literal=0;physical=0
    for n in range(1,6):
        for mask in range(1<<(n*(n-1)//2)):
            a=tournament(n,mask);p=bin_polynomial(a);cs=coordinates(p,n%2)
            require(all(v>=0 for v in cs),'negative physical basis coefficient')
            require(cs[0]==h_subsets(a)[-1], 'readout')
            require(all(not p[j] for j in range(len(p)) if j%2!=n%2),'parity')
            require(all(v>=0 for v in p),'negative bin monomial')
            census+=1
            if n<=4:
                require(cover_polynomial(a)==literal_cover(a),'literal cover mismatch');literal+=1
                require(bin_polynomial(wrapper(a))==D(p),'physical D mismatch')
                require(bin_polynomial(edge_wrapper(a))==E(p),'physical E mismatch');physical+=2
                h=cs[0];c1=cs[1] if len(cs)>1 else 0
                require(value(D(p),1)==(3*h+6*c1 if n%2 else 5*h+12*c1),'D scalar')
                require(value(E(p),1)==(5*h+18*c1 if n%2 else 11*h+36*c1),'E scalar')
    report.update(physical_profiles=census,literal_covers=literal,wrapper_graph_identities=physical)
    identities=0
    for k in range(21):
        q=basis(k,1);r=basis(k,0)
        for par in (0,1):
            a=k*k+6*k+(3 if par else 5);b=2*k*k*(2*k+1) if par else 2*k*(k+1)*(2*k+1)
            rhs=add(basis(k+1,par),scale(basis(k,par),a),scale(basis(k-1,par),b) if k else [0])
            require(D(basis(k,par))==rhs,('D basis',k,par));identities+=1
        require(E(q)==add(basis(k+1,0),scale(r,k*k+10*k+5),scale(basis(k-1,0),2*k*(2*k+1)**2) if k else [0]),('Eq',k))
        rhs=add(basis(k+2,1),scale(basis(k+1,1),2*k*k+14*k+13),scale(q,(k+1)*(k**3+19*k*k+27*k+11)),scale(basis(k-1,1),2*k*k*(k+1)*(2*k+1)**2) if k else [0])
        require(E(r)==rhs,('Er',k));identities+=2
        require(E(q)==add(D(mul(mon(1),q)),mul(mon(1),D(q)),scale(mul(mon(3),q),-1),scale(mul(mon(1),q),-2)),'operator E decomposition')
    report['operator_basis_controls']=identities
    joins=0
    for i,j,pi,pj in product(range(7),range(7),(0,1),(0,1)):
        require(all(c>=0 for c in coordinates(mul(basis(i,pi),basis(j,pj)),(pi+pj)%2)),'join sign');joins+=1
    report['positive_join_products']=joins
    moments={}
    for par,seed in [(1,(0,)),(0,transitive(2))]:
        p=bin_polynomial(seed);seq=[];a=seed
        for t in range(13):
            seq.append(value(p,1))
            if t<=5:require(h_subsets(a)[-1]==seq[-1],('physical vacuum',par,t));a=wrapper(a)
            p=D(p)
        for N in range(1,7):
            got=det([[seq[i+j] for j in range(N)] for i in range(N)]);want=1
            for k in range(1,N):
                b=2*k*k*(2*k+1) if par else 2*k*(k+1)*(2*k+1);want*=b**(N-k)
            require(got==want,('Hankel product',par,N))
        moments[str(par)]=seq[:8]
    report['vacuum_moments']=moments
    hostile={}
    for n in (4,5):
        a=transitive(n);hs=[]
        for t in range(3):hs.append(h_subsets(a)[-1]);a=wrapper(a)
        hdet=hs[0]*hs[2]-hs[1]**2;require(hdet<0,'missing negative Hankel');hostile[str(n)]=dict(counts=hs,determinant=hdet)
    report['negative_hankel']=hostile
    for m,want in [(6,-8752),(7,-224218032)]:
        p=bin_polynomial(wrapper(transitive(m)));formula=add(mon(m+2),mul(mon(1),add(power([1,1],m),scale(power([-1,1],m),-1))))
        require(p==formula,('clone formula',m));q=p[(m%2)::2]
        while q[0]==0:q=q[1:]
        require(discriminant(q)==want,('discriminant',m,q));hostile['clone_'+str(m)]={'reduced_polynomial':q,'discriminant':want}
    triangle=(2,4,1);body=triangle
    for _ in range(2):body=join(body,triangle)
    physical_cycles=wrapper(body);pp=packing_polynomial(physical_cycles)
    require(pp==[1,507,147,19],('cycle packing hostile',pp))
    require(value(pp,2)==h_subsets(physical_cycles)[-1]==1755,'odd cycle Hamilton readout')
    report['cycle_packing_hostile']={'rows':physical_cycles,'polynomial':pp,'derivative_discriminant':(2*147)**2-4*(3*19)*507}
    h11=edge_wrapper(transitive(2));require(h_subsets(h11)[-1]==11,'physical escape')
    report['edge_wrapper_escape']={'rows':h11,'H':11}
    for target in (13,19):
        sol=[]
        for h in range(1,target+1,2):
            for c1 in range(target+1):
                for form in (3*h+6*c1,5*h+12*c1,5*h+18*c1,11*h+36*c1):
                    if form==target:sol.append((h,c1,form))
        require(not sol,('prime grammar obstruction',target,sol))
    print(json.dumps(report,sort_keys=True,indent=2))
if __name__=='__main__':run()
