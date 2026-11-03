#!/usr/bin/env python3
"""Independent finite controls for the S164 profile/guard algebras."""
from collections import deque
from itertools import product
from math import comb, factorial


def require(c, message):
    if not c: raise RuntimeError(message)


def tournament(n, bits):
    rows = [0] * n
    k = 0
    for i in range(n):
        for j in range(i+1,n):
            if (bits >> k) & 1: rows[i] |= 1 << j
            else: rows[j] |= 1 << i
            k += 1
    return tuple(rows)


def counts(rows):
    """All induced-subgraph H counts and full-graph counts by start."""
    n=len(rows); allmask=(1<<n)-1
    dp=[[0]*n for _ in range(1<<n)]
    for i in range(n): dp[1<<i][i]=1
    for mask in range(1,1<<n):
        for first in range(n):
            if mask>>first&1 and mask!=(1<<first):
                rest=mask^(1<<first)
                dp[mask][first]=sum(dp[rest][j] for j in range(n) if rows[first]>>j&1)
    h=[sum(x) for x in dp];h[0]=1
    return h,dp[allmask]


def induced(rows, vertices):
    return tuple(sum(1<<j for j,v in enumerate(vertices) if rows[u]>>v&1) for u in vertices)


def covers(rows):
    """Unordered covers by choosing the block of the least remaining vertex."""
    n=len(rows);h,_=counts(rows);size=1<<n
    p=[[0]*(n+1) for _ in range(size)];p[0][0]=1
    for mask in range(1,size):
        low=mask&-mask;sub=mask
        while sub:
            if sub&low:
                rest=mask^sub
                for k in range(1,n+1):p[mask][k]+=h[sub]*p[rest][k-1]
            sub=(sub-1)&mask
    return tuple(p[-1])


def starts_cover_poly(rows,q):
    n=len(rows);out=[0]*(n+1);full=(1<<n)-1;sub=full
    while sub:
        if sub>>q&1:
            vv=[i for i in range(n) if sub>>i&1]
            _,hs=counts(induced(rows,vv))
            a=hs[vv.index(q)]
            rest=covers(induced(rows,[i for i in range(n) if not(sub>>i&1)]))
            for k,b in enumerate(rest):out[k+1]+=a*b
        sub=(sub-1)&full
    return tuple(out)


def pointed_poly(rows,q):
    p=covers(rows);s=starts_cover_poly(rows,q)
    incoming=[i for i in range(len(rows)) if rows[i]>>q&1]
    require(len(incoming)==1,'pointed guard indegree')
    vv=[i for i in range(len(rows)) if i!=q]
    d=starts_cover_poly(induced(rows,vv),vv.index(incoming[0]))+(0,)
    return s,tuple(a-b for a,b in zip(p,s)),d


def polyadd(*polys):
    n=max(map(len,polys));return tuple(sum(p[i] if i<len(p) else 0 for p in polys) for i in range(n))


def scale(p,c):return tuple(c*a for a in p)


def number(p):
    out=[0]*(len(p)+1)
    for k,a in enumerate(p):out[k]+=k*a;out[k+1]+=a
    return tuple(out)


def Cpoly(v):
    s,u,d=v
    return polyadd(number(s),scale(s,-1),number(u)),polyadd(s,u,d)+(0,),s+(0,)


def Dpoly(v):
    s,u,d=v
    return polyadd(number(s),d),polyadd(number(u),u),polyadd(number(d),scale(d,-1))


def Gpoly(p):
    n=len(p)+2;s=[0]*n;u=[0]*n;d=[0]*n
    for j,a in enumerate(p):
        # S=z²P+(z+2z²)P'+z²P''; U=z(P+2P'+P''); d=zP.
        s[j+2]+=a;s[j]+=j*a+j*(j-1)*a;s[j+1]+=2*j*a
        u[j+1]+=a;u[j]+=2*j*a
        if j:u[j-1]+=j*(j-1)*a
        d[j+1]+=a
    return tuple(s),tuple(u),tuple(d)


def ordinal(x,y):
    nx,ny=len(x),len(y)
    return tuple([r|(((1<<ny)-1)<<nx) for r in x]+[r<<nx for r in y])


def triangle(x):
    n=len(x);a=n;b=n+1
    return tuple([r|(1<<a) for r in x]+[1<<b,(1<<n)-1])


def star(p,q):
    out=[0]*(len(p)+len(q)-1)
    for k,a in enumerate(p):
        for l,b in enumerate(q):
            for t in range(min(k,l)+1):out[k+l-t]+=a*b*comb(k,t)*comb(l,t)*factorial(t)
    return tuple(out)


def L(p):
    out=[0]*(len(p)+2)
    for j,a in enumerate(p):
        if j:out[j-1]+=a*j*(j-1)
        out[j]+=a*j*(j+2)
        out[j+1]+=a*(2*j+1)
        out[j+2]+=a
    return tuple(out)


def slots(p,r):
    f=1;total=0
    for k,a in enumerate(p):
        if k:f*=r-k+1
        total+=a*f
    return total


def profile_mod(p,m):return tuple(slots(p,r)%m for r in range(m))


def Lmod(x,m):
    return tuple((r*x[(r+1)%m]+r*r*x[r]-r*x[(r-1)%m])%m for r in range(m))


def grammar(m):
    root=tuple(range(m));seen={root:('v',)};queue=deque([root])
    while queue:
        x=queue.popleft();gx=Lmod(x,m)
        if gx not in seen:seen[gx]=('G',x);queue.append(gx)
        for y in list(seen):
            z=tuple(a*b%m for a,b in zip(x,y))
            if z not in seen:seen[z]=('*',x,y);queue.append(z)
    return list(seen),seen


def guarded(rows,q):
    incoming=[v for v in range(len(rows)) if rows[v]>>q&1]
    require(len(incoming)==1,'guard indegree')
    r=incoming[0];h,start=counts(rows)
    vv=[v for v in range(len(rows)) if v!=q]
    _,starts=counts(induced(rows,vv))
    return h[-1],h[-1]-start[q],starts[vv.index(r)]


def attach(rows,q):
    n=len(rows);a=n;b=n+1;out=list(rows)+[0,0]
    out[q]|=1<<a
    for v in range(n):
        if v!=q:out[a]|=1<<v
        out[b]|=1<<v
    out[a]|=1<<b
    return tuple(out),a,b


def one_attach(rows,q):
    n=len(rows);out=list(rows)+[0];out[q]|=1<<n
    for v in range(n):
        if v!=q:out[n]|=1<<v
    return tuple(out),n


def C(v):
    s,u,d=v;return u,s+u+d,s


def D(v):
    s,u,d=v;return s+d,2*u,0


def endpoint_state(rows,q):
    h,u,d=guarded(rows,q);return h-u,u,d


def interval(m):
    n=2*m+2;rows=[sum(1<<j for j in range(i+1,n)) for i in range(n)]
    for i in range(m):
        a,b=2*i,2*i+3;rows[a]^=1<<b;rows[b]|=1<<a
    return tuple(rows)


def A(v):
    h,u,d=v;return 3*h+2*d,2*h+2*d,0


def B(v):
    h,u,d=v;return 3*h+2*d,2*h+d,u


def strong(rows):
    n=len(rows)
    for root in range(n):
        seen=1<<root
        while 1:
            new=seen
            for i in range(n):
                if seen>>i&1:new|=rows[i]
            if new==seen:break
            seen=new
        if seen!=(1<<n)-1:return False
    return True


def main():
    ts=[tournament(n,b) for n in range(1,5) for b in range(1<<(n*(n-1)//2))]
    for x in ts:
        p=covers(x)
        require(covers(triangle(x))==L(p),'triangle cover identity')
        require(pointed_poly(triangle(x),len(x)+1)==Gpoly(p),'positive guarded triangle lift')
        for m in [3,4,6,9,12]:require(profile_mod(L(p),m)==Lmod(profile_mod(p,m),m),'finite profile transfer')
    small=[x for x in ts if len(x)<=3]
    for x,y in product(small,repeat=2):require(covers(ordinal(x,y))==star(covers(x),covers(y)),'ordinal Wick product')
    # All source graphs, including non-strong graphs and two-vertex boundaries.
    guarded_cases=0
    for x in ts:
        for q in range(len(x)):
            if sum(bool(r>>q&1) for r in x)!=1:continue
            v=guarded(x,q);out,a,b=attach(x,q)
            require(guarded(out,a)==A(v),'a guard transfer')
            require(guarded(out,b)==B(v),'b guard transfer')
            if strong(x):require(strong(out),'strength preservation')
            one,newq=one_attach(x,q)
            raw=(v[0]-v[1],v[1],v[2])
            require(endpoint_state(one,newq)==C(raw),'one-vertex C transfer')
            require(endpoint_state(one,q)==D(raw),'one-vertex D transfer')
            pp=pointed_poly(x,q)
            require(pointed_poly(one,newq)==Cpoly(pp),'full pointed C cover transfer')
            require(pointed_poly(one,q)==Dpoly(pp),'full pointed D cover transfer')
            for y in [tournament(1,0),tournament(2,0),tournament(2,1)]:
                py=covers(y)
                require(pointed_poly(ordinal(x,y),q)==tuple(star(a,py) for a in pp),'pointed ordinal module')
            if strong(x):require(strong(one),'one-vertex strength')
            two,secondq=one_attach(one,newq)
            require(two==out,'CC literal F identification')
            guarded_cases+=1
    states,witness=grammar(12);require(len(states)==9,'mod 12 state count')
    require({x[1] for x in states}=={1,3,5,9},'mod 12 scalar residues')
    # Explicit closure certificate, independent of stopping logic of saturation.
    S=set(states)
    require(all(Lmod(x,12) in S for x in S),'unary closure')
    require(all(tuple(a*b%12 for a,b in zip(x,y)) in S for x in S for y in S),'binary closure')
    print('MOD12_STATES =',states)
    print('MOD12_G_TABLE =',[states.index(Lmod(x,12)) for x in states])
    print('MOD12_PRODUCT_TABLE =',[[states.index(tuple(a*b%12 for a,b in zip(x,y))) for y in states] for x in states])
    x=(2,4,1);q=0;v=(3,2,0);seq=[]
    for k in range(7):
        require(guarded(x,q)==v,'iterated b literal profile')
        require(strong(x),'iterated b strength')
        seq.append(v[0])
        if k<6:x,a,q=attach(x,q);v=B(v)
    require(seq==[3,9,31,105,355,1201,4063],'b sequence')
    for i in range(3,len(seq)):require(seq[i]==3*seq[i-1]+seq[i-2]+seq[i-3],'b recurrence')
    for t in range(1,50):
        v=(3*t,2*t,0)
        require(A(v)==tuple(3*c for c in v),'canonical factor 3')
        require(A(A(B(v)))==tuple(31*c for c in v),'canonical factor 31')
    # The S159 interval family is exactly the odd C subsequence, by permutation.
    x=(2,4,1);q=0;guards=[0]
    for j in range(1,20):
        x,q=one_attach(x,q);guards.append(q)
        if j%2:
            m=(j+1)//2
            ordering=[v for i in range(2*m-2,-1,-2) for v in [guards[i],guards[i+1]]]+[1,2]
            require(induced(x,ordering)==interval(m),'general bridge permutation control')
    interval_profiles=[covers(interval(m)) for m in range(5)]
    for r in range(1,5):
        f=[slots(p,r) for p in interval_profiles]
        require(f[:3]==[r*r,r**4+4*r*r,r**6+8*r**4+8*r*r],'interval initial cover profiles')
        for m in range(3,5):
            require(f[m]==(r*r+2)*f[m-1]+(2*r*r-1)*f[m-2]+r*r*f[m-3],'full interval-profile recurrence')
    for t in range(1,30):
        v=(t,2*t,0)
        for n in range(8):
            w=v
            for _ in range(n):w=D(w)
            w=D(C(w))
            require(w==tuple((1+(1<<(n+1)))*c for c in v),'canonical Mersenne-neighbor multiplier')
    for k in range(1,13):
        m=1<<k;seen={1};queue=deque([1])
        while queue:
            a=queue.popleft()
            for b in (3,5):
                c=a*b%m
                if c not in seen:seen.add(c);queue.append(c)
        require(seen==set(range(1,m,2)),'2-adic unit generators 3 and 5')
    # Modular guard reachability is a finite control only, not a universal theorem.
    for m in [3,5,7,9,12,15,16,21,25,27,32]:
        root=tuple(c%m for c in (3,2,0));seen={root};queue=deque([root])
        while queue:
            v=queue.popleft()
            for f in (A,B):
                w=tuple(c%m for c in f(v))
                if w not in seen:seen.add(w);queue.append(w)
        expected=set(range(m)) if m%2 else set(range(1,m,2))
        require({x[0] for x in seen}==expected,'guard finite residue control')
        print('guard modulus',m,'states',len(seen))
    print('PASS:',len(ts),'triangle sources;',len(small)**2,'ordinal pairs;',guarded_cases,'guarded sources with full pointed C/D controls; 330 pointed ordinal grafts; 75 guarded triangle lifts; seven literal b-family controls; mod12 closure; power-of-two multiplier controls.')

if __name__=='__main__':main()
