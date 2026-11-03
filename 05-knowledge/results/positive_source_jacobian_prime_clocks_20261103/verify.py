"""Exact finite controls for the accompanying all-order proofs. No asserts.

Run: python -B verify.py > verification.json (also run with -O).
Universe: all labelled tournaments n<=5; all recursive-tree birth words n<=8;
all even moduli 4<=m<=20000. No random or numerical-rank tests.
"""
from collections import Counter
from functools import lru_cache
from itertools import combinations, permutations
from fractions import Fraction as F
from math import gcd, factorial
from pathlib import Path
import bisect
import json
import sympy as s


def require(ok, msg):
    if not ok:
        raise RuntimeError(msg)


def rows(n, code):
    out = [0]*n
    for k, (i,j) in enumerate(combinations(range(n),2)):
        a,b = (i,j) if code>>k&1 else (j,i)
        out[a] |= 1<<b
    return tuple(out)


def birth(a, x):
    n = len(a)
    return tuple(b | ((1<<n) if i != x else 0) for i,b in enumerate(a))+(1<<x,)


def sink(a):
    return tuple(b | (1<<len(a)) for b in a)+(0,)


def profile(a):
    n=len(a); full=(1<<n)-1
    dp=[[0]*n for _ in range(1<<n)]
    hp=[0]*(1<<n)
    for i in range(n): dp[1<<i][i]=1
    for mask in range(1,full+1):
        for i in range(n):
            if dp[mask][i]:
                avail=a[i]&~mask
                while avail:
                    bit=avail&-avail; avail-=bit; j=bit.bit_length()-1
                    dp[mask|bit][j]+=dp[mask][i]
        hp[mask]=sum(dp[mask])
    @lru_cache(None)
    def covers(mask):
        if not mask: return (1,)+(0,)*n
        z=[0]*(n+1); least=mask&-mask; sub=mask
        while sub:
            if sub&least:
                rest=covers(mask^sub)
                for r in range(1,n+1): z[r]+=hp[sub]*rest[r-1]
            sub=(sub-1)&mask
        return tuple(z)
    return covers(full)


def independent_profile(a):
    """Permute, cut, divide by r! to forget the order of the path blocks."""
    n=len(a); counts=[0]*(n+1)
    for perm in permutations(range(n)):
        bad=sum(1<<i for i in range(n-1) if not (a[perm[i]]>>perm[i+1]&1))
        for cuts in range(1<<(n-1)):
            if cuts&bad==bad: counts[1+cuts.bit_count()]+=1
    return tuple(counts[r]//factorial(r) for r in range(n+1))


def transfer(c,n):
    at=lambda r:c[r] if 0<=r<len(c) else 0
    return tuple(F(at(r-1))+r*at(r)+F(r*(r+1),n)*at(r+1) for r in range(n+2))


def tournament_checks():
    checked=0; children=0
    for n in range(1,6):
        for code in range(1<<(n*(n-1)//2)):
            a=rows(n,code); c=profile(a)
            require(c==independent_profile(a), ('profile',n,code))
            sums=[0]*(n+2)
            for x in range(n):
                b=profile(birth(a,x)); children+=1
                sums=[i+j for i,j in zip(sums,b)]
            require(tuple(F(z,n) for z in sums)==transfer(c,n),('birth',n,code))
            sc=profile(sink(a))
            require(all(sc[r]==(c[r-1] if r else 0)+(r*c[r] if r<=n else 0)
                        for r in range(n+2)),('sink',n,code))
            if n==4: require(c==(0,c[1],c[1]+6,6,1),('four',code))
            checked+=1
    family=[(0,)]; moments=[]; aggregate=(0,1)
    for n in range(1,9):
        dist=Counter(profile(a)[1] for a in family)
        mass=sum(dist.values()); mean=F(sum(h*c for h,c in dist.items()),mass)
        second=F(sum(h*h*c for h,c in dist.items()),mass)
        require(mass==factorial(n-1),('word mass',n))
        require(mean==F(aggregate[1],mass),('first moment recurrence',n))
        moments.append(dict(n=n,words=mass,mean=str(mean),variance=str(second-mean*mean)))
        if n==4: require(dist==Counter({1:1,3:1,5:4}),'prefix law')
        if n<8: family=[birth(a,x) for a in family for x in range(n)]
        aggregate=tuple(n*v for v in transfer(aggregate,n))
    t=s.symbols('t'); A=(0,1); hs=[]; sturm=[]
    for n in range(1,21):
        require(all(z.denominator==1 if isinstance(z,F) else True for z in A),'integrality')
        hs.append(int(A[1]))
        if n<=12:
            poly=s.Poly(sum(int(a)*t**r for r,a in enumerate(A)),t)
            require(poly.count_roots(-s.oo,0)==n,('real roots',n))
            sturm.append(n)
        A=tuple(n*v for v in transfer(A,n))
    d=(0,F(1),F(1),F(0),F(0)); alphas=[]
    for n in range(4,21):
        alphas.append(dict(n=n,alpha=str(d[1]),variance_lower=str(F(7,3)*d[1]**2)))
        require(n==4 or F(7,3)*d[1]**2>=F(21,4),'variance lower bound')
        d=transfer(d,n)
    return dict(all_tournaments=checked,birth_children=children,profile_algorithms=2,
                tree_moments=moments,aggregate_H=hs,sturm_orders=sturm,amplification=alphas)


def jacobian_checks():
    x,y,z,v,t,u=s.symbols('x y z v t u')
    R=2*x-3*x*x*y-x**3*z
    Q=y+3*x*(1+x*y)**2*z+3*x*y*y*(4+3*x*y)
    P=(1+x*y)**3*z+y*y*(1+x*y)*(4+3*x*y)
    K=s.Matrix([R,Q,P])
    require(s.expand(K.jacobian([x,y,z]).det())==2,'K determinant')
    ell=2-3*v-t
    A=v+3*(1+v)**2*t+3*v*v*(4+3*v)
    B=(1+v)**3*t+v*v*(1+v)*(4+3*v)
    V=s.expand(ell*A); T=s.expand(ell**2*B)
    require(s.expand(s.det(s.Matrix([V,T]).jacobian([v,t]))-2*ell**2)==0,'quotient det')
    sub={v:V,t:T}
    l2=s.expand(ell*ell.subs(sub,simultaneous=True))
    A2=s.cancel(A.subs(sub,simultaneous=True)/ell)
    B2=s.cancel(B.subs(sub,simultaneous=True)/ell**2)
    a2=s.Poly(A2,v,t); b2=s.Poly(B2,v,t)
    require(all(i+2*j>=1 for (i,j),c in a2.terms()),'A2 integrality')
    require(all(i+2*j>=2 for (i,j),c in b2.terms()),'B2 integrality')
    bracket=lambda f,g:s.diff(f,v)*s.diff(g,t)-s.diff(f,t)*s.diff(g,v)
    require(s.expand(l2*bracket(A2,B2)+A2*bracket(l2,B2)+2*B2*bracket(A2,l2))==4,
            'composed twisted equation')
    require(s.Poly(l2,v,t).total_degree()==7,'nonlinear orbit')
    q=(0,0,-s.Rational(1,4))
    pts=[q,(1,-s.Rational(3,2),s.Rational(13,2)),(-1,s.Rational(3,2),s.Rational(13,2))]
    for p in pts: require(tuple(K.subs(dict(zip((x,y,z),p))))==q,('collision',p))
    leading=lambda f:sum(c*v**i*t**j for (i,j),c in s.Poly(f,v,t).terms()
                          if i+j==s.Poly(f,v,t).total_degree())
    require(s.expand(leading(V)+3*(v*(3*v+t))**2)==0,'V leading')
    require(s.expand(leading(T)-(v*(3*v+t))**3)==0,'T leading')
    # Independent exact specialization over a finite field avoids huge coefficients.
    # Its degree is a lower bound; symbolic degree propagation supplies the upper bound.
    mod=1009
    vv=s.Poly(u,u,modulus=mod); tt=vv; L=s.Poly(1,u,modulus=mod)
    one=s.Poly(1,u,modulus=mod); specialized=[]
    for k in range(1,6):
        ll=2-3*vv-tt
        L=L*ll
        delta=(3*5**(k-1)-1)//2
        require(L.degree()==delta,('degree specialization',k))
        specialized.append(dict(k=k,orbit_degree=delta,generic_degree=3**k,determinant=2**k))
        if k<5:
            aa=vv+3*(one+vv)**2*tt+3*vv**2*(4+3*vv)
            bb=(one+vv)**3*tt+vv**2*(one+vv)*(4+3*vv)
            vv,tt=ll*aa,ll**2*bb
    return dict(det_K=2,det_quotient='2*(2-3v-t)^2',l2_degree=7,
                l2_terms=len(s.Poly(l2,v,t).terms()),A2_degree=a2.total_degree(),
                B2_degree=b2.total_degree(),twisted_constant=4,collision=[list(map(str,p)) for p in pts],
                iterate_degree_controls=specialized)


def prime_checks():
    limit=2000000; isprime=bytearray(b'\1')*(limit+1); isprime[:2]=b'\0\0'
    for p in range(2,int(limit**0.5)+1):
        if isprime[p]: isprime[p*p::p]=b'\0'*(((limit-p*p)//p)+1)
    primes=[p for p in range(2,limit+1) if isprime[p]]
    def onset(m):
        k=1
        while (k+1)*m<=limit:
            for sign in [-1,1]:
                p=k*m+sign
                if isprime[p]:
                    factors=s.factorint(m); count=bisect.bisect_right(primes,p)-len(factors)
                    return dict(m=m,k=k,sign=sign,p=p,eligible=count,phi=int(s.totient(m)),
                                normalized=str(F(2*count,int(s.totient(m)))))
            k+=1
        raise RuntimeError(('sieve too short',m))
    expected=[(262,263,54),(278,277,57),(292,293,60),(120,239,49),
              (144,431,81),(208,1039,173),(290,1451,227)]
    table=[]
    for m,p,count in expected:
        row=onset(m); require((row['p'],row['eligible'])==(p,count),('table',m))
        explicit=[q for q in primes if q<=p and gcd(q,m)==1]
        require(len(explicit)==count and explicit[-1]==p,('independent eligible',m))
        require(all(q%m not in [1,m-1] for q in explicit[:-1]),('first',m))
        table.append(row)
    records=[]; best=F(0)
    for m in range(4,20001,2):
        row=onset(m); ratio=F(row['normalized'])
        if ratio>best: records.append(row); best=ratio
    factors={str(n):{str(p):int(e) for p,e in s.factorint(n).items()}
             for n in [289,291,579,581,869,871,1159,1161,1449,1451]}
    roots={str(sign):[a for a in range(290) if a*a%290==sign%290] for sign in [-1,1]}
    require(roots['-1']==[17,133,157,273],'roots -1')
    gaps=[]
    for g in [(4,4,4),(1,1,10)]:
        N=sum(g); targets=[0]; run=0
        for a in g[:-1]: run+=a; targets.append(run)
        waits=[next(j+1 for j in range(N) if (start+j)%N in targets) for start in range(N)]
        mean=F(sum(waits),N); formula=F(sum(a*a for a in g)+N,2*N)
        require(mean==formula,('waiting',g))
        gaps.append(dict(gaps=g,return_mean=str(F(N,len(g))),wait_mean=str(mean)))
    from sympy.ntheory.modular import crt
    L=6; qs=list(s.primerange(7,60))[:2*L]
    residues=[]; covers=[]
    for k in range(1,L+1):
        for sign in [-1,1]:
            q=qs[len(residues)]; residues.append((-sign*pow(k,-1,q))%q)
            covers.append((k,sign,q))
    base,period=map(int,crt([2]+qs,[0]+residues)); m=base+period
    for k,sign,q in covers: require((k*m+sign)%q==0 and k*m+sign>q,'CRT composite')
    return dict(table=table,universe_even_moduli=[4,20000],sieve_limit=limit,
                normalized_records=records,factors_290=factors,roots_290=roots,
                cyclic_wait_controls=gaps,CRT_control=dict(L=L,m=m,period=period,divisors=covers))


if __name__=='__main__':
    result=dict(tournaments=tournament_checks(),jacobian=jacobian_checks(),primes=prime_checks())
    print(json.dumps(result,indent=2,sort_keys=True))
