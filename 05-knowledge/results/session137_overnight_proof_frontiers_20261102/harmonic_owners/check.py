#!/usr/bin/env python3
"""Exact missing-speed chamber and harmonic-boundary verifier; dependency-free."""
from fractions import Fraction as F
from itertools import product
from math import gcd
import hashlib
import json


def require(ok,msg):
    if not ok: raise RuntimeError(msg)


def setup(q,A):
    if type(q) is not int or q<3: raise ValueError('q must be integer >=3')
    A=frozenset(A)
    if any(type(b) is not int or not 1<=b<q for b in A): raise ValueError('invalid attachment')
    d=next((b for b in range(1,q) if b not in A),q)
    U=tuple(u for u in range(1,q) if gcd(u,q)==1)
    E=tuple(u+q*(u in A) for u in U)
    C=tuple(range(1,q))+tuple(q+b for b in sorted(A))
    return A,d,U,E,C


def intersect(X,Y):
    out=[];i=j=0
    while i<len(X) and j<len(Y):
        lo=max(X[i][0],Y[j][0]);hi=min(X[i][1],Y[j][1])
        if lo<=hi: out.append((lo,hi))
        if X[i][1]<Y[j][1]:i+=1
        elif Y[j][1]<X[i][1]:j+=1
        else:i+=1;j+=1
    return out


def windows(V,r):
    out=[(F(0),F(1))]
    for v in sorted(set(V)):
        out=intersect(out,[(F(j,v)+r/v,F(j+1,v)-r/v) for j in range(v)])
        if not out:break
    return out


def safe(V,t,r):
    a,b=t.numerator,t.denominator
    return all(r.denominator*min(a*v%b,-a*v%b)>=r.numerator*b for v in V)


def owners(V,t,r):
    out=[]
    for v in V:
        z=(t*v)%1
        if z==r:out.append((v,'+'))
        if z==1-r:out.append((v,'-'))
    return sorted(out)


def merge(X):
    out=[]
    for l,r in sorted(X):
        if out and l<=out[-1][1]:out[-1]=(out[-1][0],max(out[-1][1],r))
        else:out.append((l,r))
    return out


def oracle(V,r):
    """Independent wall/midpoint decision, not repeated interval intersection."""
    walls=sorted({F(0),F(1)}|{F(j,v)+r/v for v in V for j in range(v)}|
                 {F(j+1,v)-r/v for v in V for j in range(v)})
    atoms=[(x,x) for x in walls if safe(V,x,r)]
    for x,y in zip(walls,walls[1:]):
        if safe(V,(x+y)/2,r):atoms.append((x,y))
    return merge(atoms)


def atlas(q,A,r):
    A,d,U,E,C=setup(q,A)
    eps=F(1,q)-r
    require(F(1,q+d)<=r<=F(1,q),'atlas radius outside closed chamber')
    out=[]
    for a in U:
        p=pow(a,-1,q); ep=p+q*(p in A);em=q-p+q*(q-p in A)
        out.append((F(a,q)-eps/ep,F(a,q)+eps/em))
    return out


def predicted(q,A,r):
    A,d,U,E,C=setup(q,A);n=q+d
    out=atlas(q,A,r)
    if r==F(1,n):out += [(F(b,n),F(b,n)) for b in range(1,n) if gcd(b,n)==1]
    return merge(out)


def endpoint_owners(q,A,r,a,left=True):
    A,d,U,E,C=setup(q,A);p=pow(a,-1,q)
    if r==F(1,q):
        out=[(p,'+'),(q-p,'-')]
        if p in A:out.append((q+p,'+'))
        if q-p in A:out.append((2*q-p,'-'))
    else:
        u=p if left else q-p;sgn='+' if left else '-'
        out=[(u+q*(u in A),sgn)]
        if r==F(1,q+d) and u==d and 2*d in A:out.append((q+2*d,sgn))
    return sorted(out)


def grid_owners(q,A,b):
    A,d,U,E,C=setup(q,A);n=q+d;u=pow(b,-1,n)
    return sorted([(v,'+') for v in (u,u+n) if v in C]+
                  [(v,'-') for v in (n-u,2*n-u) if v in C])


def ceiling_points(q,A,m,h):
    A,d,U,E,C=setup(q,A);n=q+h
    require(1<=h<=d,'h range')
    points=[]
    for e in E:
        a=pow(e,-1,q)
        for j in range(1,m*h//e+1):
            if gcd(j,n)==1:
                t=F(a,q)-F(j,q*m*n)
                points.extend([t,1-t])
    require(len(points)==len(set(points)),'q-window maximum duplicates')
    if h==d and gcd(q*m,n)==1:
        points += [F(b,n) for b in range(1,n) if gcd(b,n)==1]
    return sorted(set(points))


def weight(n,j):
    ell=pow(j,-1,n)
    return F(1,ell)+F(1,n-ell)


def ceiling_mass(q,A,m,h):
    A,d,U,E,C=setup(q,A);n=q+h
    B=F(2,q*m)*sum((weight(n,j) for e in E for j in range(1,m*h//e+1) if gcd(j,n)==1),F(0))
    if h==d and gcd(q*m,n)==1:
        H=sum((F(1,u) for u in range(1,n) if gcd(u,n)==1),F(0))
        B+=F(2,q*m)*(H-weight(n,m))
    return B


def exact_at(V,r):
    """Positive witness + negative certificate and denominator-separation proof."""
    M=max(V);eps=F(1,8*M*M)
    X=windows(V,r)
    require(X and not windows(V,r+eps),('not exact maximum',V,r))
    require(r.denominator<=2*M,'predicted maximum denominator')
    require(safe(V,X[0][0],r),'direct residue witness')
    return X


def main():
    core_cases=owner_cases=oracle_cases=harmonic_cases=ceiling_cases=mass_cases=0
    digest=hashlib.sha256()
    # Exhaustive universe: every A, q=3,...,12; four closed-chamber radii.
    for q in range(3,13):
        for mask in range(1<<(q-1)):
            A={b for b in range(1,q) if mask>>(b-1)&1}
            A,d,U,E,C=setup(q,A);n=q+d;lo=F(1,n);hi=F(1,q)
            for r in (lo,(2*lo+hi)/3,(lo+2*hi)/3,hi):
                X=windows(C,r);P=predicted(q,A,r)
                require(X==P,('closed atlas',q,sorted(A),r,X,P))
                digest.update(repr((q,mask,str(r),[(str(x),str(y)) for x,y in X])).encode())
                for a,(l,u) in zip(U,atlas(q,A,r)):
                    require(owners(C,l,r)==endpoint_owners(q,A,r,a,True),('left owners',q,A,r,a))
                    require(owners(C,u,r)==endpoint_owners(q,A,r,a,False),('right owners',q,A,r,a))
                    owner_cases+=2
                core_cases+=1
                if q<=7 and r in (lo,hi):
                    require(oracle(C,r)==X,('independent wall oracle',q,A,r))
                    oracle_cases+=1
            X=windows(C,lo)
            phi_n=sum(gcd(b,n)==1 for b in range(1,n))
            require(sum(l<u for l,u in X)==len(U),'positive component count')
            require(sum(l==u for l,u in X)==phi_n-2*(gcd(q,d)==1),'isolated count')
            for b in range(1,n):
                if gcd(b,n)==1:
                    require(owners(C,F(b,n),lo)==grid_owners(q,A,b),'all boundary-grid owners')
            # Complete harmonic finite universe, independent exact maxima certificates.
            if q<=7:
                s=min(E)
                for m in range(1,4):
                    N=n//gcd(q*m,n)
                    for k in range(1,n):
                        V=C+tuple(q*m*j for j in range(1,k+1))
                        if d*m>=s:
                            r=min(F(m,q*m+s),F(1,k+1))
                            exact_at(V,r)
                        elif k<N:
                            exact_at(V,lo)
                        else:
                            require(not windows(V,lo),('strictly below floor',q,A,m,k))
                        harmonic_cases+=1
                    for h in range(1,d+1):
                        V=C+tuple(q*m*j for j in range(1,q+h))
                        r=F(1,q+h);X=windows(V,r)
                        P=ceiling_points(q,A,m,h)
                        require(X==[(t,t) for t in P],('ceiling maximizers',q,A,m,h,X,P))
                        require(bool(P)==(m*h>=s),'ceiling preservation threshold')
                        B=ceiling_mass(q,A,m,h)
                        direct=F(0)
                        for t in P:
                            O=owners(V,t,r)
                            plus=max(v for v,sgn in O if sgn=='+')
                            minus=max(v for v,sgn in O if sgn=='-')
                            direct+=F(1,plus)+F(1,minus)
                        require(B==direct,'all-owner harmonic terminal mass')
                        if P and q<=6 and m<=2:
                            eta=F(1,64*max(V)**3)
                            mass=sum((u-l for l,u in windows(V,r-eta)),F(0))
                            require(mass==eta*B,'independent below-maximum measure slope')
                            mass_cases+=1
                        if h==d:
                            base=2*sum(sum(gcd(j,n)==1 for j in range(1,m*d//e+1)) for e in E)
                            count=base+(phi_n-2)*(gcd(q*m,n)==1)
                            require(len(P)==count,'boundary maximum count')
                        ceiling_cases+=1
    # Larger q: all prefix promotions, boundary plus two interior probes.
    prefix_cases=0
    for q in range(13,41):
        for h in range(q):
            A=set(range(1,h+1));A,d,U,E,C=setup(q,A)
            for r in (F(1,q+d),(F(1,q+d)+F(1,q))/2,F(1,q)):
                require(windows(C,r)==predicted(q,A,r),'larger prefix atlas')
                prefix_cases+=1
    # Named hostile controls and serialized all-owner exact-max certificates.
    controls=[]
    for name,q,A,m,k in [
        ('isolated_floor_maxima',3,{1},1,4),
        ('same_sign_boundary_tie',5,{1,4},1,6),
        ('deeper_prefix_gain',6,{1,2,3,4},2,6),
        ('deeper_prefix_floor',6,{1,2,3,4},1,7),
        ('fully_promoted_below_old_floor',3,{1,2},1,3),
        ('nonunit_first_gap_strict_below',6,{1},1,4),
    ]:
        A,d,U,E,C=setup(q,A);V=C+tuple(q*m*j for j in range(1,k+1))
        # Exact bisection + denominator recovery only for named controls.
        lo=F(0);hi=F(1,2);M=max(V)
        for _ in range((32*M*M).bit_length()+2):
            mid=(lo+hi)/2
            if windows(V,mid):lo=mid
            else:hi=mid
        r=((lo+hi)/2).limit_denominator(2*M)
        X=exact_at(V,r)
        controls.append({'name':name,'q':q,'A':sorted(A),'d':d,'m':m,'k':k,
            'maximum':str(r),'maximum_denominator_bound':2*M,
            'empty_above':str(r+F(1,8*M*M)),
            'components':[{'left':str(l),'right':str(u),
                'left_owners':owners(V,l,r),'right_owners':owners(V,u,r)} for l,u in X]})
    q=4;A={1};B={1,2};r=F(1,6)
    require(setup(q,A)[3]==setup(q,B)[3],'same owner spectrum hostile')
    require(safe(setup(q,A)[4],r,r) and not safe(setup(q,B)[4],r,r),'nonunit blindness hostile')
    require(not any(l<=F(1,5)<=u for l,u in atlas(3,{1},F(1,5))), 'old atlas excludes genuine isolated point')
    require(safe(setup(3,{1})[4],F(1,5),F(1,5)),'isolated witness survives')
    out={'status':'PASS / FINITE-EXACT; universal claims have separate proofs; LRC(14) OPEN',
         'exhaustive_q':[3,12],'exhaustive_attachment_sets':4092,
         'closed_atlas_cases':core_cases,'endpoint_owner_checks':owner_cases,
         'independent_wall_oracle_cases':oracle_cases,'exhaustive_digest_sha256':digest.hexdigest(),
         'larger_prefix_q':[13,40],'larger_prefix_cases':prefix_cases,
         'harmonic_universe':'q=3..7, every A, m=1..3, every k=1..q+d-1',
         'harmonic_exact_or_below_floor_cases':harmonic_cases,
         'complete_ceiling_maximizer_cases':ceiling_cases,
         'all_owner_terminal_mass_cases':ceiling_cases,'independent_mass_slope_cases':mass_cases,
         'controls':controls,
         'same_owner_nonunit_hostile':{'q':4,'A':[1],'B':[1,2],'owners':list(setup(q,A)[3]),
                                      'point':'1/6','radius':'1/6','A_safe':True,'B_safe':False}}
    print(json.dumps(out,indent=2,sort_keys=True))

if __name__=='__main__':main()
