#!/usr/bin/env python3
"""Exact response controls for extreme triangles and two-root readiness."""
from functools import lru_cache
from hashlib import sha256
from itertools import permutations
from pathlib import Path
from fractions import Fraction
import json

HERE = Path(__file__).resolve().parent
ROOT = HERE.parents[3]
C = (2,4,1)

def require(ok, message):
    if not ok:
        raise RuntimeError(message)

def rows_from_mask(n, mask):
    rows = [0]*n
    bit = 0
    for i in range(n):
        for j in range(i+1,n):
            if mask >> bit & 1:
                rows[i] |= 1 << j
            else:
                rows[j] |= 1 << i
            bit += 1
    return tuple(rows)

def transitive(n):
    return tuple(sum(1 << j for j in range(i+1,n)) for i in range(n))

def ordinal(A,B):
    a,b = len(A),len(B)
    return tuple(x | (((1 << b)-1) << a) for x in A)+tuple(x << a for x in B)

def birth(T,r):
    n=len(T)
    return tuple(x | ((1 << n) if i==r else 0) for i,x in enumerate(T))+(((1<<n)-1)^(1<<r),)

def induced(T,keep):
    return tuple(sum(((T[u] >> v)&1) << j for j,v in enumerate(keep)) for u in keep)

def skew_det(T):
    n=len(T)
    A=[[Fraction(0 if i==j else 1 if T[i] >> j & 1 else -1)
        for j in range(n)] for i in range(n)]
    result=Fraction(1)
    for j in range(n):
        pivot=next((i for i in range(j,n) if A[i][j]),None)
        if pivot is None: return 0
        if pivot!=j:
            A[j],A[pivot]=A[pivot],A[j]
            result=-result
        v=A[j][j];result*=v
        for i in range(j+1,n):
            q=A[i][j]/v
            for k in range(j,n):
                A[i][k]-=q*A[j][k]
    require(result.denominator==1,'integral determinant')
    return int(result)

def reverse(T):
    n=len(T)
    return tuple(((1<<n)-1)^x^(1<<i) for i,x in enumerate(T))

@lru_cache(None)
def data(T):
    n=len(T)
    ends=[[0]*n for _ in range(1<<n)]
    for v in range(n):
        ends[1<<v][v]=1
    for mask in range(1,1<<n):
        for v in range(n):
            if not (mask >> v & 1) or mask == 1 << v:
                continue
            rest=mask^(1<<v)
            ends[mask][v]=sum(ends[rest][u] for u in range(n) if T[u] >> v & 1)
    counts=[sum(row) for row in ends]
    counts[0]=1
    return counts,ends

def hp(T):
    return data(T)[0][-1]

def profile(T,r):
    counts,ends=data(T)
    starts=data(reverse(T))[1]
    full=(1<<len(T))-1
    F=sum(ends[A][r]*counts[full^A] for A in range(full+1) if A >> r & 1)
    return (counts[full],F,starts[full][r],counts[full^(1<<r)])

def polynomials(T,mmax):
    counts,ends=data(T)
    starts=data(reverse(T))[1]
    full=(1<<len(T))-1
    @lru_cache(None)
    def P(mask,m):
        if not m:
            return int(mask==0)
        total=0
        sub=mask
        while True:
            total += counts[sub]*P(mask^sub,m-1)
            if sub==0: return total
            sub=(sub-1)&mask
    Pfull=[P(full,m) for m in range(mmax+1)]
    E=[];S=[]
    for m in range(1,mmax+1):
        E.append([sum(ends[A][r]*P(full^A,m-1) for A in range(full+1) if A >> r & 1) for r in range(len(T))])
        S.append([sum(starts[A][r]*P(full^A,m-1) for A in range(full+1) if A >> r & 1) for r in range(len(T))])
    return Pfull,E,S

def literal_hp(T):
    return sum(all(T[p[i]] >> p[i+1] & 1 for i in range(len(T)-1)) for p in permutations(range(len(T))))

def gadget(m,I):
    # a=0,b=1,c=2, r_i=i+3; a,b dominate R.
    T=list(ordinal(C,transitive(m)))
    for j in range(m):
        r=j+3
        if I >> j & 1:
            T[2] &= ~(1<<r)
            T[r] |= 1<<2
    return tuple(T)

def main():
    controls={'residual_inputs':0,'polynomial_equalities':0,'triangle_births':0,
              'first_birth_dichotomies':0,'literal_path_counts':0,'two_root_cuts':0}
    examples=[]
    for n in range(5):
        for mask in range(1 << (n*(n-1)//2)):
            K=rows_from_mask(n,mask)
            h=hp(K)
            front=ordinal(C,K)
            back=ordinal(K,C)
            pk,_,_=polynomials(K,3)
            pf,ef,sf=polynomials(front,3)
            pb,eb,sb=polynomials(back,3)
            controls['residual_inputs']+=1
            for m in range(4):
                require(pf[m]==pb[m]==(m**3+2*m)*pk[m],('P identity',n,mask,m))
                controls['polynomial_equalities']+=2
            for m in range(1,4):
                u=m*m-m+1
                for r in range(3):
                    require(sf[m-1][r]==eb[m-1][n+r]==u*pk[m],('entry identity',n,mask,m,r))
                    require(ef[m-1][r]==sb[m-1][n+r]==u*pk[m-1],('exit identity',n,mask,m,r))
                    controls['polynomial_equalities']+=4
            for r in range(3):
                require(profile(front,r)==(3*h,3*h,h,h),('front profile',n,mask,r))
                first=birth(front,r)
                second=birth(first,len(front))
                require(profile(first,len(front))[:3]==(5*h,6*h,2*h),'first response')
                require(profile(second,len(first))[:3]==(9*h,11*h,3*h),'second response')
                controls['triangle_births']+=2
            for r in range(len(front)):
                value=hp(birth(front,r))
                require(value==5*h if r<3 else value%3==0,('first-birth dichotomy',n,mask,r))
                controls['first_birth_dichotomies']+=1
            if n<=2:
                for T in (front,back,birth(front,0)):
                    require(literal_hp(T)==hp(T),'literal paths')
                    controls['literal_path_counts']+=1
            # Keep caches bounded rather than accumulating every labelled atom.
            data.cache_clear()
    ready=[]
    for m in range(1,8):
        for I in range(1,1<<m):
            T=gadget(m,I)
            Ha,Fa,sa,qa=profile(T,0)
            H,F,s,q=profile(T,1)
            require(Ha==H and Fa==H and sa==qa,'root a zero')
            b=F-H
            initial = (I & (I+1))==0
            j=I.bit_length()
            require((b==1)==initial,('ready1 classification',m,I,b))
            require((b==2)==(I==2),('ready2 classification',m,I,b))
            if b==1:
                require(H==1+(1<<(j+1)) and H%7!=0,'ready1 count')
            if b==2:
                require(H==9,'ready2 count')
            if b<=2:
                ready.append({'m':m,'I':I,'H':H,'defect':b,'start':s,'deletion':q})
            controls['two_root_cuts']+=1
            if m<=3:
                require(literal_hp(T)==H,'literal gadget')
                controls['literal_path_counts']+=1
            data.cache_clear()
    for m,I in ((1,1),(2,2),(3,7)):
        T=gadget(m,I)
        H,F,s,q=profile(T,1)
        child=birth(birth(T,1),len(T))
        value=hp(child)
        require(value==3*H+2*(F-H),'actual nonzero branch')
        indegree_one=[r for r in range(len(child))
                      if sum(bool(child[v] >> r & 1) for v in range(len(child)))==1]
        require(indegree_one==[len(T),len(T)+1],'two newborn indegree-one roots')
        require(profile(child,len(T))[1]-value==0,'first newborn zero defect')
        require(profile(child,len(T)+1)[1]-value==H-s>=3,'latest newborn defect')
        examples.append({'m':m,'I':I,'rows':T,'profile_H_F_s_q':[H,F,s,q],
                         'child_rows':child,'child_H':value,
                         'child_defects':[profile(child,r)[1]-value for r in range(len(child))]})
        data.cache_clear()
    staircase=[]
    for j in range(2,9):
        T=gadget(j,(1<<j)-1)
        keep=[0,3,2]+list(range(4,j+3))
        child=induced(T,keep)
        require(child==gadget(j-1,(1<<(j-1))-1),'literal deletion staircase')
        require(hp(T)==1+(1<<(j+1)),'staircase source count')
        require(hp(child)==(hp(T)+1)//2,'guarded even branch')
        require(hp(T)%8==1,'staircase branch guard')
        staircase.append({'j':j,'rows':T,'retained_order':keep,'child_rows':child,
                          'H':hp(T),'child_H':hp(child)})
        data.cache_clear()
    base=gadget(1,1)
    nine=birth(base,1)
    seventeen=birth(nine,len(base))
    require(base==(10,12,1,4) and nine==(10,28,1,4,13)
            and seventeen==(10,28,1,4,45,15),'cycle literal rows')
    require((hp(base),hp(nine),hp(seventeen))==(5,9,17),'cycle counts')
    require(induced(seventeen,range(5))==nine and induced(nine,range(4))==base,'exact cycle deletions')
    canonical=gadget(3,7)
    require(hp(canonical)==17 and skew_det(canonical)==1 and skew_det(seventeen)==49,'H17 chart determinant contrast')
    cycle={'base_rows':base,'first_birth_rows':nine,'second_birth_rows':seventeen,
           'arithmetic_count_cycle':[5,17,9,5],'retained_base_names':[0,1,2,3],
           'canonical_H17_rows':canonical,'canonical_H17_determinant':1,
           'transient_H17_determinant':49}
    source_names=[
        '05-knowledge/results/triples_ports_sieve_20261112/ports/FOUR_RETURNS.md',
        '05-knowledge/results/polarized_path_calculus_20261103/FUNNELS.md',
        '05-knowledge/results/two_entry_funnels_20261103/READY.md',
        '05-knowledge/results/core_descent_endpoints_20261112/README.md',
        '01-canon/theorems/THM-1370-h-spectrum-omits-7-21-all-n.md',
    ]
    print(json.dumps({'status':'PASS; FINITE-EXACT controls',
                      'checks':controls,'ready_cuts':ready,'actual_branches':examples,
                      'collatz_staircase':staircase,'collatz_cycle':cycle,
                      'source_sha256':{p:sha256((ROOT/p).read_bytes()).hexdigest() for p in source_names},
                      'script_sha256':sha256(Path(__file__).read_bytes()).hexdigest(),
                      'scope':'Exact declared small universes; no unrestricted return or spectrum classification'},
                     indent=2,sort_keys=True))

if __name__=='__main__':
    main()
