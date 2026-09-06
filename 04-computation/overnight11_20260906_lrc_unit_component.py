"""Exact controls for the actual-entry primitive-unit-component closure.

Standard library only. No producer imports and no large-height census.
"""
from fractions import Fraction as F
from itertools import combinations
from functools import reduce
from math import gcd, prod
import sys

sys.stdout.reconfigure(newline="\n")
Q=91**6
gates=0


def need(ok,label):
    global gates
    gates+=1
    if not ok:
        raise ArithmeticError(label)


def norm(w,x):
    r=(w*x.numerator)%x.denominator
    return F(min(r,x.denominator-r),x.denominator)


def ceil(x):
    return -((-x.numerator)//x.denominator)


def atlas(p,q):
    if not 1<=p<q or gcd(p,q)!=1 or p+q>356:
        return False
    m=p+q
    prime=2
    while prime*prime<=m:
        if m%prime==0:
            exponent=0
            while m%prime==0:
                exponent+=1
                m//=prime
            if prime%3!=2 or exponent>2:
                return False
        prime+=1
    return m==1 or m%3==2


def bounded_division(x,K,bound):
    a=min(x//K,bound)
    r=x-a*K
    return a,r


def pair_phase(p,q):
    m=p+q
    z=F((pow(p,-1,m)*(m//2))%m,m)
    need(min(norm(p,z),norm(q,z))==F(m//2,m)>=F(1,3), "explicit safe pair phase")
    return z


def grid_witness(V,t,g,p,q,y):
    K=max(V)
    need(len(V)==11 and len(set(V))==11 and min(V)==1, "eleven distinct primitive unit-core speeds")
    need(min(norm(v,y) for v in V)>=F(1,12), "literal supplier control at core threshold")
    need(gcd(t,g)==1 and atlas(p,q), "coprime scales and actual pair atlas")
    S=tuple(t*v for v in V)+(g*p,g*q)
    need(len(set(S))==13 and reduce(gcd,S)==1, "actual distinct primitive physical row")
    z=pair_phase(p,q)
    if 11*t>=21*q:
        left=z-F(11,42*q)
        right=z+F(11,42*q)
        need(0<left<right<1 and right-left>=F(1,t), "pair safe arc covers one t-grid spacing")
        k=ceil(t*left-g*y)
        j=(pow(g,-1,t)*k)%t
        x=(y+j)/t
        target=(g*y+k)/t
        need(left<=target<=right and g*x-target==int(g*x-target), "literal transported pair clock")
        branch="large core scale"
    else:
        delta=gcd(t,p)
        need(t<=677<Q and g*p//delta>Q*(K+1), "small-scale survivor inequality")
        need(g>42*K, "scale forced above complete core-grid threshold")
        left=y-F(1,84*K)
        right=y+F(1,84*K)
        k=ceil(g*left-t*z)
        if (t*z+k)/g==left:
            k+=1
        j=(pow(t,-1,g)*k)%g
        x=(z+j)/g
        target=(t*z+k)/g
        need(left<target<right and t*x-target==int(t*x-target), "literal transported primitive core clock")
        branch="forced large pair scale"
    clearance=min(norm(w,x) for w in S)
    need(clearance>=F(1,14), "literal simultaneous physical witness")
    return (K,t,g,p,q,branch,str(x),str(clearance))


def main():
    toy=0
    for bound in range(2,14):
        for K in range(2,bound+1):
            for x in range(1,bound*(K+1)+1):
                a,r=bounded_division(x,K,bound)
                need(x==a*K+r and 0<=a<=bound and 0<=r<=bound, "inclusive bounded integer compiler")
                toy+=1
            x=bound*(K+1)+1
            # Even allowing signed coefficients and a larger external
            # multiplier cannot exceed the exact box radius Q(K+1).
            need(all(abs(a*K+r)<x for a in range(-bound,bound+1)
                     for r in range(-bound,bound+1)), "sharp first gap of the two-core-coordinate box")
    generalized=0
    for bound in (4,7,11):
        for K in range(2,bound+1):
            for t in sorted({1,2,3,bound}):
                for p in (1,2,3,5,7):
                    for g in range(1,18):
                        if gcd(t,g)!=1:
                            continue
                        delta=gcd(t,p)
                        x=g*p//delta
                        if x>bound*(K+1):
                            continue
                        a,r=bounded_division(x,K,bound)
                        c=t//delta
                        need(c*g*p-a*t*K-r*t==0 and 1<=c<=bound
                             and max(a,r)<=bound, "cleared primitive-scale crossing relation")
                        generalized+=1
    pairs=[(p,q) for q in range(2,356) for p in range(1,min(q,357-q)) if atlas(p,q)]
    need(len(pairs)==5855, "inherited all-scale inert pair atlas")
    need(max(p for p,q in pairs)==177 and max(q for p,q in pairs)==355, "sharp smaller/larger coordinate caps")
    for p,q in pairs:
        need(Q>42*p, "uniform pair-to-core gluing inequality")
        need((21*q-1)//11<=677<Q, "strict residual core-scale coefficient bound")
    # Boundary and huge-integer controls retain the inclusive endpoint.
    for K in (11,22,Q):
        a,r=bounded_division(Q*(K+1),K,Q)
        need((a,r)==(Q,Q), "inclusive endpoint requires non-Euclidean final remainder Q")
    # A genuinely typed inherited nonvacuity control, THM4049/4117.
    U=(1,4,6,8,10,12,14,15,16,18,22)
    g=2**45
    S=U+(g,3*g)
    edges=[(i,j) for i,j in combinations(range(13),2)
           if atlas(min(S[i],S[j])//gcd(S[i],S[j]),max(S[i],S[j])//gcd(S[i],S[j]))]
    parent=list(range(13))
    def find(i):
        while parent[i]!=i:
            i=parent[i]
        return i
    for i,j in edges:
        parent[find(i)]=find(j)
    parts={}
    for i in range(13):
        parts.setdefault(find(i),[]).append(i)
    need(sorted(parts.values(),key=len)==[[11,12],list(range(11))], "literal decoder 11+2 graph of the canonical unit-core row")
    need(sum(S)<=Q**2 and g>2*Q*max(U), "canonical finite box and support-three crossing exclusion by dominance")
    # If two core coordinates occur, their total magnitude is <=2QK<g.
    # If both pair coordinates occur, their weighted sum is a multiple of g;
    # a nonzero core coefficient has magnitude <=QK<g. Thus no crossing row.
    canonical=grid_witness(U,1,g,1,3,F(9,19))
    # Primitive gcd one does not imply a literal unit or K<=Q.
    primes=(37,43,61,67,73,79,97,103,127)
    P=15*prod(primes)
    hostile=tuple(2*P//r for r in primes)+(P//3,P//5)
    need(reduce(gcd,hostile)==1 and min(hostile)>1 and max(hostile)>Q,
         "canonical primitive nonunit shape defeats unlicensed normalization")
    need(max(max(a,b)//gcd(a,b) for a,b in combinations(hostile,2))==127,
         "large primitive nonunit core still has all small internal pair heights")
    controls=[]
    cores=[(tuple(range(1,12)),F(1,12)),(U,F(9,19)),(tuple(range(1,22,2)),F(1,2))]
    for V,y in cores:
        for p,q in ((1,3),(1,4),(1,355),(177,178)):
            t=(21*q+10)//11
            controls.append(grid_witness(V,t,1,p,q,y))
            for t in (1,2,17,677):
                if 11*t>=21*q:
                    continue
                delta=gcd(t,p)
                g=delta*Q*(max(V)+1)//p+1
                while gcd(t,g)!=1:
                    g+=1
                controls.append(grid_witness(V,t,g,p,q,y))
    # d=3 consumer conversion: K=max(3M,E), g=3h.
    for K in (11,42*95,Q):
        for p in (1,17,177):
            h=Q*(K+1)//(3*p)+1
            need(14*h>=29*K, "unit-core crossing survivor enters both THM4448 d3 inequalities")
    print("STATUS: PASS; scoped actual-entry primitive-unit eleven-component closure")
    print("BOUNDED DIVISION:",toy,"inclusive controls;",generalized,"scale-cleared crossing controls")
    print("ATLAS:",len(pairs),"pairs; p<=177,q<=355; residual t<=677")
    print("NONVACUOUS ACTUAL ENTRY CONTROL:",canonical)
    print("GENERAL CLOCK WITNESSES:",len(controls),"exact physical controls in both branches")
    print("NONUNIT HOSTILE: primitive K=",max(hostile),"but all internal pair heights<=127")
    print("CONCLUSION: literal unit in primitive eleven-component shape closes the actual W=Vdec branch")
    print("SCOPE: no arbitrary gcd-one core closure; no decoder entry inferred from this finite bank")
    print("ACTIVE GATES:",gates)


if __name__=="__main__":
    main()
