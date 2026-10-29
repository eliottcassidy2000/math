#!/usr/bin/env python3
"""Exact, bounded checks for the affine first-extinction recurrence packet.

No assertions: checks also run under python -O. Standard library only.
"""
from collections import defaultdict
from fractions import Fraction as Q
from itertools import combinations, permutations
from decimal import Decimal, localcontext
from math import comb, factorial
import json

LIMIT = 16
CHECKS = 0


def eq(a, b, label):
    global CHECKS
    CHECKS += 1
    if a != b:
        raise RuntimeError((label, a, b))


def clean(p):
    return {k: v for k, v in p.items() if v}


def add(*ps):
    out = defaultdict(int)
    for p in ps:
        for k, v in p.items():
            out[k] += v
    return clean(out)


def scale(p, c):
    return clean({k: c*v for k, v in p.items()})


def shift(p, a=1, b=1):
    return { (j,d): v*a**j*b**d for (j,d),v in p.items() }


def mul(p, q):
    out = defaultdict(int)
    for (j,d), v in p.items():
        for (k,e), w in q.items():
            if d+e <= LIMIT:
                out[j+k,d+e] += v*w
    return clean(out)


ONE = {(0,0): 1}


def power(p, n):
    out = ONE
    for _ in range(n):
        out = mul(out, p)
    return out


def resolve(numerator, coefficient):
    out, term = {}, numerator
    while term:
        out = add(out, term)
        term = mul(term, coefficient)
    return out


def concat(A, B, r):
    return add(*(scale(mul(shift(A[i], b=2**(r-i)),
                           shift(B[r-i], a=3**i)), comb(r,i))
                 for i in range(r+1)))


def triple(F, r, mixed=False):
    out = {}
    for i in range(r+1):
        for j in range(r-i+1):
            k = r-i-j
            if mixed and max(i,j,k) == r:
                continue
            term = mul(mul(shift(F[i], b=2**(j+k)),
                           shift(F[j], a=3**i, b=2**k)),
                       shift(F[k], a=3**(i+j)))
            out = add(out, scale(term, factorial(r)//factorial(i)//factorial(j)//factorial(k)))
    return out


def trees(d=3, r=2, selected=None):
    """All chronological first-hit words with valuation clock <= LIMIT."""
    def visit(word, queue, clock, carry):
        if queue == 0:
            yield word, clock, carry
            return
        for a in range(1, LIMIT-clock+1):
            branch=(a%r==0) if selected is None else ((a-1)%r+1 in selected)
            nextq = queue-1+(d if branch else 0)
            # Every pending vertex costs at least one more unit of clock.
            if nextq > LIMIT-clock-a:
                continue
            yield from visit(word+(a,), nextq, clock+a, 3*carry+2**clock)
    yield from visit((), 1, 0, 0)


def parse(word, d=3, r=2, selected=None):
    def take(i):
        start = i
        a = word[i]
        i += 1
        children = []
        branch=(a%r==0) if selected is None else ((a-1)%r+1 in selected)
        if branch:
            for _ in range(d):
                i, child = take(i)
                children.append(child)
        return i, (word[start:i], children)
    end, tree = take(0)
    eq(end, len(word), 'unique full parse')
    return tree


def affine(word):
    c, clock = 0, 0
    for a in word:
        c, clock = 3*c+2**clock, clock+a
    return len(word), clock, c


def native_guard(word, clock, carry):
    modulus = 2**(clock+1)
    n = ((2**clock-carry)*pow(3**len(word), -1, modulus)) % modulus
    for source in (n, n+modulus):
        x = source
        for a in word:
            x = 3*x+1
            got = (x & -x).bit_length()-1
            eq(got, a, 'native valuation')
            x >>= a
        eq(x, (3**len(word)*source+carry)//2**clock, 'native affine endpoint')


def univariate_mul(a,b,n):
    out = [Q(0)]*(n+1)
    for i,x in enumerate(a):
        for j,y in enumerate(b[:n+1-i]):
            out[i+j] += x*y
    return out


def univariate_pow(a,d,n):
    out = [Q(1)]+[Q(0)]*n
    for _ in range(d):
        out = univariate_mul(out,a,n)
    return out


def coefficients(d,r,n=12,probs=None):
    p, q = (Q(1,2**r-1), Q(3,4**r-1)) if probs is None else probs
    M, L = [1-p], [1-q]
    for m in range(1,n+1):
        M.append(p*univariate_pow(M,d,m-1)[m-1])
        L.append(q*univariate_pow(L,d,m-1)[m-1])
    for m in range(n+1):
        catalan = Q(comb(d*m+1,m),d*m+1)
        eq(M[m], catalan*(1-p)**((d-1)*m+1)*p**m, 'tree mass closed form')
        eq(L[m], catalan*(1-q)**((d-1)*m+1)*q**m, 'slope mass closed form')
    S = [Q(0)]*(n+1)
    for i in range(d):
        term = univariate_mul(univariate_pow(L,i,n),univariate_pow(M,d-1-i,n),n)
        S = [x+y for x,y in zip(S,term)]
    B = []
    for m in range(n+1):
        B.append(L[m]/3 + (p*sum(S[j]*B[m-1-j] for j in range(m)) if m else 0))
    return M,L,B


def shape_moments(d,r,m,probs=None):
    """Independent ordered-shape sum; each infinite letter-gap sum is exact."""
    p,q = (Q(1,2**r-1),Q(3,4**r-1)) if probs is None else probs
    mass,slope,carry = Q(0),Q(0),Q(0)
    size=d*m+1
    for chosen in combinations(range(size),m):
        branch=set(chosen); queue=1
        valid=True
        for i in range(size):
            queue+=d*(i in branch)-1
            if queue==0 and i+1<size:
                valid=False; break
        if not valid:
            continue
        eq(queue,0,'shape final queue')
        ps=[p if i in branch else 1-p for i in range(size)]
        ls=[q if i in branch else 1-q for i in range(size)]
        suffix=[Q(1)]*(size+1)
        for i in range(size-1,-1,-1):
            suffix[i]=suffix[i+1]*ls[i]
        prefix=Q(1)
        for i in range(size):
            carry+=prefix*suffix[i]/3
            prefix*=ps[i]
        mass+=prefix; slope+=suffix[0]
    return mass,slope,carry


def abelian_probes():
    out={}
    with localcontext() as ctx:
        ctx.prec=70
        for r in (2,3,4):
            d=2**r-1
            p,q=Decimal(1)/d,Decimal(3)/(4**r-1)
            values=[]
            for exponent in (4,8,12):
                delta=Decimal(10)**(-exponent); t=1-delta
                roots=[]
                for b in (p,q):
                    lo,hi=Decimal(0),Decimal(1)
                    for _ in range(230):
                        mid=(lo+hi)/2
                        if 1-b+t*b*mid**d-mid>0:
                            lo=mid
                        else:
                            hi=mid
                    roots.append((lo+hi)/2)
                M,L=roots
                S=sum(L**i*M**(d-1-i) for i in range(d))
                scaled=delta.sqrt()*(L/3)/(1-t*p*S)
                target=(Decimal(2*d)/(d-1)).sqrt()/3
                values.append({'delta':str(delta),'scaled_carry':str(scaled),
                               'target':str(target)})
                eq(abs(scaled-target)<10*delta.sqrt(),True,'Abelian numeric control')
            out[str(d)]=values
    return out


def main():
    # Rational supersolutions, independent of any floating approximation.
    t,ms,ls=Q(4,3),Q(6,5),Q(21,20)
    eq(Q(6,7)+t*ms**3/7<=ms,True,'mass supersolution')
    eq(Q(20,21)+t*ls**3/21<=ls,True,'slope supersolution')
    eq(1-t*(ls*ls+ls*ms+ms*ms)/7,Q(193,700),'carry denominator')
    eq((ls/3)/Q(193,700),Q(245,193),'carry exponential budget')
    eq(2*Q(3,4)**3,Q(27,32),'source refinement factor')
    eq(2*243**2-343**2,449,'two-branch coarse budget hostile')
    odd_arities={}
    for d in range(3,102,2):
        h=1
        while pow(2,h,d)!=1:
            h+=1
        numerator=(2**h-1)//d
        selected={j for j in range(1,h+1) if numerator & (1<<(h-j))}
        p=Q(sum(2**(h-j) for j in selected),2**h-1)
        q=Q(3*sum(4**(h-j) for j in selected),4**h-1)
        eq(p,Q(1,d),'critical selector source mass')
        eq(p*p<q<3*p*p,True,'base-four tilt bounds')
        eq(d*q<Q(3,4),True,'all odd tilted subcriticality')
        odd_arities[str(d)]={'period':h,'selector':sorted(selected),'q':str(q)}
        if d in (5,9):
            M,L,B=coefficients(d,h,3,(p,q))
            for m in range(4):
                eq(shape_moments(d,h,m,(p,q)),(M[m],L[m],B[m]),'multiresidue carry recurrence')
            counts=defaultdict(int)
            for word,clock,carry in trees(d,h,selected):
                parse(word,d,h,selected); native_guard(word,clock,carry)
                counts[len(word),clock]+=1
            lp={(1,a):1 for a in range(1,LIMIT+1) if (a-1)%h+1 not in selected}
            ep={(1,a):1 for a in range(1,LIMIT+1) if (a-1)%h+1 in selected}
            z={}
            while True:
                nz=add(lp,mul(ep,power(z,d)))
                if nz==z:
                    break
                z=nz
            eq(z,dict(counts),'multiresidue full clock grammar')
    for r in range(1,9):
        d=2**r+1; h=2*r
        q=Q(3*sum(4**(h-j) for j in range(r+1,h+1)),4**h-1)
        eq(q,Q(1,4**r+1),'plus-one arity tilt')
    b=Q(9*(8+1),8**4-1); a=Q(9,7)-b
    eq((a,b),(Q(576,455),Q(9,455)),'five-child square-slope weights')
    eq(a+b*2**5,Q(864,455),'five-child square-slope supersolution')
    eq(a+b*2**5<2,True,'finite slope second moment control')
    leaf = {(1,a):1 for a in range(1,LIMIT+1,2)}
    branch = {(1,a):1 for a in range(2,LIMIT+1,2)}
    f0 = {}
    while True:
        nxt = add(leaf,mul(branch,power(f0,3)))
        if nxt == f0:
            break
        f0 = nxt
    F = [f0]
    for r in range(1,5):
        P, Qp = shift(f0,a=3**r), shift(f0,b=2**r)
        S = add(power(P,2),mul(P,Qp),power(Qp,2))
        mixed = triple(F,r,True)
        numerator = add(leaf,mul(shift(branch,b=2**r),mixed))
        for i in range(1,r+1):
            numerator = add(numerator,scale(mul(shift(branch,b=2**(r-i)),
                                        shift(triple(F,r-i),a=3**i)),comb(r,i)))
        F.append(resolve(numerator,mul(shift(branch,b=2**r),S)))
    direct = [defaultdict(int) for _ in range(5)]
    word_count = 0
    collision = defaultdict(list)
    for word, clock, carry in trees():
        word_count += 1
        parse(word)
        native_guard(word,clock,carry)
        for r in range(5):
            direct[r][len(word),clock] += carry**r
        collision[len(word),clock].append((word,carry))
        for cut in range(len(word)+1):
            ju,du,cu = affine(word[:cut])
            jv,dv,cv = affine(word[cut:])
            eq(3**jv*cu+2**du*cv,carry,'concatenation carry')
    for r in range(5):
        eq(F[r],dict(direct[r]),'jet '+str(r))
        H = [triple(F,j) for j in range(r+1)]
        eq(F[r],add(leaf,concat([branch]*(r+1),H,r)),'grammar '+str(r))
    A=[{(1,1):1,(1,4):1}]*5
    V=[{(1,2):1,(1,3):1}]*5
    W=[{(1,1):1,(1,2):1}]*5
    AV=[concat(A,V,r) for r in range(5)]
    VW=[concat(V,W,r) for r in range(5)]
    for r in range(5):
        eq(concat(AV,W,r),concat(A,VW,r),'shifted jet associativity')
    P,Qp = shift(f0,a=3),shift(f0,b=2)
    eq(F[1],resolve(scale(P,Q(1,3)),mul(shift(branch,b=2),
              add(power(P,2),mul(P,Qp),power(Qp,2)))),'compact first jet')

    # Different independent mechanism: sum over tree size after all gap values.
    coefficient_rows = {}
    for d,r in ((3,2),(3,3),(7,3),(2,2),(1,2),(15,4)):
        M,L,B = coefficients(d,r)
        for m in range(4 if d>3 else 6):
            eq(shape_moments(d,r,m),(M[m],L[m],B[m]),'independent ordered shape moments')
        coefficient_rows[f'{d},{r}'] = {'mass':list(map(str,M[:5])),
            'slope':list(map(str,L[:5])), 'intercept':list(map(str,B[:5]))}
        # A second bivariate grammar, checked against chronological enumeration.
        lp = {(1,a):1 for a in range(1,LIMIT+1) if a%r}
        ep = {(1,a):1 for a in range(1,LIMIT+1) if a%r==0}
        z = {}
        while True:
            nz = add(lp,mul(ep,power(z,d)))
            if nz == z:
                break
            z = nz
        counts = defaultdict(int)
        for word,clock,carry in trees(d,r):
            counts[len(word),clock] += 1
            parse(word,d,r)
        eq(z,dict(counts),'general arity grammar')
    # Direct conditional source-cylinder sums for the cheaper rule.
    # Each prefix guard intersection is evaluated exactly, keeping D and C.
    cylinders={k:defaultdict(Q) for k in range(7)}
    for word,clock,carry in trees(3,3):
        m=(len(word)-1)//3
        guard=((2**clock-carry)*pow(3**len(word),-1,2**(clock+1)))%2**(clock+1)
        value=Q(carry,2**clock)*Q(4,3)**m
        for k,rows in cylinders.items():
            if clock>=k:
                rows[guard%2**(k+1)]+=Q(2**k,2**clock)*value
            else:
                for residue in range(guard,2**(k+1),2**(clock+1)):
                    rows[residue]+=value
    for k,rows in cylinders.items():
        for residue in range(1,2**(k+1),2):
            eq(rows[residue]<=2**k*Q(245,193),True,'fixed source-cylinder carry bound')
    for m in range(13):
        M,L,_ = coefficients(3,2)
        eq(L[m]/M[m],Q(6,5)*Q(108,125)**m,'conditional slope')

    # Exact historical hostile controls and the surviving corrected statements.
    for n in range(1,1001):
        a,carry,count = n,0,0
        while a or carry:
            digit=a%3
            carry=(2*digit+carry)//3
            count += carry
            a//=3
        c=comb(2*n,n); val=0
        while c%3==0:
            val+=1; c//=3
        eq(count,val,'full ternary carry chain')
    eq(comb(10,5)%9,0,'n=5 has two carries')
    tournament_tests=0
    for n in range(1,8):
        for perm in permutations(range(n)):
            order=sorted(range(n),key=perm.__getitem__)
            for i in range(n):
                for j in range(i+1,n):
                    eq(perm[order[i]]<perm[order[j]],True,'value-order tournament')
            tournament_tests+=1
    examples = []
    for jd, vals in sorted(collision.items(),key=lambda x:(x[0][1],x[0][0])):
        if len({c for _,c in vals})>1:
            first=vals[0]
            other=next(v for v in vals if v[1]!=first[1])
            examples.append({'J,D':jd,'first':first,'second':other})
            break
    abelian=abelian_probes()
    print(json.dumps({'status':'FINITE-EXACT (Decimal Abelian probes are VERIFIED)','clock_limit':LIMIT,'jet_orders':list(range(5)),
        'critical_words':word_count,'checks':CHECKS,'permutation_controls':tournament_tests,
        'carry_collision':examples,'size_coefficients':coefficient_rows,'abelian_probes':abelian,
        'odd_arity_selectors':odd_arities},indent=2))


if __name__ == '__main__':
    main()
