"""Exact finite controls for reciprocal frames, prime clocks, and return Grams.

Decimal evaluations are explicitly numerical, not minimal-polynomial proofs.
All mathematical assertions remain active under python -O.
"""
from pathlib import Path
from fractions import Fraction as Q
from itertools import product
from math import isqrt,gcd,comb,lcm,prod
from decimal import Decimal as D,localcontext
import json
HERE=Path(__file__).resolve().parent
CHECKS=0
def need(v,m):
    global CHECKS
    CHECKS+=1
    if not v:raise RuntimeError(m)

def lucas(n):
    a,b=2,1
    for _ in range(n):a,b=b,a+b
    return a
def fib(n):
    a,b=0,1
    for _ in range(n):a,b=b,a+b
    return a
def sieve(n):
    a=bytearray(b'\1')*(n+1);a[:2]=b'\0\0'
    for p in range(2,isqrt(n)+1):
        if a[p]:a[p*p:n+1:p]=b'\0'*((n-p*p)//p+1)
    return [i for i in range(2,n+1) if a[i]]
def leg(a,p):
    v=pow(a%p,(p-1)//2,p)
    return -1 if v==p-1 else v
def matmul(a,b,p=None):
    x=(a[0]*b[0]+a[1]*b[2],a[0]*b[1]+a[1]*b[3],
       a[2]*b[0]+a[3]*b[2],a[2]*b[1]+a[3]*b[3])
    return tuple(v%p for v in x) if p else x
def matpow(a,k,p=None):
    out=(1,0,0,1)
    while k:
        if k&1:out=matmul(out,a,p)
        a=matmul(a,a,p);k//=2
    return out
def mobius(a,x,p):
    if x==p:num,den=a[0],a[2]
    else:num,den=(a[0]*x+a[1])%p,(a[2]*x+a[3])%p
    return num*pow(den,-1,p)%p if den else p
def cycles(perm):
    unseen=set(range(len(perm)));out=[]
    while unseen:
        start=min(unseen);c=[];v=start
        while v in unseen:unseen.remove(v);c.append(v);v=perm[v]
        need(v==start,'permutation cycle closes');out.append(c)
    return out
def order_scalar_matrix(a,p):
    b=(1,0,0,1)
    for k in range(1,2*(p+1)+1):
        b=matmul(b,a,p)
        if b[1]==b[2]==0 and b[0]==b[3]:return k
    raise RuntimeError(('projective order bound',p))
def golden_mul(x,y,p):
    a,b=x;c,d=y
    return ((a*c+b*d)%p,(a*d+b*c+b*d)%p)
def golden_pow(x,k,p):
    out=(1,0)
    while k:
        if k&1:out=golden_mul(out,x,p)
        x=golden_mul(x,x,p);k//=2
    return out

def prime_clocks():
    records=[];incidences=0;bound=1000;index_bound=300
    M=(1,1,1,0)
    for p in sieve(bound):
        if p in (2,5):continue
        eps=leg(5,p);d=order_scalar_matrix(M,p)
        need((p-eps)%d==0,('torus order',p,d))
        # alpha=phi/phi'=-phi^2=-1-phi, represented in the full etale algebra.
        alpha=(p-1,p-1)
        need(golden_pow(alpha,d,p)==(1,0),'projective/torus order')
        for q in sieve(d):
            if d%q==0:need(golden_pow(alpha,d//q,p)!=(1,0),'order exact')
        perm=[mobius(M,x,p) for x in range(p+1)];cs=cycles(perm)
        fixed=[c[0] for c in cs if len(c)==1]
        need(len(fixed)==1+eps,'split fixed-point count')
        need(all(len(c) in (1,d) for c in cs),'uniform nonfixed clock cycles')
        a,b=2,1
        for n in range(index_bound+1):
            divides=a%p==0
            prediction=d%2==0 and n%d==d//2
            need(divides==prediction,('Lucas progression',p,n,d))
            if divides:
                incidences+=1
                A=matpow(M,n,p)
                need((A[0]+A[3])%p==0,'trace-zero half turn')
                need(order_scalar_matrix(A,p)==2,'half-turn order')
                if n%2:need(eps==1,'odd Lucas divisor splits')
            a,b=b,a+b
        records.append({'prime':p,'split':eps,'projective_order':d,
                        'first_Lucas_index':d//2 if d%2==0 else None,
                        'fixed_points':fixed,'cycles':len(cs)})
    # Exact complete trial factorization, rather than probable-prime labels.
    primes=sieve(isqrt(lucas(60))+1);factor_rows=[]
    for n in range(2,61):
        value=lucas(n);rem=value;factors=[]
        for p in primes:
            if p*p>rem:break
            if rem%p:continue
            e=0
            while rem%p==0:rem//=p;e+=1
            factors.append([p,e])
        if rem>1:factors.append([rem,1])
        factor_product=1
        for p,e in factors:factor_product*=p**e
        need(factor_product==value,'full factor product')
        isprime=len(factors)==1 and factors[0][1]==1
        if isprime:need(n in sieve(n) or n&(n-1)==0,'prime index necessary gate')
        if n%2 and n in sieve(n):
            for p,e in factors:
                if p!=2:need(p%(2*n)==1 and leg(5,p)==1,'prime-index factor restriction')
        if n&(n-1)==0:
            for p,e in factors:
                eps=leg(5,p)
                need(p%(4*n)==(1 if eps==1 else 2*n-1),'dyadic split factor classes')
        factor_rows.append({'index':n,'value':value,'factors':factors,'prime':isprime})
    # The power/valuation law, checked only where its odd-multiple hypothesis holds.
    valuation_cases=0
    for r in records:
        p=r['prime'];base=r['first_Lucas_index']
        if base is None or base>40:continue
        def valuation(x):
            e=0
            while x%p==0:x//=p;e+=1
            return e
        for multiplier in (1,3,5,7,9,p):
            need(valuation(lucas(base*multiplier))==valuation(lucas(base))+valuation(multiplier),
                 ('odd multiplier valuation',p,base,multiplier))
            valuation_cases+=1
    r6=(0,6,1,0);half=matpow(M,5,11)
    need(half==(8,5,5,3),'literal fifth power')
    fixed_r6=[x for x in range(12) if mobius(r6,x,11)==x]
    fixed_half=[x for x in range(12) if mobius(half,x,11)==x]
    need(fixed_r6==[] and fixed_half==[4,8],'nonconjugate eleven involutions')
    reciprocity=[]
    for p in sieve(200):
        if p in (2,3):continue
        cs=cycles([mobius(r6,x,p) for x in range(p+1)])
        need(sum(len(c)==1 for c in cs)==1+leg(6,p),'reciprocal involution fixed points')
        need(len(cs)==(p+2+leg(6,p))//2,'reciprocal orbit count')
        reciprocity.append([p,leg(6,p),leg(-6,p),len(cs)])
    # A small finite-wheel hostile which also passes the prime-index gate.
    need(lucas(23)==139*461 and all(lucas(23)%p for p in sieve(100)),
         'L23 evades every prime <=100 but is composite')
    clocks={r['prime']:(r['projective_order'] if r['projective_order']%2==0 else 1)
            for r in records};clocks.update({2:3,4:6,5:1})
    banks=[]
    for bank in ((2,3,11),(2,5,7),(4,2,5,3,7,13),tuple(sieve(29))):
        periods=[clocks[p] for p in bank];period=lcm(*periods)
        labels=[];a,b=2,1
        for n in range(2*period):
            labels.append(tuple(a%p==0 for p in bank));a,b=b,a+b
        need(labels[:period]==labels[period:],'full joint output repeats')
        candidates=[d for d in range(1,period+1) if period%d==0]
        least=next(d for d in candidates if all(labels[i]==labels[i%d] for i in range(period)))
        need(least==period,'minimal shared phase from direct recurrence')
        banks.append({'bank':bank,'individual_periods':periods,'minimal_states':period,
                      'independent_product_states':prod(periods),'compatible_fraction_denominator':prod(periods)//period})
    # Distinct admissible indices have disjoint prime support. A fixed bank
    # can therefore hit only finitely many of these candidates.
    candidate_indices=sorted(set(sieve(index_bound))|{2**k for k in range(1,10) if 2**k<=index_bound})
    candidate_values={n:lucas(n) for n in candidate_indices}
    for i,n in enumerate(candidate_indices):
        for m in candidate_indices[:i]:
            need(gcd(candidate_values[n],candidate_values[m])==1,'coprime Lucas candidates')
    bank_candidate_hits=[]
    for bank in ((2,3,11),(2,5,7),(4,2,5,3,7,13),tuple(sieve(29))):
        hits={p:[n for n in candidate_indices if candidate_values[n]%p==0] for p in bank}
        need(all(len(indices)<=1 for indices in hits.values()),'one candidate per fixed divisor')
        bank_candidate_hits.append({'bank':bank,'hits':hits})
    gcd_cases=0
    values={n:lucas(n) for n in range(1,121)}
    for n,value in values.items():
        e=0
        while value%2==0:value//=2;e+=1
        need(e==(0 if n%3 else (2 if n%2 else 1)),'two-adic Lucas valuation')
    for n in range(1,121):
        for m in range(1,n+1):
            g=gcd(m,n)
            prediction=values[g] if (m//g)%2 and (n//g)%2 else (2 if g%3==0 else 1)
            need(gcd(values[m],values[n])==prediction,'exact Lucas gcd composition')
            gcd_cases+=1
    return {'prime_bound':bound,'index_bound':index_bound,'records':records,
            'incidences':incidences,'full_factorizations':factor_rows,
            'valuation_cases':valuation_cases,'reciprocal_prime_table':reciprocity,
            'eleven_halfturn':list(half),'eleven_fixed_sets':[fixed_r6,fixed_half],
            'minimal_joint_banks':banks,'gcd_cases':gcd_cases,
            'candidate_indices':candidate_indices,'bank_candidate_hits':bank_candidate_hits}

def trace_frames():
    def trace(u):return u+1/u
    counts=[]
    for k in range(1,9):
        # Path and heap trees are two different LCA charts on the same exact fiber.
        values=list(map(Q,range(2,k+2)));edge_sets=[[(i-1,i) for i in range(1,k)],
                                                [((i-1)//2,i) for i in range(1,k)]]
        for edges in edge_sets:
            signatures={}
            for signs in product((-1,1),repeat=k):
                lifted=[u**e for u,e in zip(values,signs)]
                sig=tuple(trace(lifted[i]*lifted[j]) for i,j in edges)
                orbit=min(signs,tuple(-e for e in signs))
                if sig in signatures:need(signatures[sig]==orbit,'tree traces recover relative branches')
                signatures[sig]=orbit
            need(len(signatures)==2**(k-1),'exact global inversion quotient')
        counts.append([k,2**(k-1)])
    cut_codes=[]
    for k in range(2,7):
        edges=[(i,j) for i in range(k) for j in range(i+1,k)]
        words={tuple(s[i]^s[j] for i,j in edges) for s in product((0,1),repeat=k)}
        need(len(words)==2**(k-1),'cut-code dimension')
        for word in words:
            e=dict(zip(edges,word))
            need(all(e[i,j]^e[j,l]^e[i,l]==0 for i in range(k) for j in range(i+1,k)
                     for l in range(j+1,k)),'cycle realization guard')
        cut_codes.append([k,len(edges),len(words)])
    for a in range(2,9):
        for b in range(2,9):
            u,v=Q(a),Q(b);x,y=trace(u),trace(v)
            plus,minus=trace(u*v),trace(u/v)
            need(plus+minus==x*y and plus*minus==x*x+y*y-4,'Fricke pair polynomial')
            need((plus-minus)**2==(x*x-4)*(y*y-4),'branch discriminant')
    need(18+3==3*7 and 18*3==3**2+7**2-4,'golden trace hostile')
    gaussian65=[]
    for s,t in product((-1,1),repeat=2):
        x,y=6-2*s*t,4*t+3*s
        legs=sorted((abs(x*x-y*y),abs(2*x*y)))
        need(x*x+y*y==65 and sum(v*v for v in legs)==65**2,'Gaussian 65 allocation')
        need(legs==([33,56] if s==t else [16,63]),'relative conjugation branch preserved')
        gaussian65.append({'signs':[s,t],'Gaussian_product':[x,y],'PPT':legs+[65]})
    return {'tree_fibers':counts,'binary_cut_codes':cut_codes,
            'golden_same_marginals':[3,7],'possible_joint_traces':[18,3],
            'Gaussian65':gaussian65}

def reversal_holonomy():
    def compose(a,b):return tuple(a[b[i]] for i in range(9))
    def inverse(a):return tuple(a.index(i) for i in range(9))
    ident=tuple(range(9));t0=tuple(3*((-i)%3)+(-j)%3 for i in range(3) for j in range(3))
    group=[tuple(3*((i+a)%3)+(j+bs[i])%3 for i in range(3) for j in range(3))
           for a in range(3) for bs in product(range(3),repeat=3)]
    def edge(x,y):
        i,j=divmod(x,3);k,l=divmod(y,3)
        return (l-j)%3==1 if i==k else (k-i)%3==1
    for g in group:
        need(all(edge(x,y)==edge(g[x],g[y]) for x in range(9) for y in range(9) if x!=y),'actual wreath automorphism')
    marks=[compose(g,t0) for g in group if compose(compose(g,t0),compose(g,t0))==ident]
    central=[g for g in group if compose(g,t0)==compose(t0,g)]
    need(len(group)==81 and len(marks)==27 and len(central)==3,'wreath reversal control')
    def sqrt_odd(a):
        p=ident;powers=[]
        while True:
            powers.append(p);p=compose(a,p)
            if p==ident:break
        need(len(powers)%2==1,'odd-order square root')
        return powers[(len(powers)+1)//2%len(powers)]
    def transport(t,s):
        g=sqrt_odd(compose(t,s))
        need(compose(compose(g,s),inverse(g))==t,'transport realizes marked relation')
        return g
    frequencies={};witness=None
    for t in marks:
        for u in marks:
            h=compose(inverse(transport(u,t0)),compose(transport(u,t),transport(t,t0)))
            need(h in central,'loop remains in actual centralizer')
            frequencies[h]=frequencies.get(h,0)+1
            if h!=ident and witness is None:witness={'middle':t,'last':u,'holonomy':h}
    need(set(frequencies)==set(central),'triangles realize the full C3 residual')
    orbitset={tuple(sorted({g[x] for g in central})) for x in range(9)}
    need(sorted(map(len,orbitset))==[1,1,1,3,3],'maximal flat vertex quotient')
    return {'automorphisms':len(group),'reversal_markings':len(marks),
            'centralizer':central,'based_triangle_counts':[{'permutation':h,'count':v} for h,v in sorted(frequencies.items())],
            'maximal_path_independent_vertex_orbits':sorted(orbitset),'witness':witness}

def sparse_apply(A,v):
    out={}
    for i,row in enumerate(A):
        s=sum(a*v.get(j,0) for j,a in row.items())
        if s:out[i]=s
    return out
def rank_mod(A,p=1000003):
    rows=[[r.get(j,0)%p for j in range(len(A))] for r in A];rank=0
    for j in range(len(A)):
        pos=next((i for i in range(rank,len(rows)) if rows[i][j]),None)
        if pos is None:continue
        rows[rank],rows[pos]=rows[pos],rows[rank]
        inv=pow(rows[rank][j],-1,p);rows[rank]=[(v*inv)%p for v in rows[rank]]
        for i in range(rank+1,len(rows)):
            c=rows[i][j]
            if c:rows[i]=[(x-c*y)%p for x,y in zip(rows[i],rows[rank])]
        rank+=1
    return rank
def stars_and_lucas_graphs():
    stars=[]
    for m in range(1,10):
        states=list(range(1<<m))+[1<<m];ids={s:i for i,s in enumerate(states)}
        forbidden={}
        for s in states:
            for j in range(m+1):
                t=s^(1<<j)
                if t not in ids:forbidden.setdefault(t,[]).append(ids[s])
        A=[{} for s in states]
        for pred in forbidden.values():
            for i in pred:
                for j in pred:A[i][j]=A[i].get(j,0)+1
        # Check the full polynomial identity on every coordinate vector.
        for j in range(len(states)):
            v={j:1};Av=sparse_apply(A,v);A2v=sparse_apply(A,Av);A3v=sparse_apply(A,A2v)
            need(all(A3v.get(i,0)-(m+2)*A2v.get(i,0)+(m+1)*Av.get(i,0)==0
                     for i in set(A3v)|set(A2v)|set(Av)),'native Gram cubic identity')
        need(sum(row.get(i,0) for i,row in enumerate(A))==2**m+m-1,'native Gram trace')
        if m<=6:need(rank_mod(A)==2**m-1,'independent modular rank lower bound')
        stars.append({'leaves':m,'dimension':len(states),'rank':2**m-1,
                      'eigenvalue_multiplicities':[[0,2],[1,2**m-2],[m+1,1]]})
    graphs=[]
    for n in range(3,14):
        states=[s for s in range(1<<n) if all(not((s>>i&1) and(s>>((i+1)%n)&1)) for i in range(n))]
        bank=set(states);edges=[(s,s|(1<<i)) for s in states for i in range(n)
                              if not(s>>i&1) and s|(1<<i) in bank]
        need(len(states)==lucas(n),'Lucas word count')
        need(len(edges)==n*fib(n-1),'all-size edge formula control')
        degree={s:0 for s in states}
        for a,b in edges:degree[a]+=1;degree[b]+=1
        need([s for s in states if degree[s]==n]==[0],'empty vertex intrinsically fixed')
        need((len(edges)%2==1)==(n%6 in (3,5)),'six-class signed cancellation')
        graphs.append({'period':n,'vertices':len(states),'edges':len(edges),
                       'forced_self_converse_zero':len(edges)%2==1})
    return {'stars':stars,'lucas_graphs':graphs}

def forms(discriminant):
    need(discriminant<0 and discriminant%4 in (0,1),'negative quadratic discriminant')
    out=[]
    for a in range(1,isqrt((-discriminant)//3)+1):
        for b in range(-a,a+1):
            if (b*b-discriminant)%(4*a):continue
            c=(b*b-discriminant)//(4*a)
            if a>c or gcd(gcd(a,abs(b)),c)>1:continue
            if (abs(b)==a or a==c) and b<0:continue
            out.append([a,b,c])
    return out
def trace_recurrence(t,norm,k):
    a,b=2,t
    for _ in range(k):a,b=b,t*b-norm*a
    return a
def series_controls():
    degree=6
    def mul(a,b):
        out=[0]*(degree+1)
        for i,x in enumerate(a):
            for j,y in enumerate(b):
                if i+j<=degree:out[i+j]+=x*y
        return out
    def inv(a):
        need(a[0]==1,'unit power series');out=[1]+[0]*degree
        for k in range(1,degree+1):out[k]=-sum(a[j]*out[k-j] for j in range(1,min(k,len(a)-1)+1))
        return out
    records={}
    for sign in (-1,1):
        a=[1]+[0]*degree
        for step in range(1,degree+1,2):
            factor=[0]*(degree+1)
            for k in range(degree//step+1):factor[k*step]=comb(24,k)*sign**k
            a=mul(a,factor)
        b=inv(a)
        # q*(64 g^24+64 g^-24)=a+4096q^2/a, likewise for G.
        total=[a[k]+(4096*b[k-2] if k>=2 else 0) for k in range(degree+1)]
        need(total[:4]==[1,sign*24,4372,-sign*96256],'exact reciprocal product coefficients')
        records[str(sign)]=total
    need(fib(3)==2 and fib(3)%3!=0,'historical Fp divisibility hostile')
    need(lucas(23)==139*461,'prime value versus prime index')
    return records
def pi_decimal():
    def atan_inv(n):
        x=D(1)/n;term=x;s=term;k=1
        while abs(term)>D('1e-125'):
            term*=-(x*x);s+=term/(2*k+1);k+=1
        return s
    return 16*atan_inv(5)-4*atan_inv(239)
def cm_numerics():
    rows=[]
    for n in (23,29,31,47,71,79):
        d=-4*n
        f=forms(d);g=forms(-8*n)
        need(len(g)&(len(g)-1)==0,'listed doubled class numbers are powers of two')
        involutions=sum(b==0 or abs(b)==a or a==c for a,b,c in g)
        need(involutions==2,'doubled class group cyclic by unique involution')
        rows.append({'n':n,'D_original':d,'class_number_original':len(f),
                     'D_doubled':-8*n,'class_number_doubled':len(g),
                     'doubled_two_torsion':involutions,'genus_fiber':len(g)//involutions,
                     'forms_original':f,'forms_doubled':g})
    numeric={}
    with localcontext() as ctx:
        ctx.prec=120;pi=pi_decimal()
        def value(n,sign):
            q=(-pi*D(n).sqrt()).exp();prod=D(1)
            for j in range(80):prod*=1+sign*q**(2*j+1)
            return prod/(D(2)**(D(1)/4)*q**(D(1)/24))
        g22=(1+D(2).sqrt()).sqrt();G37=(6+D(37).sqrt())**(D(1)/4)
        g130=((2+D(5).sqrt())*(3+D(13).sqrt())/2).sqrt()
        g18=(D(2).sqrt()+D(3).sqrt())**(D(1)/3)
        for name,n,sign,v in [('g22',22,-1,g22),('G37',37,1,G37),('g130',130,-1,g130),('g18',18,-1,g18)]:
            discrepancy=abs(value(n,sign)-v)
            need(discrepancy<D('1e-100'),'independent infinite-product truncation check')
            numeric[name]={'value':str(v),'product_error':str(discrepancy),'status':'NUMERICAL-VERIFIED'}
        for n,N in [(22,64*trace_recurrence(2,-1,12)+24),
                    (37,64*trace_recurrence(12,-1,6)-24)]:
            exponential=(pi*D(n).sqrt()).exp();q=1/exponential;error=D(N)-exponential
            need(0<error<1,'displayed near integer below')
            need(abs(error/q-D(4372))<D('0.1'),'leading error constant')
            numeric[str(n)]={'nearest_integer':N,'exponential':str(exponential),
                             'error':str(error),'scaled_error':str(error/q)}
        q=(-pi*D(130).sqrt()).exp();N=64*trace_recurrence(4,-1,12)*trace_recurrence(3,-1,12)+24
        a=2+D(5).sqrt();b=(3+D(13).sqrt())/2
        numeric['130_trace']={'integer_trace_approximation':N,'exponential':str(1/q),
                              'trace_error':str(D(N)-1/q),'cross_conjugate_term':str(64*((a/b)**12+(b/a)**12))}
    return {'class_number_pairs':rows,'g18_class_number':len(forms(-72)),
            'g18_polynomial':[1,0,0,0,0,0,-10,0,0,0,0,0,1],'numerical':numeric}

def main(write=True):
    global CHECKS
    CHECKS=0
    out={'status':'FINITE-EXACT except labelled numerical evaluations',
         'frames':trace_frames(),'holonomy':reversal_holonomy(),'prime_clocks':prime_clocks(),
         'graph_controls':stars_and_lucas_graphs(),'series':series_controls(),'cm':cm_numerics()}
    out['active_checks']=CHECKS
    if write:(HERE/'receipt.json').write_text(json.dumps(out,indent=2)+'\n',encoding='utf-8')
    print(json.dumps({'checks':CHECKS,'clock_primes':len(out['prime_clocks']['records']),
                      'clock_incidences':out['prime_clocks']['incidences'],
                      'class_number_pairs':[[r['n'],r['class_number_original'],r['class_number_doubled']]
                                           for r in out['cm']['class_number_pairs']],
                      'near_integers':{n:out['cm']['numerical'][n] for n in ('22','37','130_trace')}},indent=2))
    return out
if __name__=='__main__':main()
