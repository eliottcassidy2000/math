#!/usr/bin/env python3
"""Independent exact controls; all checks stay active under python -O.

Does not import or modify either author's verifier. Finite tests support, rather
than replace, the unbounded proof reviewed in INDEPENDENT_REVIEW.md.
"""
from __future__ import annotations
import argparse
from collections import Counter
import hashlib
import json
import math
from pathlib import Path
import sys
import sympy as S

SOURCE = Path(__file__).resolve().parent.parent
E = {2,4,6,8,10,12,16,24,28,32,40,44,48}
CYCLES = ((16,78,46,30,22,18),(24,52,54,34),(28,36,42),
          (32,158,86,50),(40,66),(48,108,96,220,180,150,82))
CHAINS = ((20,30,22,18),(56,80,138,76,72,94,54),
          (88,108,96,220,180,150,82),(184,192,444,348,276,222,118),
          (192,444,348,276,222,118),(224,344,332,264,262,138),
          (384,892,684,528,558,286),(448,696,640,1146,580,450,232),
          (752,768,1788,1356,1032,934,474),(768,1788,1356,1032,934,474),
          (1536,3580,2700,2040,1816,1620,1230),(1792,2808,2488,2208,2266,1140))
FINITE = {row[0]:len(row)-1 for row in CHAINS}

def require(condition, detail):
    if not condition:
        raise RuntimeError(detail)

def core(n):
    require(isinstance(n,int) and n>0, ('invalid core argument',n))
    return n & -n

def F(n,t=2):
    require(isinstance(n,int) and n>0 and n%2==0 and t>=1, ('invalid map input',n,t))
    c=n & -n
    return n-n//c+2*t*c-1

def path(n,length):
    out=[n]
    for _ in range(length): out.append(F(out[-1]))
    return out

def algebra():
    c,d,e,h,p,t,k,u,v,w,x,y=S.symbols('c d e h p t k u v w x y')
    identities=[]; guards=[]; inequalities=[]
    def eq(name,left,right):
        require(S.cancel(left-right)==0, ('identity',name,left,right))
        identities.append(name)
    def step(n,C,t=t): return n-n/C+2*t*C-1
    def G(n,C): return step(n,C,S.Integer(2))
    def route(name,n,cores,expected):
        for i,(C,want) in enumerate(zip(cores,expected),1):
            n=G(n,C); eq(f'{name}, step {i}',n,want)
        return S.cancel(n)
    def odd(name,quotient,subs,variables):
        poly=S.Poly(S.cancel(quotient.subs(subs,simultaneous=True)),*variables)
        constant=poly.coeff_monomial(1)
        require(constant.q==1 and int(constant)%2==1, ('odd constant',name,poly))
        for powers,coefficient in poly.terms():
            require(coefficient.q==1, ('integral quotient',name,poly))
            if any(powers): require(int(coefficient)%2==0, ('even nonconstant',name,poly))
        guards.append({'name':name,'odd_quotient':str(poly.as_expr())})
    def negative(name,expression,subs,variables):
        poly=S.Poly(S.cancel(expression.subs(subs,simultaneous=True)),*variables)
        require(poly.coeff_monomial(1)<0, ('strict negative constant',name,poly))
        require(all(coefficient<=0 for coefficient in poly.coeffs()), ('negative coefficients',name,poly))
        inequalities.append({'name':name,'shifted_polynomial':str(poly.as_expr()),
                             'domain':'All displayed variables are nonnegative.'})
    n0=c*p; n1=step(n0,c); n2=step(n1,d)
    eq('all-t two-step identity',d*(n2-n0),(c+d)*(2*t*(d-1)-p)+2*(t-1)*d+p+1)
    eq('fixed point',step(c*(2*t*c-1),c),c*(2*t*c-1))
    negative('all-t residual upper bound',c*(1-d)-d*d-d,{c:x+2,d:y+2},(x,y))
    odd('first image core d when d<c',(c*(d*k+2*t-1)-d*k)/d,
        {c:4*u*v,d:2*v,k:2*w+1},(u,v,w,t))

    n0=e*d*(d*k-1); n1=G(n0,e*d); n2=G(n1,d)
    eq('common first iterate',n1,d*(e*(d*k+3)-k))
    eq('common second iterate',n2,e*(d-1)*(d*k+3)+d*(4-k)+k-1)
    N3=S.cancel(n2.subs(k,3)); B3=n0.subs(k,3)
    eq('k3 common second iterate',N3,3*e*(d*d-1)+d+2)
    expr=G(N3,2)-B3
    eq('k3 large-e descent',expr,-e*(3*d*d-2*d+3)/2+d/2+8)
    negative('k3 d>=4,e>=4',expr,{d:x+4,e:y+4},(x,y))
    odd('k3 large-e n2 core 2',N3/2,{d:4*u,e:4*v},(u,v))
    expr=(G(N3,4)-B3).subs(e,2)
    eq('k3 e=2 descent',expr,-3*d*d/2+11*d/4+12)
    negative('k3 e=2,d>=8',expr,{d:x+8},(x,))
    odd('k3 e=2 n2 core 4',N3/4,{e:2,d:8*u},(u,))
    eq('k3 d=2 two-step descent',(N3-B3).subs(d,2),4-e)
    negative('k3 d=2,e>=8',4-e,{e:x+8},(x,))
    diag=[3*c*c,7*c*c-4,21*c*c/4+12,63*c*c/16+24,441*c*c/128+52,1323*c*c/512+54]
    z=route('k3 d=c',c*(3*c-1),[c,c*c,4,4,8,4],diag)
    eq('k3 diagonal descent',z-c*(3*c-1),-213*c*c/512+c+54)
    negative('k3 diagonal c>=32',z-c*(3*c-1),{c:x+32},(x,))
    for index,C,threshold in [(0,c*c,2),(1,4,4),(2,4,8),(3,8,16),(4,4,32)]:
        odd(f'k3 diagonal n{index+1} core {C}',diag[index]/C,{c:threshold*u},(u,))
    lane=[c*(2*c+1),2*c*c+3*c-2,c*c+3*c/2+6]
    z=route('k1 d=2c',2*c*c-c,[c,c,2],lane)
    negative('k1 d=2c,c>=8',z-(2*c*c-c),{c:x+8},(x,))
    odd('k1 d=2c n1 core c',lane[0]/c,{c:2*u},(u,))
    odd('k1 d=2c n2 core 2',lane[1]/2,{c:4*u},(u,))
    lane=[c*(c+2),c*c+19*c/2-2,c*c/2+19*c/4+6]
    z=route('k1 d=c',c*(c-1),[c,2*c,2],lane)
    negative('k1 d=c,c>=16',z-c*(c-1),{c:x+16},(x,))
    odd('k1 d=c n1 core 2c',lane[0]/(2*c),{c:4*u},(u,))
    odd('k1 d=c n2 core 2',lane[1]/2,{c:8*u},(u,))
    N1=S.cancel(n2.subs(k,1)); B1=n0.subs(k,1)
    eq('k1 common second iterate',N1,e*(d*d+2*d-3)+3*d)
    small=G(N1,e)
    eq('k1 e<d third iterate',small,e*(d+1)**2-d*d+d+2-3*d/e)
    eq('k1 e<d descent',small-B1,e*(3*d+1)-d*d+d+2-3*d/e)
    odd('k1 e<d n2 core e',N1/e,{e:2*u,d:4*u*v},(u,v))
    negative('k1 e<=d/4 upper bound,d>=8',-d*d/4+5*d/4+2,{d:x+8},(x,))
    half=S.cancel(small.subs(e,d/2))
    eq('k1 e=d/2 third iterate',half,d**3/2+3*d/2-4)
    odd('k1 e=d/2 n3 core 4',half/4,{d:16*u},(u,))
    negative('k1 e=d/2,d>=16',G(half,4)-B1.subs(e,d/2),{d:x+16},(x,))
    diagonal=S.cancel(N1.subs(e,d)); third=G(diagonal,2*d*d)
    eq('k1 e=d second iterate',diagonal,d*d*(d+2))
    eq('k1 e=d third iterate',third,d**3+10*d*d-d/2-2)
    odd('k1 e=d n2 core 2d^2',diagonal/(2*d*d),{d:4*u},(u,))
    odd('k1 e=d n3 core 2',third/2,{d:8*u},(u,))
    negative('k1 e=d,d>=16',G(third,2)-B1.subs(e,d),{d:x+16},(x,))
    large=G(N1.subs(e,h*d),d); base=B1.subs(e,h*d)
    eq('k1 e>d third iterate',large,h*(d-1)**2*(d+3)+7*d-4)
    odd('k1 e>d n2 core d',N1.subs(e,h*d)/d,{d:2*u,h:2*v},(u,v))
    odd('k1 h=2 n3 core 2',large/2,{h:2,d:4*u},(u,))
    odd('k1 h>=8,d>=8 n3 core 4',large/4,{h:8*v,d:8*u},(u,v))
    odd('k1 h=4 n3 core 8',large/8,{h:4,d:16*u},(u,))
    negative('k1 h=2,d>=4',(G(large,2)-base).subs(h,2),{d:x+4},(x,))
    negative('k1 h>=8,d>=8',G(large,4)-base,{h:y+8,d:x+8},(x,y))
    negative('k1 h=4,d>=16',(G(large,8)-base).subs(h,4),{d:x+16},(x,))
    fourth=G(large.subs(d,4),8)
    eq('k1 d=4 fourth iterate',fourth,441*h/8+52)
    odd('k1 d=4,h>=64 n3 core 8',large.subs(d,4)/8,{h:64*u},(u,))
    odd('k1 d=4,h>=64 n4 core 4',fourth/4,{h:64*u},(u,))
    negative('k1 d=4,h>=64',G(fourth,4)-base.subs(d,4),{h:x+64},(x,))
    odd('k1 d=2,h>=16 n3 core 2',large.subs(d,2)/2,{h:16*u},(u,))
    negative('k1 d=2,h>=16',(G(large,2)-base).subs(d,2),{h:x+16},(x,))

    D=c+d-1; A=2*t-1
    P=2*t*(d-1)+A*(2*d-1)/D
    Q=2*t*(c-1)+A*(2*c-1)/D
    eq('two-cycle first equation',(c-1)*P+2*t*c-1,d*Q)
    eq('two-cycle second equation',(d-1)*Q+2*t*d-1,c*P)
    eq('two-cycle fractions sum',A*(2*d-1)/D+A*(2*c-1)/D,4*t-2)
    eq('Mersenne gcd linear combination',2*D-(2*d-1),2*c-1)
    eq('larger core bound simplification',(4*t-3)*(2*t-3)+(2*t-1),8*(t-1)**2)
    eq('sharp first odd part',P.subs({t:c/2+1,d:2*c*c}),2*c**3+4*c*c+c-1)
    eq('sharp second odd part',Q.subs({t:c/2+1,d:2*c*c}),c*c+c-1)
    return {'identities_checked':len(identities),'identities':identities,
            'unbounded_exact_core_certificates':guards,'unbounded_strict_inequality_certificates':inequalities}

def branch(n):
    c=core(n); q=n//c+1; d=core(q); k=q//d
    if q>4*c: return 'immediate',1,[]
    if q==4*c: return 'fixed',0,[]
    if k>=5: return 'general odd part',2,[(1,d)]
    def boundary():
        if n in E:return 'exception',0,[]
        require(n in FINITE, ('uncovered exponent region',n,c,d,k))
        return 'finite chain',FINITE[n],[]
    if k==3:
        require(d<=c,('k3 growth range',n))
        if d==c:
            if c>=32:return 'k3 diagonal',6,[(1,c*c),(2,4),(3,4),(4,8),(5,4)]
            return boundary()
        e=c//d
        if d>=4 and e>=4:return 'k3 large e',3,[(1,d),(2,2)]
        if e==2 and d>=8:return 'k3 e2',3,[(1,d),(2,4)]
        if d==2 and e>=8:return 'k3 d2',2,[(1,d)]
        return boundary()
    require(k==1 and d<=2*c,('residual coefficient',n))
    if d==2*c:
        return ('k1 d2c',3,[(1,c),(2,2)]) if c>=8 else boundary()
    if d==c:
        return ('k1 dc',3,[(1,2*c),(2,2)]) if c>=16 else boundary()
    e=c//d
    if e<d:
        if e<=d//4:return 'k1 small e',3,[(1,d),(2,e)]
        require(2*e==d,('adjacent powers',n))
        return ('k1 half e',4,[(1,d),(2,e),(3,4)]) if d>=16 else boundary()
    if e==d:
        return ('k1 diagonal e',4,[(1,d),(2,2*d*d),(3,2)]) if d>=16 else boundary()
    h=e//d
    require(e==h*d and h>=2,('large e ratio',n))
    if h==2 and d>=4:return 'k1 h2',4,[(1,d),(2,d),(3,2)]
    if h>=8 and d>=8:return 'k1 large h',4,[(1,d),(2,d),(3,4)]
    if h==4 and d>=16:return 'k1 h4',4,[(1,d),(2,d),(3,8)]
    if d==4 and h>=64:return 'k1 d4',5,[(1,d),(2,d),(3,8),(4,4)]
    if d==2 and h>=16:return 'k1 d2',4,[(1,d),(2,d),(3,2)]
    return boundary()

def finite_and_grid(max_exponent,max_n):
    members=set()
    for row in CHAINS:
        require(path(row[0],len(row)-1)==list(row),('finite chain',row))
        require(row[-1]<row[0] and len(row)<=7,('chain descent',row))
    for cycle in CYCLES:
        require(len(set(cycle))==len(cycle) and not (members&set(cycle)),('cycle simplicity',cycle))
        require(path(cycle[0],len(cycle))==list(cycle)+[cycle[0]],('cycle edges',cycle))
        members.update(cycle)
    for n in E:
        z=n; visited=set()
        while z not in members:
            require(z not in visited,('exception failed to reach listed cycle',n))
            visited.add(z); z=F(z)
    counts=Counter(); coefficients=(1,3,5,7,9,11,17,257)
    for a in range(1,max_exponent+1):
        for b in range(1,max_exponent+1):
            for k in coefficients:
                n=(1<<a)*((1<<b)*k-1)
                label,j,guards=branch(n); counts[label]+=1
                if j==0:
                    require(F(n)==n if label=='fixed' else n in E,('zero-step case',n,label))
                    continue
                orbit=path(n,j)
                require(1<=j<=6 and orbit[-1]<n,('predicted descent',n,label,orbit))
                for index,C in guards:require(core(orbit[index])==C,('guard',n,label,index,C))
    first_descent=Counter(); actual_exceptions=[]; fixed_count=0
    for n in range(2,max_n+1,2):
        z=n
        for j in range(1,7):
            z=F(z)
            if z<n:first_descent[j]+=1;break
        else:
            if F(n)==n:fixed_count+=1
            else:actual_exceptions.append(n)
    require(actual_exceptions==sorted(n for n in E if n<=max_n),('six-step exceptions',actual_exceptions))
    sharp=path(56,6)
    require(all(z>=56 for z in sharp[1:6]) and sharp[6]<56,('six-step sharpness',sharp))
    return {'exponent_grid':{'a_max':max_exponent,'b_max':max_exponent,'k':coefficients,
                            'inputs':sum(counts.values()),'branch_counts':dict(sorted(counts.items()))},
            'literal_short_block_interval':{'max_n':max_n,'even_inputs':max_n//2,
               'first_descent_counts':dict(sorted(first_descent.items())),
               'fixed_starts':fixed_count,'nonfixed_no_descent_starts':actual_exceptions},
            'finite_chains':len(CHAINS),'listed_cycles':len(CYCLES),'sharp_six_step_witness':sharp}

def general_reduction():
    count=0
    for t in range(1,65):
        for a in range(1,41):
            c=1<<a
            for b in range(1,41):
                d=1<<b
                for k in (2*t+1,2*t+3,2*t+17,2*t+255):
                    n=c*(d*k-1)
                    if d*k>=2*t*c:continue
                    first=F(n,t)
                    require(first>n and d<c and core(first)==d and F(first,t)<n,
                            ('all-t two-step growth reduction',t,a,b,k))
                    count+=1
    return {'parameters':64,'a_max':40,'b_max':40,'qualifying_growth_inputs':count}

def predicted_cycles(t):
    out=set(); c=2
    while c<=2*t-2:
        d=2*c
        while d<=8*(t-1)**2:
            D=c+d-1; gcd=math.gcd(D,2*d-1)
            if (2*t-1)%(D//gcd)==0:
                j=(2*t-1)*(2*d-1)//D
                p=2*t*(d-1)+j; q=2*t*(c-1)+(4*t-2-j)
                n=c*p; m=d*q
                require(p>0 and q>0 and p%2==q%2==1,('positive exact odd parts',t,c,d))
                require(n!=m and F(n,t)==m and F(m,t)==n,('predicted cycle',t,c,d))
                out.add(tuple(sorted((n,m))))
            d*=2
        c*=2
    return out

def two_cycles():
    total=0; pair_tests=0
    for t in range(1,1025):
        total+=len(predicted_cycles(t))
        for a in range(1,25):
            c=1<<a
            for b in range(a+1,33):
                d=1<<b; D=c+d-1; g=(1<<math.gcd(a+1,b+1))-1
                require(math.gcd(D,2*d-1)==g,('Mersenne gcd',a,b))
                integral=(2*t-1)*(2*d-1)%D==0
                require(integral==((2*t-1)%(D//g)==0),('gcd divisibility criterion',t,a,b))
                if integral:require(c<=2*t-2 and d<=8*(t-1)**2,('core bounds',t,a,b))
                pair_tests+=1
    for t in range(1,33):
        observed=set()
        for n in range(2,200001,2):
            m=F(n,t)
            if n!=m and F(m,t)==n:observed.add(tuple(sorted((n,m))))
        require(observed=={cycle for cycle in predicted_cycles(t) if cycle[0]<=200000},('literal two-cycle census',t))
    for a in range(1,257):
        c=1<<a; t=c//2+1; d=2*c*c
        n=c*(2*c**3+4*c*c+c-1); m=d*(c*c+c-1)
        require(core(n)==c and core(m)==d and F(n,t)==m and F(m,t)==n,
                ('sharp family',a))
        require(c==2*t-2 and d==8*(t-1)**2,('simultaneous sharpness',a))
    require(predicted_cycles(1)==set() and predicted_cycles(2)=={(40,66)},'t=1,2 specialization')
    return {'classified_parameters':1024,'cycles_total':total,'gcd_and_bound_tests':pair_tests,
            'literal_parameters':32,'literal_even_starts_per_parameter':100000,
            'sharp_family_exponents':256,'t2_two_cycles':[[40,66]]}

def main():
    parser=argparse.ArgumentParser()
    parser.add_argument('--max-exponent',type=int,default=256)
    parser.add_argument('--max-n',type=int,default=2000000)
    parser.add_argument('--output',type=Path)
    parser.add_argument('--negative-control',action='store_true')
    args=parser.parse_args()
    if args.negative_control:require(False,'Intentional negative control: checks remain active under -O.')
    names=('FINITE_ODD_PART_REDUCTION.md','T2_GLOBAL_CLASSIFICATION.md','TWO_CYCLE_CLASSIFICATION.md',
           'verify_t2_global.py','t2-verification.json','verify_two_cycles.py',
           'verify_two_cycles.out','verify_two_cycles.optimized.out')
    hashes={name:hashlib.sha256((SOURCE/name).read_bytes()).hexdigest() for name in names}
    result={'status':'PASS','check_implementation':'Explicit require checks; no Python assert statements.',
            'scope':'Independent exact algebra, universal polynomial guard/inequality certificates, and finite controls.',
            'source_sha256':hashes,'symbolic':algebra(),
            't2':finite_and_grid(args.max_exponent,args.max_n),
            'all_t_reduction':general_reduction(),'all_t_two_cycles':two_cycles()}
    require(hashes=={name:hashlib.sha256((SOURCE/name).read_bytes()).hexdigest() for name in names},'Source files changed during audit.')
    rendered=json.dumps(result,indent=2)+'\n'
    if args.output:args.output.write_text(rendered)
    print(rendered,end='')

if __name__=='__main__':main()
