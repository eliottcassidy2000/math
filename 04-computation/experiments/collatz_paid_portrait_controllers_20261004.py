"""Paid cyclic controllers, anchor-switch shells, and finite quadratic portraits.

All arithmetic is exact. No search hit is treated as a universal home certificate.
The inherited sixteen-debt/ternary bank is the only coverage comparison claimed.
"""
from fractions import Fraction as F
from itertools import product
from pathlib import Path
import json
import math

ROOT = Path(__file__).resolve().parents[2]
CHECKS = 0
CYCLES = ((1, (1,)), (5, (1, 2)), (17, (1, 1, 1, 2, 1, 1, 4)))


def check(ok, why):
    global CHECKS
    CHECKS += 1
    if not ok:
        raise ArithmeticError(why)


def v2(x):
    x = F(x)
    if not x:
        raise ValueError('zero has infinite valuation')
    def val(z):
        z = abs(z)
        return (z & -z).bit_length()-1
    return val(x.numerator)-val(x.denominator)


def coefficients(word):
    P, Q, B = 1, 1, 0
    for a in word:
        P, Q, B = 3*P, Q*2**a, 3*B+Q
    return P, Q, B


def anchor(word):
    P, Q, B = coefficients(word)
    return F(B, Q-P)


def step(n):
    z = 3*n+1
    a = v2(z)
    return z//2**a, a


def replay(n, word):
    for a in word:
        n, actual = step(n)
        check(a == actual, 'literal exact integer valuation')
    return n


def rank(n):
    if n == 1:
        return 0, 0
    k = v2(n-1)
    return 3**k*((n-1)//2**k)**2, k


def critical(n):
    return n > 1 and n % 4 == 3 and n % 3 != 2 and n % 32 != 11


def anchored_runs():
    cases, longest = 0, None
    words = [w for length in range(1, 5) for w in product(range(1, 5), repeat=length)
             if sum(w) <= 9]
    for w in words:
        rho, A = anchor(w), sum(w)
        check(rho.denominator % 2 == 1 and v2(rho) == 0, 'odd periodic anchor')
        P, Q, B = coefficients(w)
        check(F(P*rho+B, Q) == rho, 'anchor fixed by formal word')
        lam=F(P,Q);fold=lam+1/lam
        check(coefficients(w+w)==(P*P,Q*Q,B*(P+Q)), 'word repetition transports multiplier and carry')
        check(lam*lam+1/(lam*lam)==fold*fold-2 and anchor(w+w)==rho,
              'Chebyshev operation quotient with retained anchor')
        for n in range(-101, 102, 2):
            if n == rho:
                continue
            K = v2(n-rho)
            q = (K-1)//A
            x = n
            for _ in range(q):
                x = replay(x, w)
            check(x == rho+F(P, Q)**q*(n-rho), 'exact repeated affine endpoint')
            check(v2(x-rho) == K-q*A, 'exact fuel remaining')
            y, still_matches = x, True
            for a in w:
                y, actual = step(y)
                if actual != a:
                    still_matches = False
                    break
            check(not still_matches, 'maximal whole-word run ends within the next word')
            cases += 1
            if longest is None or q > longest['repeats']:
                longest = dict(n=n, word=w, anchor=str(rho), repeats=q)
    return dict(words=len(words), signed_inputs_per_word=102, cases=cases, longest=longest)


def switch_shells():
    words = ((1,), (2,), (1,2), (1,1,2), (1,1,1,2,1,1,4))
    anchors = [anchor(w) for w in words]
    count = cancellation = recharge = 0
    for w in words:
        P, Q, B = coefficients(w)
        A = sum(w)
        for rho, sigma in product(anchors, repeat=2):
            delta = F(P*rho+B, Q)-sigma
            for n in range(-101, 102, 2):
                y = F(P*n+B, Q)
                if n == rho or y == sigma:
                    continue
                K, L = v2(n-rho), v2(y-sigma)
                if not delta:
                    check(L == K-A, 'transported anchor consumes its exact cost')
                else:
                    h = v2(delta)
                    if K-A != h:
                        check(L == min(K-A, h), 'ultrametric switch outside cancellation shell')
                    else:
                        check(L > h, 'cancellation of two odd normalized terms')
                        cancellation += 1
                    if L > K:
                        check(K == A+h, 'all fuel increases lie on the stated shell')
                        recharge += 1
                count += 1
    examples = []
    for h in range(3, 41):
        modulus = 2**(h-1)
        # n=32t-5, K(n+5)=5; require v2(9t-1)=h-2 exactly.
        t = ((1+2**(h-2))*pow(9,-1,modulus)) % modulus
        n, y = 32*t-5, 36*t-5
        check(replay(n,(1,2)) == y, 'recharge actual word')
        check(v2(n+5) == 5 and v2(y+1) == h, 'arbitrarily large recharge from fixed old fuel')
        if h in (3, 5, 6, 20, 40):
            examples.append(dict(h=h, t=t, n=n, child=y))
    separation = [[str(v2(a-b)) if a != b else 'infinity' for b in anchors] for a in anchors]
    return dict(cases=count, cancellation_cases=cancellation, strict_recharges=recharge,
                anchors=list(map(str,anchors)), separation=separation, examples=examples)


def symbolic_word(n, period, word):
    for a in word:
        check((3*n+1) % 2**(a+1) == 2**a and 3*period % 2**(a+1) == 0,
              'all-height exact word guard')
        n, period = (3*n+1)//2**a, 3*period//2**a
    return n, period


def threshold(r, A, q):
    a = 2
    while 2**(A*q+a) <= 3**(r*q+1):
        a += 1
    return a


def paid_families():
    certificates, density_rows = [], []
    for c, w in CYCLES:
        r, A = len(w), sum(w)
        check(anchor(w) == -c and replay(-c,w) == -c, 'inherited signed cycle')
        for q in range(1, 65):
            amin = threshold(r,A,q)
            L = amin+1
            power, d = 3**(r*q+1), (3*c-1)//2
            check(F(power,2**(A*q+L)) < F(1,2), 'safe exit pays at least half the multiplier')
            check(2**(A*q+1)-c >= c, 'smallest source exceeds anchor magnitude')
            t0 = (d*pow(power,-1,2**(L-1))) % 2**(L-1)
            residue, modulus = 2**(A*q+1)*t0-c, 2**(A*q+L)
            density_rows.append(dict(c=c,q=q,L=L,residue=residue,modulus=modulus,
                                     odd_density=str(F(2,modulus))))
            for a in range(amin,amin+3):
                t = ((d+2**(a-1))*pow(power,-1,2**a)) % 2**a
                N, P = 2**(A*q+1)*t-c, 2**(A*q+1+a)
                M, Q = (power*t-d)//2**(a-1), 2*power
                check(P > Q, 'strict all-height paid slope')
                shift = 0 if N > M else (M-N)//(P-Q)+1
                N, M = N+shift*P, M+shift*Q
                word = w*q+(a,)
                check(symbolic_word(N,P,word) == (M,Q), 'symbolic all-height actual forward exit')
                check(0 < M < N and rank(M) < rank(N), 'paid family seed')
                for j in (0,1,17,10**6):
                    x,y=N+P*j,M+Q*j
                    check(0 < y < x and rank(y) < rank(x), 'all-height rank payment control')
                certificates.append(dict(c=c,q=q,a=a,n=N,P=P,m=M,Q=Q,tail_shift=shift))
    return certificates,density_rows


def rational_cycle_bank(H=24, Qmax=24):
    rows=[]
    for h in range(1,H+1):
        w=(1,)*h+(2,)
        P,Q,B=coefficients(w)
        c=-anchor(w)
        check(c == F(3**(h+1)-2**(h+1),3**(h+1)-2**(h+2)),
              'entire one-run reset-two anchor atlas')
        check(v2(c-1)==h+1, 'anchor atlas approaches minus one at a known binary depth')
        for q in range(2,Qmax+1):
            a0=threshold(h+1,h+2,q);L=a0+1
            p,den=P**q,Q**q
            carry=c*(p-den)
            check(carry.denominator==1,'repeated word carry integral')
            totalP,totalB=3*p,3*int(carry)+den
            modulus=den*2**L
            residue=(-totalB*pow(totalP,-1,modulus)) % modulus
            check(2**((h+2)*q)>=c.numerator,'all-q safe source-height premise')
            check(F(totalP,modulus)<F(1,2),'rational-anchor exit strict half multiplier')
            for j in (0,1,17):
                n=residue+modulus*j
                x=replay(n,w*q)
                m,a=step(x)
                check(a>=L and 0<m<n and rank(m)<rank(n),'rational-controller literal payment')
                check(v2(n-1)==1 and n%32!=11,'rational bank lies in binary critical domain')
            rows.append(dict(h=h,q=q,anchor=str(-c),L=L,residue=residue,modulus=modulus,
                             odd_density=str(F(2,modulus))))
    inner=sum(F(1,2**((h+2)*(Qmax+1)+2))/(1-F(1,2**(h+2))) for h in range(1,H+1))
    outer=F(8,7)*F(1,2**(2*(H+1)+6))/(1-F(1,4))
    return dict(H=H,Q=Qmax,rows=rows,tail_bound=str(inner+outer))


def powers_of_three(bank):
    rows=[]
    for row in bank['rows']:
        if row['h']!=1 or row['q']>12:
            continue
        residue,modulus=row['residue'],row['modulus']
        bits=modulus.bit_length()-1
        exponent=1
        check(residue%8==3,'power-three subgroup at the base precision')
        for m in range(3,bits):
            candidates=[exponent,exponent+2**(m-2)]
            valid=[a for a in candidates if pow(3,a,2**(m+1))==residue%2**(m+1)]
            check(len(valid)==1,'unique exponent lift in the cyclic group generated by three')
            exponent=valid[0]
        period=2**(bits-2)
        check(exponent%8==3,'new families lie in the incoming residual exponent class')
        check(pow(3,exponent,modulus)==residue and pow(3,period,modulus)==1,
              'whole exponent progression enters paid source cylinder')
        literal=None
        if exponent<=10000:
            n=3**exponent
            x=replay(n,(1,2)*row['q'])
            child,a=step(x)
            check(a>=row['L'] and 0<child<n and rank(child)<rank(n),'literal power-three paid controller')
            literal=dict(source_bits=n.bit_length(),child_bits=child.bit_length(),exit=a)
        rows.append(dict(q=row['q'],L=row['L'],exponent=exponent,period=period,
                         source_modulus=modulus,literal=literal))
    return rows


def switched_bank(Qmax=24,Rmax=64):
    rows=[]
    for q in range(1,Qmax+1):
        for r in range(1,Rmax+1):
            L=2
            while 2**(3*q+r+L)<=2*3**(2*q+r+1):L+=1
            prefix=(1,2)*q+(1,)*r
            P,Q,B=coefficients(prefix+(L,))
            residue=-B*pow(P,-1,Q)%Q
            density=F(2,Q)
            holes=[]
            if (q,r)==(1,1):
                p,z,b=coefficients(prefix+(6,))
                hole=((z-b)*pow(p,-1,2*z))%(2*z)
                holes=[dict(residue=hole,modulus=2*z)]
                check((hole,2*z)==(315,2048),'exact inherited debt-row overlap')
                density-=F(1,1024)
            for j in (0,1,17):
                n=residue+Q*j
                x=replay(n,prefix)
                m,a=step(x)
                check(a>=L and 0<m<n and rank(m)<rank(n),'paid switch retains original source rank')
                check(n%32==27,'switched bank binary critical class')
            rows.append(dict(q=q,r=r,L=L,residue=residue,modulus=Q,
                             holes=holes,odd_density=str(density)))
    tail=F(4,7)*F(1,2**(3*(Qmax+1)))+F(1,14*2**Rmax)
    powers=[]
    for row in rows:
        if row['q']!=1 or row['r'] not in (1,2,3,4):continue
        residue,modulus=row['residue'],row['modulus']
        bits=modulus.bit_length()-1;exponent=1
        for b in range(3,bits):
            valid=[a for a in (exponent,exponent+2**(b-2))
                   if pow(3,a,2**(b+1))==residue%2**(b+1)]
            check(len(valid)==1,'switched-controller exponent lift')
            exponent=valid[0]
        powers.append(dict(q=1,r=row['r'],L=row['L'],exponent=exponent,
                           period=2**(bits-2),excluded_source_cells=row['holes']))
    checkpoint=dict(residue=155,modulus=2048,odd_density=str(F(1,1024)),
                    child=111,child_period=1458)
    check(symbolic_word(155,2048,(1,2,1,1,1,2))==(445,5832),
          'incoming checkpoint identity verified at every height')
    for j in (0,1,17,10**6):
        n,m=155+2048*j,111+1458*j
        x=replay(n,(1,2,1,1,1,2))
        check(x==4*m+1 and step(x)[0]==step(m)[0] and rank(m)<rank(n),
              'incoming supplied-child common future and original rank payment')
    check(pow(3,483,2048)==155 and pow(3,512,2048)==1,'incoming exponent progression')
    return dict(Q=Qmax,R=Rmax,rows=rows,tail_bound=str(tail),power_rows=powers,checkpoint=checkpoint)


def coverage(density_rows, rational_bank, switched):
    old=json.loads((ROOT/'05-knowledge/results/collatz_binary_ternary_guard_fusion_20261004.json').read_text())
    debt1=[row for row in old['debt_rows'] if row['e']==1]
    check([(row['s'],row['source_word'][:-1]) for row in debt1] == [(1,[1,6]),(2,[6])],
          'frozen binary bank has exactly the two inherited debt-one prefixes')
    selected=[row for row in density_rows if row['c']==17]
    for row in selected:
        check(row['residue'] % 32 == 15, 'minus-seventeen bank binary critical class')
    one_lower=sum(F(row['odd_density']) for row in selected+rational_bank['rows'])
    switch_lower=sum(F(row['odd_density']) for row in switched['rows'])
    lower=one_lower+switch_lower+F(1,1024)
    tail=F(rational_bank['tail_bound'])+F(1,2**(11*65+2))/(1-F(1,2**11))+F(switched['tail_bound'])
    upper=lower+tail
    d2=F(old['binary_necessary_density'])
    lo3,hi3=F(old['ternary']['lower']),F(old['ternary']['upper'])
    newlo,newhi=(d2-upper)*(1-hi3),(d2-lower)*(1-lo3)
    gainlo,gainhi=lower*(1-hi3),upper*(1-lo3)
    check(newhi < F(old['fused_necessary_density_interval'][0]), 'strict improvement of named inherited bank')
    origin_rows=[(1,2)]
    for k in range(1,32):
        d=1
        while 3**d<=2**(d+2*k-1):d+=1
        residue=-((4**k+2)//3)*pow(4**k,-1,3**d) % 3**d
        if not any(residue%3**j==r for j,r in origin_rows):
            origin_rows.append((d,residue))
    originlo=sum(F(1,3**d) for d,_ in origin_rows)
    d=1
    while 3**d<=2**(d+63):d+=1
    originhi=originlo+F(27,26*3**d)
    check(F('0.4458180914749391')<=originlo<originhi<=F('0.4458180914749392'),
          'independent recovery of incoming strengthened origin bank')
    origin_residual=((d2-upper)*(1-originhi),(d2-lower)*(1-originlo))
    origin_gain=(lower*(1-originhi),upper*(1-originlo))
    c5lower=sum(F(row['odd_density']) for row in density_rows if row['c']==5 and row['q']>=2)
    c5upper=c5lower+F(1,2**(3*65+2))/(1-F(1,8))
    # Independent finite census uses the live map, not the density formula.
    counts=dict(odd_sources=50000,critical=0,integer_cycle_safe_union=0,integer_cycle_safe_critical=0,new_union=0,new_critical=0)
    examples=[]
    allcells=selected+rational_bank['rows']+switched['rows']+[switched['checkpoint']]
    active_cells=[row for row in allcells if row['residue']<100000]
    def contains(row,n):
        return n%row['modulus']==row['residue'] and not any(
            n%hole['modulus']==hole['residue'] for hole in row.get('holes',[]))
    def old_binary_applies(n):
        x,a=step(n)
        if a!=1:
            return False
        while a==1:
            y,a=step(x)
            if a==1:x=y
        if a!=2:
            return False
        e=0
        while True:
            y,a=step(x)
            if a!=2:break
            e+=1;x=y
        for row in old['debt_rows']:
            if row['e']!=e:continue
            y=x
            for a in row['source_word'][:-1]:
                y,b=step(y)
                if a!=b:break
            else:return True
        return False
    check(old_binary_applies(2363),'positive control for inherited binary bank')
    for n in range(1,100000,2):
        crit=critical(n)
        counts['critical']+=crit
        admitted=[]
        for c,w in CYCLES:
            K=v2(n+c);r,A=len(w),sum(w)
            if (K-1)%A or K<=1:
                continue
            q=(K-1)//A
            L=threshold(r,A,q)+1
            x=F(-c)+F(3**r,2**A)**q*(n+c)
            check(x.denominator==1, 'formula endpoint integral in live census')
            y,a=step(int(x))
            if a>=L:
                check(0<y<n and rank(y)<rank(n), 'live census paid endpoint')
                admitted.append(c)
                isnew=c==17 or c==5 and q>=2
                if isnew and crit and len(examples)<12:
                    examples.append(dict(n=n,m=y,c=c,q=q,a=a))
        check(len(admitted)<=1, 'three safe cycle banks are disjoint')
        counts['integer_cycle_safe_union']+=bool(admitted)
        counts['integer_cycle_safe_critical']+=bool(admitted) and crit
        new=17 in admitted
        h=v2(n+1)-1
        if h>=1:
            w=(1,)*h+(2,);rho=anchor(w);A=h+2
            K=v2(n-rho)
            if K>=2*A+1 and (K-1)%A==0:
                q=(K-1)//A;L=threshold(h+1,A,q)+1
                x=rho+F(3**(h+1),2**A)**q*(n-rho)
                check(x.denominator==1,'rational controller census endpoint integral')
                m,a=step(int(x))
                if a>=L:
                    check(not new,'rational atlas disjoint from minus-seventeen bank')
                    check(replay(n,w*q)==x and rank(m)<rank(n),'independent rational census payment')
                    new=True
        K=v2(n+5);q=(K-1)//3
        if q>=1:
            x=F(-5)+F(9,8)**q*(n+5)
            check(x.denominator==1,'switch controller first phase integral')
            x=int(x);r=v2(x+1)-1
            if r>=1:
                L=2
                while 2**(3*q+r+L)<=2*3**(2*q+r+1):L+=1
                z=F(3,2)**r*(x+1)-1
                check(z.denominator==1,'switch controller second phase integral')
                m,a=step(int(z))
                if a>=L and not (q==r==1 and a==6):
                    check(not new,'switched bank disjoint from single-anchor banks')
                    check(rank(m)<rank(n),'switched selector pays original rank')
                    new=True
        if n%2048==155:
            check(not new,'incoming checkpoint cell disjoint from both new controller banks')
            new=True
        counts['new_union']+=new
        counts['new_critical']+=new and crit
        check(new==any(contains(row,n) for row in active_cells),
              'independent finite cylinder census matches actual dynamic selector')
        if new:
            check(not old_binary_applies(n),'new selector disjoint from inherited binary bank')
    hostiles=[n for n in (7,27,703) if not any(contains(row,n) for row in active_cells)]
    check(hostiles==[7,27,703],'canonical sources remain outside this safe-exit bank')
    return dict(compared_bank='binary16 and ternary sibling bank from collatz_binary_ternary_guard_fusion_20261004.json',
                single_anchor_binary_decimal=float(one_lower),switched_binary_decimal=float(switch_lower),
                new_binary_density_interval=list(map(str,(lower,upper))),new_binary_decimal=float(lower),
                uncovered_density_interval=list(map(str,(newlo,newhi))),uncovered_decimal=float(newhi),
                added_outside_old_bank_interval=list(map(str,(gainlo,gainhi))),added_decimal=float(gainlo),
                strengthened_origin_density_interval=list(map(str,(originlo,originhi))),
                strengthened_fused_residual_interval=list(map(str,origin_residual)),
                strengthened_fused_residual_decimal=float(origin_residual[1]),
                added_outside_strengthened_origin_interval=list(map(str,origin_gain)),
                added_outside_strengthened_origin_decimal=float(origin_gain[0]),
                single_anchor_repaired_fraction_of_exponents_3_mod8_decimal=float(16*c5lower),
                repaired_fraction_of_exponents_3_mod8_interval=list(map(str,(
                    16*(c5lower+switch_lower+F(1,1024)),
                    16*(c5upper+switch_lower+F(switched['tail_bound'])+F(1,1024))))),
                repaired_fraction_of_exponents_3_mod8_decimal=float(16*(c5lower+switch_lower+F(1,1024))),
                census=counts,examples=examples,hostiles=hostiles)


def group_and_portrait():
    def compose(f,g):
        a,b=f;c,d=g
        return a*c,a*d+b
    def inv(f):
        a,b=f
        return 1/a,-b/a
    D,T=(F(2),F(0)),(F(3),F(1))
    check(compose(compose(compose(D,T),inv(D)),inv(T)) == (1,1), 'affine commutator is translation one')
    finite=[]
    for modulus in (5,7,11,13,19,25,35,95):
        H={1};queue=[1]
        while queue:
            x=queue.pop()
            for g in (2,3):
                y=x*g%modulus
                if y not in H:H.add(y);queue.append(y)
        finite.append(dict(modulus=modulus,multiplier_group_order=len(H),affine_group_order=modulus*len(H)))
    portraits={0:(-1,0,1),-1:(-1,0,1),-2:(-2,-1,0,1,2)}
    # Explicit pullback quotient identities, checked by direct evaluation.
    quotients={0:lambda x:x*(x*x+1),-1:lambda x:x*(x*x-2),
               -2:lambda x:x*(x*x-2)*(x*x-3)}
    rows=[]
    for c,vertices in portraits.items():
        def portrait(x):
            return math.prod(x-r for r in vertices)
        edges=[(x,x*x+c) for x in vertices]
        check(all(y in vertices for _,y in edges),'finite portrait invariant')
        for x in range(-30,31):
            check(portrait(x*x+c)==portrait(x)*quotients[c](x),'portrait divisor pullback')
        rows.append(dict(c=c,edges=edges))
    return dict(finite_affine_quotients=finite,quadratic_portraits=rows)


def main():
    report=dict(status='PROVED scoped controllers and shell laws; FINITE-EXACT audits; Collatz OPEN')
    report['runs']=anchored_runs()
    report['switches']=switch_shells()
    report['groups_and_portraits']=group_and_portrait()
    certificates,density_rows=paid_families()
    report['paid_certificates']=certificates
    report['safe_cylinders']=density_rows
    report['rational_bank']=rational_cycle_bank()
    report['power_three_repairs']=powers_of_three(report['rational_bank'])
    report['switched_bank']=switched_bank()
    report['coverage']=coverage(density_rows,report['rational_bank'],report['switched_bank'])
    report['checks']=CHECKS
    stem=ROOT/'05-knowledge/results/collatz_paid_portrait_controllers_20261004'
    stem.with_suffix('.json').write_text(json.dumps(report,indent=2)+'\n')
    lines=['PAID PORTRAIT CONTROLLERS: universal coverage OPEN',
           'RUNS '+json.dumps(report['runs']),
           'SWITCHES '+json.dumps({k:v for k,v in report['switches'].items() if k not in ('separation','examples')}),
           'PAID FAMILIES '+str(len(certificates)),
           'POWER THREE REPAIRS '+json.dumps(report['power_three_repairs']),
           'COVERAGE '+json.dumps(report['coverage']),
           'GROUPS AND PORTRAITS '+json.dumps(report['groups_and_portraits']),
           f'PASS: {CHECKS} explicit checks.']
    output='\n'.join(lines)+'\n'
    stem.with_suffix('.out').write_text(output)
    print(output,end='')


if __name__=='__main__':
    main()
