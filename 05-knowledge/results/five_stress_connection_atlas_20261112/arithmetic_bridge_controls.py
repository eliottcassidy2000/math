"""Independent small controls on recovered arithmetic bridges, not new theorems.

Sources: endpoint_coverage_compiler_20261111/CONNECTIONS.md E7/E9;
catalan_elliptic_20260921_catalan.md; Session260 universal congruence.
"""
from pathlib import Path
from itertools import combinations,permutations
from math import gcd,lcm,isqrt
import json
HERE=Path(__file__).resolve().parent
def need(x,msg):
    if not x:raise RuntimeError(msg)
def prime(n):return n>=2 and all(n%d for d in range(2,isqrt(n)+1))
def order(a,p):return next(k for k in range(1,p) if pow(a,k,p)==1)
def crt(a,m,b,n):
    g=gcd(m,n)
    if (b-a)%g:return None
    q=n//g
    t=((b-a)//g*pow(m//g,-1,q))%q if q>1 else 0
    return ((a+m*t)%lcm(m,n),lcm(m,n))
def endpoint(rows):
    n=len(rows);H=0;starts=[0]*n;B=[[0]*n for _ in rows]
    for p in permutations(range(n)):
        legal=[bool(rows[p[i]]>>p[i+1]&1) for i in range(n-1)]
        if all(legal):H+=1;starts[p[0]]+=1
        for cut in range(1,n):
            if all(legal[:cut-1]+legal[cut:]):B[p[cut-1]][p[cut]]+=1
    return H,B[0][3]+B[3][0],[starts[j]+sum(B[i][j] for i in range(n)) for j in [0,3]]
def run():
    primes=[p for p in range(2,100) if prime(p)];need(len(primes)==25,'prime cutoff')
    classes=[];absent=[]
    for p in primes:
        if p==2 or 99%p==0:absent.append(p);continue
        m=order(2,p);hits=[k for k in range(m) if (99*pow(2,k,p)+23)%p==0]
        need(len(hits)<=1,'single discrete-log class')
        if hits:classes.append((p,m,hits[0]))
        else:absent.append(p)
    period=lcm(*(m for p,m,a in classes));survivors=period
    for size in range(1,len(classes)+1):
        for subset in combinations(classes,size):
            pair=(0,1)
            for p,m,a in subset:
                pair=crt(*pair,a,m)
                if pair is None:break
            if pair:survivors+=(-1)**size*(period//pair[1])
    need((period,survivors)==(428468040,226472946),'recovered exact wheel')
    for k in range(1,5001):
        need(all((99*pow(2,k,p)+23)%p for p in primes)==all(k%m!=a for p,m,a in classes),'direct exponent control')
    need(prime(419) and order(2,419)==418,'composite progression divisor')
    M=lcm(period,418);need(M==8140892760,'progression period')
    need(all((99*pow(2,2,p)+23)%p for p in primes),'progression survives wheel')
    rows=[40,17,3,6,13,30];h,w,u=endpoint(rows)
    rev=[sum(1<<j for j in range(6) if rows[j]>>i&1) for i in range(6)]
    hr,wr,ur=endpoint(rev)
    need((h,w,u)==(37,14,[46,53]) and (hr,wr,sum(ur))==(37,14,117),'literal endpoint enumeration')
    def carry(exponents):
        K=0;B=0
        for a in exponents:B=3*B+2**K;K+=a
        return {'exponents':exponents,'K':K,'L':len(exponents),'B':B,'gap':2**K-3**len(exponents)}
    ca=carry([1,3]);cb=carry([2,2]);neg=carry([1,1,1,2,1,1,4])
    need((ca['B'],cb['B'],ca['gap'],cb['gap'])==(5,7,7,7),'same clock different carry')
    need((neg['B'],neg['gap'])==(2363,-139),'negative cycle nonunit gap')
    n=-17;orbit=[n]
    for a in neg['exponents']:
        numerator=3*n+1;need(numerator%(2**a)==0 and numerator%(2**(a+1))!=0,'exact valuation')
        n=numerator//2**a;orbit.append(n)
    need(n==-17,'literal signed cycle')
    congruences={m:all(pow(a,m,m)==a for a in range(m)) for m in [7,9,561]}
    need(congruences=={7:True,9:False,561:True},'prime/composite congruence distinction')
    out={'scope':'Independent finite controls of explicitly credited prior mechanisms; no new wheel, Catalan theorem, Carmichael classification or tournament ray claimed.','prime_wheel':{'primes':primes,'classes_prime_period_residue':classes,'absent':absent,'period':period,'survivors':survivors,'direct_exponents':5000,'composite_progression':{'start':2,'step':M,'proper_divisor':419}},'base_endpoint_state':{'rows':rows,'selected':[0,3],'H':h,'pair_tax':w,'selected_U':u,'reverse_U_sum':sum(ur)},'carry_controls':[ca,cb,neg],'signed_cycle':orbit,'universal_power_congruence_controls':congruences}
    tag='.optimized' if not __debug__ else ''
    (HERE/('arithmetic_bridge_controls'+tag+'.json')).write_text(json.dumps(out,indent=2)+'\n')
    print(json.dumps({'status':'PASS','period':period,'survivors':survivors,'endpoint_state':[h,w,u],'signed_cycle':orbit}),flush=True)
if __name__=='__main__':run()
