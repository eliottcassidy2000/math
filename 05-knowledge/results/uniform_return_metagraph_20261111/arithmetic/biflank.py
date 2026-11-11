"""Exact two-growing-flank family and a bounded balanced-offset hostile probe."""
from collections import defaultdict
from math import isqrt
from pathlib import Path
import importlib.util,json
import chart_inverse as inv

HERE=Path(__file__).resolve().parent


def phi(coeff,A,B):
    h,d,L,i,o,s,e=coeff
    return h+d*(inv.g(A)+inv.g(B)-2)+L*A*B+i*A+o*B+s*(A&1)+e*(B&1)


def main():
    source=inv.RESULTS/'session251_creative_completion_20261111/ternary/compiler.py'
    spec=importlib.util.spec_from_file_location('signed',source)
    oracle=importlib.util.module_from_spec(spec);spec.loader.exec_module(oracle)
    positive=[];checks=0
    for record in inv.raw_records():
        if record['cut']:continue
        h,d,L,i,o,s,e=record['positive_observer']
        inv.require((h,i,o,s,e)==(d,L-d,0,d,0),'uniform-cut identity')
        sig=oracle.signature(inv.CORES[record['core_index']],0)
        for alpha in (1,2,3,5,17,51):
            for beta in (1,2,3,5,17,51):
                C=d*(inv.g(alpha)-alpha+(alpha&1)+inv.g(beta-1)-1)
                for a in range(1,6):
                    for b in range(1,6):
                        A=alpha*2**a;B=beta*2**b-1
                        target=L*alpha*beta*2**(a+b)+C
                        inv.require(phi(record['positive_observer'],A,B)==target,'uniform product law')
                        inv.require(oracle.count(sig,A,B)==target,'independent signed product law')
                        checks+=1
        positive.append({'core':record['core_index'],'d':d,'L':L})
    C102=[{'alpha':a,'beta':51//a,
           'constant':inv.g(a)-a+(a&1)+inv.g(51//a-1)-1} for a in (1,3,17,51)]
    inv.require([x['constant'] for x in C102]==[37,15,-5,-39],'102 coefficient constants')
    # In the balanced probe A=alpha*2^m+t,B=beta*2^m+u, the fixed-offset
    # lemma gives g(A)=gamma_A*2^m+kappa_A for every m>=M. These exact
    # coefficients are computed by differences; a third point is a control.
    offset_bound=4096;M=14;D=2**M;cache={}
    profile_checks=0
    def profiles(a):
        nonlocal profile_checks
        if a not in cache:
            data=[]
            for t in range(-offset_bound,offset_bound+1):
                g0=inv.g(a*D+t);g1=inv.g(a*2*D+t)
                gammaD=g1-g0;kappa=2*g0-g1
                inv.require(inv.g(a*4*D+t)==4*gammaD+kappa,'fixed-offset affine control')
                data.append((t,gammaD,kappa));profile_checks+=1
            cache[a]=data
        return cache[a]
    balanced=[]
    for j in range(10):
        lead=102*2**j;hits=[];zero_drift=0;charts=0;factor_pairs=0
        for record in inv.raw_records():
            h,d,L,i,o,s,e=record['positive_observer']
            if lead%L:continue
            charts+=1;product=lead//L
            small=[a for a in range(1,isqrt(product)+1) if product%a==0]
            divisors=sorted(set(small+[product//a for a in small]))
            for alpha in divisors:
                beta=product//alpha;factor_pairs+=1
                bucket=defaultdict(list)
                for u,gb,kb in profiles(beta):bucket[(alpha*L*u)*D+d*gb].append((u,kb))
                for t,ga,ka in profiles(alpha):
                    key=-(alpha*i+beta*(L*t+o))*D-d*ga
                    for u,kb in bucket.get(key,[]):
                        zero_drift+=1
                        C=h+d*(ka+kb-2)+L*t*u+i*t+o*u+s*(t&1)+e*(u&1)
                        if C==19:hits.append({'core':record['core_index'],'cut':record['cut'],
                                             'alpha':alpha,'beta':beta,'t':t,'u':u})
        balanced.append({'j':j,'leading_coefficient':lead,'eligible_raw_charts':charts,
                         'ordered_factor_pairs':factor_pairs,'zero_linear_drift':zero_drift,'hits':hits})
    report={'status':'PROVED independent-product family; FINITE-EXACT balanced-offset probe',
        'independent_signed_product_checks':checks,'uniform_core_coefficients':positive,
        'coefficient102_odd_factor_pairs':C102,
        'balanced_probe':{'scope':'All90 raw cuts; alpha,beta positive with L*alpha*beta=102*2^j; j0..9; offsets t,u in[-4096,4096]; eventual identity for every m>=14.',
                          'offset_bound':offset_bound,'first_m':M,'profile_controls':profile_checks,'rows':balanced},
        'stopping_reason':'No +19 balanced fixed-offset identity in this explicit finite box. A growing offset, larger shift, different core grammar or source-changing operation remains possible.'}
    out=json.dumps(report,indent=2,sort_keys=True)+'\n'
    (HERE/'biflank.json').write_text(out,encoding='utf-8')
    print(json.dumps({'status':report['status'],'product_checks':checks,'balanced_hits':sum(len(r['hits']) for r in balanced)}))


if __name__=='__main__':main()
