"""Complete near-square intersection of branching root defects at most six.

Source-language completeness is inherited from the capped-DFA classification.
Coefficients are extracted from literal kernel path tables, not fitted rays.
Modular certificates and finite factor-pair equations cover all exponents.
"""
from hashlib import sha256
from math import isqrt
from pathlib import Path
import json
import audit_defect as audit
import maxima
import quadratic_primes as qp

HERE=Path(__file__).resolve().parent
SOURCE=HERE.parent/'defect/two_inport.json'
MODULI=(16,3,5,7,13,17,31)
PERIOD=120


def coefficients(tail):
    rows=audit.graph(tail);n=len(rows);r=n-1;full=(1<<n)-1
    starts=audit.starts(rows)
    ends=audit.starts([sum(((rows[w]>>v)&1)<<w for w in range(n)) for v in range(n)])
    H=[1]+[sum(x) for x in starts[1:]]
    F=sum(ends[S][r]*H[full^S] for S in range(1,full+1) if S>>r&1)
    q=H[full^(1<<r)];s=starts[full][r]
    E=F-q;B=H[full]+q-s-E
    for k in (1,2):
        actual=sum(audit.starts(audit.graph((3,)*k+tail))[-1])
        audit.check(actual==E*2**k+B,'independent literal prefix transfer')
    return E,B,{'H':H[full],'F':F,'q':q,'s':s}


def square_difference(C,c):
    """All k>=1,n>0 with n²-c²=C*2^k, for c=2 or4, by finite factor pairs."""
    audit.check(c in (2,4),'power-two offset domain')
    t=c.bit_length()-1;odd=C;v=0
    while odd%2==0:odd//=2;v+=1
    solutions=set();cert=[]
    for k in range(1,max(0,2*t+2-v)+1):
        n=isqrt(C*2**k+c*c)
        if n*n==C*2**k+c*c:solutions.add((k,n))
    for a in range(1,odd+1):
        if odd%a:continue
        b=odd//a
        for sign in (-1,1):
            numerator=a+sign
            if numerator<=0 or numerator%b:continue
            power=numerator//b
            if power<2 or power&(power-1):continue
            K=2*t+2+(power.bit_length()-1);k=K-v
            if k<1:continue
            first=2**(t+1)*a;second=2**(K-t-1)*b
            audit.check(abs(first-second)==2*c,'factor difference')
            n=(first+second)//2
            audit.check(n*n-c*c==C*2**k,'factor solution')
            solutions.add((k,n));cert.append({'a':a,'b':b,'sign':sign,'power':power,'K':K,'k':k,'n':n})
    return sorted(solutions),cert


def main():
    for p in MODULI[1:]:audit.check(pow(2,PERIOD,p)==1,'certified odd-modulus period')
    residues={p:{n*n%p for n in range(p)} for p in MODULI}
    data=json.loads(SOURCE.read_text());rays=[];isolated=[];finite=set();literal_checks=0
    for auto in data['automata']:
        budget=auto['budget']
        if budget>6:continue
        for regex in auto['exact_disjoint_path_regex']:
            if not any(x.endswith('*') for x in regex):
                word=tuple(map(int,regex));H=sum(audit.starts(audit.graph(word))[-1]);n=isqrt(H-1)
                hit=n*n==H-1
                isolated.append({'budget':budget,'word':word,'H':H,'square_plus_one':hit})
                if hit:finite.add(H)
                continue
            audit.check(regex[:2]==['3','3*'] and not any(x.endswith('*') for x in regex[2:]),'one uniform prefix')
            tail=tuple(map(int,regex[2:]));E,B,kernel=coefficients(tail);literal_checks+=2
            blockers=[];allowed=[]
            for residue in range(PERIOD):
                k=residue+PERIOD
                bad=next((p for p in MODULI if (E*pow(2,k,p)+B-1)%p not in residues[p]),0)
                blockers.append(bad)
                if not bad:allowed.append(residue)
            small=[]
            for k in range(1,4):
                H=E*2**k+B;n=isqrt(H-1)
                if n*n==H-1:small.append((k,n))
            factor=[];infinite=None
            if not allowed:
                solutions=small
            elif B==1:
                odd=E;v=0
                while odd%2==0:odd//=2;v+=1
                audit.check(isqrt(odd)**2==odd,'remaining zero-offset square odd part')
                infinite={'k_mod2':v%2,'first_k':1 if v%2 else 2,'odd_square_root':isqrt(odd),'v2_E':v}
                solutions=[]
            else:
                c=isqrt(B-1)
                audit.check(c*c==B-1 and c in (2,4),'remaining finite factor shape')
                solutions,factor=square_difference(E,c)
            hits=[{'k':k,'n':n,'H':n*n+1,'word':(3,)*k+tail} for k,n in solutions]
            for hit in hits:finite.add(hit['H'])
            rays.append({'budget':budget,'tail':tail,'E':E,'B':B,'kernel':kernel,
                         'allowed_k_mod120_for_k_ge4':allowed,'first_blocker_by_k_mod120':blockers,
                         'finite_hits':hits,'infinite_square_family':infinite,'factor_certificate':factor})
    audit.check(len(rays)==52 and len(isolated)==8,'complete low-defect regex inventory')
    infinite=[r for r in rays if r['infinite_square_family']]
    audit.check(len(infinite)==1 and (infinite[0]['budget'],infinite[0]['tail'],infinite[0]['E'],infinite[0]['B'])==(4,(1,3),16,1),'unique infinite near-square tail')
    # The ray contributes4^m+1,m>=3; finite5,17 complete m1,2.
    exceptions=sorted(finite-{5,17,65})
    audit.check(exceptions==[37,101,145,485,3845,4357,10001,15877,35345],'exact near-square exception spectrum')
    primes=[];spf=maxima.spf_sieve(126)
    for H in (37,101,4357,15877):
        n=isqrt(H-1);base,factors=qp.lucas_certificate(n,spf)
        primes.append({'H':H,'n':n,'Lucas_base':base,'prime_divisors_H_minus_1':factors})
    composite={145:5,485:5,3845:5,10001:73,35345:5}
    for H,p in composite.items():audit.check(1<p<H and H%p==0,'sporadic composite certificate')
    report={'status':'PROVED complete classified-spectrum intersection; exact finite modular and factor certificates',
        'scope':'Strong rooted tournaments with indegree>=2 and positive defect<=6. d1 is excluded. Completeness of source families is inherited from defect/PROOF.md.',
        'source_normalized_sha256':sha256(SOURCE.read_bytes().replace(b'\r\n',b'\n')).hexdigest(),
        'moduli':MODULI,'period':PERIOD,'literal_transfer_checks':literal_checks,
        'rays':rays,'isolated':isolated,
        'full_spectrum':{'infinite_family':'4^m+1, every m>=1','additional_values':exceptions},
        'prime_spectrum':'All prime members of4^m+1,m>=1 (equivalently Fermat primes index>=1), plus37,101,4357,15877.',
        'sporadic_prime_certificates':primes,'sporadic_composite_factors':composite,
        'excluded_other_sectors':'d3 minimum6 collars and d2/C3-residual6 have no square-plus-one values by MINIMUM_SQUARES.md.'}
    out=json.dumps(report,indent=2,sort_keys=True)+'\n'
    (HERE/'tail_squares.json').write_text(out,encoding='utf-8')
    print(json.dumps({'status':report['status'],'rays':len(rays),'isolates':len(isolated),'additional_values':exceptions}))


if __name__=='__main__':main()
