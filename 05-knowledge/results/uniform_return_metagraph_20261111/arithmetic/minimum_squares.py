"""Exact minimum-root-defect / square-plus-one / Fermat-candidate interface."""
from math import isqrt
from pathlib import Path
import json
import audit_defect as audit
import maxima
import quadratic_primes as qp

HERE=Path(__file__).resolve().parent


def collar(d,k):
    # Own gauge: chain first, then d transitive inports, then root.
    root=k+d;rows=[0]*(root+1)
    for v in range(k+d):
        rows[v]=sum(1<<w for w in range(v+1,k+d))
    for v in range(k,k+d):rows[v]|=1<<root
    rows[root]=(1<<k)-1
    return rows


def main():
    arithmetic=0
    for d in range(3,129):
        for k in range(1,129):
            H=(2**d-1)*(2**k-1)+2
            r=isqrt(H-1)
            audit.check((r*r==H-1)==(k==1 and d%2==0),'all-degree square intersection')
            arithmetic+=1
    for k in range(1,129):
        H=9*2**k-3
        audit.check(H%3==0 and isqrt(H-1)**2!=H-1,'triangle exclusion')
    literal=[]
    for d in range(3,9):
        for k in range(1,4):
            rows=collar(d,k)
            H=sum(audit.starts(rows)[-1])
            audit.check(H==(2**d-1)*(2**k-1)+2,'independent literal collar count')
            literal.append({'d':d,'k':k,'order':len(rows),'H':H})
    spf=maxima.spf_sieve(256)
    examples=[]
    for j in (2,3,4):
        d=2**j;r=2**(d//2);H=r*r+1
        a,factors=qp.lucas_certificate(r,spf)
        examples.append({'Fermat_index':j,'d':d,'H':H,'square_root':r,'defect':H-3,
                         'source_rows':collar(d,1),'root':d+1,
                         'Lucas_base':a,'complete_prime_divisors_of_H_minus_1':factors})
    residual=[]
    for E,B in ((3,-1),(6,1),(5,-1),(13,3),(14,3)):
        for k in range(129):
            H=12*E*2**k+3*B
            audit.check(H%3==0 and isqrt(H-1)**2!=H-1,'residual-C3 square exclusion')
        residual.append({'E':E,'B':B,'leading':12*E,'constant':3*B})
    out={'status':'PROVED classified-stratum iff; FINITE-EXACT controls',
         'arithmetic_transitive_collar_cases':arithmetic,'triangle_cases':128,
         'literal_controls':literal,'Fermat_prime_examples':examples,
         'residual_C3_excluded_rays':residual,
         'scope':'All-degree classification is inherited; congruence proof gives iff. Prime members are exactly Fermat primes of index>=2, without asserting their infinitude or that the examples exhaust them.'}
    text=json.dumps(out,indent=2,sort_keys=True)+'\n'
    (HERE/'minimum_squares.json').write_text(text,encoding='utf-8')
    print(json.dumps({'status':out['status'],'arithmetic_cases':arithmetic,'literal_cases':len(literal)}))


if __name__=='__main__':main()
