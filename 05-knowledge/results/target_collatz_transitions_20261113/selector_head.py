"""Complete all-order head of the polarized constructor below39971."""
from pathlib import Path
from math import isqrt
import gzip,json
from polarized import primitive_guard,need,arithmetic_candidates

HERE=Path(__file__).resolve().parent

def prime(n):
    return n>=2 and all(n%d for d in range(2,isqrt(n)+1))

def main():
    bank=json.loads(gzip.decompress((HERE/'response_bank.json.gz').read_bytes()))
    small=[v for v in bank if v['p']<=512]
    targets={}
    for M in small:
        h,p=M['h'],M['p']
        for N in small:
            j,q=N['h'],N['p'];Y=36*p*q+6*h*q+18*p*j-h*j
            if Y>=39971:continue
            item=dict(M=M,N=N,primitive=primitive_guard(h,p,j,q))
            targets.setdefault(Y,[]).append(item)
    need(len(targets)==205,'complete inherited-bank head')
    need(881 not in targets and not arithmetic_candidates(881),'arithmetic exclusion881')
    for Y in (191,347,371,5299):
        candidates=arithmetic_candidates(Y)
        need(any(all(item['M'][k]==c[v] for k,v in (('h','h'),('p','p'))) and
                 all(item['N'][k]==c[v] for k,v in (('h','j'),('p','q')))
                 for item in targets[Y] for c in candidates),'factor inverse reconstructs witness')
    records=[dict(Y=Y,prime=prime(Y),seven_prime=(Y%7==0 and prime(Y//7)),witnesses=items)
             for Y,items in sorted(targets.items())]
    result=dict(status='FINITE-EXACT with proved all-order cutoff',scope='polarized nonempty-module constructor; allY<39971, inherited completeordinarymodules≤8; no global spectrum exclusion',targets=len(targets),prime_targets=sum(x['prime'] for x in records),seven_prime_targets=sum(x['seven_prime'] for x in records),records=records)
    (HERE/'selector_head.json').write_text(json.dumps(result,indent=2)+'\n')
    print({k:v for k,v in result.items() if k!='records'})

if __name__=='__main__':main()
