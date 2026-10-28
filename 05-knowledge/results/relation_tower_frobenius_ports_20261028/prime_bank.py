"""Finite rank-11 square-sigma realization, with recursive Lucas certificates.

SymPy discovers factors/primes; verify_certificate uses only exact integers.
No claim of prime supply beyond this finite bank.
"""
from math import gcd, prod, isqrt
from pathlib import Path
import json
from sympy import factorint, isprime

BASE=[2,3,5,7,11,13,17,19,23,29,31]
CERT={2:{'factors':{},'witness':1}}
def certificate(p):
    if p in CERT: return
    factors={int(q):int(e) for q,e in factorint(p-1).items()}
    for q in factors: certificate(q)
    a=2
    while not (pow(a,p-1,p)==1 and all(gcd(pow(a,(p-1)//q,p)-1,p)==1 for q in factors)):
        a+=1
        if a>=p: raise RuntimeError('prime certificate failure')
    CERT[p]={'factors':factors,'witness':a}

def verify_certificate(cert):
    verified={2}
    if cert[2]!={'factors':{},'witness':1}: raise RuntimeError('base certificate')
    for p in sorted(set(cert)-{2}):
        row=cert[p]; factors={int(q):e for q,e in row['factors'].items()}
        if not all(q in verified and isinstance(e,int) and e>0 for q,e in factors.items()): raise RuntimeError('unproved factor')
        if prod(q**e for q,e in factors.items())!=p-1: raise RuntimeError('factorization')
        a=row['witness']
        if pow(a,p-1,p)!=1 or any(gcd(pow(a,(p-1)//q,p)-1,p)!=1 for q in factors): raise RuntimeError('Lucas witness')
        verified.add(p)
    return verified

def main():
    rows=[]
    for mask in range(1,1<<len(BASE)):
        d=prod(p for i,p in enumerate(BASE) if mask>>i&1)
        a=1
        while not isprime(d*a*a-1): a+=1
        p=d*a*a-1
        certificate(p)
        rows.append({'mask':mask,'squarefree_part':d,'a':a,'prime':p})
    proved=verify_certificate(CERT)
    if not all(row['prime'] in proved for row in rows): raise RuntimeError('uncertified bank member')
    if len({row['prime'] for row in rows})!=2047: raise RuntimeError('repeated prime')
    for row in rows:
        if row['prime']+1!=row['squarefree_part']*row['a']**2: raise RuntimeError('realization')
    out={'status':'FINITE-EXACT','coordinate_primes':BASE,'rank':11,'length':2047,'kernel_dimension':2036,'largest_prime':max(row['prime'] for row in rows),'largest_multiplier':max(row['a'] for row in rows),'rows':rows,'prime_certificates':CERT}
    Path(__file__).with_name('prime_bank.json').write_text(json.dumps(out,indent=2)+'\n',encoding='utf-8')
    print(json.dumps({k:v for k,v in out.items() if k not in ('rows','prime_certificates')}))
    print('Recursive exact primality certificates:',len(CERT))

if __name__=='__main__': main()
