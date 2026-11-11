"""Exact square-plus-one sieve with a short Lucas certificate for every prime.

The indexed variable is the square root r, not the prime's sequence index.
No probabilistic primality test and no prime-infinitude assumption is used.
"""
from array import array
from collections import Counter
from hashlib import sha256
from math import gcd,isqrt
from pathlib import Path
import argparse
import json
import maxima


def need(ok,message):
    if not ok:
        raise RuntimeError(message)


def prime_factors(n,spf):
    result=[]
    while n>1:
        p=spf[n]
        result.append(p)
        while n%p==0:
            n//=p
    return result


def quadratic_sieve(bound):
    prime=maxima.prime_sieve(bound)
    least=array('I',[0])*(bound+1)
    active_primes=0
    for p in range(5,bound+1,4):
        if not prime[p]:
            continue
        active_primes+=1
        # If a is a nonresidue, r=a^((p-1)/4) is a square root of -1.
        a=2
        while pow(a,(p-1)//2,p)!=p-1:
            a+=1
        r=pow(a,(p-1)//4,p)
        need((r*r+1)%p==0,'modular root')
        for residue in (r,p-r):
            first=residue if residue%2==0 else residue+p
            if first*first+1==p:
                first+=2*p
            for n in range(first,bound+1,2*p):
                if least[n]==0:
                    least[n]=p
    return least,active_primes


def lucas_certificate(r,spf):
    p=r*r+1
    factors=prime_factors(r,spf)
    a=2
    while a<p:
        if pow(a,p-1,p)==1 and all(gcd(pow(a,(p-1)//q,p)-1,p)==1 for q in factors):
            return a,factors
        a+=1
    raise RuntimeError(f'No Lucas witness for {p}')


def main():
    ap=argparse.ArgumentParser()
    ap.add_argument('--root-bound',type=int,default=1_000_000)
    ap.add_argument('--output',type=Path)
    ap.add_argument('--table',type=Path)
    args=ap.parse_args()
    N=args.root_bound
    need(N>=2,'root bound >=2')
    least,active=quadratic_sieve(N)
    spf=maxima.spf_sieve(N)
    values=[2]
    rows=[]
    lucas_hash=sha256(b'1 2 special\n')
    factor_hash=sha256()
    bases=Counter()
    max_base=(0,0)
    for r in range(2,N+1,2):
        p=r*r+1
        if least[r]:
            q=int(least[r])
            need(1<q<p and p%q==0,'proper divisor certificate')
            factor_hash.update(f'{r} {q}\n'.encode())
        else:
            a,qs=lucas_certificate(r,spf)
            # Full factorization p-1=r^2; every prime q is independently
            # certified by the small Eratosthenes/SPF table, no big MR test.
            need(pow(a,p-1,p)==1,'Fermat equality')
            need(all(gcd(pow(a,(p-1)//q,p)-1,p)==1 for q in qs),'full order')
            values.append(p)
            bases[a]+=1
            if a>max_base[0]:max_base=(a,p)
            lucas_hash.update(f'{r} {p} {a} {qs}\n'.encode())
            if len(rows)<15:rows.append((r,p,a,qs))
    # Independent small primality path uses literal division, not root sieving.
    trial_limit=min(N,10_000)
    trial=[]
    for r in range(1,trial_limit+1):
        p=r*r+1
        if all(p%d for d in range(2,isqrt(p)+1)):
            trial.append(p)
    need(trial==[p for p in values if p<=trial_limit*trial_limit+1],'trial-division prefix')
    # Cross-sequence extremal equivalence, over the actual computed A135244 range.
    maximum,_,_=maxima.enumerate_maxima(1_000_000)
    bridge=[]
    for r in range(1,isqrt(999_998)+1):
        p=r*r+1
        is_quadratic_prime=p in values
        need((maximum[p+1]==p*p)==is_quadratic_prime,'aliquot extremal equivalence')
        if is_quadratic_prime:bridge.append((r,p,int(maximum[p+1])))
    text=''.join(f'{i} {p}\n' for i,p in enumerate(values,1))
    if args.table:args.table.write_text(text,encoding='utf-8')
    receipt={
        'status':'FINITE-EXACT; A002496 infinitude remains OPEN',
        'root_bound':N,'value_bound':N*N+1,'prime_count':len(values),
        'active_modulus_primes':active,
        'certified_composites':N//2-(len(values)-1),
        'prime_table_sha256':sha256(text.encode()).hexdigest(),
        'lucas_certificate_stream_sha256':lucas_hash.hexdigest(),
        'composite_factor_stream_sha256':factor_hash.hexdigest(),
        'largest_lucas_base_and_prime':max_base,
        'first15_lucas_certificates':rows,
        'last10':values[-10:],
        'trial_division_root_bound':trial_limit,
        'trial_division_prime_count':len(trial),
        'aliquot_bridge_root_bound':isqrt(999_998),
        'aliquot_bridge_prime_count':len(bridge),
        'aliquot_bridge_first10':bridge[:10],
    }
    encoded=json.dumps(receipt,indent=2,sort_keys=True)+'\n'
    if args.output:args.output.write_text(encoded,encoding='utf-8')
    print(encoded)


if __name__=='__main__':main()
