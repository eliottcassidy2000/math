"""Certify attainable vectors of finite-field path-cover readings with ordinal products."""
from math import comb, factorial
from pathlib import Path
import sys

ROOT = Path(__file__).resolve().parents[3]
sys.path.insert(0,str(ROOT / '05-knowledge/results/chanpak_gap_star_ports_20261122'))
from replay import subpaths


def profile(rows):
    """Counts of unordered partitions into directed paths, indexed by their size."""
    n = len(rows)
    H = subpaths(rows)
    dp = [[0]*(n+1) for _ in H]
    dp[0][0] = 1
    for mask in range(1,1<<n):
        fixed = mask & -mask
        rest = mask ^ fixed
        part = rest
        while True:
            sub = part | fixed
            for k in range((mask ^ sub).bit_count()+1):
                dp[mask][k+1] += H[sub]*dp[mask ^ sub][k]
            if not part:
                break
            part = (part-1)&rest
    return tuple(dp[-1][k]*factorial(k) for k in range(n+1))


def evaluations(runs, prime):
    n = len(runs)-1
    # Odd-order seeds are padded by one ordinal singleton to put all factors in even parity.
    return tuple(pow(j, n % 2, prime)*sum(comb(j,k)*runs[k] for k in range(1,min(n,j)+1)) % prime
                 for j in range(1,(prime+1)//2))


def rank_insert(basis, row, prime):
    row = [v % prime for v in row]
    for i in sorted(basis):
        if row[i]:
            row = [(v-row[i]*w)%prime for v,w in zip(row,basis[i])]
    if not any(row):
        return
    i = next(i for i,v in enumerate(row) if v)
    inv = pow(row[i],-1,prime)
    basis[i] = [(v*inv)%prime for v in row]


def logs(p):
    for gen in range(2,p):
        table={pow(gen,j,p):j for j in range(p-1)}
        if len(table)==p-1:
            return table
    raise ValueError('no generator')


def certify(primes, bank):
    columns = [p for p in primes for _ in range((p-1)//2)]
    columns_by_prime = {q:[j for j,p in enumerate(columns) if (p-1)%q==0]
                        for q in range(2,max(primes)) if all(q%d for d in range(2,int(q**.5)+1))
                        and any((p-1)%q==0 for p in primes)}
    bases = {q:{} for q in columns_by_prime}
    tables = {p:logs(p) for p in primes}
    zeros=set()
    for n, run in bank:
        row=[a for p in primes for a in evaluations(run,p)]
        missing=[i for i,a in enumerate(row) if not a]
        if len(missing)==1:
            zeros.add(missing[0])
        elif not missing:
            vec=[tables[p][a] for a,p in zip(row,columns)]
            for q,ix in columns_by_prime.items():
                if len(bases[q])<len(ix):
                    rank_insert(bases[q],[vec[j] for j in ix],q)
        if len(zeros)==len(columns) and all(len(bases[q])==len(v) for q,v in columns_by_prime.items()):
            return n, len(columns)
    return None,len(zeros),len(columns),[(q,len(bases[q]),len(ix)) for q,ix in columns_by_prime.items()]


def profiles(max_order):
    pool,seen=[],set()
    for n in range(1,max_order+1):
        p=ROOT/'session212'/'cycle-levels'/f'classes{n}.rows'
        for line in p.read_text().splitlines():
            coeff=profile(tuple(map(int,line.split())))
            if coeff not in seen:
                seen.add(coeff)
                pool.append((n,coeff))
    return pool


if __name__=='__main__':
    order = int(sys.argv[1]) if len(sys.argv)>1 else 8
    bank=profiles(order)
    for ps in ((3,), (5,), (7,), (11,), (13,), (17,), (19,), (31,), (41,),
               (3,5,7,13), (3,5,7,11,13), (3,5,7,11,13,17,19)):
        result = certify(ps,bank)
        if order == 8:
            assert len(result) == 2 and result[0] is not None
        print(ps,result,flush=True)
