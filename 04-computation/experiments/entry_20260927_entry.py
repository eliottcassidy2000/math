"""Exact entry reductions and a terminating recursive certificate interpreter."""
from collections import Counter
from functools import lru_cache
from pathlib import Path

ROOT=Path(__file__).resolve().parents[2]


def need(test,label):
    if not test:
        raise RuntimeError(label)


def v2(n):
    return (n & -n).bit_length()-1


def U(n):
    z=3*n+1
    return z >> v2(z)


def iterate(n,j):
    for _ in range(j):
        n=U(n)
    return n


def read_bank():
    text=(ROOT/'05-knowledge/results/reset_20260926_swaplift.out').read_text(encoding='utf-8-sig')
    rows=[tuple(map(int,line.split())) for line in text.splitlines()
          if len(line.split())==9 and all(s.isdecimal() for s in line.split())]
    need(len(rows)==171,'inherited bank size')
    return rows


ROWS=read_bank()


@lru_cache(None)
def rule(n):
    if n==1:
        return 1,0,'terminal'
    if n%4==1:
        target=U(n)
        need(target<n,('quarter descent',n))
        return target,1,'quarter'
    for q,J,A,P,b,clip,R,residue,K in ROWS:
        if n%(1 << K)==residue:
            target=iterate(n,J+1)
            need(target<n,('bank descent',n,q))
            return target,J+1,'bank'
    return None,0,'unmatched'


def certify(n,fuel):
    """Total for each input: every move strictly lowers the integer n*2**fuel."""
    while n!=1:
        old_rank=n << fuel
        target,steps,tag=rule(n)
        if target is not None:
            n=target
        else:
            if fuel==0:
                return False
            fuel-=1
            n=U(n)
        need((n << fuel)<old_rank,'recursive search rank')
    return True


def main():
    entry_counts=Counter()
    for source in range(3,32768,2):
        n=source
        while n>1 and n%4==1:
            n=U(n)
        if n==1:
            entry_counts['terminal_before_rise']+=1
        else:
            k=v2(n+1)
            need(k>=2,'rise begins')
            landing=iterate(n,k-2)
            need(landing==3**(k-2)*((n+1) >> (k-2))-1 and landing%8==3,
                 ('entry3mod8',source))
            need(v2(3*(landing-2)+1)==2, ('fixedq1R1',source))
            entry_counts['entry3mod8']+=1
        k=v2(source+1)
        landing=iterate(source,k-1)
        need(landing%4==1,('entry local descending region',source))
        if landing>1:
            need(U(landing)<landing,('local cert',source))

    # Exact policy cost on a frozen finite universe; cap is not a conjecture.
    debt={1:0}
    for source in range(3,100000,2):
        n=source
        path=[]
        seen=set()
        for _ in range(1024):
            if n in debt:
                break
            need(n not in seen,('finite policy cycle',source,n))
            seen.add(n)
            target,steps,tag=rule(n)
            charge=int(target is None)
            path.append((n,charge))
            n=U(n) if charge else target
        else:
            raise RuntimeError(('finite policy cap',source))
        cost=debt[n]
        for point,charge in reversed(path):
            cost+=charge
            debt[point]=cost
    levels=(0,1,2,4,8,16,32,64)
    counts={d:sum(debt[n]<=d for n in range(3,100000,2)) for d in levels}
    max_debt=max(debt[n] for n in range(3,100000,2))
    first_max=next(n for n in range(3,100000,2) if debt[n]==max_debt)
    for n in range(3,4096,2):
        for d in levels:
            need(certify(n,d)==(debt[n]<=d),('independent bounded interpreter',n,d))

    Kmax=max(r[-1] for r in ROWS)
    need(Kmax==65 and all(r[-2]!=(1 << r[-1])-1 for r in ROWS),
         'bank excludes negative-one class')
    for fuel in range(65):
        H=fuel+Kmax
        n=(1 << H)-1
        need(not certify(n,fuel),('unbounded required fuel',fuel))
        for j in range(fuel+1):
            current=3**j*(1 << (H-j))-1
            need(rule(current)[0] is None,('Mersenne unmatched prefix',fuel,j))

    # A countermodel to the general entry-implies-termination inference.
    n=2
    for k in range(100):
        need(n==2+2*k,'countermodel even source')
        n=n+3
        need(n%2==1,'countermodel enters local descent region')
        n=n-1
    print('entry reduction: odd sources3..32767:',dict(sorted(entry_counts.items())))
    print('recursive policy universe:49999 odd sources3..99999, cap1024 decisions; all finite walks completed')
    print('certified by exception fuel:',counts)
    print('maximum required fuel in finite universe:',max_debt,'; first source attaining it:',first_max)
    print('q1R1 source27 required exception fuel:',debt[27])
    print('bounded interpreter cross-check:2047 sources times8 fuel levels')
    print('Mersenne nonuniformity:65 controls H=65+d; each C_d rejects its source')
    print('abstract local-entry countermodel:100 divergent grow/drop pairs')
    print('PASS: finite counts are not a universal basin-coverage claim')


if __name__=='__main__':
    main()
