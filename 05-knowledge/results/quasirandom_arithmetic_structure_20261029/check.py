"""Exact finite guards for the proofs; asymptotics use cited analytic inputs."""
from collections import Counter
from itertools import combinations, product
from math import gcd, prod
from pathlib import Path
import json

HERE = Path(__file__).resolve().parent
def need(ok, msg):
    if not ok: raise RuntimeError(msg)

def factors(n):
    out = {}; p = 2
    while p*p <= n:
        while n % p == 0:
            out[p] = out.get(p, 0)+1; n //= p
        p += 1
    if n > 1: out[n] = out.get(n, 0)+1
    return out

def jacobi(a, n):
    need(n > 0 and n % 2 == 1, 'odd Jacobi modulus')
    a %= n; sign = 1
    while a:
        while a % 2 == 0:
            a //= 2
            if n % 8 in (3, 5): sign = -sign
        a, n = n, a
        if a % 4 == n % 4 == 3: sign = -sign
        a %= n
    return sign if n == 1 else 0

def primes(limit):
    s = bytearray(b'\1')*(limit+1); s[:2] = b'\0\0'
    for p in range(2, int(limit**.5)+1):
        if s[p]: s[p*p::p] = b'\0'*(((limit-p*p)//p)+1)
    return [i for i in range(2, limit+1) if s[i]]

def order(a, p, k=1):
    m = p**k; d = (p-1)*p**(k-1)
    for r in factors(d):
        while d % r == 0 and pow(a, d//r, m) == 1: d //= r
    return d

def sigma(n):
    return prod(sum(p**j for j in range(e+1)) for p,e in factors(n).items())

def parity_columns(values):
    fs = [factors(n) for n in values]
    ps = sorted(set().union(*(f.keys() for f in fs)))
    cols = [sum((f.get(p,0)%2)<<i for i,p in enumerate(ps)) for f in fs]
    return ps, cols

def pattern(v, cols):
    return tuple((v & c).bit_count()%2 for c in cols)

def code_check(bank, use_sigma):
    vals = [sigma(b) if use_sigma else b for b in bank]
    ps, cols = parity_columns(vals)
    patterns = Counter(pattern(t, cols) for t in range(1<<len(ps)))
    modulus = 8*prod(p for p in ps if p != 2)
    direct = Counter()
    for a in range(3, modulus, 4):
        if gcd(a, modulus) == 1:
            direct[tuple((1-jacobi(x,a))//2 for x in vals)] += 1
    need(set(direct) == set(patterns), 'CRT code image')
    need(len(set(direct.values())) == 1, 'uniform CRT fibres')
    kernel = []
    for mask in range(1<<len(cols)):
        z = 0
        for i,c in enumerate(cols):
            if mask>>i & 1: z ^= c
        if z == 0: kernel.append(mask)
    need(all(sum(v[i] for i in range(len(cols)) if m>>i&1)%2 == 0
             for m in kernel for v in patterns), 'dual orthogonality')
    return dict(bank=bank, values=vals, coordinates=ps, columns=cols,
                crt_modulus=modulus, crt_each=next(iter(direct.values())),
                patterns=[''.join(map(str,v)) for v in sorted(patterns)],
                relation_weight_enumerator=dict(sorted(Counter(m.bit_count() for m in kernel).items())))

def phi(x, m=None):
    y = x*x-x+1
    return y if m is None else y % m

def orbit(x, m):
    seen = {}; vals = []
    while x not in seen:
        seen[x] = len(vals); vals.append(x); x = phi(x,m)
    return vals[:seen[x]], vals[seen[x]:]

def main():
    pp = primes(100000)
    for p in pp[:100]:
        if p == 2: continue
        for a in range(p):
            z = pow(a,(p-1)//2,p); z = -1 if z == p-1 else z
            need(jacobi(a,p) == z, 'Euler/Jacobi comparison')
    syl = [2]
    for _ in range(11): syl.append(phi(syl[-1]))
    for i,j in combinations(range(1,len(syl)),2):
        need(syl[j] % syl[i] == 1, 'Sylvester future reset')
        need(jacobi(syl[j],syl[i]) == 1 and jacobi(syl[i],syl[j]) == -1,
             'transitive Sylvester tournament')
    need(syl[4] == 13*139 and syl[5] == 3263443, 'specified terms')
    for p in (13,139):
        need(order(43,p) == 6, 'cyclotomic factors')
    need(phi(phi(phi(phi(2))))-2 == 5*19**2, '19 return depth')
    lifts = set(a+19*t for a in (2,3,7,5) for t in range(19))
    hist = Counter()
    while lifts:
        tail, cyc = orbit(min(lifts),19**2)
        need(not tail and set(cyc) <= lifts, 'cycle lift')
        hist[len(cyc)] += 1; lifts.difference_update(cyc)
    need(hist == {4:1,12:6}, '19 cycle histogram')
    g = 5417
    for _ in range(4): g = phi(g,19**3)
    need(g == 5417, 'unique next periodic lift')
    clocks = []
    for k in range(1,7):
        tail, cyc = orbit(2,19**k)
        need(not tail, '19 seed lies on cycle')
        clocks.append(len(cyc))
    exact_return, derivative = 2, 1
    for _ in range(12):
        derivative *= 2*exact_return-1
        exact_return = phi(exact_return)
    need((exact_return-2) % 19**4 == 17*19**3 and derivative % 19**2 == 172,
         'all-size clock seed certificate')
    need(clocks == [4,4]+[12*19**(k-3) for k in range(3,7)], 'proved period formula')
    tower = []
    for p in (3,5,7,13,17,19,139,257,1093,3511,65537):
        d = order(2,p); s = 1
        while pow(2,d,p**(s+1)) == 1: s += 1
        orders = [order(2,p,k) for k in range(1,5)]
        need(orders == [d*p**max(0,k-s) for k in range(1,5)], 'split grow law')
        tower.append(dict(p=p,d=d,s=s,orders=orders))
    codes = [code_check([4,2,5,3,7,13],False),
             code_check([4,2,5,3,7,13],True),
             code_check([2,5,7,19,29,59,89],True)]
    need(codes[2]['relation_weight_enumerator'] == {0:1,3:7,4:7,7:1}, 'Hamming code')
    need(Counter(x.count('1') for x in codes[2]['patterns']) == {0:1,4:7}, 'simplex code')
    for i,j in combinations(range(7),2):
        need(set(Counter((x[i],x[j]) for x in codes[2]['patterns']).values()) == {2},
             'pairwise independent simplex coordinates')
    fs = [2**(2**j)+1 for j in range(8)]
    for i,j in combinations(range(8),2):
        need(fs[j]%fs[i] == 2 and gcd(fs[i],fs[j]) == 1, 'Fermat residue reset')
        if i >= 2: need(jacobi(fs[j],fs[i]) == jacobi(fs[i],fs[j]) == 1,
                        'symmetric Fermat character relation')
    pp3 = [p for p in pp if p%4 == 3 and p not in (3,7)]
    vals = codes[1]['values']
    empirical = Counter(''.join(str((1-jacobi(a,p))//2) for a in vals) for p in pp3)
    need(set(empirical) == set(codes[1]['patterns']), 'observed code fibres')
    gram = []
    for p in (3,7,19,43,139):
        rows = [[jacobi(y-x,p) for y in range(p)] for x in range(p)]
        for x in range(p):
            for y in range(p):
                need(sum(a*b for a,b in zip(rows[x],rows[y])) == (p-1 if x == y else -1),
                     'Paley Gram')
        gram.append(p)
    need(order(2,3511) == 1755 and pow(2,182,1093**2) == 1093**2-1,
         'doubling symmetry hostile')
    aa,bb,dd = 1,0,1
    for bit in '11101000100':
        if bit == '1': aa,bb,dd = 3*aa,3*bb+dd,2*dd
        else: dd *= 2
    need((dd-aa,bb,gcd(dd-aa,bb)) == (1805,475,95), 'shared modulus carry hostile')
    monster = {2:46,3:20,5:9,7:6,11:2,13:3,17:1,19:1,23:1,29:1,31:1,41:1,47:1,59:1,71:1}
    baby = {2:41,3:13,5:6,7:2,11:1,13:1,17:1,19:1,23:1,31:1,47:1}
    ratio = prod(p**(e-baby.get(p,0)) for p,e in monster.items())
    def baseval(n,b):
        e = 0
        while n % b == 0: n //= b; e += 1
        return e
    need([baseval(ratio,b) for b in [4,2,5,3,7,13]] == [2,5,3,7,4,2], 'ratio functional graph')
    need([baseval(ratio//2,b) for b in [4,2,5,3,7,13]] == [2,4,3,7,4,2], 'actual index graph')
    # Same final multiplier but distinct Collatz carries for two words.
    def affine(word):
        a,b,d = 1,0,1
        for v in word: a,b,d = 3*a,3*b+d,d*2**v
        return a,b,d
    need(affine([1,3]) == (9,5,16) and affine([3,1]) == (9,11,16), 'ordered carry hostile')
    out = dict(status='FINITE-EXACT guards; analytic claims proved separately',
               sylvester_terms=[str(x) for x in syl], sylvester_pairs=55,
               phi6_mod19_lift_cycles=dict(hist), seed2_periods_mod19powers=clocks,
               twelve_return_certificate=dict(valuation=3, first_digit=17, derivative_mod361=172),
               doubling_towers=tower, codes=codes,
               primes3mod4_below100000=len(pp3), sigma_bank_prime_frequencies=dict(sorted(empirical.items())),
               paley_gram_checked_orders=gram, monster_ratio=str(ratio),
               fermat_pair_checks=28, collatz_shared_modulus=dict(D=dd-aa,B=bb,gcd=gcd(dd-aa,bb)),
               same_clock_carries=[affine([1,3]),affine([3,1])])
    (HERE/'checks.json').write_text(json.dumps(out,indent=2)+'\n',encoding='utf8')
    print(json.dumps({k:v for k,v in out.items() if k not in ('sylvester_terms','doubling_towers','codes')},indent=2))

if __name__ == '__main__': main()
