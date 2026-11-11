"""A135244 by prime-factor strata; exact sieve of k^2+1; saved primary controls."""
from pathlib import Path
from array import array
from math import isqrt
from collections import Counter
import csv, hashlib, json, time, urllib.request

HERE = Path(__file__).resolve().parent
N = 10000
K = 100000

def require(ok, msg):
    if not ok:
        raise RuntimeError(msg)

def primes_upto(n):
    flags = bytearray(b'\x01')*(n+1)
    flags[:2] = b'\x00\x00'
    for p in range(2,isqrt(n)+1):
        if flags[p]:
            flags[p*p::p] = b'\x00'*((n-p*p)//p+1)
    return flags, [p for p in range(2,n+1) if flags[p]]

def primary(seq):
    p = HERE / ('b'+seq[1:]+'.txt')
    url = 'https://oeis.org/'+seq+'/b'+seq[1:]+'.txt'
    if not p.exists():
        req = urllib.request.Request(url,headers={'User-Agent':'Mathematical research reproducibility; Python urllib'})
        p.write_bytes(urllib.request.urlopen(req,timeout=30).read())
    pairs = []
    for line in p.read_text(encoding='utf-8-sig').splitlines():
        line = line.strip()
        if line and not line.startswith('#'):
            i,v = line.split()[:2]
            pairs.append((int(i),int(v)))
    canonical=p.read_bytes().replace(b'\r\n',b'\n')
    return dict(pairs), dict(url=url,sha256_canonical_lf=hashlib.sha256(canonical).hexdigest(),rows=len(pairs))

def trial_sum(n):
    if n == 1:
        return 0
    total = 1
    for d in range(2,isqrt(n)+1):
        if n%d == 0:
            total += d
            if d*d != n:
                total += n//d
    return total

def mr(n):
    # Deterministic on unsigned 64-bit integers; trial division is a separate small control.
    if n < 2:
        return False
    for p in [2,3,5,7,11,13,17,19,23,29,31,37]:
        if n%p == 0:
            return n == p
    d = n-1
    s = 0
    while d%2 == 0:
        d//=2
        s+=1
    for a in [2,325,9375,28178,450775,9780504,1795265022]:
        if a%n == 0:
            continue
        x = pow(a,d,n)
        if x in (1,n-1):
            continue
        for _ in range(s-1):
            x = x*x%n
            if x == n-1:
                break
        else:
            return False
    return True

def run():
    start = time.perf_counter()
    # Omega(m)>=3 forces m<=n^(3/2); the sieve need not scan the quadratic tail.
    bound = isqrt(N**3)
    sums = array('I',[0])*(bound+1)
    for d in range(1,bound//2+1):
        for m in range(2*d,bound+1,d):
            sums[m] += d
    omega = bytearray(bound+1)
    for p in range(2,bound+1):
        if omega[p] == 0:
            power = p
            while power <= bound:
                for m in range(power,bound+1,power):
                    omega[m] += 1
                power *= p
    flags, primes = primes_upto(max(N,K))
    high = [0]*(N+1)
    for m in range(2,bound+1):
        n = sums[m]
        if omega[m]>=3 and 2<=n<=N:
            high[n] = m
    rows = []
    max_inverse = [0]*(N+1)
    for n in range(2,N+1):
        square = (n-1)**2 if flags[n-1] else 0
        semiprime = 0
        if n%2 == 0:
            q = n-3
            if q>2 and flags[q]:
                semiprime = 2*q
        else:
            # Distinct primes: equality p=q has divisor sum p+1, not 2p+1.
            p = (n-2)//2
            while p>=2:
                q = n-1-p
                if flags[p] and flags[q]:
                    semiprime = p*q
                    break
                p-=1
        a = max(square,semiprime,high[n])
        max_inverse[n] = a
        family = 'none' if not a else ('prime_square' if a==square else
                 'distinct_semiprime' if a==semiprime else 'Omega_ge_3')
        rows.append([n,a,family,square,semiprime,high[n],
                     int(semiprime>isqrt(n**3))])
    with (HERE/'A135244_10000.tsv').open('w',newline='',encoding='utf-8') as f:
        w = csv.writer(f,delimiter='\t')
        w.writerow(['n','largest_predecessor','winning_family','prime_square_max','distinct_semiprime_max',
                    'Omega_ge_3_max','semiprime_exceeds_high_Omega_bound'])
        w.writerows(rows)
    reference, src1 = primary('A135244')
    checked = 0
    for n,a in reference.items():
        if 2<=n<=N:
            require(max_inverse[n]==a, 'OEIS A135244 at '+str(n))
            checked += 1
    brute = [0]*151
    for m in range(2,149**2+1):
        n = trial_sum(m)
        if 2<=n<=150:
            brute[n] = m
    require(brute[2:] == max_inverse[2:151], 'independent bounded inverse fibres')
    # A complete congruence sieve: every composite k^2+1 has a prime divisor <=k.
    candidate = bytearray(b'\x01')*(K+1)
    candidate[0] = 0
    candidate[3::2] = b'\x00'*len(candidate[3::2])
    roots = []
    for p in primes:
        if p%4 != 1:
            continue
        g = 2
        while pow(g,(p-1)//2,p) != p-1:
            g += 1
        r = pow(g,(p-1)//4,p)
        require(r*r%p == p-1, 'root of minus one')
        roots.append((p,r,p-r))
        for a in [r,p-r]:
            for k in range(a,K+1,p):
                if k*k+1 != p:
                    candidate[k] = 0
    values = [(k,k*k+1) for k in range(1,K+1) if candidate[k]]
    # Independent exact primality path on the entire input interval.
    for k in range(1,K+1):
        require(bool(candidate[k]) == mr(k*k+1), 'Miller-Rabin quadratic agreement')
        if k<=1000:
            v = k*k+1
            require(bool(candidate[k]) == all(v%d for d in range(2,isqrt(v)+1)), 'trial quadratic agreement')
    with (HERE/'A002496_inputs_100000.tsv').open('w',newline='',encoding='utf-8') as f:
        w = csv.writer(f,delimiter='\t')
        w.writerow(['index','k','prime'])
        for i,(k,p) in enumerate(values,1):
            w.writerow([i,k,p])
    ref2, src2 = primary('A002496')
    checked2 = 0
    for i,v in ref2.items():
        if 1<=i<=len(values):
            require(values[i-1][1]==v, 'OEIS A002496 at '+str(i))
            checked2 += 1
    # Exact bridge: A135244(k^2+2)=(k^2+1)^2 iff k^2+1 is prime.
    bridge = 0
    for k in range(1,isqrt(N-2)+1):
        p = k*k+1
        require((max_inverse[p+1]==p*p)==bool(flags[p]), 'quadratic equality bridge')
        bridge += 1
    controlA, src0 = primary('A000568')
    birth = list(csv.DictReader((HERE/'birth_sequences.tsv').open(encoding='utf-8'),delimiter='\t'))
    checkedA = 0
    for row in birth:
        n = int(row['source_n'])
        if n in controlA:
            require(int(row['A000568'])==controlA[n], 'OEIS A000568 at '+str(n))
            checkedA += 1
    result = dict(status='PASS',aliquot_target_bound=N,aliquot_sieved_source_bound=bound,
        original_quadratic_bound=(N-1)**2,aliquot_OEIS_rows=checked,independent_inverse_targets=149,
        winning_families=dict(Counter(r[2] for r in rows)),quadratic_input_bound=K,
        quadratic_prime_count=len(values),quadratic_OEIS_rows=checked2,A000568_OEIS_rows=checkedA,
        quadratic_inverse_bridge_checks=bridge,
        named_numbers={str(n):max_inverse[n] for n in [2,3,5,7,11,13,19,27,33,1039]},
        primary_sources=[src0,src1,src2])
    (HERE/'arithmetic_checks.json').write_text(json.dumps(result,indent=2)+'\n',encoding='utf-8')
    print(json.dumps(result))

if __name__ == '__main__':
    run()
