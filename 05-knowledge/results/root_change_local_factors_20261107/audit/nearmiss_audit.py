"""nearmiss_audit.py - auditor's check of the near-miss claims (PROOF section 2 'Near-misses are rough', HYP-9506 (D)).
Reads the lists from the S14 packet files (rc_1e9.txt: misses above 10^5; explore/rc_low_2p28.txt: M_2 <= 12 above
2^22), factors them with a self-contained trial division, and compares the part below 2^24 with the auditor's own census
output (census_2p24.txt, produced by census_lf_audit.c).
usage: python3 nearmiss_audit.py <S14 packet dir> <census_2p24.txt>
"""
import re
import sys

base, census = sys.argv[1], sys.argv[2]


def factor(n):
    f = []; p = 2
    while p * p <= n:
        while n % p == 0:
            f.append(p); n //= p
        p += 1 if p == 2 else 2
    if n > 1:
        f.append(n)
    return f


txt = open(base + '/rc_1e9.txt').read()
misses = [int(x) for x in re.search(r'misses \(<=2\) above 10\^5: ([\d ]+)', txt).group(1).split()]
low = open(base + '/explore/rc_low_2p28.txt').read()
lowl = [(int(a), int(b)) for a, b in re.findall(r'(\d+):(\d+)', low.split('multiplicity <= 12:')[1])]
le6 = [(v, m) for v, m in lowl if m <= 6]


def describe(name, vals):
    print('%s: %d values' % (name, len(vals)))
    res = {0: 0, 1: 0, 2: 0}; primes = 0; small = {}
    for v in vals:
        f = factor(v); res[v % 3] += 1; primes += len(f) == 1
        for p in set(f):
            if p < 20:
                small[p] = small.get(p, 0) + 1
        print('   %11d  mod3=%d  %s%s' % (v, v % 3, ' * '.join(map(str, f)), '  (prime)' if len(f) == 1 else ''))
    print('   residues mod 3 (0,1,2): %d, %d, %d; primes: %d; prime factors below 20 (factor: #values): %s'
          % (res[0], res[1], res[2], primes, small))


describe('misses of F_<=2 above 10^5 (rc_1e9.txt)', misses)
describe('M_2 <= 6 in (2^22, 2^28) (rc_low_2p28.txt)', [v for v, m in le6])
print('277438727 factors: %s' % factor(277438727))
print('all 130 values with M_2 <= 12 in (2^22,2^28): residues mod 3 = %s' % {r: sum(1 for v, m in lowl if v % 3 == r) for r in range(3)})
# cross-check with the auditor's own census below 2^24
c = open(census).read()
mine_m = [int(x) for x in re.search(r'misses \(odd N <= X with M_2 = 0\):([\d ]+)', c).group(1).split()]
mine_low = [(int(a), int(b)) for a, b in re.findall(r'(\d+):(\d+)', re.search(r'low values \(2\^22 < N <= X, M_2 <= 12\):(.*)', c).group(1))]
print('own census: %d misses, largest %d; misses above 10^5 equal to the S14 list: %s'
      % (len(mine_m), mine_m[-1], [v for v in mine_m if v > 10 ** 5] == misses))
print('own census: values with M_2 <= 12 in (2^22, 2^24]: %d; identical (value and count) to the S14 list restricted to <= 2^24: %s'
      % (len(mine_low), mine_low == [(v, m) for v, m in lowl if v <= 2 ** 24]))

# small-prime content of the 130 values with M_2 <= 12 in (2^22, 2^28) and of the 15 misses, against random odd
# integers coprime to 3 (expected share divisible by p: 1/p; with a prime factor in 5..19: 1 - prod(1 - 1/p))
import math
small = [5, 7, 11, 13, 17, 19]
exp_any = 1 - math.prod(1 - 1 / p for p in small)
for name, vals in (('130 values M_2 <= 12', [v for v, m in lowl]), ('15 misses > 10^5', misses)):
    cnt = {p: sum(1 for v in vals if v % p == 0) for p in small}
    anyp = sum(1 for v in vals if any(v % p == 0 for p in small))
    print('%s: divisible by p (observed/expected): %s; with a prime factor in 5..19: %d/%d = %.3f (random odd, 3-free: %.3f)'
          % (name, ', '.join('%d: %d/%.1f' % (p, cnt[p], len(vals) / p) for p in small), anyp, len(vals), anyp / len(vals), exp_any))
