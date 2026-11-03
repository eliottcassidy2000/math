from math import gcd, isqrt

def check(ok, message):
    if not ok:
        raise RuntimeError(message)

def is_prime(n):
    return n >= 2 and all(n % d for d in range(2, isqrt(n) + 1))

def squarefree(n):
    return all(n % (p * p) for p in range(2, isqrt(n) + 1) if is_prime(p))

def divisor_counts(n):
    proper = [d for d in range(2, n) if n % d == 0]
    return len(proper), sum(squarefree(d) for d in proper), sum(is_prime(d) for d in proper)

for n, want in [(125, (2, 1, 1)), (5525, (10, 7, 3)), (36725, (10, 7, 3)), (145, (2, 2, 2))]:
    check(divisor_counts(n) == want, ('divisor profile', n))

# Independent Euclid-parameter enumeration, with odd leg first.
def ppt_fiber(c):
    out = set()
    for n in range(2, isqrt(c) + 1):
        m = isqrt(c - n * n)
        if 0 < m < n and n*n + m*m == c and gcd(n, m) == 1 and (n-m) % 2:
            out.add((n*n - m*m, 2*n*m, c))
    return out

check(ppt_fiber(125) == {(117,44,125)}, '125 fiber')
check(ppt_fiber(5525) == {(5133,2044,5525),(5427,1036,5525),(2163,5084,5525),(4557,3124,5525)}, '5525 fiber')
check((813,36716,36725) in ppt_fiber(36725), '36725 fiber')

def plus(s, t):
    check(s > 0 and t > 0 and s%2 == 1 and t%2 == 1 and gcd(s,t) == 1, 'input guard')
    x, k = 3*s+t, 0
    while x % 2 == 0:
        x //= 2
        k += 1
    g = gcd(x,t)
    return x//g, t//g, k, g

check(plus(271,3) == (17,1,4,3), '36725 return')
check((17, (17*17-1)//2, (17*17+1)//2) == (17,144,145), '145 target')
for s in range(1, 150, 2):
    for t in range(1, 150, 2):
        if gcd(s,t) != 1:
            continue
        sn,tn,k,g = plus(s,t)
        check(g == (3 if t%3 == 0 else 1), ('split',s,t))
        check(gcd(sn,tn) == 1 and sn%2 and tn%2, ('reduced',s,t))
        if t%3 == 0:
            check(3*(sn*sn+tn*tn) < s*s+t*t, ('contraction',s,t))
check(plus(3,271)[1] == 271, 'reversed signed B')
print('PASS independent codec, divisor counts, guard, contraction, sign hostile')
