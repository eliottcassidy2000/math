#!/usr/bin/env python3
"""Base-phi (Bergman) expansions of primes <= 500, read as binary under three conventions:
   I = integer part only; F = full string without the radix point; R = reverse of F.
Exact arithmetic in Z[phi] (pairs (a,b) = a + b phi)."""
from sympy import isprime, primerange
import random

def fib(k):
    a, b = 0, 1
    for _ in range(k): a, b = b, a + b
    return a

def phipow(k):
    if k >= 0:
        return (1, 0) if k == 0 else (fib(k - 1), fib(k))
    m = -k
    s = -1 if m % 2 else 1
    return (s * fib(m + 1), -s * fib(m))

def sign(x):
    a, b = x
    # a + b phi = ((2a+b) + b sqrt5)/2
    p, q = 2 * a + b, b
    if p >= 0 and q >= 0: return 0 if p == 0 and q == 0 else 1
    if p <= 0 and q <= 0: return -1
    # opposite signs: compare p^2 vs 5 q^2
    if p > 0: return 1 if p * p > 5 * q * q else -1
    return 1 if 5 * q * q > p * p else -1

def sub(x, y): return (x[0] - y[0], x[1] - y[1])

def expand(n):
    """Greedy base-phi expansion of positive integer n; returns dict exponent->1."""
    x = (n, 0); digits = {}
    k = 0
    while sign(sub(x, phipow(k + 1))) >= 0: k += 1
    while sign(x) != 0:
        while sign(sub(x, phipow(k))) < 0: k -= 1
        digits[k] = 1; x = sub(x, phipow(k)); k -= 1
        assert k > -200
    return digits

def strings(n):
    d = expand(n)
    hi = max(d); lo = min(d)
    ip = ''.join('1' if d.get(i) else '0' for i in range(hi, -1, -1))
    fp = ''.join('1' if d.get(i) else '0' for i in range(-1, lo - 1, -1)) if lo < 0 else ''
    return ip, fp

def check_value(n):
    ip, fp = strings(n)
    # re-evaluate exactly
    tot = (0, 0)
    for i, c in enumerate(reversed(ip)):
        if c == '1': p = phipow(i); tot = (tot[0] + p[0], tot[1] + p[1])
    for j, c in enumerate(fp):
        if c == '1': p = phipow(-(j + 1)); tot = (tot[0] + p[0], tot[1] + p[1])
    assert tot == (n, 0), (n, ip, fp, tot)
    assert '11' not in ip + fp
    return ip, fp

def main():
    targets = [105, 223, 233, 332, 425]
    print('binary of targets:', {t: bin(t)[2:] for t in targets}, '-> every one contains "11":',
          all('11' in bin(t)[2:] for t in targets))
    print('base-phi of targets:')
    for t in targets + [47, 322, 121, 242, 377, 27]:
        ip, fp = check_value(t)
        print(f'   {t:4d} = {ip}.{fp}_phi   I={int(ip,2)}  F={int(ip+fp,2)}  R={int((ip+fp)[::-1],2)}')
    rows = []
    for p in primerange(2, 501):
        ip, fp = check_value(p)
        I = int(ip, 2); F = int(ip + fp, 2); R = int((ip + fp)[::-1], 2)
        rows.append((p, ip, fp, I, F, R))
    print('\nfirst primes:')
    for p, ip, fp, I, F, R in rows[:16]:
        print(f'   {p:3d} = {ip}.{fp}  I={I} F={F} R={R}')
    allvals = set()
    for r in rows: allvals |= {r[3], r[4], r[5]}
    print('targets among readings:', [t for t in targets if t in allvals])
    # length law: number of fractional digits vs integer digits
    lens = [(len(ip), len(fp)) for _, ip, fp, *_ in rows]
    print('integer-part length minus fractional length, over primes <= 500:', sorted(set(a - b for a, b in lens)))
    # H: fractional part = reverse of integer part with some fixed rule?
    # Test palindrome of the full word with radix centered
    pal = sum(1 for _, ip, fp, *_ in rows if (ip + fp) == (ip + fp)[::-1])
    print('primes whose full digit word is a palindrome:', pal, 'of', len(rows),
          [p for p, ip, fp, *_ in rows if (ip + fp) == (ip + fp)[::-1]])
    # Primality of readings vs composites (base rate)
    comps = [n for n in range(4, 501) if not isprime(n)]
    def rate(nums, idx):
        c = 0
        for n in nums:
            ip, fp = check_value(n)
            v = [int(ip, 2), int(ip + fp, 2), int((ip + fp)[::-1], 2)][idx]
            c += isprime(v)
        return c / len(nums)
    for idx, name in enumerate('IFR'):
        print(f'   P(reading {name} is prime): primes {rate([r[0] for r in rows], idx):.3f}  composites {rate(comps, idx):.3f}')
    # Lucas L_{2k} -> 2^{4k}+1 law
    print('L_{2k} full readings:', [(L, int(''.join(strings(L)), 2)) for L in (3, 7, 18, 47, 123, 322, 843, 2207)])
    # Collatz: T-parity word of p versus its base-phi word (both golden words)
    def tword(n, L):
        w = []
        for _ in range(L):
            w.append(n % 2); n = n // 2 if n % 2 == 0 else 3 * n + 1
        return ''.join(map(str, w))
    agree = []; agree_ctrl = []
    random.seed(2)
    for p, ip, fp, *_ in rows:
        word = ip + fp
        L = len(word)
        tw = tword(p, L)
        agree.append(sum(a == b for a, b in zip(word, tw)) / L)
        q = random.randrange(3, 500)
        tw2 = tword(q, L)
        agree_ctrl.append(sum(a == b for a, b in zip(word, tw2)) / L)
    print(f'agreement of base-phi word of p with T-parity word of p: {sum(agree)/len(agree):.3f}; '
          f'with T-parity word of random q: {sum(agree_ctrl)/len(agree_ctrl):.3f}')
    # stopping time of F-readings
    def stop(n):
        m = n; t = 0
        while m >= n and t < 10000:
            m = m // 2 if m % 2 == 0 else 3 * m + 1; t += 1
            if m < n: break
        return t
    st = [stop(r[4]) for r in rows if r[4] > 2]
    print('stopping times of F-readings (all F end in ...01, so F = 1 mod 4 -> stopping time 3):', sorted(set(st)))
    # reading (ii): binary digits of p evaluated at phi
    print('\n(ii) binary of p evaluated at phi, V(p) = A + B phi; norm N = A^2 + AB - B^2')
    hit = []
    for p in primerange(2, 501):
        bits = bin(p)[2:]
        A, B = 0, 0
        for i, c in enumerate(reversed(bits)):
            if c == '1':
                q = phipow(i); A += q[0]; B += q[1]
        N = A * A + A * B - B * B
        for t in targets:
            if t in (A, B, abs(N), A + B):
                hit.append((p, t, (A, B, N)))
    print('   target hits among A, B, |N|, A+B for primes <= 500:', hit)

if __name__ == '__main__':
    main()
