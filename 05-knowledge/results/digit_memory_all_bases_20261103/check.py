#!/usr/bin/env python3
"""Exact finite controls; all quantifiers and proofs are in PROOFS.md.

No repository imports, floats, random sampling or optimization-sensitive asserts.
"""
import itertools
import json
import math
from collections import Counter


def require(ok, payload):
    if not ok:
        raise RuntimeError(payload)


def valuation(n, p):
    if n == 0:
        return math.inf
    a = 0
    while n % p == 0:
        a += 1
        n //= p
    return a


def split_base(B, m):
    u = m
    while (g := math.gcd(u, B)) > 1:
        u //= g
    s = m // u
    L, power = 0, 1
    while power % s:
        L += 1
        power *= B
    return s, u, L


def signature(B, m, r):
    s, u, L = split_base(B, m)
    K, power = L, B**L
    while power < m:
        K += 1
        power *= B
    short = []
    for k in range(K):
        v = (-pow(B, k, m) * r) % m
        short.append(v if v < B**k else -1)
    return (r % u, *short)


def dfa_partition(B, m):
    # Independent refinement by all one-digit successors, until stable.
    labels = [int(r == 0) for r in range(m)]
    while True:
        bank = {}
        new = []
        for r in range(m):
            key = (r == 0, *(labels[(B*r+d) % m] for d in range(B)))
            new.append(bank.setdefault(key, len(bank)))
        if all((labels[i] == labels[j]) == (new[i] == new[j])
               for i in range(m) for j in range(i)):
            return new
        labels = new


def run_collatz(n, k):
    word, ordinary = 0, 0
    for j in range(k):
        b = n % 2
        word |= b << j
        ordinary += 1+b
        n = (3*n+1)//2 if b else n//2
    return word, n, ordinary


def carrier(word, k):
    A, C, D = 1, 0, 1
    for j in range(k):
        b = (word >> j) & 1
        A, C, D = (1+2*b)*A, (1+2*b)*C+b*D, 2*D
    r = (-C*pow(A, -1, D)) % D if D > 1 else 0
    return A, C, D, r


def main():
    counts = Counter()
    samples = []
    # Exact complete finite remainder universes, all bases and moduli declared.
    for B in range(2, 17):
        for m in range(1, 101):
            sig = [signature(B, m, r) for r in range(m)]
            labels = dfa_partition(B, m)
            for i in range(m):
                for j in range(i):
                    require((sig[i] == sig[j]) == (labels[i] == labels[j]),
                            ('DFA', B, m, i, j))
            counts['base_modulus_automata'] += 1
            if (B, m) in [(2,4), (2,8), (10,9), (10,16), (12,18), (16,15)]:
                samples.append({'base':B, 'modulus':m, 'minimal_states':len(set(labels))})
            s, u, L = split_base(B, m)
            for n in range(2*m+1):
                t = n % (B**L)
                require((n % s == 0) == (t % s == 0), ('suffix', B,m,n))
                if u > 1:
                    a,d = divmod(n,B)
                    require((n % u == 0) == ((a+pow(B,-1,u)*d) % u == 0),
                            ('trim', B,m,n))
                counts['divisibility_sources'] += 1

    # Translation periods preserving legality AND output modulo M.
    for a in range(-6, 7):
        for b in range(-4, 5):
            for d in range(1, 9):
                for M in range(1, 10):
                    def f(n):
                        return ((a*n+b)//d) % M if (a*n+b) % d == 0 else None
                    expected = (d*M // math.gcd(a,d*M)
                                if b % math.gcd(a,d) == 0 else 1)
                    values = [f(n) for n in range(d*M)]
                    actual = next(P for P in range(1,d*M+1)
                                  if all(values[n] == values[(n+P) % (d*M)]
                                         for n in range(d*M)))
                    require(actual == expected, ('quotient',a,b,d,M,actual,expected))
                    counts['affine_quotient_parameters'] += 1

    separators = []
    for B in range(2, 21):
        for p in (2,3,5,7):
            a = valuation(B,p)
            C = B//p**a
            for n in range(1, 65):
                for m in range(n+1, 65):
                    if C == 1:
                        if valuation(n,p) == valuation(m,p):
                            for digit in range(B):
                                require(valuation(B*n+digit,p) == valuation(B*m+digit,p),
                                        ('pure-base', B,p,n,m,digit))
                            counts['pure_base_equal_valuation_pairs'] += 1
                    else:
                        r = valuation(m-n,p)
                        k = 1
                        while C**k < p**(r+1):
                            k += 1
                        t = (-C**k*n) % p**(r+1)
                        suffix = p**(a*k)*t
                        x,y = B**k*n+suffix, B**k*m+suffix
                        require(0 <= suffix < B**k, ('suffix-domain',B,p,n,m))
                        require(valuation(x,p) >= a*k+r+1 and valuation(y,p) == a*k+r,
                                ('separator',B,p,n,m,k,suffix))
                        counts['mixed_base_distinct_pairs'] += 1
                        if (B,p,n,m) in [(10,2,1,3),(3,2,1,3),(12,3,1,4)]:
                            separators.append([B,p,n,m,k,suffix,x,y])

    memory = []
    for k in range(0, 10):
        for M in (1,2,3,4,5,6,7,8,9,11,18,25,27,54,81):
            response = Counter()
            # One complete period. n=0 represents the positive congruence class.
            for n in range((2**k)*M):
                w,y,clock = run_collatz(n,k)
                require(clock == k+w.bit_count(), ('clock',n,k))
                response[w,y % M] += 1
            predicted = sum(math.comb(k,o)*M//math.gcd(3**o,M) for o in range(k+1))
            require(len(response) == predicted, ('count',k,M,len(response),predicted))
            for (w,phase),multiplicity in response.items():
                A,C,D,r = carrier(w,k)
                g = math.gcd(A,M)
                d0 = (A*r+C)//D
                require((phase-d0) % g == 0 and multiplicity == g,
                        ('fibre',k,M,w,phase,multiplicity,g))
            counts['collatz_complete_periods'] += 1
            counts['collatz_sources'] += (2**k)*M
            if M in (3,9,27) and k in (0,1,2,3,6,9):
                memory.append([k,M,predicted])

    # Direct ordinary-clock controls, CRT precision, and all-odd hostile family.
    for k in range(1, 9):
        for w in range(2**k):
            A,C,D,r = carrier(w,k)
            n = r+D
            ww,y,clock = run_collatz(n,k)
            require(ww == w and D*y == A*n+C, ('native',k,w))
            z=n
            for _ in range(clock):
                z = 3*z+1 if z%2 else z//2
            require(z == y, ('ordinary-clock',k,w))
            for extra in range(3):
                for ternary in range(3):
                    delta = 2**(k+extra)*3**ternary
                    w2,y2,_ = run_collatz(n+delta,k)
                    require(w2 == w and y2-y == 2**extra*3**(ternary+w.bit_count()),
                            ('precision',k,w,extra,ternary))
                    counts['precision_cells'] += 1
            for M0 in (1,2,3,5,6,9,12,27):
                o = w.bit_count()
                M = M0*3**o if M0 % 3 == 0 else M0
                phases = set()
                for t in range(M0):
                    yy = (A*r+C)//D + A*t
                    phases.add(yy % M)
                require(len(phases) == M0, ('adaptive',k,w,M0,M))
                # Necessary and sufficient quotient period; minimal modulus.
                require(M//math.gcd(A,M) == M0, ('adaptive-period',k,w,M0,M))
                counts['adaptive_source_cells'] += 1
        for q in (1,2,7):
            n=q*2**k-1
            w,y,_=run_collatz(n,k)
            require(w == 2**k-1 and y == q*3**k-1, ('odd-run',k,q))
    print(json.dumps({'status':'FINITE-EXACT; proofs in PROOFS.md',
                      'checks':dict(counts), 'minimal_divisibility_states':samples,
                      'mixed_base_separator_examples':separators,
                      'collatz_response_counts_k_M_N':memory}, indent=2, sort_keys=True))


if __name__ == '__main__':
    main()
