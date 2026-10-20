#!/usr/bin/env python3
"""Exact, standalone controls. No repository mathematical-code imports."""
from itertools import product, combinations
from math import gcd, comb
from hashlib import sha256
import json

checks = 0

def require(condition, message):
    global checks
    checks += 1
    if not condition:
        raise RuntimeError(message)

def carry(w):
    b, scale = 0, 1
    for bit in w:
        if bit:
            b = 3*b + scale
        scale *= 2
    return b

def positional(w):
    p = sum(w)
    return sum(2**i * 3**(p-j-1) for j, i in enumerate(i for i,b in enumerate(w) if b))

def rotate(w, r):
    return w[r:] + w[:r]

def core_difference(a, b):
    n = len(a)
    left = 0
    while left < n and a[left] == b[left]:
        left += 1
    right = n
    while right > left and a[right-1] == b[right-1]:
        right -= 1
    return carry(a[left:right]) - carry(b[left:right])

neighbor_cases = 0
for q in range(2, 13):
    for w in product((0,1), repeat=q):
        p = sum(w)
        if p in (0,q):
            continue
        D = 2**q - 3**p
        B = carry(w)
        require(B == positional(w), 'independent carry sum')
        g = gcd(abs(D), B)
        for i in range(q):
            a, b = rotate(w,i), rotate(w,(i+1)%q)
            delta = core_difference(a,b)
            require(delta != 0, 'neighbor core nonzero')
            require(gcd(abs(D),abs(delta)) == g, 'single chronological edge')
            neighbor_cases += 1

repetition_cases = 0
for q in range(2,9):
    for w in product((0,1), repeat=q):
        p = sum(w)
        if p in (0,q):
            continue
        D, B = 2**q-3**p, carry(w)
        g = gcd(abs(D), B)
        for k in range(1,6):
            Dk, Bk = 2**(k*q)-3**(k*p), carry(w*k)
            H = Dk//D
            require(Dk == H*D and Bk == H*B, 'power cocycle')
            require(gcd(abs(Dk),Bk) == H*g, 'power content')
            require(abs(Dk)//gcd(abs(Dk),Bk) == abs(D)//g, 'denominator invariant')
            repetition_cases += 1

gap_cases = 0
for k in range(2,201):
    D = 2**(2*k+1)-9
    for s in range(k):
        d = k-s
        B = 3 + 2**d
        M = 2**(2*s+1)-1
        g = gcd(D,B)
        require(g == gcd(M,B), 'Mersenne content formula')
        require(g < D, 'positive two-one exclusion')
        next_k = k+2*s+1
        require(gcd(M,3+2**(next_k-s)) == g, 'content periodicity')
        if s == 1:
            require(g == (7 if k%3 == 0 else 1), 'one-swap exact classification')
        gap_cases += 1

hostiles = []
for s,k,expected in [(10,27,49),(36,99,439),(58,102,6559)]:
    D = 2**(2*k+1)-9
    B = 3+2**(k-s)
    require(gcd(D,B) == expected, 'higher displacement content')
    hostiles.append({'s':s,'k':k,'q':2*k+1,'short_gap':k-s,'g':expected})

prime_power_word = tuple(map(int,'11101000100'))
require(carry(prime_power_word) == 475, 'prime support hostile carry')
require(2**11-3**5 == 1805 and gcd(1805,475) == 95, 'prime support hostile content')

spectra = {}
for q,p in [(7,2),(11,5),(15,8)]:
    hist = {}
    necklace_seen = set()
    necklace_hist = {}
    for positions in combinations(range(q),p):
        w = tuple(int(i in positions) for i in range(q))
        g = gcd(abs(2**q-3**p),carry(w))
        hist[g] = hist.get(g,0)+1
        canonical = min(rotate(w,r) for r in range(q))
        if canonical not in necklace_seen:
            necklace_seen.add(canonical)
            necklace_hist[g] = necklace_hist.get(g,0)+1
    require(sum(hist.values()) == comb(q,p), 'fixed shape census size')
    require(all(count%q == 0 for count in hist.values()), 'free rotations at coprime shape')
    require(sum(necklace_hist.values())*q == comb(q,p), 'Burnside free count')
    spectra[f'{q},{p}'] = {'word_content_histogram':hist,'necklace_content_histogram':necklace_hist}
require(set(map(int,spectra['15,8']['word_content_histogram'])) == {1,73,359}, 'mixed-prime incompatible witnesses')
require(set(map(int,spectra['11,5']['word_content_histogram'])) == {1,5,19,95}, 'full prime powers required')

# Enumerate directed cyclic swaps of zero-intercept coprime mechanical words.
# Repeated output necklaces are retained. This is not a distinct-necklace census.
digest = sha256()
swap_cases = nonunit = positive_cases = negative_cases = 0
maximum = {'g':0}
integer_candidates = []
for q in range(3,401):
    powers2 = [2**i for i in range(q)]
    for p in range(1,q):
        if gcd(q,p) != 1:
            continue
        c = tuple((i+1)*p//q-i*p//q for i in range(q))
        D = 2**q-3**p
        B = carry(c)
        powers3 = [3**i for i in range(p+1)]
        suffix = [0]*(q+1)
        for i in range(q-1,-1,-1):
            suffix[i] = suffix[i+1] + c[i]
        for i in range(q):
            j = (i+1)%q
            if c[i] == c[j]:
                continue
            if j:
                Bw = B + (1 if c[i] else -1)*powers2[i]*powers3[suffix[i+2]]
            else:
                w = list(c)
                w[i],w[j] = w[j],w[i]
                Bw = carry(w)
            if q <= 30:
                w = list(c)
                w[i],w[j] = w[j],w[i]
                require(Bw == positional(w), 'independent one-swap carry')
            g = gcd(abs(D),Bw)
            swap_cases += 1
            nonunit += g > 1
            positive_cases += D > 0
            negative_cases += D < 0
            if g > maximum['g']:
                maximum = {'q':q,'p':p,'swap_start':i,'g':g,'D':D,'B':Bw}
            if abs(D)>1 and g == abs(D):
                integer_candidates.append({'q':q,'p':p,'swap_start':i,'D':D,'B':Bw})
            digest.update(f'{q}:{p}:{i}:{g}\n'.encode())

print(json.dumps({
    'status':'PASS',
    'exact_checks':checks,
    'neighbor_edge_universe':'every binary word containing both bits, 2<=q<=12, every chronological neighboring rotation',
    'neighbor_edge_cases':neighbor_cases,
    'word_power_universe':'every binary word containing both bits, 2<=q<=8, repetitions 1<=k<=5',
    'word_power_cases':repetition_cases,
    'two_one_universe':'2<=k<=200, 0<=s<k, q=2k+1 and short cyclic one-gap k-s',
    'two_one_cases':gap_cases,
    'higher_displacement_hostiles':hostiles,
    'fixed_shape_spectra':spectra,
    'swap_universe':'all directed unequal cyclic adjacent swaps of zero-intercept lower mechanical words, gcd(q,p)=1, 3<=q<=400; duplicate outputs retained',
    'swap_occurrences':swap_cases,
    'positive_clock_occurrences':positive_cases,
    'negative_clock_occurrences':negative_cases,
    'nonunit_occurrences':nonunit,
    'maximum_content':maximum,
    'integer_candidates_abs_D_gt_1':integer_candidates,
    'swap_digest_sha256':digest.hexdigest(),
    'boundary':'The q<=400 absence is finite evidence only. The all-q one-swap conjecture remains OPEN; the two-one formula and the other algebraic results have all-length proofs.'
}, indent=2, sort_keys=True))
