#!/usr/bin/env python3
"""Independent exact controls. Does not import or write the producer packet."""
from fractions import Fraction as Q
from itertools import product
import json

checks = 0

def require(condition, message):
    global checks
    checks += 1
    if not condition:
        raise RuntimeError(message)

def odd_trace(n, word, stop=False):
    values = [n]
    for exponent in word:
        if stop and n == 1:
            return None
        value = 3*n + 1
        a = 0
        while value % 2 == 0:
            a += 1
            value //= 2
        if a != exponent:
            return None
        n = value
        values.append(n)
    return values

def literal_trace(n, length):
    out = [n]
    for _ in range(length):
        n = n//2 if n % 2 == 0 else (3*n+1)//2
        out.append(n)
    return out

def affine(word):
    # Ordinary slope/intercept multiplication, independent of integer carry
    # recurrence in the producer implementation.
    slope, intercept = Q(1), Q(0)
    for a in word:
        slope = Q(3, 2**a)*slope
        intercept = Q(3, 2**a)*intercept + Q(1, 2**a)
    m = 2**sum(word)
    s, c = slope*m, intercept*m
    require(s.denominator == c.denominator == 1, 'integer coefficients')
    s, c = int(s), int(c)
    r = (pow(s, -1, 2*m)*(m-c)) % (2*m)
    return s, m, c, r

parameters = list(product([Q(1,2), Q(8,9), Q(1), Q(4,3), Q(2), Q(3)],
                          [Q(-1), Q(0), Q(1,3), Q(5,9), Q(1)]))
affine_words = 0
for length in range(4):
    for word in product(parameters, repeat=length):
        a, b, multiplier, debt, slack = Q(1), Q(0), Q(1), Q(0), Q(0)
        for l, d in word:
            a, b = a/l, (b+d)/l
            slack += multiplier*(l-1-d)
            debt += multiplier*d
            multiplier *= l
        require(a == 1/multiplier and b == debt/multiplier, 'forward normalization')
        require(multiplier - 1 - debt == slack, 'weighted slack telescoping')
        if all(0 <= d <= l-1 for l, d in word):
            require(0 <= debt <= multiplier-1, 'nonnegative budget cone closure')
        if word and all(l > 1 for l, _ in word):
            prefix = Q(1)
            weights, barycentre = Q(0), Q(0)
            for l, d in word:
                weight = prefix*(l-1)/(multiplier-1)
                require(weight > 0, 'positive barycentric coefficient')
                weights += weight
                barycentre += weight*d/(l-1)
                prefix *= l
            require(weights == 1 and barycentre == debt/(multiplier-1), 'convex debt')
        affine_words += 1

# Nonnegative debt ALONE does not control the sign in equation (2.6).
l1, d1, l2, d2 = Q(1,2), Q(0), Q(2), Q(1)
defect_difference = (1+d1+l1*d2)/(l1*l2) - (1+d1)*(1+d2)/(l1*l2)
require(defect_difference == Q(-1,2), 'sign needs eta_1 nonnegative')
# Negative lambda reverses substitution, even for nonnegative resources.
require(Q(0) <= Q(1) and Q(0) <= -Q(0) and not Q(0) <= -Q(1),
        'negative-lambda hostile')

native_words = 0
for word in product(product([Q(1,2), Q(7,10), Q(9,10)], [Q(0), Q(4), Q(7), Q(9)]), repeat=3):
    original, normalized = Q(13), Q(13,7)
    multiplier, debt = Q(1), Q(0)
    for rho, residual in word:
        original = rho*original+(1-rho)*residual
        l, d = 1/rho, (1/rho-1)*residual/7
        normalized = (normalized+d)/l
        debt += multiplier*d
        multiplier *= l
    require(original/7 == normalized == (Q(13,7)+debt)/multiplier, 'native conversion')
    native_words += 1

valuation_words = residues = inverse_candidates = root_trimmed_words = 0
for length in range(1,5):
    for word in product(range(1,5), repeat=length):
        s, m, c, r = affine(word)
        good = []
        for n in range(1, 2*m, 2):
            if odd_trace(n, word) is not None:
                good.append(n)
            residues += 1
        require(good == [r], 'independent exhaustive residue cylinder')
        endpoint = (s*r+c)//m
        for t in range(4):
            n = r+2*m*t
            trace = odd_trace(n, word)
            stopped = odd_trace(n, word, stop=True)
            require(trace[-1] == endpoint+2*s*t, 'exact image progression')
            require(literal_trace(n, sum(word))[-1] == trace[-1], 'literal versus odd clocks')
            require((stopped is not None) == all(v > 1 for v in trace[:-1]), 'ROOT proper prefix')
            if t:
                require(stopped is not None, 'large lifts never hit ROOT in proper prefixes')
        root_trimmed = odd_trace(r, word, stop=True) is None
        root_trimmed_words += root_trimmed
        nmin = r+2*m*int(root_trimmed)
        criterion = m > s and (m-s)*nmin > c
        for t in range(4):
            n = nmin+2*m*t
            descend = odd_trace(n, word, stop=True)[-1] < n
            require(not criterion or descend, 'cylinder descent sufficient')
            if t == 0:
                require(criterion == descend, 'least-source descent necessary')
        for y in range(1, endpoint+4*s+1, 2):
            preimage = Q(m*y-c, s)
            lawfully_inverts = (preimage.denominator == 1 and preimage > 0
                               and preimage.numerator % 2 == 1
                               and odd_trace(int(preimage), word) is not None)
            guarded = y >= endpoint and (y-endpoint) % (2*s) == 0
            require(lawfully_inverts == guarded, 'inverse congruence and lower bound iff')
            inverse_candidates += 1
        valuation_words += 1

bank = {'A': (1,2,3), 'B': (1,3), 'E': (1,1,3)}
macro_data = {}
for name, word in bank.items():
    s,m,c,r = affine(word)
    # Each literal value is affine in t; derive its coefficients at t=0,1.
    base, lift = literal_trace(r,sum(word)), literal_trace(r+2*m,sum(word))
    differences = [(base[i]-r, lift[i]-base[i]-2*m) for i in range(1,len(base))]
    first = next(i for i,(a,b) in enumerate(differences,1) if a < 0 and b <= 0)
    require(all(a >= 0 and b >= 0 for a,b in differences[:first-1]), 'uniform first-exit prefix')
    require(c <= 4*(m-s) and Q(s,m) <= Q(27,32) and r > 4, 'common barrier and guard')
    macro_data[name] = {'word': list(word), 'coefficients': [s,m,c],
                        'source_residue': r, 'source_modulus': 2*m,
                        'image_residue': base[-1], 'image_modulus': 2*s,
                        'fixed_point': str(Q(c,m-s)), 'first_exit': first,
                        'literal_values_at_first_lift': base}
require(odd_trace(11,bank['A']) == [11,17,13,5], 'growing reset repaid')
require(odd_trace(55,bank['E']) == [55,83,125,47], 'unit budget not necessary')
require(all(odd_trace(5,w) is None for w in bank.values()), 'source 11 leaves bank at 5')
require(Q(64*7-23,27).denominator != 1, 'inverse ray alone not a legal inverse branch')
require(odd_trace(5,(4,2)) == [5,1,1] and odd_trace(5,(4,2),True) is None,
        'ROOT hostile rejects continued root')
require(odd_trace(133,(4,2),True) == [133,25,19], 'ROOT head removal survivor')

first_reset_words = 0
for h in range(151):
    a = 2
    while 2**(h+a) < 2*3**(h+1)-2**(h+1):
        a += 1
    word = (1,)*h+(a,)
    s,m,c,r = affine(word)
    require(c == 3**(h+1)-2**(h+1), 'first-reset carry closed form')
    require(c <= m-s and Q(s,m) <= Q(3,4), 'first-reset budget and uniform slope')
    n = r+2*m if r == 1 else r
    actual = literal_trace(n,sum(word))
    first = next(i for i,y in enumerate(actual) if y < n)
    N = 3**(h+1)*(n+1)-2**(h+1)
    j = next(j for j in range(1,a+1) if N < 2**(h+j)*n)
    require(first == h+j and odd_trace(n,word,True) is not None, 'first-reset literal first exit')
    first_reset_words += 1

schedule_count = 0
for length in range(1,9):
    for names in product(bank, repeat=length):
        slope, intercept = Q(1), Q(0)
        for name in names:
            s,m,c = macro_data[name]['coefficients']
            slope,intercept = Q(s,m)*slope, Q(s,m)*intercept+Q(c,m)
        require(intercept <= 4*(1-slope), 'arbitrary bank finite semigroup control')
        schedule_count += 1

for depth in range(1,201):
    s,m,c,r = affine(bank['A']*depth)
    require(37*c == 23*(m-s), 'repeated-A closed intercept')
    require((37*r-23) % (2*m) == 0, 'unique nested two-adic residue')
    s,m,c,r = affine((2,)*depth)
    require(r == 1 and 2*m == 2*4**depth, 'repeated-two exact cylinder')
require(Q(23,37) != int(Q(23,37)), 'no integer equals the repeated-A address')
fixed_trace = [Q(23,37)]
for a in bank['A']:
    z = 3*fixed_trace[-1]+1
    numerator = z.numerator
    valuation = 0
    while numerator % 2 == 0:
        valuation += 1
        numerator //= 2
    require(z.denominator % 2 == 1 and valuation == a, 'two-adic rational exact valuation')
    fixed_trace.append(z/2**a)
require(fixed_trace == [Q(23,37), Q(53,37), Q(49,37), Q(23,37)], 'rational macro cycle')

result = {
    'status': 'PASS', 'checks': checks,
    'independent_affine_words': affine_words,
    'native_recurrence_words': native_words,
    'valuation_words': valuation_words,
    'positive_odd_residue_representatives': residues,
    'inverse_candidates': inverse_candidates,
    'ROOT_trimmed_words': root_trimmed_words,
    'macro_schedules': schedule_count,
    'first_reset_words': first_reset_words,
    'repeated_word_depth': 200,
    'macros': macro_data,
    'wording_controls': {
        'nonnegative_debt_only_sign_hostile': str(defect_difference),
        'repeated_A_integer_intersection': 'empty',
        'repeated_A_Z2_intersection': '23/37',
        'positive_source_bank_exit': [11,17,13,5],
    },
    'global_Collatz_coverage_or_convergence_claim': False,
}
print(json.dumps(result, indent=2, sort_keys=True))
