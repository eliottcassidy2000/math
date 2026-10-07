"""Exact dyadic joint cylinders, retained information, and finite-coset flushing.

No statistical independence of the full paired streams is assumed. ROOT
is not absorbing here: this is the ordinary odd 2-adic dynamical map.
"""
from fractions import Fraction as F
from itertools import product
from collections import Counter
from math import comb
import json


def need(ok, message):
    if not ok:
        raise ValueError(message)


def v2(n):
    need(type(n) is int and n != 0, "nonzero exact integer required")
    n = abs(n)
    return (n & -n).bit_length()-1


def step(n, sign=1):
    need(type(n) is int and n % 2 == 1 and type(sign) is int and sign in (-1, 1), "odd source and sheet required")
    value = 3*n+sign
    a = v2(value)
    return value//(1 << a), a


def word(n, length, sign=1):
    need(type(length) is int and length >= 0, "nonnegative exact length required")
    need(type(n) is int and n % 2 == 1 and type(sign) is int and sign in (-1, 1),
         "odd source and sheet required")
    result = []
    for _ in range(length):
        n, a = step(n, sign)
        result.append(a)
    return tuple(result)


def cylinder(w):
    need(type(w) is tuple and all(type(a) is int and a >= 1 for a in w), "valuation tuple required")
    p, carry, total = 1, 0, 0
    for a in w:
        p, carry, total = 3*p, 3*carry+(1 << total), total+a
    modulus = 1 << (total+1)
    r = ((1 << total)-carry)*pow(p, -1, modulus) % modulus
    return r, total


def pullback(w, alpha, v, beta):
    need(type(alpha) is int and alpha % 2 and type(v) is int and v >= 0,
         "odd affine unit and nonnegative scale required")
    need(type(beta) is int and beta % 2 == (1 if v else 0), "odd-output affine map required")
    r, total = cylinder(w)
    power = total+1
    if power <= v:
        return (1, 0) if (r-beta) % (1 << power) == 0 else None
    if (r-beta) % (1 << v):
        return None
    modulus = 1 << (power-v)
    residue = ((r-beta)//(1 << v))*pow(alpha, -1, modulus) % modulus
    return (residue, total-v) if residue % 2 else None


def joint(wy, wx, alpha, v, beta):
    ry, ay = cylinder(wy)
    pulled = pullback(wx, alpha, v, beta)
    if pulled is None:
        return None
    rx, ex = pulled
    if (ry-rx) % (1 << (min(ay, ex)+1)):
        return None
    return {'mass': F(1, 1 << max(ay, ex)), 'y_mass': F(1, 1 << ay),
            'x_mass': F(1, 1 << ex), 'shared_bits': min(ay, ex),
            'precision_bits': max(ay, ex)}


def flush_coset(residue, bits):
    """Push forward Haar on residue mod2^(bits+1) until full odd Haar.

    Only the coset law is returned; endpoint samples and word histories
    are not discarded when used in an actual source certificate.
    """
    need(type(bits) is int and bits >= 0 and type(residue) is int and residue % 2,
         "odd dyadic coset required")
    residue %= 1 << (bits+1)
    steps = 0
    while bits:
        z = (3*residue+1) % (1 << (bits+1))
        if z == 0:
            return steps+1
        a = v2(z)
        bits = max(0, bits-a)
        residue = ((3*residue+1)//(1 << a)) % (1 << (bits+1))
        steps += 1
    return steps


def nb_tail(s, bits):
    need(type(s) is int and s >= 1 and type(bits) is int and bits >= 0,
         "positive step count and nonnegative bit budget required")
    return F(sum(comb(bits, k) for k in range(min(s-1, bits)+1)), 1 << bits)


def covariance_squared_bound(s, t, v):
    need(type(v) is int and v >= 0 and type(t) is int and t >= v+1,
         "second observation must follow the initial coset flush")
    return min(F(4), 38*nb_tail(s, t-1-v))


def main():
    checks = 0
    def check(ok, message):
        nonlocal checks
        need(ok, message); checks += 1

    max_bits = 6
    words = [()]
    for length in range(1, 4):
        words += [w for w in product(range(1, max_bits+1), repeat=length) if sum(w) <= max_bits]
    maps = [(1, 0, 0), (3, 0, 0), (-1, 0, 0), (3, 1, 1), (3, 2, 1), (3, 5, 1)]
    for alpha, v, beta in maps:
        observed = Counter()
        for y in range(1, 1 << (max_bits+1), 2):
            x = alpha*(1 << v)*y+beta
            yw = [word(y, k) for k in range(4)]
            xw = [word(x, k) for k in range(4)]
            for wy in yw:
                for wx in xw:
                    if sum(wy) <= max_bits and sum(wx) <= max_bits:
                        observed[(wy, wx)] += 1
        for wy in words:
            for wx in words:
                result = joint(wy, wx, alpha, v, beta)
                actual = F(observed[(wy, wx)], 1 << max_bits)
                check(actual == (result['mass'] if result else 0), "independent joint cylinder census")
                if result:
                    check(result['mass']/(result['x_mass']*result['y_mass']) ==
                          2**result['shared_bits'], "pointwise mutual-information density")

    histogram = Counter()
    for bits in range(0, 12):
        longest = 0
        for r in range(1, 1 << (bits+1), 2):
            count = flush_coset(r, bits)
            check(count <= bits, "deterministic finite-information forgetting time")
            longest = max(longest, count)
        check(longest == bits, "all-valuation-one coset attains cutoff")
        histogram[bits] = longest

    for x in range(1, 100, 2):
        for y in range(1, 100, 2):
            if x != y:
                gap = v2(y-x)
                wx, wy = word(x, gap+1), word(y, gap+1)
                common = next(k for k,(a,b) in enumerate(zip(wx,wy)) if a != b)
                check(sum(wx[:common]) < gap <= sum(wx[:common+1]), "same-sheet valuation prefix cutoff")
            gap = v2(x+y)
            wx, wy = word(x, gap+1), word(y, gap+1, -1)
            common = next(k for k,(a,b) in enumerate(zip(wx,wy)) if a != b)
            check(sum(wx[:common]) < gap <= sum(wx[:common+1]), "opposite-sheet reflection cutoff")

    for s in range(1, 8):
        for bits in range(0, 16):
            counted = sum(sum(w) < s for w in product((0,1), repeat=bits))
            check(nb_tail(s,bits) == F(counted, 1 << bits), "independent binomial tail enumeration")

    table = [{'s':s, 't':t, 'v':2, 'covariance_squared_upper':str(covariance_squared_bound(s,t,2))}
             for s in (1,4,16) for t in (16,32,64,128)]
    hostiles = (lambda:cylinder((True,)), lambda:pullback((1,),2,1,1),
                lambda:pullback((1,),3,2,0), lambda:flush_coset(2,4),
                lambda:nb_tail(0,2), lambda:covariance_squared_bound(1,2,2),
                lambda:step(1,True), lambda:word(1,True), lambda:word(1,-1),
                lambda:word(True,0), lambda:word(2,0), lambda:word(1,0,True))
    for bad in hostiles:
        try: bad()
        except ValueError: check(True,"malformed input rejected")
        else: check(False,"malformed input accepted")
    print(json.dumps({'status':'PASS; exact prefix laws and long-lag bounds, recurrence OPEN',
                      'checks':checks, 'joint_word_max_sum':max_bits, 'affine_maps':maps,
                      'flush_worst_steps_by_known_bits':histogram,
                      'long_gap_bounds':table,
                      'scope':'Whole-history information coupling does not determine near-diagonal recurrence'},
                     indent=2,sort_keys=True))


if __name__ == '__main__': main()
