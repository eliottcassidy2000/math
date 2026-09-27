"""Exact exponent-lifted families with a recursive Collatz certificate.

The 27-containing family is constructed, not a coverage theorem for all
integers or for every 4*8**k-5. All checks survive python -O.
"""
from collections import Counter
import json

checks = Counter()


def need(ok, label):
    checks[label] += 1
    if not ok:
        raise RuntimeError(label)


def valuation(n, p):
    if n <= 0:
        raise ValueError('positive valuation input required')
    a = 0
    while n % p == 0:
        a += 1
        n //= p
    return a


def U(n):
    z = 3*n+1
    return z // (z & -z)


def twoadic(n):
    return (n & -n).bit_length()-1


def exponent_lifts(c, d, depth):
    """Unique h_m in [0,3**m) with c*4**h_m+d=0 mod3**(m+1)."""
    need(c % 3 != 0 and (c+d) % 3 == 0, 'lift initial unit')
    h = 0
    out = [h]
    digits = []
    for m in range(depth):
        modulus = 3**(m+2)
        residue = (c*pow(4, h, modulus)+d) % modulus
        need(residue % 3**(m+1) == 0, 'lift prior precision')
        error = residue // 3**(m+1)
        digit = (-error*pow(c, -1, 3)) % 3
        h += digit*3**m
        need((c*pow(4, h, modulus)+d) % modulus == 0, 'lift new precision')
        need(0 <= h < 3**(m+1), 'least exponent representative')
        out.append(h)
        digits.append(digit)
    return out, digits


def family47(k, h):
    numerator = 2*(47*(1 << (2*h))+7)
    denominator = 3**(2*k+1)
    if numerator % denominator:
        raise ValueError('inadmissible family parameter')
    coefficient = numerator // denominator
    return (coefficient << (3*k))-5


def decode47(n):
    """Recognize the k>=1 family from the actual integer, without orbit search."""
    if n == 27:
        return 1, 0
    if n <= 1 or n % 2 == 0:
        return None
    e = twoadic(n+5)
    if e < 4 or e % 3 != 1:
        return None
    k = (e-1)//3
    a = (n+5) >> (3*k)
    numerator = 3**(2*k+1)*a-14
    if numerator <= 0 or numerator % 94:
        return None
    power = numerator // 94
    if power & (power-1):
        return None
    exponent = power.bit_length()-1
    if exponent % 2:
        return None
    h = exponent//2
    return k, h


def first_descent(n, cap):
    current = n
    for j in range(1, cap+1):
        current = U(current)
        if current < n:
            return j, current
    return None


def main():
    lifts, digits = exponent_lifts(47, 7, 60)
    for m in range(1, 9):
        modulus = 3**(m+1)
        residue = 1
        solutions = []
        for h in range(3**m):
            if (47*residue+7) % modulus == 0:
                solutions.append(h)
            residue = 4*residue % modulus
        need(solutions == [lifts[m]], 'independent exponent census')
    for m in range(61):
        need(pow(4, 3**m, 3**(m+2)) == 1+3**(m+1), 'elementary lifting coefficient')
    for h in range(41):
        for d in range(1, 101):
            difference = 47*4**h*(4**d-1)
            need(valuation(difference, 3) == 1+valuation(d, 3), 'ternary index isometry')
    # One complete exponent period: count depth min(K(h),4) without big powers.
    depth_census = Counter()
    modulus = 3**9
    power = 1
    for h in range(9**4):
        residue = (47*power+7) % modulus
        depth = 4 if residue == 0 else (valuation(residue, 3)-1)//2
        depth_census[depth] += 1
        power = 4*power % modulus
    need(dict(depth_census) == {0:5832, 1:648, 2:72, 3:8, 4:1}, 'complete exponent depth distribution')

    suffix = [47]
    while suffix[-1] != 1:
        need(len(suffix) < 100, 'fixed suffix cap')
        suffix.append(U(suffix[-1]))
    need(len(suffix)-1 == 38, 'fixed47 suffix length')
    need(first_descent(27, 50) == (37, 23), '27 finite first descent')

    family_rows = []
    for k in range(1, 7):
        for t in range(3):
            h = lifts[2*k]+t*9**k
            n = family47(k, h)
            need(decode47(n) == (k, h), 'direct membership decoding')
            need(n % 8 == 3 and twoadic(3*(n-2)+1) == 2, 'fixed q1 R1 interface')
            current = n
            for j in range(k):
                need(current == family47(k-j, h), 'recursive actual phase boundary')
                need(twoadic(3*current+1) == 1, 'first half valuation')
                middle = U(current)
                need(twoadic(3*middle+1) == 2, 'second half valuation')
                need(middle > n, 'all first phase odd steps above source')
                current = U(middle)
                need(current > n, 'all first phase even steps above source')
                need((8*current-5) % 9 == 0 and (8*current-5)//9 == family47(k-j, h), 'actual inverse block')
            expected = (94*(1 << (2*h))-1)//3
            need(current == expected, 'exact preterminal endpoint')
            need(twoadic(3*current+1) == 2*h+1, 'terminal division clock')
            need(U(current) == 47, 'terminal target47')
            need(current > 1 and all(s > 1 for s in suffix[:-1]), 'no earlier terminal1')
            tau, lower = first_descent(n, 2*k+40)
            need(tau == (37 if n == 27 else 2*k+1), 'exact first descent time')
            if n != 27:
                need(n > 47 and lower == 47 and twoadic(n+5) == 3*k+1, 'ordinary family decoder type')
            family_rows.append(dict(k=k, t=t, h=h, terminal_division=2*h+1,
                                    source_bits=n.bit_length(), first_descent=tau,
                                    first_source=n if n.bit_length() <= 128 else None))

    # Independent finite membership census from direct exponent enumeration.
    generated = {family47(1, h) for h in range(28) if (47*4**h+7) % 27 == 0}
    for k in range(2, 6):
        need(family47(k, lifts[2*k]) > 100000, 'small decoder higher depth excluded')
    need(2*8**6-5 > 100000 and family47(1, 36) > 100000, 'small decoder census bounds')
    accepted = {n for n in range(3, 100000, 2) if decode47(n) is not None}
    need(accepted == {n for n in generated if n < 100000} == {27}, 'small-source decoder census')
    for k in range(2, 129):
        need(decode47(4*8**k-5) is None, 'fixed coefficient family hostile')

    # Target1 control uses the same clock and closes directly after the long rise.
    ones, _ = exponent_lifts(1, 14, 12)
    target1_rows = []
    for k in range(1, 7):
        h = ones[2*k]
        a = ((1 << (2*h))+14)//3**(2*k+1)
        n = (a << (3*k))-5
        need(a % 2 == 0, 'target1 positive even coefficient')
        need(first_descent(n, 2*k+1) == (2*k+1, 1), 'target1 exact first descent')
        target1_rows.append(dict(k=k, h=h, source_bits=n.bit_length(), first_descent=2*k+1))

    # Mersenne endpoint valuations alone do not pay exponential growth.
    mersenne_descending = []
    for H in range(1, 1025):
        source = (1 << H)-1
        raw = 3**H-1
        exponent = 1 if H % 2 else 2+twoadic(H)
        need(twoadic(raw) == exponent, 'Mersenne exact clipping valuation')
        endpoint = raw >> exponent
        if endpoint < source:
            mersenne_descending.append(H)
        if H >= 9:
            need(raw > 4*H*source and endpoint > source, 'Mersenne all-large endpoint hostile')
        if H <= 128:
            current = source
            for _ in range(H):
                current = U(current)
            need(current == endpoint, 'Mersenne actual clock')
    need(mersenne_descending == [2, 4, 8], 'Mersenne endpoint finite boundary')

    print(json.dumps(dict(
        status='PROVED scoped family, FINITE-EXACT controls, OPEN universal coverage',
        target47_minimal_exponents=[lifts[2*k] for k in range(1, 13)],
        target47_initial_ternary_digits=digits[:24],
        exponent_depth_census_through3power8=dict(sorted(depth_census.items())),
        target47_suffix=suffix,
        family47_controls=family_rows,
        target1_controls=target1_rows,
        mersenne_phase_end_descending_H=mersenne_descending,
        universes=dict(exponent_lift_depth=60, independent_exponent_depth=8,
                       actual_family_k=6, h_period_offsets=3, membership_source_bound=100000,
                       exponent_depth_census_period=6561,
                       mersenne_valuation_H=1024, mersenne_orbit_H=128),
        checks=dict(sorted(checks.items())), total=sum(checks.values()), result='PASS'
    ), indent=2))


if __name__ == '__main__':
    main()
