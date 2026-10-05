"""Rational mixtures of rooted sibling flows and finite backward approximants.

All exported source weights require an actual finite ROOT certificate.
Backward box construction does not search unknown forward orbits.
"""
from dataclasses import dataclass
from fractions import Fraction as F
from math import comb, factorial
import json


def need(ok, message):
    if not ok:
        raise ValueError(message)


def natural(n, name='parameter'):
    need(type(n) is int and n >= 0, name+' must be an exact nonnegative integer')


def odd(n):
    need(type(n) is int and n > 0 and n % 2, 'exact positive odd source')


def sibling(n, j):
    odd(n)
    natural(j)
    return 4**j*n+(4**j-1)//3


def base(n):
    odd(n)
    j = 0
    while n % 8 == 5:
        n = (n-1)//4
        j += 1
    return n, j


def valuation(n):
    need(type(n) is int and n > 0, 'positive exact valuation argument')
    return (n & -n).bit_length()-1


def predecessor(y):
    odd(y)
    if y % 3 == 0:
        return None
    return ((4 if y % 3 == 1 else 2)*y-1)//3


def fraction_rate(rho):
    need(type(rho) in (int, F) and 0 < rho < 1, 'exact rational rate in(0,1)')
    return F(rho)


def fixed_r_weight(T, L, r, rho=F(1, 2)):
    natural(T); natural(L)
    r, rho = fraction_rate(r), fraction_rate(rho)
    return rho**T*(1-r)**T*r**L


def mixed_weight(T, L, rho=F(1, 2)):
    natural(T); natural(L)
    rho = fraction_rate(rho)
    return 2*rho**T*F(factorial(L)*factorial(T+1), factorial(L+T+2))


def double_weight(T, L):
    natural(T); natural(L)
    return F(4*factorial(L)*factorial(T), (T+2)*factorial(L+T+2))


def harmonic(n):
    natural(n)
    return sum((F(1, j) for j in range(1, n+1)), F(0))


def box_errors(N, M, rho=F(1, 2)):
    need(type(N) is int and N >= 1 and type(M) is int and M >= 2, 'positive box, M>=2')
    rho = fraction_rate(rho)
    single = 2*rho**N/(1-rho)+F(2, M+1)/(1-rho)**2
    double = F(4, N+1)+4*(harmonic(M)-1)/(M-1)
    return single, double


def replay(source, word):
    odd(source)
    need(type(word) is tuple and all(type(a) is int and a >= 1 for a in word),
         'exact positive valuation tuple')
    n = source
    for a in word:
        need(n != 1, 'no outgoing ROOT padding')
        raw = 3*n+1
        need(a == valuation(raw), 'actual valuation mismatch')
        n = raw >> a
    need(n == 1, 'certificate must finish at ROOT')


def shifted_word(b, word, j):
    """Transport a supplied ROOT word to its j-th sibling, with the ROOT exception."""
    replay(b, word)
    natural(j)
    if j == 0:
        return word
    if b == 1:
        return (2+2*j,)
    return (word[0]+2*j,)+word[1:]


@dataclass(frozen=True)
class Receipt:
    source: int
    T: int
    K: int
    j: int
    word: tuple

    @property
    def L(self):
        return self.K+self.j


def decode_receipt(source, word):
    """Authenticate source and derive all weight parameters; no guessed rank accepted."""
    replay(source, word)
    b, j = base(source)
    T = K = 0
    # Each base move decreases a supplied strict ROOT rank by at least one.
    for _ in range(len(word)+1):
        if b == 1:
            return Receipt(source, T, K, j, word)
        y = (3*b+1) >> valuation(3*b+1)
        b, k = base(y)
        K += k
        T += 1
    raise ValueError('quotient path exceeds authenticated ROOT rank')


def receipt_weight(source, word, rho=None):
    receipt = decode_receipt(source, word)
    return (double_weight(receipt.T, receipt.L) if rho is None else
            mixed_weight(receipt.T, receipt.L, rho))


def backward_box(N, M):
    """All rooted sources with T<N and K+j<M; finite inverse generation only."""
    need(type(N) is int and N >= 1 and type(M) is int and M >= 1, 'positive exact box')
    nodes = {1: (0, 0, ())}
    frontier = [1]
    for T in range(N-1):
        following = []
        for parent in frontier:
            oldT, K, word = nodes[parent]
            need(oldT == T, 'breadth-first base depth')
            for k in range(M-K):
                target = sibling(parent, k)
                child = predecessor(target)
                if child is None or child == 1:
                    continue
                need(child not in nodes, 'inverse base address must be unique')
                a = valuation(3*child+1)
                tail = shifted_word(parent, word, k)
                actual = (a,)+tail
                replay(child, actual)
                nodes[child] = (T+1, K+k, actual)
                following.append(child)
        frontier = following
    sources = {}
    for b, (T, K, word) in nodes.items():
        for j in range(M-K):
            n = sibling(b, j)
            need(n not in sources, 'unique primitive sibling base')
            actual = shifted_word(b, word, j)
            replay(n, actual)
            sources[n] = Receipt(n, T, K, j, actual)
    return sources


def literal_word(n, cap):
    """Independent finite-test helper; not used by the backward constructor."""
    odd(n)
    out = []
    while n != 1:
        need(len(out) < cap, 'finite control horizon exceeded')
        a = valuation(3*n+1)
        out.append(a)
        n = (3*n+1) >> a
    return tuple(out)


def rising(a, count):
    value = 1
    for i in range(count):
        value *= a+i
    return value


def main():
    checks = 0

    def check(ok, message):
        nonlocal checks
        checks += 1
        need(ok, message)

    rho = F(1, 2)
    for T in range(9):
        for L in range(13):
            integral = sum((F((-1)**i*comb(T+1, i), L+i+1)
                            for i in range(T+2)), F(0))
            check(mixed_weight(T, L, rho) == 2*rho**T*integral,
                  'independent exact beta integral')
            check(double_weight(T, L) == mixed_weight(T, L, rho)/rho**T
                  *F(2, (T+1)*(T+2)), 'independent rho integration')
            count = comb(L+T, T)
            check(count*mixed_weight(T, L, rho)
                  == F(2*(T+1), (L+T+1)*(L+T+2))*rho**T,
                  'weak-composition single-mixture majorant')
            check(count*double_weight(T, L)
                  == F(4, (T+2)*(L+T+1)*(L+T+2)), 'double majorant')
            for r in (F(1, 16), F(1, 4), F(1, 2), F(3, 4), F(15, 16)):
                for rate in (F(1, 3), F(1, 2), F(3, 4)):
                    fixed = fixed_r_weight(T, L, r, rate)
                    total = T+L
                    check(mixed_weight(T, L, rate) >= fixed*F(2*(T+1), (total+1)*(total+2)),
                          'single mixture polynomial comparison')
                    check(double_weight(T, L) >= fixed*F(4, (T+2)*(total+1)*(total+2)),
                          'double mixture polynomial comparison')
    for T in range(7):
        for L in range(10):
            check(mixed_weight(T, L+1)/mixed_weight(T, L) == F(L+1, L+T+3),
                  'adaptive sibling ratio')
            for k in range(7):
                factor = F((T+2)*rising(L+1, k), rising(L+T+3, k+1))
                check(mixed_weight(T+1, L+k)/mixed_weight(T, L) == rho*factor,
                      'adaptive base-edge ratio')
                check(double_weight(T+1, L+k)/double_weight(T, L)
                      == F(T+1, T+3)*factor, 'double posterior update')
    for J in range(1, 41):
        rootmass = sum((double_weight(0, j) for j in range(J)), F(0))
        check(rootmass+F(2, J+1) == 2, 'exact infinite rooted-ray tail')
        uniform = sum((F(1, j+1) for j in range(J)), F(0))
        check(uniform == harmonic(J), 'uniform refuel mixture harmonic hostile')

    box = backward_box(6, 12)
    larger = backward_box(7, 14)
    histogram = {}
    for n, receipt in box.items():
        check(decode_receipt(n, receipt.word) == receipt, 'source and weight fields independently decoded')
        check(receipt_weight(n, receipt.word) == double_weight(receipt.T, receipt.L),
              'authenticated source weight')
        check(larger[n] == receipt, 'nested backward approximants')
        key = (receipt.T, receipt.L)
        histogram[key] = histogram.get(key, 0)+1
    for (T, L), count in histogram.items():
        check(count <= comb(T+L, T), 'actual finite addresses below composition bound')
    for n in range(1, 4096, 2):
        word = literal_word(n, 1000)
        decoded = decode_receipt(n, word)
        expected = decoded.T < 6 and decoded.L < 12
        check((n in box) == expected, 'complete box membership on independent literal controls')
    normF = sum((mixed_weight(x.T, x.L) for x in box.values()), F(0))
    normD = sum((double_weight(x.T, x.L) for x in box.values()), F(0))
    check(normF < 3 and normD < 3, 'finite lower norms obey the sharpened global bound')
    # Exact infinite incoming sum for every target in a finite authenticated universe.
    incoming_checks = 0
    for n, rec in box.items():
        if n == 1:
            continue
        b = predecessor(n)
        if b is None:
            check(n % 3 == 0, 'empty incoming row remains strict, not equality')
            continue
        a = valuation(3*b+1)
        inverse = decode_receipt(b, (a,)+rec.word)
        T, L = inverse.T, inverse.K
        incomingF = 2*rho**T*F(factorial(L)*factorial(T), factorial(L+T+1))
        incomingD = F(4*factorial(L)*factorial(T), (T+1)*(T+2)*factorial(L+T+1))
        check(incomingF == rho*mixed_weight(rec.T, rec.L), 'single-mixture full incoming sum')
        check(incomingD == F(rec.T+1, rec.T+3)*double_weight(rec.T, rec.L),
              'double-mixture exact strict row ratio')
        incoming_checks += 1
    hostiles = [lambda: decode_receipt(True, ()), lambda: decode_receipt(1.0, ()),
                lambda: decode_receipt(2, ()), lambda: decode_receipt(1, (2,)),
                lambda: decode_receipt(3, (1, 2)), lambda: mixed_weight(True, 0),
                lambda: double_weight(0, -1), lambda: backward_box(0, 5),
                lambda: receipt_weight(3, (1.0, 4)), lambda: fixed_r_weight(0, 0, 1)]
    for action in hostiles:
        try:
            action()
        except ValueError:
            check(True, 'invalid certificate or parameter rejected')
        else:
            raise ValueError('hostile accepted')
    eF, eD = box_errors(64, 4096)
    simple_double_bound = F(4, 65)+F(48, 4095)
    check(eD <= simple_double_bound, 'dyadic harmonic bound gives compact rational error')
    print(json.dumps({
        'checks': checks, 'rho': str(rho),
        'backward_box_T_lt6_L_lt12': {'sources': len(box), 'base_sources': sum(x.j == 0 for x in box.values()),
            'single_weight_sum': str(normF), 'double_weight_sum': str(normD)},
        'larger_box_T_lt7_L_lt14_sources': len(larger),
        'independent_literal_membership_sources': 2048, 'exact_incoming_rows': incoming_checks,
        'rooted_sibling_family_total_weight_both_mixtures': 2,
        'certified_parameter_box_64_4096_error_bounds': {'single': str(eF),
            'double_exact_expression': '4/65 + 4*(H_4096-1)/4095',
            'double_rational_upper_bound': str(simple_double_bound)},
        'hostiles': len(hostiles)}, indent=2, sort_keys=True))
    print('PASS: finite backward receipts and unconditional weight tails; all-source positivity remains OPEN.')


if __name__ == '__main__':
    main()
