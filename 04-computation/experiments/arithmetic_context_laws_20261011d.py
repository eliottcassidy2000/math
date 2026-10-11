"""Exact sampling identities and two finite arithmetic context switches.

No asymptotic theorem or prime-representation conjecture is inferred from tests.
The wheel transfer preserves divisibility data, not primality or intervals.
"""

from collections import Counter
from fractions import Fraction
from itertools import product
from math import isqrt, prod

from abc_difference_families_20261011c import factors, prime


def integer(n, lower=0):
    if type(n) is not int or n < lower:
        raise ValueError('exact integer above stated bound required')


def prime_set(primes):
    if type(primes) is not tuple or any(type(p) is not int or not prime(p) for p in primes):
        raise ValueError('prime moduli required')
    if primes != tuple(sorted(set(primes))):
        raise ValueError('sorted tuple of distinct primes required')


def mobius(n):
    integer(n, 1)
    ff = factors(n)
    return 0 if any(a > 1 for _, a in ff) else (-1)**len(ff)


def point_statistics(n):
    integer(n, 1)
    ff = factors(n)
    return (prod(a + 1 for _, a in ff), len(ff), sum(a for _, a in ff),
            int(prime(n)), int(all(a == 1 for _, a in ff)))


def sample_sums(n):
    """Totals of tau, omega, Omega, prime indicator, squarefree indicator."""
    integer(n, 1)
    primes = tuple(p for p in range(2, n + 1) if prime(p))
    omega = sum(n // p for p in primes)
    bigomega = 0
    for p in primes:
        power = p
        while power <= n:
            bigomega += n // power
            power *= p
    return (sum(n // d for d in range(1, n + 1)), omega, bigomega,
            len(primes), sum(mobius(d) * (n // (d*d)) for d in range(1, isqrt(n) + 1)))


def sample_means(n):
    integer(n, 1)
    return tuple(Fraction(x, n) for x in sample_sums(n))


def squarefree_flag_law(primes, flags):
    """Exact conditional law on Z/(product p^2), given no tested p^2 divides."""
    prime_set(primes)
    if type(flags) is not tuple or len(flags) != len(primes) or any(type(x) is not int or x not in (0, 1) for x in flags):
        raise ValueError('one exact divisibility bit per prime required')
    return prod(Fraction(1 if bit else p, p + 1) for p, bit in zip(primes, flags))


def forms_input(forms):
    if type(forms) is not tuple or not forms:
        raise ValueError('nonempty tuple of affine forms required')
    if any(type(f) is not tuple or len(f) != 2 or any(type(x) is not int for x in f) for f in forms):
        raise ValueError('exact affine integer coefficient pairs required')


def local_bad(forms, p):
    forms_input(forms)
    if type(p) is not int or not prime(p):
        raise ValueError('prime required')
    return tuple(x for x in range(p) if any((a*x+b) % p == 0 for a, b in forms))


def local_factor(forms, p):
    """Normalization is against k generic unit densities (1-1/p)^k."""
    return Fraction(p - len(local_bad(forms, p)), p) / Fraction(p - 1, p)**len(forms)


def wheel(forms, primes):
    forms_input(forms)
    prime_set(primes)
    modulus = prod(primes)
    return modulus, tuple(x for x in range(modulus)
                          if all(all((a*x+b) % p for a, b in forms) for p in primes))


def pair_forms(n, kind):
    integer(n, 1)
    if type(kind) is not str or kind not in ('goldbach', 'lemoine'):
        raise ValueError('exact pair context required')
    return ((1, 0), (-1 if kind == 'goldbach' else -2, n))


def transfer_multiplier(primes):
    prime_set(primes)
    if not primes or primes[0] != 2:
        raise ValueError('wheel must include 2')
    q = prod(primes[1:])
    return 2*q, q + 4


def transfer_residue(n, primes, residue, inverse=False):
    """For odd n: Lemoine n residue <-> Goldbach 2n residue on one wheel."""
    integer(n, 1)
    if n % 2 != 1 or type(inverse) is not bool:
        raise ValueError('odd target and exact direction flag required')
    modulus, multiplier = transfer_multiplier(primes)
    integer(residue)
    if residue >= modulus:
        raise ValueError('canonical input residue required')
    if inverse:
        multiplier = pow(multiplier, -1, modulus)
    return multiplier * residue % modulus


def prime_representations(n, coefficient):
    """Actual ordered p+coefficient*q=n receipts, not a probability model."""
    integer(n, 1)
    if type(coefficient) is not int or coefficient not in (1, 2):
        raise ValueError('coefficient 1 or 2 required')
    return tuple((n - coefficient*q, q)
                 for q in range(2, (n - 2)//coefficient + 1)
                 if prime(q) and prime(n - coefficient*q))


def mm(a, b, modulus=None):
    result = tuple(tuple(sum(a[i][k]*b[k][j] for k in range(2)) for j in range(2)) for i in range(2))
    return result if modulus is None else tuple(tuple(x % modulus for x in row) for row in result)


def golden_product(x, y, p):
    """Internal exact coordinate product in F_p[t]/(t^2-t-1)."""
    a, b = x
    c, d = y
    return ((a*c+b*d) % p, (a*d+b*c+b*d) % p)


def golden_quotient(p):
    if type(p) is not int or not prime(p):
        raise ValueError('exact prime required')
    roots = tuple(x for x in range(p) if (x*x-x-1) % p == 0)
    kind = {0: 'inert', 1: 'ramified', 2: 'split'}[len(roots)]
    state = (1, 0)
    for order in range(1, p*p + 1):
        state = golden_product(state, (0, 1), p)
        if state == (1, 0):
            break
    else:
        raise RuntimeError('unit order not found in finite quotient')
    scalar_orders = tuple(next(j for j in range(1, p) if pow(x, j, p) == 1) for x in roots)
    return kind, roots, order, scalar_orders


def main():
    checks = 0
    def check(ok, why):
        nonlocal checks
        checks += 1
        if not ok:
            raise RuntimeError(why)
    def reject(call):
        nonlocal checks
        checks += 1
        try:
            call()
        except ValueError:
            return
        raise RuntimeError('malformed input accepted')

    cumulative = [0]*5
    for n in range(1, 257):
        row = point_statistics(n)
        cumulative = [a+b for a, b in zip(cumulative, row)]
        check(sample_sums(n) == tuple(cumulative), 'independent summatory identities')
        d = isqrt(n)
        check(cumulative[0] == 2*sum(n//j for j in range(1, d+1)) - d*d, 'hyperbola exact identity')
    samples = []
    for n in (100, 1000, 10000):
        direct = tuple(sum(point_statistics(i)[j] for i in range(1, n+1)) for j in range(5))
        check(sample_sums(n) == direct, 'larger exact sample')
        samples.append((n, sample_means(n)))
    check(point_statistics(1) == (1, 0, 0, 0, 1), 'unit boundary')
    for p in (2, 3, 5, 7, 11, 101):
        check(point_statistics(p) == (2, 1, 1, 1, 1), 'prime-conditioned context')
    check(point_statistics(12) == (6, 2, 3, 0, 0), 'distinct versus repeated factors')

    squarefree_tables = []
    for primes in ((2,), (2, 3), (2, 3, 5), (2, 3, 5, 7)):
        modulus = prod(p*p for p in primes)
        counts = Counter(tuple(int(x % p == 0) for p in primes)
                         for x in range(modulus) if all(x % (p*p) != 0 for p in primes))
        size = sum(counts.values())
        check(size == prod(p*p-1 for p in primes), 'squarefree wheel size')
        for flags in product((0, 1), repeat=len(primes)):
            check(Fraction(counts[flags], size) == squarefree_flag_law(primes, flags), 'exact conditional product law')
        squarefree_tables.append((primes, modulus, size))

    pairs_checked = transfer_points = 0
    for primes in ((2,), (2, 3), (2, 3, 5), (2, 3, 5, 7), (2, 3, 5, 7, 11)):
        for n in range(1, 66, 2):
            lf, gf = pair_forms(n, 'lemoine'), pair_forms(2*n, 'goldbach')
            modulus, left = wheel(lf, primes)
            _, right = wheel(gf, primes)
            check(len(left) == len(right) == prod(p-len(local_bad(lf,p)) for p in primes), 'CRT pair counts')
            check({transfer_residue(n, primes, x) for x in left} == set(right), 'source-bound finite wheel bijection')
            for p in primes:
                check(local_factor(lf, p) == local_factor(gf, p), 'identical normalized local factors')
                bad = 1 if p == 2 or n % p == 0 else 2
                check(len(local_bad(lf, p)) == bad, 'exact exceptional-prime split')
            for x in range(modulus):
                y = transfer_residue(n, primes, x)
                check(transfer_residue(n, primes, y, inverse=True) == x, 'full residue inverse')
                transfer_points += 1
            pairs_checked += 1
    check(len(wheel(pair_forms(30, 'goldbach'), (2,3,5))[1]) == 8, 'correlated target30')
    check(len(wheel(pair_forms(32, 'goldbach'), (2,3,5))[1]) == 3, 'same unit marginals different joint count')
    check(len(wheel(pair_forms(31, 'goldbach'), (2,))[1]) == 0, 'odd Goldbach parity obstruction beyond prime2 exception')
    check(local_bad(((0, 0), (1, 1)), 3) == (0, 1, 2), 'zero form has full local obstruction')
    check(local_factor(((0, 0), (1, 1)), 3) == 0, 'zero local factor')
    check(transfer_residue(13, (2,3,5), 7) == 13, 'wheel-specific multiplier')
    check(transfer_residue(13, (2,3,5,7), 7) == 133, 'refining wheel need not keep representative')

    representations = []
    for n in (7, 9, 11, 13, 15, 21, 35, 105):
        g, le = prime_representations(2*n, 1), prime_representations(n, 2)
        check(all(p+q == 2*n and prime(p) and prime(q) for p,q in g), 'actual Goldbach receipts')
        check(all(p+2*q == n and prime(p) and prime(q) for p,q in le), 'actual Lemoine receipts')
        representations.append((n, len(g), len(le)))
    check(len(prime_representations(26, 1)) == 5 and len(prime_representations(13, 2)) == 2, 'same local law is not same global count')
    check((3,2) in prime_representations(7,2) and 2 not in wheel(pair_forms(7,'lemoine'),(2,))[1], 'sieving prime itself is a finite boundary exception')

    # Recovered roles {2,3,11}: prime2 operator on conductor11 oldforms,
    # then reduce its lattice modulo3. This is not a prime-triple classification.
    a = ((-2, 1), (-2, 0)); identity = ((1,0),(0,1))
    power = identity; orders = []
    for j in range(1, 9):
        power = mm(power, a)
        if tuple(tuple(x % 3 for x in row) for row in power) == identity:
            orders.append(j)
        if j == 4:
            check(power == ((-4,0),(0,-4)), 'signed oldform fourth power')
    check(power == ((16,0),(0,16)) and orders == [8], 'mod3 golden order8, actual eighth power16')
    golden_rows = tuple((p, golden_quotient(p)) for p in (2,3,5,7,11,13))
    check(tuple(row[1][2] for row in golden_rows) == (3,8,20,16,10,28), 'exact golden unit orders')
    check(golden_quotient(11) == ('split',(4,8),10,(5,10)), 'split components retain distinct orders')
    nilpotent = (2,1)  # t-3 modulo5
    check(nilpotent != (0,0) and golden_product(nilpotent,nilpotent,5) == (0,0), 'ramified quotient is not a field')
    check(golden_quotient(5) == ('ramified',(3,),20,(4,)), 'residue field loses nilpotent period factor5')
    reject(lambda: sample_sums(True))
    reject(lambda: sample_sums(0))
    reject(lambda: squarefree_flag_law((2,3),(True,0)))
    reject(lambda: squarefree_flag_law((2,3),(1,)))
    reject(lambda: wheel(((1,0),),(2,2)))
    reject(lambda: wheel(((1,0),),(2,9)))
    reject(lambda: wheel(((1.0,0),),(2,3)))
    reject(lambda: local_factor(((1,0),),3.0))
    reject(lambda: transfer_residue(12,(2,3),1))
    reject(lambda: transfer_residue(13,(3,5),1))
    reject(lambda: transfer_residue(13,(2,3),6))
    reject(lambda: transfer_residue(13,(2,3),1,inverse=1))
    reject(lambda: prime_representations(15,True))
    reject(lambda: golden_quotient(5.0))
    reject(lambda: golden_quotient(9))
    reject(lambda: wheel(((1,0),),([2],)))

    print('Arithmetic context laws: exact sampling, squarefree conditioning, and pair-wheel transfer')
    print('Sample columns: E[tau], E[omega], E[Omega], P(prime), P(squarefree); uniform1..N')
    for n, row in samples:
        print('N', n, ':', ', '.join(str(x) for x in row))
    print('Squarefree finite-wheel (primes,modulus,admissible):', squarefree_tables)
    print('Conditional divisibility probabilities: p2=1/3,p3=1/4,p5=1/6,p7=1/8; fixed-set product law')
    print('Lemoine N / Goldbach2N: oddN1..65, five wheels;', pairs_checked, 'contexts;', transfer_points, 'full residue inverse controls')
    print('Same unit marginals on wheel30: Goldbach30 joint8, Goldbach32 joint3')
    print('Actual comparison (odd N, ordered Goldbach2N count, LemoineN count):', representations)
    print('Inherited roles2/3/11: U2 matrix reduces modulo3 to golden matrix of order8; conductor11 retained')
    print('Golden quotient (prime,(type,roots,phi-order,root-orders)):', golden_rows)
    print('Checks:', checks)
    print('Scope: exact finite laws only; PNT/mean asymptotics cited, Goldbach/Lemoine/Hardy-Littlewood still conjectural, Chen is prime+P2')


if __name__ == '__main__':
    main()
