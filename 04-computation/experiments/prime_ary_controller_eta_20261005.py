"""Exact controls for ternary controller codes, the eta11 local2 reader,
and the height boundary of a rooted ternary tree. Standard library only.

Run from the repository root with python3 -B (also with -O).
The JSON on stdout is the reproducible artifact; no files are mutated.
"""
from collections import Counter
from fractions import Fraction as F
from itertools import product
from math import comb
import json

# p, q, b, half-ternary-depth, binary-cost, credit, source-residue, period
OPS = {
    'H': (729, 1024, 669, 3, 10, 2, 155, 2048),
    'G': (9, 8, 5, 1, 3, -1, 11, 16),
    'A': (81, 128, 85, 2, 7, 1, 187, 256),
    'B': (81, 128, 73, 2, 7, 1, 7, 256),
    'L': (9, 16, -3, 1, 4, 4, 219, 256),
}
CHECKS = Counter()


def check(ok, category):
    CHECKS[category] += 1
    if not ok:
        raise ValueError(category)


def expand(word):
    return word.replace('H', 'LB')


def compress(word):
    return word.replace('LB', 'H')


def fuse(left, right):
    return compress(expand(left)+expand(right))


def stats(word):
    return tuple(sum(OPS[c][i] for c in word) for i in (3, 4, 5))


def sign(word):
    return (-1)**(len(word)-word.count('H'))


def words(limit, alphabet, native=False):
    stack = [('', 0)]
    while stack:
        word, grade = stack.pop()
        yield word, grade
        for c in alphabet:
            new = grade+OPS[c][3]
            if new <= limit and not (native and word.endswith('L') and c == 'B'):
                stack.append((word+c, new))


def affine(word):
    p, q, b = 1, 1, 0
    for c in word:
        r, s, t = OPS[c][:3]
        p, q, b = r*p, s*q, r*b+t*q
    return p, q, b


def native_at(word, value):
    for c in word:
        p, q, b, _, _, _, residue, modulus = OPS[c]
        if value % modulus != residue:
            return False
        numerator = p*value+b
        if numerator % q:
            raise ValueError('native map nonintegral')
        value = numerator//q
    return True


def native_port(word):
    if 'LB' in word:
        return None
    p, q, b = affine(word)
    s, t = (16, 11) if word.endswith('L') else (2, 1)
    modulus = s*q
    return ((t*q-b)*pow(p, -1, modulus)) % modulus, modulus


def native_dp(limit, signed=False, source_weight=False):
    states = [[F(0), F(0)] for _ in range(limit+1)]
    states[0][0] = F(1)
    for j in range(limit+1):
        for state in (0, 1):
            for c, data in OPS.items():
                if state == 1 and c == 'B':
                    continue
                new = j+data[3]
                if new <= limit:
                    weight = sign(c) if signed else 1
                    if source_weight:
                        weight = F(weight, data[1])
                    states[new][c == 'L'] += states[j][state]*weight
    return [a+b*(F(1, 8) if source_weight else 1) for a, b in states]


def rational_series(numerator, denominator, limit):
    out = []
    for j in range(limit+1):
        value = numerator[j] if j < len(numerator) else 0
        value -= sum(denominator[k]*out[j-k]
                     for k in range(1, min(j, len(denominator)-1)+1))
        out.append(F(value, denominator[0]))
    return out


def eta_product(limit):
    # Direct multiplication, including the initial q shift.
    out = [1]+[0]*(limit-1)
    for n in range(1, limit):
        for _ in range(1+(n % 11 == 0)):
            for j in range(limit-1, n-1, -1):
                out[j] -= 2*out[j-n]
                if j >= 2*n:
                    out[j] += out[j-2*n]
    return [0]+out


def eta_logarithm(limit):
    # Independent divisor convolution / logarithmic derivative.
    divisor = [0]*limit
    for d in range(1, limit):
        for j in range(d, limit, d):
            divisor[j] += d*(2+2*(d % 11 == 0))
    out = [1]
    for j in range(1, limit):
        numerator = -sum(divisor[k]*out[j-k] for k in range(1, j+1))
        check(numerator % j == 0, 'Euler-logarithm integrality')
        out.append(numerator//j)
    return [0]+out


def matmul(a, b, modulus=None):
    size = len(a)
    c = tuple(tuple(sum(a[i][k]*b[k][j] for k in range(size))
                    for j in range(size)) for i in range(size))
    return tuple(tuple(v % modulus for v in row) for row in c) if modulus else c


def matpow(a, n, modulus=None):
    out = tuple(tuple(int(i == j) for j in range(len(a))) for i in range(len(a)))
    for _ in range(n):
        out = matmul(out, a, modulus)
    return out


def kappa(parent):
    # Independent exhaustive discrete logarithm at the one fixed small modulus.
    power = 1
    for k in range(1, 730):
        power = power*4 % 2187
        if power*(3*parent+1) % 2187 == 334:
            return k
    raise ValueError('no parent phase')


def child_mod(parent, branch, modulus):
    k = kappa(parent)+729*branch
    value = (1024*pow(4, k, 2187*modulus)*(3*parent+1)-3031) % (2187*modulus)
    check(value % 2187 == 0, 'tree division precision')
    return value//2187


def branch_digits(parent, target, depth):
    t, period, rows = 0, 1, []
    for _ in range(depth):
        new_modulus = 3*period
        old = child_mod(parent, t, new_modulus)
        check((target-old) % period == 0, 'tree compatible prefix')
        digit = ((target-old)//period) % 3
        t += digit*period
        period = new_modulus
        check(child_mod(parent, t, period) == target % period, 'tree exact residue')
        rows.append((period, t, digit))
    return rows


def child(parent, branch):
    return (1024*4**(kappa(parent)+729*branch)*(3*parent+1)-3031)//2187


def height_bound(parent, target):
    """Largest t>=0 with child(parent,t)<=target, or -1 if none."""
    if child(parent, 0) > target:
        return -1
    low, high = 0, 1
    while child(parent, high) <= target:
        low, high = high, 2*high
    while high-low > 1:
        middle = (low+high)//2
        if child(parent, middle) <= target:
            low = middle
        else:
            high = middle
    return low


def bounded_branch_decode(parent, target):
    """Complete fixed-parent membership: ternary prefix plus real height."""
    bound = height_bound(parent, target)
    if bound < 0:
        return None, bound, 0
    period, depth = 1, 0
    while period <= bound:
        period *= 3
        depth += 1
    candidate = branch_digits(parent, target, depth)[-1][1] if depth else 0
    if candidate > bound or child(parent, candidate) != target:
        return None, bound, depth
    return candidate, bound, depth


def valuation(n, prime):
    if n == 0:
        raise ValueError('finite valuation expected')
    n, out = abs(n), 0
    while n % prime == 0:
        n //= prime
        out += 1
    return out


def main():
    native = list(words(10, 'HGABL', native=True))
    unsigned, signed = Counter(), Counter()
    for w, grade in native:
        v = expand(w)
        check(compress(v) == w, 'native-shadow inverse')
        check('LB' not in w and 'H' not in v, 'native-shadow alphabets')
        j, a, credit = stats(w)
        jj, aa, cc = stats(v)
        check(j == jj == grade, 'ternary depth preserved')
        check(a == aa-w.count('H'), 'binary fusion defect')
        check(credit == cc-3*w.count('H'), 'credit fusion defect')
        check(sign(w) == (-1)**len(v), 'eta sign preserved')
        unsigned[grade] += 1
        signed[grade] += sign(w)
    shadows = list(words(10, 'GLAB'))
    check(len(native) == len(shadows), 'equal enumeration sizes')
    for v, _ in shadows:
        check(expand(compress(v)) == v, 'shadow-native inverse')

    small = [w for w, _ in words(4, 'HGABL', native=True)]
    for u, v in product(small, repeat=2):
        expected = u[:-1]+'H'+v[1:] if u.endswith('L') and v.startswith('B') else u+v
        check(fuse(u, v) == expected, 'boundary-only fusion')
    tiny = [w for w, _ in words(2, 'HGABL', native=True)]
    for u, v, w in product(tiny, repeat=3):
        check(fuse(fuse(u, v), w) == fuse(u, fuse(v, w)), 'fusion associativity')

    u = native_dp(120)
    s = native_dp(120, signed=True)
    source = native_dp(120, signed=True, source_weight=True)
    for j in range(121):
        explicit = sum(comb(j-k, k)*2**(j-k) for k in range(j//2+1))
        signed_explicit = sum(comb(j-k, k)*(-2)**(j-k) for k in range(j//2+1))
        check(u[j] == explicit, 'two-colour tiling count')
        check(s[j] == signed_explicit, 'signed tiling count')
        if j <= 10:
            check(u[j] == unsigned[j] and s[j] == signed[j], 'DP versus literal words')
    for actual, expected in (
        (u, rational_series([1], [1, -2, -2], 120)),
        (s, rational_series([1], [1, 2, 2], 120)),
        (source, rational_series([1, F(7, 128)],
                                [1, F(3, 16), F(1, 64), -F(1, 2048)], 120))):
        check(actual == expected, 'generating function')

    p, q, b = affine('LB')
    ph, qh, bh = affine('H')
    check(F(ph, qh) == 2*F(p, q) and F(bh, qh) == 2*F(b, q)-F(1, 4),
          'fusion affine correction')
    check(native_port('LB') is None and native_at('H', 155), 'fusion semantic hostile')

    # Independent literal guard census, every odd input modulo a common period.
    grade3 = [w for w, j in native if j == 3]
    census_modulus = max(native_port(w)[1] for w in grade3)
    literal_sum = 0
    legal_at = {}
    for w in grade3:
        residue, period = native_port(w)
        count = 0
        for n in range(1, census_modulus, 2):
            legal = native_at(w, n)
            check(legal == (n % period == residue), 'literal native guard census')
            count += legal
        check(F(count, census_modulus//2) == F(2, period), 'source cylinder density')
        literal_sum += sign(w)*count
    for n in (155, 219):
        legal_at[n] = [w for w in grade3 if native_at(w, n)]
    check(sum(sign(w) for w in legal_at[155]) == 1, 'fixed-source signed hostile')
    check(F(literal_sum, census_modulus//2) == source[3] == F(27, 32768),
          'source-resolved cancellation defect')

    coeff = eta_product(1024)
    check(coeff == eta_logarithm(1024), 'independent eta expansion')
    eta_powers = {}
    for p in (2, 3, 5, 7):
        hp, last, pp, rows = 1, coeff[p], p, [1, coeff[p]]
        while pp*p <= 1024:
            hp, last = last, coeff[p]*last-p*hp
            pp *= p
            rows.append(last)
            check(coeff[pp] == last, 'eta good-prime recurrence')
        eta_powers[p] = rows
    for j in range(11):
        check(coeff[2**j] == s[j], 'eta controller prime2 bridge')
    check(coeff[11] == coeff[121] == 1 and coeff[22] == -2, 'eta level11 controls')
    companion = ((-2, -2), (1, 0))
    golden = ((0, 1), (1, 1))
    swap = ((0, 1), (1, 0))
    check(matpow(companion, 4) == ((-4, 0), (0, -4)), 'integral signed clock')
    check(matpow(companion, 1, 3) == matmul(matmul(swap, golden, 3), swap, 3),
          'golden conjugacy modulo3')
    check(matpow(companion, 8, 3) == ((1, 0), (0, 1)), 'golden order upper bound')
    for i in range(1, 8):
        check(matpow(companion, i, 3) != ((1, 0), (0, 1)), 'golden exact order')

    source_clock = ((-3, -4, 2), (1, 0, 0), (0, 1, 0))
    basis = ((1, 1, 0), (1, 1, 2), (1, 2, 2))
    split = ((1, 0, 0), (0, 0, 2), (0, 2, 2))  # 1 plus minus golden M
    check(matmul(source_clock, basis, 3) == matmul(basis, split, 3),
          'source-clock fixed-line golden-plane intertwiner')
    rvalues = []
    for j in range(121):
        r = 8*16**j*source[j]
        check(r.denominator == 1, 'integral source clock')
        rvalues.append(int(r))
        check((r-2*(-1)**j*s[j]) % 3 == 0, 'source eta congruence modulo3')
    check((rvalues[1]-2*(-1)*s[1]) % 9 != 0, 'modulo9 congruence hostile')
    identity3 = ((1, 0, 0), (0, 1, 0), (0, 0, 1))
    source_periods = []
    for k in range(1, 7):
        modulus = 3**k
        state = identity3
        for period in range(1, 8*3**(k-1)+1):
            state = matmul(state, source_clock, modulus)
            if state == identity3:
                break
        check(period == 8*3**(k-1) and state == identity3, 'exact lifted source-clock order')
        source_periods.append((modulus, period))
    source_eighth = matpow(source_clock, 8)
    deviations = [source_eighth[i][j]-identity3[i][j] for i in range(3) for j in range(3)]
    check(min(valuation(x, 3) for x in deviations if x) == 1, 'matrix lifting valuation')

    # The translation alphabet forms disjoint 3-adic balls.
    residues = [(b*pow(q, -1, 9)) % 9 for p, q, b, *_ in OPS.values()]
    check(len(set(residues)) == 5 and 0 not in residues, 'disjoint ternary tags')
    kraft = sum(F(1, data[0]) for data in OPS.values())
    check(kraft == F(181, 729), 'formal ternary cylinder mass')
    for w, j in words(6, 'HGABL', native=True):
        p, q, b = affine(w)
        check(p == 3**(2*j), 'ternary grade is valuation')
        check(F(p, q)*q*F(1, p) == 1, 'two-prime slope product')
    for k in range(1, 15):
        p, q, b = affine('L'*k)
        check(F(b, q) == -F(3, 7)*(1-F(9, 16)**k), 'infinite-address hostile')

    parent = 3
    first = (1024*4**kappa(parent)*(3*parent+1)-3031)//2187
    check(first.bit_length() == 1353 and first > 155, 'height escaping tree control')
    targets = (0, 1, 155, 223, 233)
    tree_rows = {n: branch_digits(parent, n, 24) for n in targets}
    for depth in range(1, 7):
        period = 3**depth
        images = [child_mod(parent, t, period) for t in range(period)]
        check(len(set(images)) == period, 'complete ternary branch level')
        for target in range(period):
            rows = branch_digits(parent, target, depth)
            check(images[rows[-1][1]] == target, 'ternary inverse table')
    for t, r in product(range(8), repeat=2):
        if t != r:
            modulus = 3**8
            diff = (child_mod(parent, t, modulus)-child_mod(parent, r, modulus)) % modulus
            check(valuation(diff, 3) == valuation(t-r, 3), 'tree isometry control')
    # A genuine branch has finitely supported address digits; precision 12 sees it.
    for t in (0, 1, 7, 729, 12003):
        target = child_mod(parent, t, 3**12)
        check(branch_digits(parent, target, 12)[-1][1] == t, 'positive finite-address control')

    decoded_examples = []
    for t in (0, 1, 2, 7, 19, 40):
        n = child(parent, t)
        decoded, bound, depth = bounded_branch_decode(parent, n)
        check(decoded == bound == t, 'height-completed positive decoder')
        check(child(parent, bound) <= n < child(parent, bound+1), 'exact branch height bound')
        rejected, gap_bound, gap_depth = bounded_branch_decode(parent, n+2048)
        check(rejected is None and gap_bound == t, 'same-guard finite rejection')
        decoded_examples.append((t, depth, gap_depth, n.bit_length()))
    for n in (155, 223, 233):
        check(bounded_branch_decode(parent, n) == (None, -1, 0), 'small-source finite rejection')

    # An explicit radial law on a p-ary tree, with Fourier value rho^conductor.
    # This is a controlled model, not the Syracuse law.
    moment_rows = []
    rho = F(1, 2)
    for p in (2, 3, 5, 7, 11):
        for depth in range(1, 6):
            for r in (1, 2):
                total = F(0)
                for frequency in range(p**depth):
                    conductor = 0 if frequency == 0 else depth-valuation(frequency, p)
                    total += rho**(2*r*conductor)
                total /= p**depth
                expected = F(1, p**depth)+F(p-1, p)*sum(
                    F(1, p**j)*rho**(2*r*(depth-j)) for j in range(depth))
                check(total == expected, 'prime-tree moment stratification')
                if depth == 5:
                    moment_rows.append((p, r, str(total)))

    result = {
        'status': 'Exact scoped identities and finite controls; Collatz and fixed-seed H1 OPEN',
        'universe': {'native_grade_limit': 10, 'native_words': len(native),
                     'shadow_words': len(shadows), 'DP_grade_limit': 120,
                     'eta_coefficient_limit': 1024, 'grade3_odd_source_modulus': census_modulus},
        'unsigned_counts': [int(u[j]) for j in range(13)],
        'signed_counts_equal_eta_a_2power': [int(s[j]) for j in range(13)],
        'source_averaged_signed_counts': [str(source[j]) for j in range(7)],
        'integral_source_clock_initial_values': rvalues[:13],
        'integral_source_clock_eighth_power': source_eighth,
        'source_clock_modulus_and_exact_order': source_periods,
        'grade3_literal_signed_sum': literal_sum,
        'grade3_legal_words_at_fixed_sources': legal_at,
        'eta_first_22': coeff[1:23],
        'eta_prime_power_controls': eta_powers,
        'translation_tags_mod9': residues,
        'translation_first_level_mass': str(kraft),
        'parent3_kappa': kappa(parent),
        'parent3_first_child_binary_digits': first.bit_length(),
        'tree_target155_first12_prefixes': tree_rows[155][:12],
        'tree_target155_depth24_prefix': tree_rows[155][-1],
        'height_decoder_branch_depth_gapdepth_sourcebits': decoded_examples,
        'radial_model_rho_one_half_moments_depth5': moment_rows,
        'checks': dict(sorted(CHECKS.items())),
        'total_checks': sum(CHECKS.values()),
    }
    print(json.dumps(result, indent=2, sort_keys=True))


if __name__ == '__main__':
    main()
