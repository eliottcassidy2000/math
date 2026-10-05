"""Exact valuation-code and native-guard information controls.

No probabilistic convergence inference and no orbit-stopping oracle.
Integer/Fraction controls survive -O; printed entropy values are numerical.
"""
from fractions import Fraction as F
from itertools import combinations, product
from math import comb, log2
import translation_phase_decoder_20261005 as controller

SIX = dict(controller.LETTERS, J=(5, 8, 147, 799, 1024))


def need(ok, message):
    if not ok:
        raise ValueError(message)


def v2(n):
    need(type(n) is int and n != 0, 'nonzero integer valuation')
    n = abs(n)
    return (n & -n).bit_length()-1


def step(n):
    need(type(n) is int and n % 2, 'odd integer source')
    need(3*n+1 != 0, 'finite valuation')
    a = v2(3*n+1)
    return (3*n+1)//(1 << a), a


def bits(n, length):
    need(type(n) is int and n > 0 and n % 2, 'positive odd source')
    need(type(length) is int and length >= 0, 'nonnegative bit budget')
    out = ''
    while len(out) < length:
        n, a = step(n)
        out += '0'*(a-1)+'1'
    return out[:length]


def shortcut_bits(n, length):
    out = ''
    for _ in range(length):
        n = (3*n+1)//2 if n % 2 else n//2
        out += str(n % 2)
    return out


def compositions(total):
    if total == 0:
        yield ()
    else:
        for a in range(1, total+1):
            for tail in compositions(total-a):
                yield (a,)+tail


def carrier(word):
    p, q, b = 1, 1, 0
    for a in word:
        need(type(a) is int and a >= 1, 'positive valuation')
        p, q, b = 3*p, q*(1 << a), 3*b+q
    return p, q, b


def cylinder(word):
    p, q, b = carrier(word)
    return ((q-b)*pow(p, -1, 2*q)) % (2*q), 2*q


def replay_word(n, word):
    for expected in word:
        n, actual = step(n)
        need(actual == expected, 'literal exact valuation')
    return n


def intersect_dyadic(left, right):
    if left is None or right is None:
        return None
    a, m = left
    b, n = right
    if (a-b) % min(m, n):
        return None
    return left if m >= n else right


def pulled_guard(word, letters=controller.LETTERS):
    """Independent intersection of every local native input condition."""
    answer = (1, 2)
    p, q, b = 1, 1, 0
    for symbol in word:
        r, a, carry, residue, modulus = letters[symbol]
        local = ((q*residue-b)*pow(p, -1, q*modulus) % (q*modulus), q*modulus)
        answer = intersect_dyadic(answer, local)
        p, q, b = 3**r*p, (1 << a)*q, 3**r*b+carry*q
    return answer


def native_prefix(word):
    """Exact unary/parity source prefix from inherited controller interfaces."""
    need(type(word) is str and all(c in SIX for c in word), 'six-letter word')
    prefixes = {'G': '101', 'A': '1011001', 'B': '1101001',
                'H': '1011110100'}
    out = ''
    for symbol in reversed(word):
        if symbol == 'J':
            if not out:
                out = '111101001'
            elif out.startswith('1'):
                out = '111101001'+out[1:]
            else:
                return None
        elif symbol != 'L':
            out = prefixes[symbol]+out
        elif not out:
            out = '1011100'
        elif out.startswith('101'):
            out = '1011100'+out[3:]
        else:
            return None
    return out


def union_mass(cells):
    """Normalized odd-Haar mass after removing nested duplicate cylinders."""
    retained = []
    for residue, modulus in sorted(set(cells), key=lambda c: c[1]):
        need(modulus >= 2 and modulus & (modulus-1) == 0 and residue % 2,
             'odd dyadic cylinder')
        if not any((residue-r) % m == 0 for r, m in retained):
            retained.append((residue, modulus))
    return sum((F(2, modulus) for _, modulus in retained), F(0)), tuple(retained)


def entropy(theta):
    if theta in (0, 1):
        return 0.0
    return -theta*log2(theta)-(1-theta)*log2(1-theta)


def local_dimension(z, mean):
    theta = 1/mean
    return -theta*log2(1-z)-(1-theta)*log2(z)


def main():
    prefix_count = 0
    for length in range(1, 11):
        values = range(1, 1 << (length+1), 2)
        codes = [bits(n, length) for n in values]
        need(len(set(codes)) == 1 << length, 'complete source-bit bijection')
        for n, code in zip(values, codes):
            need(code == shortcut_bits(n, length), 'unary equals dropped-initial parity')
            prefix_count += 1
    pairs = 0
    for x, y in combinations(range(1, 512, 2), 2):
        a, b = bits(x, 9), bits(y, 9)
        common = next(i for i in range(9) if a[i] != b[i])
        need(common == v2(x-y)-1, 'normalized odd metric isometry')
        pairs += 1
    print('complete binary-prefix bijections depths1..10:', prefix_count,
          '; independent isometry pairs', pairs)

    words = probability_checks = 0
    for total in range(1, 12):
        by_rank = {}
        for word in compositions(total):
            r, modulus = cylinder(word)
            need(r % 2 == 1 and modulus == 1 << (total+1), 'exact odd cylinder')
            for lift in (0, 1, 3):
                n = r+lift*modulus
                target = replay_word(n, word)
                p, q, b = carrier(word)
                need(target == (p*n+b)//q, 'ordered affine carry retained')
            code = ''.join('0'*(a-1)+'1' for a in word)
            need(bits(r, total) == code, 'cylinder selects the encoded word')
            by_rank[len(word)] = by_rank.get(len(word), 0)+1
            for z in (F(1, 4), F(1, 2), F(2, 3)):
                renewal = F(1)
                for a in word:
                    renewal *= (1-z)*z**(a-1)
                binary = (1-z)**len(word)*z**(total-len(word))
                need(renewal == binary, 'geometric renewal is Bernoulli in this code')
                if z == F(1, 2):
                    need(binary == F(1, 1 << total), 'odd-Haar cylinder mass')
                probability_checks += 1
            words += 1
        need(by_rank == {r: comb(total-1, r-1) for r in range(1, total+1)},
             'composition count by exact mean')
    print('all positive valuation words of cost1..11:', words,
          '; exact geometric probability checks', probability_checks)

    guard_words = valid_words = 0
    for length in range(1, 5):
        for letters in product(controller.LETTERS, repeat=length):
            word = ''.join(letters)
            guard = pulled_guard(word)
            need(guard == controller.native_guard(controller.encode(word)),
                 'independent local intersections equal complete native guard')
            need((guard is None) == ('LB' in word), 'native forbidden-pair boundary')
            need((native_prefix(word) is None) == (guard is None),
                 'prefix interfaces independently find the forbidden pair')
            if guard is not None:
                p, q, _ = controller.carrier(word)
                overhead = 3 if word.endswith('L') else 0
                extra_ones = 2 if word.endswith('L') else 0
                rank = sum(controller.LETTERS[c][0] for c in word)
                need(F(2, guard[1]) == F(1, q*(1 << overhead)),
                     'native information equals A plus terminal-L overhead')
                code = native_prefix(word)
                need(code == bits(guard[0], len(code)), 'actual source selects the native code')
                need(len(code) == q.bit_length()-1+overhead and
                     code.count('1') == rank+extra_ones,
                     'terminal interface retains exact bits and ones')
                likelihood = F(3**code.count('1'), 1 << len(code))
                need(likelihood == F(p, q)*(F(9, 8) if overhead else 1),
                     'native tilted likelihood equals coefficient with interface factor')
                valid_words += 1
            guard_words += 1
    cells = [pulled_guard(c) for c in controller.LETTERS]
    naive = sum((F(2, m) for _, m in cells), F(0))
    combined, retained = union_mass(cells)
    funder_mass, _ = union_mass([pulled_guard(c) for c in 'HABL'])
    need((naive, combined, funder_mass) == (F(153, 1024), F(17, 128), F(25, 1024)),
         'overlaps are removed and native is not paid')
    need(sum(any(n % m == r for r, m in cells) for n in range(1, 2048, 2)) == 136,
         'literal independent primitive-guard union')
    print('all controller words lengths1..4:', guard_words, '; native', valid_words)
    print('primitive native masses: naive sum', naive, '; true union', combined,
          '; first-letter funder union HABL', funder_mass, '; retained union cells', retained)

    six_count = six_valid = 0
    for length in range(1, 5):
        for symbols in product(SIX, repeat=length):
            word = ''.join(symbols)
            guard, code = pulled_guard(word, SIX), native_prefix(word)
            need((guard is None) == ('LB' in word or 'LJ' in word),
                 'incoming six-letter exact forbidden pairs')
            need((guard is None) == (code is None), 'six-letter prefix and congruence agree')
            if guard is not None:
                p = q = 1
                b = rank = 0
                for symbol in word:
                    r, a, c, _, _ = SIX[symbol]
                    p, q, b = 3**r*p, (1 << a)*q, 3**r*b+c*q
                    rank += r
                extra_bits, extra_ones = (3, 2) if word.endswith('L') else (1, 1) if word.endswith('J') else (0, 0)
                need(len(code) == q.bit_length()-1+extra_bits and
                     code.count('1') == rank+extra_ones, 'six-letter terminal interface counts')
                target, factor = (11, 16) if word.endswith('L') else (3, 4) if word.endswith('J') else (1, 2)
                expected = ((target*q-b)*pow(p, -1, factor*q) % (factor*q), factor*q)
                need(guard == expected and bits(guard[0], len(code)) == code,
                     'six-letter formula and actual source code')
                correction = F(9, 8) if word.endswith('L') else F(3, 2) if word.endswith('J') else F(1)
                need(F(3**code.count('1'), 1 << len(code)) == F(p, q)*correction,
                     'six-letter likelihood with terminal interface')
                six_valid += 1
            six_count += 1
    for t in range(32):
        n, child = 799+1024*t, 759+972*t
        need((243*n+147)//256 == child < n, 'incoming native J affine receipt')
        need(replay_word(n, (1, 1, 1, 1, 2, 3)) == replay_word(child, (1,)),
             'independent J actual common-future words')
    print('incoming six-letter J extension: words lengths1..4', six_count,
          '; native', six_valid, '; 32 literal J joins; terminal-J likelihood factor3/2')

    need(carrier((1, 2)) == (9, 8, 5) and carrier((2, 1)) == (9, 8, 7),
         'same mass and multiplier lose carry/order')
    need((9*27-3)//16 == 15 and 27 % 256 != 219,
         'formula-only L source is not a native receipt')
    need(native_prefix('L') == native_prefix('LG') == '1011100',
         'same guard does not decode controller history')
    need(controller.carrier('L') == (9, 16, -3) and
         controller.carrier('LG') == (81, 128, 53),
         'same guard retains different affine maps')
    need(controller.operate(219, 'L') == 123 and
         controller.operate(123, 'G') == 139,
         'same-guard programs give different paid endpoints')
    for r in range(1, 17):
        p, q, b = carrier((1,)*r)
        need(F(p*(-1)+b, q) == -1 and F(p, q) > 1,
             'negative fixed point cancels the real expanding coefficient')
        p, q, b = carrier((2,)*r)
        need(F(p+b, q) == 1 and F(p, q) < 1,
             'positive root cancels contracting coefficient by its carry')
    need(bits(1, 20) == '01'*10, 'root continuation is computable with mean2')
    for length in range(1, 20):
        weighted = sum((F(1, 1 << a)*F(3, 1 << a) for a in range(1, length+1)), F(0))
        need(weighted+F(1, 4**length) == 1, 'exact likelihood martingale tail')
    print('tilted likelihood: coefficient3^r/2^A at valuation boundaries; native terminal-L factor9/8; 19 exact martingale-tail controls')
    print('hostiles: carry5 versus7; L/LG same guard with endpoints123/139; formula L(27)=15 lacks native guard; -1 fixed despite expanding slope; root code(01)^infinity')

    rows = []
    for mean in (1, 1.5, log2(3), 2, 3, 4):
        rows.append((mean, entropy(1/mean), log2(3)-mean,
                     local_dimension(0.5, mean), local_dimension(0.25, mean)))
    print('NUMERICAL (mean,frequency-set dimension,coefficient drift,local dim z=.5,local dim z=.25)')
    for row in rows:
        print(tuple(round(x, 12) for x in row))
    print('PROVED scopes: exact code/guard information; fixed arithmetic inputs remain computable; no probabilistic Collatz conclusion')


if __name__ == '__main__':
    main()
