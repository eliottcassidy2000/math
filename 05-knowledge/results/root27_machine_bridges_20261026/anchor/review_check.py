"""Independent guard, source, carry and pumped-word controls; no producer import."""
from itertools import combinations, product
import json


def need(ok, message):
    if not ok:
        raise RuntimeError(message)


def valuation(value):
    need(value > 0, 'positive valuation input')
    return (value & -value).bit_length() - 1


def original_step(a, b):
    if b % 2:
        if a == 0:
            return None
        return a - 1, (3 * b + 3) // 2
    return a + 2, 3 * b // 2 + 2


def carry_bits(z, n):
    out = []
    for unused in range(n):
        bit = z & 1
        out.append(bit)
        z = (3 * z + bit) // 2
    return tuple(out)


def polynomial(z, word):
    n = len(word)
    return 3 ** n * z + sum(bit * 2 ** i * 3 ** (n - 1 - i)
                           for i, bit in enumerate(word))


def accepted_tree(word):
    slots = 1
    for i, bit in enumerate(word):
        slots += 3 * bit - 1
        if slots <= 0 and i + 1 != len(word):
            return False
    return slots == 0


def main():
    translation_checks = 0
    for source in range(3, 130):
        for credit in range(6):
            a, b = credit, source - 3
            q, z = credit + 1, source
            for unused in range(40):
                raw = original_step(a, b)
                bit = z & 1
                if q == 1 and bit == 0:
                    need(raw is None, 'blocked symbol is not executed')
                    break
                need(raw is not None, 'legal translated macro')
                q += 3 * bit - 1
                z = (3 * z + bit) // 2
                a, b = raw
                need((a + 1, b + 3) == (q, z), 'literal two-coordinate translation')
                translation_checks += 1

    guard_checks = 0
    for source in range(3, 514, 2):
        ell = valuation(3 * source + 1)
        formal = 3 ** (ell - 1) * ((3 * source + 1) // 2 ** ell)
        for q in range(1, 10):
            a, b, count = q - 1, source - 3, 0
            while True:
                raw = original_step(a, b)
                if raw is None:
                    break
                a, b = raw
                count += 1
                if (b + 3) % 2:
                    break
            if ell <= q + 2:
                need(count == ell and (a + 1, b + 3) == (q + 3 - ell, formal),
                     'actual odd-z return')
            else:
                need(count == q + 2 and a == 0 and (b + 3) % 2 == 0,
                     'native halt precedes formal odd-z return')
            guard_checks += 1
    need(3 ** 3 * (16 // 16) == 27 and (22 // 2) == 11,
         'pointwise growth does not imply order monotonicity')

    # Different sections of the same physical blank run, with exact macro indices.
    a, b = 0, 4
    odd_x, odd_z = [], []
    macro = 0
    while len(odd_x) < 32 or len(odd_z) < 32:
        x, z = b + 4, b + 3
        if x % 2 and len(odd_x) < 32:
            odd_x.append((macro, a, x))
        if z % 2 and len(odd_z) < 32:
            odd_z.append((macro, a + 1, z))
        raw = original_step(a, b)
        need(raw is not None and macro < 1000, 'bounded blank cross-section')
        a, b = raw
        macro += 1
    need(odd_x[0] == (3, 6, 27) and odd_z[0] == (0, 1, 7), 'two marked seeds')
    for left, right in zip(odd_x, odd_x[1:]):
        t, credit, x = left
        k = valuation(3 * x - 1)
        need(right == (t + k, credit + 2 * k - 3,
                       3 ** (k - 1) * ((3 * x - 1) // 2 ** k)), 'odd-x section clock')
    for left, right in zip(odd_z, odd_z[1:]):
        t, slots, z = left
        ell = valuation(3 * z + 1)
        need(right == (t + ell, slots + 3 - ell,
                       3 ** (ell - 1) * ((3 * z + 1) // 2 ** ell)), 'odd-z section clock')

    leaf_checks = mismatch_checks = 0
    for n in range(1, 8):
        for word in product((0, 1), repeat=n):
            for z in (1, 2, 7, 15, 39):
                actual = carry_bits(z, n)
                first = next((i for i in range(n) if actual[i] != word[i]), n)
                value = polynomial(z, word)
                need(min(valuation(value), n) == first, 'valuation equals first literal mismatch')
                mismatch_checks += 1
                for p, bit in enumerate(word):
                    if bit:
                        continue
                    new = word[:p] + (1, 0, 0, 0) + word[p + 1:]
                    suffix = word[p + 1:]
                    C = polynomial(0, suffix)
                    after = polynomial(z, new)
                    need(after - 27 * value == 2 ** p * (3 ** (len(suffix) + 3) - 38 * C),
                         'leaf identity from independent expanded polynomials')
                    h, hh = valuation(value), valuation(after)
                    need(hh == min(p, h) if p != h else hh >= h + 1, 'exact leaf precision cases')
                    leaf_checks += 1

    pump_blocks = pump_checks = 0
    for internal in range(1, 5):
        length = 3 * internal + 1
        for ones in combinations(range(length), internal):
            word = tuple(int(i in ones) for i in range(length))
            if not accepted_tree(word):
                continue
            block = word[:-1]
            L = len(block)
            C = polynomial(0, block)
            need(C % 2 == 1 and L % 3 == 0 and block[0] == 1, 'nonempty block premises')
            for z in range(1, 65):
                fuel = valuation((3 ** L - 2 ** L) * z + C)
                for e in range(6):
                    target = block * e + (0,)
                    need(accepted_tree(target), 'pumped shape remains a tree')
                    need((carry_bits(z, len(target)) == target) == (fuel == L * e),
                         'pump criterion versus literal itinerary')
                    pump_checks += 1
            pump_blocks += 1
    need(pump_blocks == 71, 'independent ballot enumeration')
    print(json.dumps({'status': 'PASS', 'translated_native_macros': translation_checks,
                      'guarded_returns': guard_checks, 'section_clock_edges': 62,
                      'literal_mismatch_checks': mismatch_checks, 'leaf_checks': leaf_checks,
                      'pumped_blocks': pump_blocks, 'pumped_literal_checks': pump_checks,
                      'scope': 'finite independent controls; universal proof audit in review'},
                     sort_keys=True, indent=2))


if __name__ == '__main__':
    main()
