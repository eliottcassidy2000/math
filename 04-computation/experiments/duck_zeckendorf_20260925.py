"""Exact Fibonacci charge, three nonzero XOR colors, and diagonal carries.

Source note: 05-knowledge/results/duck_zeckendorf_20260925.md.
This explicitly defines the color convention; no unlocated historical theorem
or conjectured Collatz/graceful equivalence is assumed.
"""

from functools import lru_cache
from itertools import product
from math import gcd, isqrt


def check(test, message):
    if not test:
        raise RuntimeError(message)


FIB = [1, 2]
CHARGE = [(1, 0), (0, 1)]
for _ in range(48):
    FIB.append(FIB[-1]+FIB[-2])
    CHARGE.append(tuple(CHARGE[-1][j]+CHARGE[-2][j] for j in (0, 1)))


def beta(n):
    """floor((n+1)/phi^2), using only integer square root."""
    check(type(n) is int and n >= 0, "nonnegative integer")
    m = n+1
    return (3*m-isqrt(5*m*m)-1)//2


def charge_formula(n):
    b = beta(n)
    return n-2*b, b


def greedy(n):
    bits = 0
    for i in range(len(FIB)-1, -1, -1):
        if FIB[i] <= n:
            n -= FIB[i]
            bits |= 1 << i
    check(n == 0 and not (bits & (bits >> 1)), "greedy normal form")
    return bits


def value(bits):
    return sum(FIB[i] for i in range(bits.bit_length()) if bits >> i & 1)


def charge(bits):
    return tuple(sum(CHARGE[i][j] for i in range(bits.bit_length()) if bits >> i & 1)
                 for j in (0, 1))


def color(n):
    a, b = charge_formula(n)
    return (a % 2, b % 2)


def successors(bits):
    for i in range(bits.bit_length()):
        if bits >> i & 7 == 3:
            yield bits-(3 << i)+(1 << (i+2))


@lru_cache(None)
def all_normal_forms(bits):
    children = list(successors(bits))
    if not children:
        return frozenset((bits,))
    return frozenset().union(*(all_normal_forms(child) for child in children))


def delta(a, b):
    return beta(a)+beta(b)-beta(a+b)


def xor(a, b):
    return a[0] ^ b[0], a[1] ^ b[1]


def matrix_color(c):
    return c[1], c[0] ^ c[1]


def color_digits(n, depth):
    palette = ((1, 0), (0, 1), (1, 1))
    out = []
    for _ in range(depth):
        out.append(palette[n % 3])
        n //= 3
    return out


def increment_colors(word):
    word = list(word)
    for i, old in enumerate(word):
        word[i] = matrix_color(old)
        if old != (1, 1):
            break
    return word


def main():
    checks = 0
    for n in range(100001):
        form = greedy(n)
        check(charge(form) == charge_formula(n), "independent greedy/floor charge")
        shifted = form << 1
        shift_n = value(shifted)
        check(shift_n == 2*n-beta(n), "exact Fibonacci shift")
        check(color(shift_n) == matrix_color(color(n)), "order3 color quotient")
        checks += 1
    print("CANONICAL FLOOR/SHIFT: all n=0..100000;", checks, "PASS")

    forms_by_value = {}
    edges = 0
    for bits in range(1 << 16):
        n = value(bits)
        check(charge(bits) == charge_formula(n), "all distinct decompositions preserve charge")
        forms_by_value.setdefault(n, 0)
        forms_by_value[n] += 1
        for child in successors(bits):
            check(value(child) == n and charge(child) == charge(bits), "binary carry preserves charge")
            check(child.bit_count()+1 == bits.bit_count(), "terminating rewrite measure")
            edges += 1
    print("DISTINCT REPRESENTATIONS: all 65536 subsets of first16 Fibonacci weights;",
          len(forms_by_value), "values;", edges, "rewrite edges PASS")
    confluent = 0
    for bits in range(1 << 12):
        check(all_normal_forms(bits) == frozenset((greedy(value(bits)),)), "all rewrite orders confluence")
        confluent += 1
    print("CONFLUENCE: all4096 subsets of first12 weights; all legal rewrite orders PASS")

    # Repeated atoms are a genuinely different universe.
    repeated = 0
    defects = set()
    for counts in product(range(4), repeat=7):
        n = sum(c*f for c, f in zip(counts, FIB))
        raw = tuple(sum(c*v[j] for c, v in zip(counts, CHARGE)) for j in (0, 1))
        defect = raw[1]-beta(n)
        check(raw == (n-2*beta(n)-2*defect, beta(n)+defect), "one dimensional carry defect")
        repaired = raw[0] % 2, (raw[1]-defect) % 2
        check(repaired == color(n), "one-bit color correction")
        repeated += 1
        defects.add(defect)
    check(xor(color(1), color(1)) != color(2), "repeated-atom minimal hostile")
    print("REPEATED REPRESENTATIONS: 4^7=", repeated, "; all one-bit repairs PASS; defect range",
          min(defects), max(defects))

    counts = {-1: 0, 0: 0, 1: 0}
    disjoint = 0
    first_split = None
    for n in range(2, 513):
        naive_colors = set()
        for a in range(1, (n-1)//2+1):
            b = n-a
            d = delta(a, b)
            check(d in counts, "two-input carry range")
            counts[d] += 1
            raw = xor(color(a), color(b))
            check(xor(raw, (0, d % 2)) == color(n), "diagonal correction")
            check(charge_formula(a)[0]+charge_formula(b)[0]-charge_formula(n)[0] == -2*d,
                  "integer charge correction")
            if not (greedy(a) & greedy(b)):
                check(d == 0, "disjoint supports need no repeated-atom carry")
                disjoint += 1
            naive_colors.add(raw)
        if len(naive_colors) > 1 and first_split is None:
            first_split = n
    check(first_split == 5, "first strict diagonal split")
    print("STRICT SUM DIAGONALS n<=512: carry counts", counts,
          "; disjoint safe pairs", disjoint, "; first split", first_split)
    print("DIAGONAL5: 1+4 has", xor(color(1), color(4)),
          ";2+3 has", xor(color(2), color(3)), ";canonical5", color(5))
    for a, b, c in product(range(65), repeat=3):
        check(delta(a, b)+delta(a+b, c) == delta(b, c)+delta(a, b+c), "carry cocycle")
    print("ASSOCIATIVE CARRY: all65^3 triples a,b,c=0..64 PASS")

    # Three cyclic colors are not additive residue classes modulo3.
    periodic_z3 = [w for w in product(range(3), repeat=3)
                   if all((w[i]+w[(i+1)%3]-w[(i+2)%3]) % 3 == 0 for i in range(3))]
    check(periodic_z3 == [(0, 0, 0)], "no nonzero period3 additive Z3 Fibonacci coloring")
    palette = ((1, 0), (0, 1), (1, 1))
    check(tuple(matrix_color(c) for c in palette) == palette[1:]+palette[:1], "three nonzero colors cycle")
    check(matrix_color((0, 0)) == (0, 0) and color(6) == (0, 0), "neutral fourth state is real")
    for depth in range(1, 8):
        modulus = 3**depth
        for n in range(modulus):
            check(increment_colors(color_digits(n, depth)) == color_digits((n+1)%modulus, depth),
                  "colored ternary odometer")
        check(color(1) != color(1+modulus), "Fibonacci color does not factor through ternary residue")
    print("TERNARY COLOR-WORD CODEC: depths1..7, all3279 residues PASS; scalar color fails every tested prefix quotient")
    clock_cases = 0
    for k0 in range(2, 6):
        scale = 2**(k0-1)
        for depth in range(1, 6):
            modulus = 3**depth
            phase = {(scale*(4**j-1)//3) % modulus: j for j in range(modulus)}
            check(len(phase) == modulus, "inherited residue-clock bijection")
            for height, j in phase.items():
                next_j = phase[(4*height+scale) % modulus]
                check(increment_colors(color_digits(j, depth)) == color_digits(next_j, depth),
                      "height-clock phase/color commuting square")
                clock_cases += 1
    print("RESIDUE CLOCK TRANSPORT: k0=2..5,depth1..5;", clock_cases,
          "height/phase/color commuting squares PASS")
    print("GROUP CONTROL: cyclic R,G,B is an order3 action on F2^2; the only period3 additive Z3 atom coloring is zero")

    for i in range(3, 50):
        m, n = CHARGE[i]
        raw = n*n-m*m, 2*m*n, n*n+m*m
        content = gcd(gcd(raw[0], raw[1]), raw[2])
        check(n*n-m*n-m*m == (-1)**(i-1), "inherited golden Euclid slice")
        check((content == 2) == (tuple(v % 2 for v in CHARGE[i]) == (1, 1)),
              "third atom color is odd-odd primitive normalization")
    print("THM3339 INTERFACE: atom charges i3..49; golden norm and color-B content2 seam PASS")
    raw = xor(color(1), color(4))
    rotated_correct = xor(matrix_color(raw), matrix_color((0, delta(1, 4) % 2)))
    rotated_wrong = xor(matrix_color(raw), (0, delta(1, 4) % 2))
    check(rotated_correct == matrix_color(color(5)) and rotated_wrong != rotated_correct,
          "carry channel must rotate with its color gauge")
    check((2-beta(1))+(8-beta(4)) != 10-beta(5), "Fibonacci shift is not additive")
    print("GAUGE HOSTILE: diagonal1+4 requires rotating carry-G toB; fixed-channel correction fails")

    examples = (((1, 0), 1, 5), ((0, 1), 4, 2), ((1, 1), 9, 3), ((0, 0), 16, 6))
    for c, square, nonsquare in examples:
        check(color(square) == color(nonsquare) == c and isqrt(square)**2 == square
              and isqrt(nonsquare)**2 != nonsquare, "color cannot decide square target")
    print("SQUARE FILTER HOSTILES: same colors at square/nonsquare pairs1/5,4/2,9/3,16/6")
    print("ALL EXACT CHECKS PASS; no square-sum/graceful/Collatz equivalence follows")


if __name__ == "__main__":
    main()
