"""Exact controls for the Kuratowski/Collatz reframe; standard library only.

Universe: multipliers 3,5; signs +/-1; odd starts <=2047; every exponent
word of length 1..4 over {1,2,3,4}. No orbit-convergence assumption.
Direct valuation iteration is independent of the affine repetition formula.
Additional controls: finite cylinders, sign reflection, actual known cycles,
Berggren counterexample, arbitrarily long macros, and unbounded budget reset.
"""
from itertools import product
from fractions import Fraction
import hashlib
from pathlib import Path


def check(ok, context):
    if not ok:
        raise RuntimeError(context)


def v2(n):
    n = abs(n)
    if not n:
        raise ValueError("valuation of zero requires separate infinity case")
    return (n & -n).bit_length() - 1


def step(n, q=3, sign=1):
    z = q * n + sign
    k = v2(z)
    return z >> k, k


def block(word, q=3, sign=1):
    a, b, c = 1, 1, 0
    for k in word:
        a, b, c = q * a, b << k, q * c + sign * b
    return a, b, c


def direct_block(n, word, q=3, sign=1):
    for k in word:
        n, actual = step(n, q, sign)
        if actual != k:
            return None
    return n


def direct_repetitions(n, word, q=3, sign=1):
    initial, t = n, 0
    while True:
        nxt = direct_block(n, word, q, sign)
        if nxt is None:
            return t
        t += 1
        if nxt == initial:
            return None  # infinity, an actual fixed point of the block
        check(t < 100, ("unexpected control cap", n, word))
        n = nxt


def residue(word, q=3, sign=1):
    a, b, c = block(word, q, sign)
    return ((b - c) * pow(a, -1, 2 * b)) % (2 * b)


def address(s, t):
    letters = []
    while (s, t) != (3, 1):
        if s > 3 * t:
            s, t, letter = s - 2 * t, t, "A"
        elif s > 2 * t:
            s, t, letter = t, s - 2 * t, "B"
        else:
            s, t, letter = t, 2 * t - s, "C"
        letters.append(letter)
        check(s > t > 0, (s, t))
    return "".join(reversed(letters))


def main():
    cases, infinite, cylinder_cases = 0, 0, 0
    for q in (3, 5):
        for sign in (1, -1):
            for length in range(1, 5):
                for word in product(range(1, 5), repeat=length):
                    a, b, c = block(word, q, sign)
                    d, s = a - b, sum(word)
                    r = residue(word, q, sign)
                    check(0 < r < 2 * b and r % 2, ("residue", word))
                    for lift in (0, 1, 3):
                        n = r + 2 * b * lift
                        end = direct_block(n, word, q, sign)
                        check(end is not None and end * b == a * n + c,
                              ("cylinder", q, sign, word, n))
                        cylinder_cases += 1
                    opposite = residue(word, q, -sign)
                    check(r + opposite == 2 * b, ("reflection", word))
                    for n in range(1, 2048, 2):
                        e = d * n + c
                        predicted = None if e == 0 else (v2(e) - 1) // s
                        actual = direct_repetitions(n, word, q, sign)
                        check(predicted == actual,
                              ("repeat", q, sign, word, n, predicted, actual))
                        infinite += actual is None
                        cases += 1
    print("FINITE-EXACT repetition controls:", cases, "cases;", infinite,
          "fixed-block cases; zero mismatches")
    print("FINITE-EXACT cylinder lift controls:", cylinder_cases,
          "; all sign complements exact")

    orbit = [23]
    for _ in range(8):
        orbit.append(step(orbit[-1], sign=-1)[0])
    check(orbit == [23, 17, 25, 37, 55, 41, 61, 91, 17], orbit)
    check(address(23, 17) == "ABCC", "first address")
    check(address(91, 17) == "ABCCAA", "second address")
    print("REFUTED pasted orbit corollary:", orbit, "; ABCC < ABCCAA")

    for sign in (1, -1):
        for N in range(1, 65):
            n = 2 * 8**N - 5 * sign
            for j in range(N):
                n = direct_block(n, (1, 2), sign=sign)
                expected = 2 * 9**(j + 1) * 8**(N - j - 1) - 5 * sign
                check(n == expected, ("macro family", sign, N, j))
            if sign == 1:
                end = direct_block(n, (3,))
                check(end == (3 * 9**N - 7) // 4, ("exit", N))
                check((end > 2 * 8**N - 5) == (N >= 9), ("exit debt", N))
                check(v2(end + 1) == 1 + v2(N), ("next rise", N))
    print("FINITE-EXACT repeated (1,2) family: N=1..64, both signs;")
    print("plus forced k=3 exit remains above start exactly at N>=9 in this range")

    reset_rows = []
    for j in range(33):
        H = 6 * j + 5
        n = (2**(H + 3) - 13) // 9
        y = direct_block(n, (1, 2))
        check(y == 2**H - 1, ("reset endpoint", H))
        check(direct_repetitions(n, (1, 2)) == 1, ("reset source", H))
        check(direct_repetitions(y, (1,)) == H - 1 if H < 100
              else v2(y + 1) - 1 == H - 1, ("reset run", H))
        check(v2(n + 5) == 5 and v2(y + 5) == 2 and v2(y + 1) == H,
              ("budget", H))
        check(Fraction(9 * n + 13, 8) == y + 1, ("pullback", H))
        if j < 3:
            reset_rows.append((H, n, y, H - 1))
    print("FINITE-EXACT unbounded-reset family: H=6j+5, j=0..32")
    print("(H, start, seam, subsequent k=1 count):", reset_rows)

    for word, sign, expected in (((2,), 1, 1), ((1,), -1, 1),
                                 ((1, 2), -1, 5),
                                 ((1, 1, 1, 2, 1, 1, 4), -1, 17)):
        reps = [residue(word * j, sign=sign) for j in range(1, 9)]
        check(reps[-1] == expected, ("cycle representatives", word, reps))
        check(all(x <= y for x, y in zip(reps, reps[1:])), reps)
        print("Representative control:", word, "sign", sign, reps)
    for q, n, word in ((5, 1, (1, 4)), (5, 13, (1, 1, 5)),
                       (5, 17, (1, 3, 3))):
        check(direct_block(n, word, q) == n, ("5x+1 cycle", n))
        a, b, c = block(word, q)
        check((a - b) * n + c == 0, ("5x+1 defect", n))
    print("Hostile 5x+1 cycles at 1,13,17 retained; formula is multiplier-general")
    print("script_sha256", hashlib.sha256(Path(__file__).read_bytes()).hexdigest())


if __name__ == "__main__":
    main()
