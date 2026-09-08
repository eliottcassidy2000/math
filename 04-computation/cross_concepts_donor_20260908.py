"""Exact checksum orientation / dyadic divided-carry controls; no repo imports.

Run from repository root with Python, and repeat with Python -O.
All checks remain active under optimization. Standard library only.
"""
from collections import Counter
from fractions import Fraction
from itertools import product
from math import comb
from pathlib import Path

GATES = 0
LINES = []


def check(ok, label):
    global GATES
    GATES += 1
    if not ok:
        raise RuntimeError(label)


def say(message):
    LINES.append(message)


def trim(a):
    a = list(a)
    while len(a) > 1 and not a[-1]:
        a.pop()
    return a


def add(*args):
    a = [0] * max(map(len, args))
    for b in args:
        for i, v in enumerate(b):
            a[i] += v
    return trim(a)


def scale(a, c):
    return trim([c * v for v in a])


def mul(a, b):
    c = [0] * (len(a) + len(b) - 1)
    for i, v in enumerate(a):
        for j, w in enumerate(b):
            c[i + j] += v * w
    return trim(c)


def monomial(n):
    return [0] * n + [1]


def qpower(n):
    return [(-1) ** k * comb(n, k) for k in range(n + 1)]


def tpower(n):
    return [0] * n + qpower(n)


def half(a):
    check(all(v % 2 == 0 for v in a), "integral division by two")
    return trim([v // 2 for v in a])


def evalpoly(a, x):
    y = 0
    for v in reversed(a):
        y = y * x + v
    return y


def native_polynomial(words, orientation):
    ans = [0]
    for word in words:
        # Native variable is the probability of the initial symbol.
        n = word.count(orientation)
        ans = add(ans, [0] * n + qpower(len(word) - n))
    return ans


def checksum(word):
    m = len(word) // 2
    if m == 1:
        return word == (0, 1)
    return sum((i + 1) * b for i, b in enumerate(word[m:])) % m < m // 2


def main():
    # Independent literal enumeration of every shell word, not all 2^(2m) words.
    total_words = 0
    for m in (1, 2, 4, 8):
        heads = {0: [], 1: []}
        layers = Counter()
        balance = Counter()
        for orientation in (0, 1):
            for tail in product((0, 1), repeat=m):
                if all(b == orientation for b in tail):
                    continue
                word = (orientation,) * m + tail
                h = checksum(word)
                total_words += 1
                layers[sum(word)] += 1
                balance[sum(word)] += 1 if h else -1
                if h:
                    heads[orientation].append(word)
                n = next(i for i, b in enumerate(word) if b != word[0])
                deadline = 2 if m == 1 else (2 * m if n == 2 * m - 1 else 2 * m - 1)
                check(deadline <= max(2, 2 * n - 1), "literal deadline")
                if m >= 2 and n < 2 * m - 1:
                    changed = word[:-1] + (1 - word[-1],)
                    check(checksum(changed) == h, "last-bit causal merge")
        check(all(v == 0 for v in balance.values()), "full Hamming shell balance")
        f = native_polynomial(heads[0], 0)
        g = native_polynomial(heads[1], 1)
        if m == 1:
            expected_f, expected_g = [0, 1, -1], [0]
        else:
            middle = add([1], scale(monomial(m), -1), scale(qpower(m), -1))
            expected_f = half(mul(monomial(m), middle))
            expected_g = add(expected_f, tpower(m))
        check(f == expected_f, "literal zero orientation formula")
        check(g == expected_g, "literal one orientation formula")
        say(f"SHELL m={m}: words={sum(layers.values())}, orientation head words={len(heads[0])}/{len(heads[1])}, all layers balanced")
    say(f"LITERAL UNIVERSE: {total_words} shell words; m=1,2,4,8; no omitted nonconstant tails")

    # Independent polynomial recurrence against direct binomial expressions.
    a, b, c = [0, 1], [0], [0]
    lacunary = [0, 1, -1]
    for r in range(1, 9):
        a2 = mul(a, a)
        a3 = mul(a2, a)
        a4 = mul(a2, a2)
        b = add(a4, scale(a3, -1), scale(mul(add(a, scale(a2, -1)), b), 2), scale(mul(b, b), 2))
        a = a2
        c = add(c, b)
        m = 2 ** r
        direct_b = half(add(tpower(m), scale(monomial(m), -1), monomial(2 * m)))
        check(b == direct_b, "nonlinear carry transition")
        lacunary = add(lacunary, tpower(m))
        direct_c = half(add(lacunary, [0, -1], monomial(2 * m)))
        check(c == direct_c, "carry accumulation")
        check(all((v - (1 if i == 1 else 0) + (1 if i == 2 * m else 0)) % 2 == 0 for i, v in enumerate(lacunary)), "finite F2 telescoping")
        f = add([0, 1, -1], scale(c, -1))
        g = add([0, 0, 1], c, scale(monomial(2 * m), -1))
        # Polynomial reflection done by fresh binomial expansion.
        reflected_g = [0]
        for degree, value in enumerate(g):
            reflected_g = add(reflected_g, scale(qpower(degree), value))
        finite_mass_twice = add([1], scale(monomial(2 * m), -1), scale(qpower(2 * m), -1))
        check(scale(add(f, reflected_g), 2) == finite_mass_twice, "actual finite fairness consequence")
        for p in (Fraction(1, 7), Fraction(1, 3), Fraction(1, 2), Fraction(2, 3)):
            check(0 <= evalpoly(f, p) <= p, "zero orientation probability")
            check(0 <= evalpoly(g, p) <= p, "one orientation probability")
        if r in (1, 2, 3, 8):
            say(f"CARRY r={r}: m={m}, polynomial degree={len(c)-1}; nonlinear/direct identities and exact finite fairness PASS")
    say("FIRST carry coefficients: C(p)=-p^3+p^4-2p^5+3p^6-2p^7+p^8+... (checked below)")
    check(c[:9] == [0, 0, 0, -1, 1, -2, 3, -2, 1], "first nontrivial lost carry")

    # Hostile: non-dyadic tail checksum no longer bisects a Hamming layer.
    hostile = Counter(checksum((0,) * 6 + tail) for tail in product((0, 1), repeat=6) if sum(tail) == 2)
    check(hostile[True] == 7 and hostile[False] == 8, "m6 weight2 hostile")
    say("HOSTILE m=6, tail weight=2: heads/tails=7/8; dyadic decoder hypothesis is necessary")

    # Exact rational finite controls for the analytic growth inequalities.
    for k in range(1, 7):
        t = 1 - Fraction(1, 2 ** k)
        partial = sum(t ** (2 ** r) for r in range(k + 4))
        check(k - 1 <= partial <= k + 3, "lacunary logarithmic bound finite control")
        for r in range(k + 1, k + 4):
            check(t ** (2 ** r) <= Fraction(1, 1 + 2 ** (r - k)), "Bernoulli upper bound")
    say("ANALYTIC CONTROLS: exact rational t=1-2^(-k), k=1..6; growth proof is all-real, not inferred from controls")
    say("SCOPE: exact inherited checksum architecture; no donor deadline improvement; no Rule30 dynamics map")
    say(f"PASS {GATES} always-active exact gates")


if __name__ == "__main__":
    main()
    transcript = "\n".join(LINES) + "\n"
    print(transcript, end="")
    destination = Path(__file__).resolve().parents[1] / "05-knowledge/results/cross_concepts_donor_20260908.out"
    destination.write_text(transcript, encoding="utf-8", newline="\n")
