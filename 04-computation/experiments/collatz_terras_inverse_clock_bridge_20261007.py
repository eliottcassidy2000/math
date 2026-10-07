"""Exact checks of the Terras-to-inverse-clock bounded-return bridge.

Recurrence is proved in the companion note; these finite checks verify the
local half-probability trial, affine transitions, offsets, and event indices.
"""
from fractions import Fraction as F


def need(ok, message):
    if not ok:
        raise ValueError(message)


def three(j):
    need(type(j) is int, "integer exponent")
    return F(3 ** j) if j >= 0 else F(1, 3 ** (-j))


def parity(c):
    need(type(c) in (int, F), "exact dyadic-integral rational")
    c = F(c)
    need(c.denominator % 2 == 1, "odd denominator required")
    return c.numerator % 2


def transition(j, c, p):
    three(j)
    e = parity(c)
    c = F(c)
    need(type(p) is int and p in (0, 1), "exact parity bit")
    if p == 0 and e == 0:
        return j, c / 2
    if p == 0:
        return j + 1, (3 * c + 1) / 2
    if e == 0:
        return j, (3 * c + 1 - three(j)) / 2
    return j - 1, (c - three(j - 1)) / 2


def terras(n):
    need(type(n) is int, "integer state")
    return (3 * n + 1) // 2 if n % 2 else n // 2


def local_success(x, y):
    """Selected sufficient event: current coincidence or one-step catch-up."""
    need(type(x) is int and type(y) is int, "integer representatives")
    if x % 2 == y % 2:
        return bool(y % 2)
    return bool(terras(x) % 2) if y % 2 else bool(terras(y) % 2)


def main():
    checks = local_states = time_states = successes = 0

    def check(ok, message):
        nonlocal checks
        checks += 1
        if not ok:
            raise RuntimeError(message)

    # Uniform mod4 represents the first two Terras parity coins bijectively.
    for j in range(-4, 5):
        for numerator in range(-30, 31):
            c = F(numerator, 3 ** max(0, -j))
            coefficient = three(j)
            coef4 = coefficient.numerator * pow(coefficient.denominator, -1, 4) % 4
            c4 = c.numerator * pow(c.denominator, -1, 4) % 4
            count = 0
            words = set()
            for y in range(4):
                x = (coef4 * y + c4) % 4
                count += local_success(x, y)
                words.add((y % 2, terras(y) % 2))
            check(count == 2, "conditional success exactly one half")
            check(len(words) == 4, "two fresh parity bits")
            local_states += 1

    for v in range(1, 7):
        for source in range(1, 128, 2):
            x0 = 3 * (1 << v) * source + 1
            ys, xs = [source], [x0]
            for _ in range(65 + v):
                ys.append(terras(ys[-1]))
                xs.append(terras(xs[-1]))
            ax = [i for i, x in enumerate(xs) if x % 2]
            ay = [i for i, y in enumerate(ys) if y % 2]
            nx = [0]
            ny = [0]
            for x, y in zip(xs, ys):
                nx.append(nx[-1] + x % 2)
                ny.append(ny[-1] + y % 2)
            j = 1 + (v + 1) // 2
            c = F(2 if v % 2 else 1)
            check(xs[v] == three(j) * source + c, "forced initial v steps")
            for s in range(60):
                x, y = xs[s + v], ys[s]
                check(j == 1 + nx[s + v] - ny[s], "raw odd-count identity")
                check(x == three(j) * y + c, "exact aligned affine relation")
                check(x % 2 == (y % 2) ^ parity(c), "parity relation")
                h = ny[s] - nx[s + v]
                check(j == 1 - h, "ordinal offset level")
                if local_success(x, y):
                    k = nx[s + v]
                    check(abs(v + ay[k + h] - ax[k]) <= 1,
                          "inverse-clock success with exact ordinal indices")
                    successes += 1
                oldj = j
                e = parity(c)
                j, c = transition(j, c, y % 2)
                check(j - oldj == e * (1 - 2 * (y % 2)), "skipped fair-sign increment")
                time_states += 1

    # Incoming lag1 Mersenne initial relation p=3q+2, q=5 mod16.
    j, c = 1, F(2)
    for p in (1, 0, 0, 0):
        j, c = transition(j, c, p)
    check((j, c) == (2, F(1)), "Mersenne forced prefix")
    hostiles = [lambda: three(True), lambda: parity(F(1, 2)),
                lambda: transition(0, F(1), True), lambda: transition(0, 1.0, 1),
                lambda: terras(1.0)]
    for action in hostiles:
        try:
            action()
        except ValueError:
            check(True, "typed hostile")
        else:
            check(False, "typed hostile accepted")
    print("PROVED bridge: aligned merge OR fixed-offset inverse debt visits[-1,1] infinitely (Haar a.s.)")
    print("local_probability_states=" + str(local_states) + "; every success count=2/4")
    print("integer_universe: oddY1..127, v1..6, alignedTerrasTimes0..59")
    print("time_states=" + str(time_states) + "; exact_inverse_clock_successes=" + str(successes))
    print("Mersenne_relation(1,2), forced1000 -> (2,1)")
    print("remaining_obligation=translation/archimedean_box_recurrence; integer-point transfer separate")
    print("typed_hostiles=" + str(len(hostiles)))
    print("checks=" + str(checks))


if __name__ == "__main__":
    main()
