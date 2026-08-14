#!/usr/bin/env python3
"""Independent PeriodicPL audit of the expanded centered continuum compiler."""
from fractions import Fraction as F
from hashlib import sha256
from bisect import bisect_right
from math import gcd
from random import Random
from subprocess import check_output


ROOT = __import__("pathlib").Path(__file__).resolve().parents[1]


def circle_indicator(x):
    x = F(x) % 1
    return int(x <= F(1, 14) or x >= F(13, 14))


def remainder_overlap(x, length):
    """Integral_0^length 1_C(x+u)du, exact and nonperiodized in u."""
    x = F(x) % 1
    length = F(length)
    answer = F(0)
    for k in range(-2, 4):
        lo = max(F(0), F(k) - F(1, 14) - x)
        hi = min(length, F(k) + F(1, 14) - x)
        if hi > lo:
            answer += hi - lo
    return answer


class GeneralBlockPL:
    """Independent exact PL construction from convolution breakpoints."""

    def __init__(self, y, alpha, d):
        self.y, self.alpha, self.d = F(y), F(alpha) % 1, d
        length = self.y / 7
        turns = length.numerator // length.denominator
        self.rem = length - turns
        self.base = F(turns, 7 * self.y)
        one_breaks = {F(0), F(1), F(1, 14), F(13, 14)}
        for h in (F(1, 14), F(13, 14)):
            one_breaks.add((h - self.rem) % 1)
        raw = {F(0), F(1)}
        for s in range(d):
            shift = s * self.alpha
            for b in one_breaks:
                raw.add((b - shift) % 1)
        self.bp = sorted(raw)
        if self.bp[-1] != 1:
            self.bp.append(F(1))
        self.value, self.slope, self.prefix_values = [], [], [F(0)]
        for left, right in zip(self.bp, self.bp[1:]):
            mid = (left + right) / 2
            value = sum((self.one_value(left + s * self.alpha) for s in range(d)), F(0))
            slope = sum((self.one_slope(mid + s * self.alpha) for s in range(d)), F(0))
            self.value.append(value)
            self.slope.append(slope)
            width = right - left
            self.prefix_values.append(self.prefix_values[-1] + value * width + slope * width * width / 2)
        self.mean = self.prefix_values[-1]
        if self.mean != F(d, 49):
            raise RuntimeError(("mean", y, alpha, d, self.mean))

    def one_value(self, x):
        return self.base + remainder_overlap(x, self.rem) / self.y

    def one_slope(self, x):
        return F(circle_indicator(x + self.rem) - circle_indicator(x), self.y)

    def at(self, x):
        x = F(x) % 1
        i = min(bisect_right(self.bp, x) - 1, len(self.value) - 1)
        return self.value[i] + self.slope[i] * (x - self.bp[i])

    def prefix(self, x):
        x = F(x)
        if x == 1:
            return self.prefix_values[-1]
        i = min(bisect_right(self.bp, x) - 1, len(self.value) - 1)
        width = x - self.bp[i]
        return self.prefix_values[i] + self.value[i] * width + self.slope[i] * width * width / 2

    def interval(self, start, length):
        start, length = F(start) % 1, F(length)
        if not 0 <= length <= 1:
            raise RuntimeError(("interval", start, length))
        if start + length <= 1:
            return self.prefix(start + length) - self.prefix(start)
        return self.mean - self.prefix(start) + self.prefix(start + length - 1)


def pos2(x):
    return x * x if x > 0 else 0


def cdf2(s, t, u):
    return pos2(s) - pos2(s - t) - pos2(s - u) + pos2(s - t - u)


def area_num(y, t, u, hperiod):
    y %= hperiod
    half = hperiod // 14
    return sum(
        cdf2(k * hperiod + half - y, t, u)
        - cdf2(k * hperiod - half - y, t, u)
        for k in range(-3, 5)
    )


def overlap_num(y, u, hperiod):
    y %= hperiod
    half = hperiod // 14
    answer = 0
    for k in range(-3, 5):
        lo = max(0, k * hperiod - half - y)
        hi = min(u, k * hperiod + half - y)
        answer += max(0, hi - lo)
    return answer


def integer_compiler_value(row):
    d, a, c, L, j, e, f, A = row
    D = d + a
    R, S = e * j % L, f * j % L
    H = 14 * L * d
    X = 14 * (R * D - S * d)
    U = 2 * L * D
    if A == 0:
        num = sum(overlap_num(X + 14 * L * a * s - U // 2, U, H) for s in range(d))
        den = D * H
    else:
        T = 14 * abs(A)
        start = X - (T if A < 0 else 0)
        num = sum(
            area_num(start + 14 * L * a * s - U // 2, T, U, H)
            for s in range(d)
        )
        den = 392 * L * d * D * abs(A)
    return F(num, den)


def periodic_pl_value(row):
    d, a, c, L, j, e, f, A = row
    D = d + a
    G = GeneralBlockPL(F(D, d), F(a, d), d)
    R, S = e * j % L, f * j % L
    # The physical convolution is centered: C++ uses X-U/2, i.e. subtracts
    # half of the normalized integration length D/(7d).
    X = (F(R * D, d) - S) / L - F(D, 14 * d)
    if A == 0:
        return G.at(X) / d
    T = F(abs(A), L * d)
    turns = T.numerator // T.denominator
    theta = T - turns
    start = X if A > 0 else X - theta
    return F(L, abs(A)) * (turns * G.mean + G.interval(start, theta))


def contexts():
    raw = check_output(
        ["git", "show", "HEAD:04-computation/lrc14_disconnected_head263_contexts_20260812.txt"],
        cwd=ROOT,
        text=True,
    )
    return tuple(tuple(map(int, line.split())) for line in raw.splitlines())


def candidate_rows(ctx):
    for d in range(1, 9):
        for a in range(0, 7 * d + 1):
            D = d + a
            g = gcd(a, d)
            for c in range(-13, 14):
                if c == 0 or c % g or (a == 0 and c < 0) or (a == 7 * d and c > 0):
                    continue
                for L, j, e, f in ctx:
                    A = L * c + e * D - d * f
                    if abs(A) * 672 <= L * d * 679:
                        yield d, a, c, L, j, e, f, A


def main():
    rng = Random(20260812)
    ordinary, superunit, zero = [], [], []
    seen = [0, 0, 0]
    targets = [256, 768, 64]
    for row in candidate_rows(contexts()):
        d, _, _, L, _, _, _, A = row
        group = 2 if A == 0 else (1 if abs(A) > L * d else 0)
        bank = (ordinary, superunit, zero)[group]
        seen[group] += 1
        if len(bank) < targets[group]:
            bank.append(row)
        else:
            k = rng.randrange(seen[group])
            if k < targets[group]:
                bank[k] = row
    rows = tuple(ordinary + superunit + zero)
    failures = []
    semantic = sha256()
    minimum = None
    for row in rows:
        left = integer_compiler_value(row)
        right = periodic_pl_value(row)
        if left != right:
            failures.append((row, left, right))
        record = (row, left)
        semantic.update((repr(record) + "\n").encode())
        if minimum is None or record[1] < minimum[1]:
            minimum = record
    if failures:
        raise RuntimeError(failures[:3])
    print("LRC14 EXPANDED CENTERED CONTINUUM INDEPENDENT PERIODICPL AUDIT")
    print("universe_seen ordinary superunit zero", tuple(seen))
    print("sampled ordinary superunit zero", (len(ordinary), len(superunit), len(zero)))
    print("exact_equalities", len(rows), "failures", len(failures))
    print("sample_minimum", minimum)
    print("semantic_sha256", semantic.hexdigest())
    print("status=PASS")


if __name__ == "__main__":
    main()
