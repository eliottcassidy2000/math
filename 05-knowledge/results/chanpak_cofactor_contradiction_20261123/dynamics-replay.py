"""Exact finite illustrations of the seven-ray prefix coverage theorem."""

from fractions import Fraction


def odd_height_in_relaxed_ray(height):
    if height < 1 or height % 2 == 0:
        return False
    q = (height - 1) // 2
    while q >= 283:
        w, leftover = divmod(q - 3, 140)
        if not leftover and w >= 2 and not (w & (w - 1)):
            return True
        q //= 3
    return False


def shell(L):
    X = 3 ** (L - 1)
    intervals = []
    s = 1
    while (q := 140 * 2 ** s + 3) < 3 * X:
        n = 1
        while q * n < X:
            n *= 3
        assert X <= q * n < 3 * X
        intervals.append((q * n, (q + 1) * n))
        s += 1
    seen = X
    covered = 0
    gap = (X, X)
    for start, end in sorted(intervals) + [(3 * X, 3 * X)]:
        if start - seen > gap[1] - gap[0]:
            gap = (seen, start)
        covered += max(0, end - max(seen, start))
        seen = max(seen, end)
    frac = Fraction(covered, 2 * X)
    assert frac < Fraction(3, 280)
    assert (gap[1] - gap[0]) * (1 + len(intervals)) >= 2 * X - covered
    if L <= 10:
        explicit = sum(odd_height_in_relaxed_ray(2 * q + 1)
                       for q in range(X, 3 * X))
        assert explicit == covered
    return frac, len(intervals), gap


if __name__ == "__main__":
    for length in range(1, 11):
        shell(length)
    for L in (8, 16, 32, 64, 128, 256):
        frac, n, (a, b) = shell(L)
        print(f"L={L}: {n} anchors, coverage={float(frac):.9%}, "
              f"largest gap fraction={float(Fraction(b-a,2*3**(L-1))):.9%}")
