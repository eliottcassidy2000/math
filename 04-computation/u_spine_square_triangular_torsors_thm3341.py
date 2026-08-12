#!/usr/bin/env python3
"""Exact companion for THM-3341; dependency-free and -O safe."""

from math import gcd, isqrt


def require(condition, message):
    if not condition:
        raise RuntimeError(message)


def c_value(t):
    return 2 * t * t + 2 * t + 1


def q_value(t):
    return 2 * c_value(t) + 1


def triangular(n):
    return n * (n + 1) // 2


def is_square(n):
    r = isqrt(n)
    return r * r == n


def triangular_index(n):
    d = 8 * n + 1
    r = isqrt(d)
    if r * r == d and r % 2 == 1:
        return (r - 1) // 2
    return None


def mmul(a, b):
    return tuple(
        tuple(sum(a[i][k] * b[k][j] for k in range(len(b))) for j in range(len(b[0])))
        for i in range(len(a))
    )


def mvec(a, v):
    return tuple(sum(a[i][j] * v[j] for j in range(len(v))) for i in range(len(a)))


def mpow(a, n):
    out = tuple(tuple(int(i == j) for j in range(len(a))) for i in range(len(a)))
    base = a
    while n:
        if n & 1:
            out = mmul(out, base)
        base = mmul(base, base)
        n //= 2
    return out


def pell_numbers(n):
    p = [0, 1]
    while len(p) <= n:
        p.append(2 * p[-1] + p[-2])
    return p


def factor(n):
    out = {}
    p = 2
    while p * p <= n:
        while n % p == 0:
            out[p] = out.get(p, 0) + 1
            n //= p
        p = 3 if p == 2 else p + 2
    if n > 1:
        out[n] = out.get(n, 0) + 1
    return out


def primitive_parents(c):
    out = set()
    for n in range(1, isqrt(c) + 1):
        m2 = c - n * n
        m = isqrt(m2)
        if m <= n or m * m != m2:
            continue
        if gcd(m, n) != 1 or (m - n) % 2 == 0:
            continue
        out.add((m * m - n * n, 2 * m * n, c))
    return out


def gamma(triple):
    a, b, c = triple
    return (abs(a * a - b * b), 2 * a * b, c * c)


def u_parent(t):
    return (2 * t + 1, 2 * t * (t + 1), 2 * t * (t + 1) + 1)


def main():
    limit = 1_000_000
    p = pell_numbers(2500)

    square_rows = []
    for k in range(30):
        m = p[2 * k + 1]
        t = (p[2 * k] + p[2 * k + 1] - 1) // 2
        require(c_value(t) == m * m, f"Pell square row failed at k={k}")
        square_rows.append((t, m))
    for k in range(28):
        t0, m0 = square_rows[k]
        t1, m1 = square_rows[k + 1]
        t2, m2 = square_rows[k + 2]
        require(m2 == 6 * m1 - m0, f"square-root recurrence failed at k={k}")
        require(t2 == 6 * t1 - t0 + 2, f"square-index recurrence failed at k={k}")

    brute_square_c = []
    brute_triangular_c = []
    brute_square_q = []
    brute_triangular_q = []
    for t in range(limit + 1):
        cv = c_value(t)
        qv = q_value(t)
        if is_square(cv):
            brute_square_c.append(t)
        if triangular_index(cv) is not None:
            brute_triangular_c.append(t)
        if is_square(qv):
            brute_square_q.append(t)
        if triangular_index(qv) is not None:
            brute_triangular_q.append(t)
        require(cv % 4 == 1 and qv % 4 == 3, f"mod-four separation failed at t={t}")

    expected_square_c = [t for t, _ in square_rows if t <= limit]
    require(brute_square_c == expected_square_c, "direct/Pell C-square census mismatch")
    require(brute_triangular_c == [0], "C triangular hostile failed")
    require(brute_square_q == [], "Q square hostile failed")

    # Square-triangular rows connect consecutive C-square roots and Markov points.
    selector_rows = []
    x, rroot = 1, 0
    for k in range(1, 16):
        x, rroot = 3 * x + 8 * rroot, x + 3 * rroot
        nindex = (x - 1) // 2
        require(triangular(nindex) == rroot * rroot, f"square-triangular row failed at k={k}")
        mminus, mplus = x - 2 * rroot, x + 2 * rroot
        tminus = 2 * rroot - nindex - 1
        tplus = 2 * rroot + nindex
        require(c_value(tminus) == mminus * mminus, f"lower C-square compiler failed at k={k}")
        require(c_value(tplus) == mplus * mplus, f"upper C-square compiler failed at k={k}")
        require(mminus * mminus + mplus * mplus + 4 == 6 * mminus * mplus, "Markov row failed")
        require(mminus * mplus - (2 * rroot) ** 2 == 1, "Pell-Cassini carry failed")
        require((tminus, mminus) == square_rows[k - 1], f"lower adjacent square row failed at k={k}")
        require((tplus, mplus) == square_rows[k], f"upper adjacent square row failed at k={k}")
        selector_rows.append((nindex, rroot, mminus, mplus, tminus, tplus))

    cannonball = [row for row in selector_rows if row[0] % 2 == 1 and row[0] >= 3]
    aligned = []
    for row in cannonball:
        nindex, rroot, _, _, _, _ = row
        h = (nindex - 1) // 2
        if h * (h + 1) * (2 * h + 1) // 6 == (2 * rroot) ** 2:
            aligned.append((h, 2 * rroot, nindex))
    require(aligned == [(24, 70, 49)], "cannonball alignment hostile failed")

    # Gaussian squaring transplants the positive/middle Berggren ray to C-square U-depths.
    U = ((1, -2, 2), (2, -1, 2), (2, -2, 3))
    Aplus = ((1, 2, 2), (2, 1, 2), (2, 2, 3))
    D = ((-1, 2, 2), (-2, 1, 2), (-2, 2, 3))
    root = (3, 4, 5)
    middle = root
    transplant_depths = []
    for j in range(9):
        t = square_rows[j + 1][0]
        image = gamma(middle)
        expected = u_parent(t)
        require(image == expected, f"Gaussian branch transplant failed at j={j}")
        require(mvec(mpow(U, t - 1), root) == expected, f"U-depth address failed at j={j}")
        transplant_depths.append(t - 1)
        middle = mvec(Aplus, middle)
    for k in range(1, 9):
        nindex = selector_rows[k - 1][0]
        require(square_rows[k][0] - square_rows[k - 1][0] == 2 * nindex + 1, "transplant drift failed")

    # Strong-divisibility controls for the unbounded Boolean-fibre proof.
    for a in range(1, 50):
        for b in range(1, 50):
            require(gcd(p[a], p[b]) == p[gcd(a, b)], f"Pell strong divisibility failed at {(a, b)}")
            if a % b == 0:
                require(p[a] % p[b] == 0, f"Pell divisibility failed at {(a, b)}")
    odd_prime_indices = (3, 5, 7, 11)
    L = 1
    for ell in odd_prime_indices:
        L *= ell
    for ell in odd_prime_indices:
        require(p[L] % p[ell] == 0, f"Pell divisor construction failed at ell={ell}")
    for i, ell in enumerate(odd_prime_indices):
        for ell2 in odd_prime_indices[i + 1 :]:
            require(gcd(p[ell], p[ell2]) == 1, "Pell divisor factors not coprime")

    first_multi = None
    for t, m in square_rows[1:]:
        if len(factor(m)) >= 2:
            first_multi = (t, m, factor(m))
            break
    require(first_multi == (696, 985, {5: 1, 197: 1}), "first square-selector fibre failed")
    c696 = c_value(696)
    parents696 = primitive_parents(c696)
    expected_parents696 = {(1393, 970224, 970225), (522767, 817344, 970225)}
    require(parents696 == expected_parents696, "c=985^2 parent fibre failed")
    other = root
    for letter in "UUUUUDADUDDU":
        other = mvec({"U": U, "A": Aplus, "D": D}[letter], other)
    require(other == (522767, 817344, 970225), "second c=985^2 ancestry word failed")
    require(mvec(mpow(U, 695), root) == (1393, 970224, 970225), "U^695 ancestry failed")

    # Q_t triangular: the two odd-X norm-17 branches.
    qtri_rows = []
    branch_rows = []
    for seed in ((0, 2), (6, 18)):
        t, nu = seed
        local = []
        for _ in range(12):
            require(q_value(t) == triangular(nu), f"norm-17 branch failed at {(t, nu)}")
            require((2 * nu + 1) ** 2 - 8 * (2 * t + 1) ** 2 == 17, "norm-17 equation failed")
            require(nu % 16 == 2, f"nu mod16 failed at {(t, nu)}")
            local.append((t, nu))
            t, nu = 17 * t + 6 * nu + 11, 48 * t + 17 * nu + 32
        branch_rows.append(local)
        qtri_rows.extend(local)
    expected_qtri = sorted(t for t, _ in qtri_rows if t <= limit)
    require(brute_triangular_q == expected_qtri, "direct/norm-17 Q-triangular census mismatch")

    first_positive_qtri = sorted((t, nu) for branch in branch_rows for t, nu in branch if t > 0)[:5]
    require(
        first_positive_qtri == [(6, 18), (23, 66), (221, 626), (798, 2258), (7524, 21282)],
        "first Q-triangular rows failed",
    )
    grade_controls = []
    for t in (6, 23, 221):
        fac = factor(c_value(t))
        grade_controls.append((t, c_value(t), fac, 2 ** (len(fac) - 1)))
    require(
        grade_controls
        == [
            (6, 85, {5: 1, 17: 1}, 2),
            (23, 1105, {5: 1, 13: 1, 17: 1}, 4),
            (221, 98125, {5: 4, 157: 1}, 2),
        ],
        "Q-triangular Boolean grades failed",
    )

    print("THM-3341 VERIFIED-EXACT")
    print("direct_search_t=0..1000000")
    print(f"C_square_t={brute_square_c}")
    print(f"C_triangular_t={brute_triangular_c}; Q_square_t={brute_square_q}")
    print(f"Q_triangular_t={brute_triangular_q}")
    print("C_square_formula=t=(P_2k+P_(2k+1)-1)/2; sqrt=P_(2k+1)")
    print("square_triangular_compiler=(M-,M+); adjacent_C_square_indices=(t-,t+)")
    print("cannonball_alignment=(h,2R,N)=(24,70,49)")
    print(f"gaussian_transplant_U_depths={transplant_depths}")
    print("first_multi_square_fibre=t696; sqrt=985=5*197; parent_count=2")
    print("c696_words=U^695 | UUUUUDADUDDU")
    print("Q_triangular_norm17_seeds=(t,nu)=(0,2),(6,18); nu=2 mod16")
    print(f"Q_triangular_grade_controls={grade_controls}")
    print("ALL CHECKS PASSED")


if __name__ == "__main__":
    main()
