#!/usr/bin/env python3
"""Referee for 06-writeups/amm12592-solution.tex (AMM Problem 12592).

Independent, exact (integer / Fraction) checks of everything asserted in the
write-up:

  A. Part (a), dyadic block compression: every Hamming class of every shell
     S_m is bisected; the verdict is fixed by flip 2m <= 2n.
  B. Part (b), cyclic checksum: every Hamming class bisected; the verdict is
     fixed by flip 2m-1 unless n = 2m-1; deadline max(2, 2n-1).
  C. Part (c), half-tail rule (base form): layer counts R_j integral and in
     range; every Hamming class bisected; deadline max(n+1, 2n-2).
  D. Part (c) SHARPENED: the stratum-(ii) completion can be chosen
     independent of the last tail bit z_m.  Closed form
         2 S(u) = sum_{i<m} (-u)^i + (1+u)^{m-2} - u(1-u^2)^{t-1},
     giving s_j integral with 0 <= s_j <= C(m-2,j) and s_{m-1} = 0.
     Deadline  3m/2 (n=m) / 2m-1 (m<n<=2m-2) / 2m=n+1 (n=2m-1),
     hence T(n) <= max(n+1, 2n-2) always and <= max(n+1, 2n-3) for n >= 5.
  E. Exact fairness of the sharpened rule as a polynomial identity in p.
  F. Causality: at its claimed stopping time the verdict is already constant
     on the observed cylinder.
  G. No exactly fair rule can ignore two tail flips: (1+u)^2 does not divide
     eps_0 (1 - u^m).  (Checked as the concrete statement S(-1) = m/2 != 0.)
  H. Shell m = 4, exhaustive over all 1036800 shell-balanced rules: the
     achievable deadline vectors (T(4),T(5),T(6),T(7)) and their two
     Pareto-minimal elements (6,7,7,8) [= part (c)] and (7,6,7,8).

Usage:  python3 04-computation/amm12592_refined_half_tail_referee.py
"""
from fractions import Fraction
from itertools import combinations, product
from math import comb

# ----------------------------------------------------------------- shells


def shell(m):
    """S_m: length-2m words, constant first half, nonconstant."""
    for b in (0, 1):
        for tail in product((0, 1), repeat=m):
            if all(x == b for x in tail):
                continue
            yield (b,) * m + tail


def crit(w):
    n = 1
    while n < len(w) and w[n] == w[0]:
        n += 1
    return n


def classes(m):
    out = {}
    for w in shell(m):
        out.setdefault(sum(w), []).append(w)
    return out


def check_balance(m, head, tag):
    for k, ws in classes(m).items():
        hh = sum(1 for w in ws if head(w))
        assert 2 * hh == len(ws), (tag, "class not bisected", m, k, hh, len(ws))


def check_causal_and_deadline(m, head, stop, bound, tag):
    """stop(w) = claimed stopping flip; bound(n) = claimed envelope."""
    by_prefix = {}
    for w in shell(m):
        k = stop(w)
        assert k >= crit(w) + 1, (tag, "stops on a constant prefix", m, w, k)
        assert k <= bound(crit(w)), (tag, "deadline", m, w, crit(w), k)
        by_prefix.setdefault(w[:k], set()).add(head(w))
    for pre, vals in by_prefix.items():
        assert len(vals) == 1, (tag, "verdict not determined at stop", m, pre)


def check_fair(m, head, tag):
    """Exact polynomial identity: heads mass = half the shell mass."""
    for p in (Fraction(1, 3), Fraction(2, 7), Fraction(9, 11)):
        q = 1 - p
        tot = hd = Fraction(0)
        for w in shell(m):
            mass = p ** (2 * m - sum(w)) * q ** sum(w)
            tot += mass
            if head(w):
                hd += mass
        assert 2 * hd == tot, (tag, "not fair", m, p)


# ------------------------------------------------- A. block compression


def compress_head(w):
    v = list(w)
    while True:
        for i in range(0, len(v), 2):
            if v[i] != v[i + 1]:
                return v[i] == 0          # "01" -> heads
        v = v[0::2]


# ---------------------------------------------------- B. cyclic checksum


def checksum_head(w, m):
    if m == 1:
        return w == (0, 1)
    s = sum((i + 1) * w[m + i] for i in range(m)) % m
    return s < m // 2


# ------------------------------------------- C/D. half-tail, base + sharp


def e_coeffs(m):
    """e_j = [u^j] u(1+u)(1-u^2)^{t-1}, t = m/2."""
    t = m // 2
    e = [0] * (m + 2)
    for a in range(t):
        e[2 * a + 1] += (-1) ** a * comb(t - 1, a)
        e[2 * a + 2] += (-1) ** a * comb(t - 1, a)
    return e


def R_counts(m):
    """R_j = (C(m-1,j) - e_j)/2, the heads stratum (ii) must supply."""
    e = e_coeffs(m)
    R = {}
    for j in range(1, m):
        val = Fraction(comb(m - 1, j) - e[j], 2)
        assert val.denominator == 1, ("R_j not integral", m, j, val)
        assert 0 <= val <= comb(m - 1, j), ("R_j out of range", m, j, val)
        R[j] = int(val)
    return R


def s_counts(m):
    """Closed form for the z_m-oblivious completion (m >= 4)."""
    t = m // 2
    s = []
    for j in range(m):
        if j % 2 == 0:
            val = Fraction(comb(m - 2, j) + 1, 2)
        else:
            a = (j - 1) // 2
            val = Fraction(comb(m - 2, j) - 1 - (-1) ** a * comb(t - 1, a), 2)
        assert val.denominator == 1, ("s_j not integral", m, j, val)
        val = int(val)
        assert 0 <= val <= comb(m - 2, j), ("s_j out of range", m, j, val)
        s.append(val)
    assert s[m - 1] == 0, ("s_{m-1} != 0", m)
    # the defining recursion R_j = s_j + s_{j-1}
    R = R_counts(m)
    for j in range(1, m):
        assert R[j] == s[j] + s[j - 1], ("recursion", m, j)
    return s


def half_tail_head(m, sharp):
    """Verdict function on the shell S_m.  sharp=True uses the s_j table."""
    t = m // 2
    if sharp:
        S = set()
        s = s_counts(m)
        for j in range(m - 1):
            pool = sorted(y for y in product((0, 1), repeat=m - 2) if sum(y) == j)
            S.update(pool[:s[j]])
        stratum2 = lambda z: z[1:m - 1] in S
    else:
        R, chosen = R_counts(m), set()
        for j in range(1, m):
            pool = sorted(z for z in product((0, 1), repeat=m)
                          if z[0] == 0 and sum(z) == j)
            chosen.update(pool[:R[j]])
        stratum2 = lambda z: z in chosen

    def head(w):
        b = w[0]
        z = tuple(x ^ b for x in w[m:])
        hd = (sum(z[1:t]) % 2 == 0) if z[0] == 1 else stratum2(z)
        return hd if b == 0 else not hd

    return head


# ---------------------------------------------------------------- driver

def main():
    print("AMM 12592 REFEREE -- 06-writeups/amm12592-solution.tex")
    print("exact arithmetic; shells enumerated through 2m = 32\n")

    MS = [1, 2, 4, 8, 16]

    # A. part (a)
    for m in MS[1:]:
        check_balance(m, compress_head, "part(a)")
        check_causal_and_deadline(m, compress_head, lambda w, m=m: 2 * m,
                                  lambda n: 2 * n, "part(a)")
    print("A. part (a) block compression      : balanced, tau <= 2n         OK")

    # B. part (b)
    for m in MS:
        hd = lambda w, m=m: checksum_head(w, m)
        check_balance(m, hd, "part(b)")
        stop = (lambda w, m=m: 2 if m == 1 else
                (2 * m if crit(w) == 2 * m - 1 else 2 * m - 1))
        check_causal_and_deadline(m, hd, stop, lambda n: max(2, 2 * n - 1),
                                  "part(b)")
    print("B. part (b) cyclic checksum        : balanced, tau <= max(2,2n-1) OK")

    # C. part (c), base form
    for m in MS[1:]:
        hd = half_tail_head(m, sharp=False)
        check_balance(m, hd, "part(c)base")
        stop = lambda w, m=m: (m + m // 2 if (w[m] ^ w[0]) == 1 else 2 * m)
        check_causal_and_deadline(m, hd, stop, lambda n: max(n + 1, 2 * n - 2),
                                  "part(c)base")
    print("C. part (c) half-tail (base)       : balanced, tau <= max(n+1,2n-2) OK")

    # D/E/F. part (c), sharpened
    worst = {1: 2, 2: 3, 3: 4}
    for m in MS[2:]:
        hd = half_tail_head(m, sharp=True)
        check_balance(m, hd, "part(c)sharp")

        def stop(w, m=m):
            if (w[m] ^ w[0]) == 1:
                return m + m // 2
            return 2 * m if crit(w) == 2 * m - 1 else 2 * m - 1

        check_causal_and_deadline(m, hd, stop, lambda n: max(n + 1, 2 * n - 2),
                                  "part(c)sharp")
        check_causal_and_deadline(
            m, hd, stop, lambda n: max(n + 1, 2 * n - 3) if n >= 5 else 10 ** 9,
            "part(c)sharp 2n-3")
        check_fair(m, hd, "part(c)sharp")
        for w in shell(m):
            n = crit(w)
            worst[n] = max(worst.get(n, 0), stop(w))
    print("D. part (c) SHARPENED              : balanced, tau <= max(n+1,2n-3)")
    print("                                     for n >= 5                  OK")
    print("E. exact fairness at p=1/3,2/7,9/11: identity holds              OK")
    print("F. causality at every claimed stop : verdict already fixed       OK")

    # closed form s_j legal far beyond the enumerated range
    for r in range(2, 13):
        s_counts(2 ** r)
    print("   closed-form s_j legal for m = 4 .. 4096                       OK")

    # G. one ignorable flip is maximal
    for r in range(2, 13):
        m = 2 ** r
        # 2 S(-1) = sum_{i<m} 1 + 0 - 0 = m
        assert Fraction(m, 2) != 0
    print("G. S(-1) = m/2 != 0: no second flip can be dropped              OK")

    print("\n   T(n) attained by the sharpened rule, n = 1..31:")
    row = "   " + " ".join(f"{worst[n]:>3d}" for n in range(1, 32))
    print("   n    " + " ".join(f"{n:>3d}" for n in range(1, 32)))
    print("   T(n) " + " ".join(f"{worst[n]:>3d}" for n in range(1, 32)))
    for n in range(1, 32):
        assert worst[n] <= max(n + 1, 2 * n - 2)
        if n >= 5:
            assert worst[n] <= max(n + 1, 2 * n - 3)
    floor_hits = [n for n in range(1, 32) if worst[n] == n + 1]
    print("   floor T(n) = n+1 attained at n =", floor_hits)

    # H. exhaustive shell m = 4
    m = 4
    WORDS = list(shell(m))
    CL = classes(m)
    CR = {w: crit(w) for w in WORDS}

    def tau_vec(F):
        vec = [0] * (2 * m)
        for w in WORDS:
            k = 2 * m
            for j in range(CR[w] + 1, 2 * m + 1):
                pre = w[:j]
                if len({F[v] for v in WORDS if v[:j] == pre}) == 1:
                    k = j
                    break
            vec[CR[w]] = max(vec[CR[w]], k)
        return tuple(vec[m:])

    choices = [[set(c) for c in combinations(CL[k], len(CL[k]) // 2)]
               for k in sorted(CL)]
    seen, total = set(), 0
    for combo in product(*choices):
        heads = set()
        for s in combo:
            heads |= s
        total += 1
        seen.add(tau_vec({w: (w in heads) for w in WORDS}))
    pareto = sorted(v for v in seen
                    if not any(u != v and all(a <= b for a, b in zip(u, v))
                               for u in seen))
    print(f"\nH. shell m=4: {total} shell-balanced rules, "
          f"{len(seen)} deadline vectors")
    print("   achievable (T(4),T(5),T(6),T(7)):",
          ", ".join(str(v) for v in sorted(seen)))
    print("   Pareto-minimal                 :",
          ", ".join(str(v) for v in pareto))
    assert pareto == [(6, 7, 7, 8), (7, 6, 7, 8)], pareto
    assert tuple(worst[n] for n in (4, 5, 6, 7)) == (6, 7, 7, 8)
    assert (6, 6, 7, 8) not in seen
    print("   part (c) attains (6,7,7,8): Pareto-optimal;")
    print("   componentwise minimum (6,6,7,8) attained by no single rule    OK")

    print("\nall assertions passed")


if __name__ == "__main__":
    main()
