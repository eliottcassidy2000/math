"""Next-level analogue of THM-4820 (C): rational n with no split of slack <= 1 (audit; bears on HYP-9426 (3)).

For Q(x) = 1 mod 3 the residue r + s of a split x = y + z is forced by beta (r + s = 1 - beta mod 3), so the slack is a
function of beta alone:  slack(beta) = floor((1 - beta)/3) for beta <= 0, and (2 beta + 1 - ((1 - beta) mod 3))/3 for
beta >= 1.  Hence slack <= 1 iff beta in {-4, ..., 2}.  For a rational n (odd, 3 not dividing n), a split with a given
beta is a representation  n^2 - 2 beta = s^2 + 3 b^2  (THM-4820 (C) translation; b != 0 for beta = 0).
Since a minimum tree's first nontrivial split has slack <= total slack, g(n) >= floor(n^2/3) + min-slack(n).

1. brute-force check of slack(beta) on all splits in a box;
2. brute-force check of the representation criterion for n <= 151 (direct enumeration of splits);
3. scan: odd n <= NMAX, 3 not dividing n, with min-slack >= 2 (so g(n) >= floor(n^2/3) + 2, window-free).
"""
import sys

import sympy

NMAX = int(sys.argv[1]) if len(sys.argv) > 1 else 3000


def Q(a, b):
    return a * a - a * b + b * b


def slack_beta(beta):
    if beta <= 0:
        return (1 - beta) // 3
    return (2 * beta + 1 - ((1 - beta) % 3)) // 3


# 1. slack is a function of beta when Q(x) = 1 mod 3
bad = 0
cnt = 0
R = 12
for xa in range(-R, R + 1):
    for xb in range(-R, R + 1):
        qx = Q(xa, xb)
        if qx % 3 != 1:
            continue
        for ya in range(-R, R + 1):
            for yb in range(-R, R + 1):
                za, zb = xa - ya, xb - yb
                if (ya, yb) == (0, 0) or (za, zb) == (0, 0):
                    continue
                qy, qz = Q(ya, yb), Q(za, zb)
                beta = qx - qy - qz
                sl = max(beta, 0) + qy // 3 + qz // 3 - qx // 3
                cnt += 1
                bad += sl != slack_beta(beta)
print(f"1. slack = slack(beta) for Q(x) = 1 mod 3: {cnt} splits checked, {bad} mismatches;  slack(beta) for beta = -6..4:",
      {b: slack_beta(b) for b in range(-6, 5)})


def rep(m, nontrivial_square=None):
    """is m = s^2 + 3 b^2 (with b != 0 if m == nontrivial_square)?"""
    if m <= 0:
        return False
    if any(e % 2 for p, e in sympy.factorint(m).items() if p % 3 == 2):
        return False
    if nontrivial_square is not None and m == nontrivial_square:
        # need b != 0: m = n^2 represented with b != 0 iff n has a prime factor = 1 mod 3 (THM-4820 (C))
        return any(p % 3 == 1 for p in sympy.factorint(int(round(m ** 0.5))))
    return True


def min_slack_arith(n, cap=3):
    for sl in range(cap):
        for beta in range(-12, 8):
            if slack_beta(beta) != sl:
                continue
            m = n * n - 2 * beta
            if rep(m, nontrivial_square=n * n if beta == 0 else None):
                # b = 0 representations of non-square m do not exist; m = s^2 exactly gives b = 0, y rational
                if beta != 0 and int(round(m ** 0.5)) ** 2 == m:
                    # check whether a b != 0 representation exists as well (brute force; only tiny n)
                    if not any((m - 3 * b * b) >= 0 and int(round((m - 3 * b * b) ** 0.5)) ** 2 == m - 3 * b * b
                               for b in range(1, int((m / 3) ** 0.5) + 1)):
                        # only b = 0: y = (s + n)/2 rational; still a genuine split unless y in {0, n}
                        s = int(round(m ** 0.5))
                        if (s + n) // 2 in (0, n) and (n - s) // 2 in (0, n):
                            continue
                return sl
    return cap


def min_slack_brute(n):
    L = 2 * n + 2
    best = 99
    for a in range(-L, L + 1):
        for b in range(-L, L + 1):
            if (a, b) in ((0, 0), (n, 0)):
                continue
            qy, qz = Q(a, b), Q(n - a, -b)
            beta = n * n - qy - qz
            best = min(best, max(beta, 0) + qy // 3 + qz // 3 - (n * n) // 3)
    return best


bad = 0
for n in range(5, 152, 2):
    if n % 3 == 0:
        continue
    a, b = min_slack_arith(n), min(min_slack_brute(n), 3)
    bad += a != b
    if a != b:
        print("   mismatch", n, a, b)
print(f"2. arithmetic min-slack == brute-force min-slack for odd n in [5, 151], 3 not dividing n: {bad} mismatches")
lvl2 = []
lvl1 = 0
for n in range(5, NMAX + 1, 2):
    if n % 3 == 0:
        continue
    ms = min_slack_arith(n)
    lvl1 += ms >= 1
    if ms >= 2:
        lvl2.append(n)
print(f"3. odd n <= {NMAX}, 3 not dividing n: {lvl1} with no tight split (THM-4820 (C) exceptions), "
      f"{len(lvl2)} with NO split of slack <= 1 (g(n) >= floor(n^2/3) + 2):")
print("  ", lvl2)
for n in lvl2[:3]:
    print(f"   n={n}:", {c: dict(sympy.factorint(n * n + c)) for c in (0, 2, -4, 6, 8)})
