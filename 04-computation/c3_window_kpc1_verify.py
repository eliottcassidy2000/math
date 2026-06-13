# c3_window_kpc1_verify.py — ADVERSARIAL VERIFICATION (session kind-pasteur-2026-06-10-S1)
# Independent check of the PROOF INGREDIENTS of claim A3 (THM-462 window step):
#   (a) exact algebra: W(n) := M(n)-M(n-1) = h(h+1)/2 for n=2h+1, h(h-1)/2 for n=2h,
#       and C(n,3) - m_nearreg(n) = M(n), and the L6 identity C(n,3)-C(n-1,2)=C(n-1,3),
#       all checked for n up to 5000 in exact integers;
#   (b) the h>=8 threshold: isqrt-exact check that sqrt(W(n)) <= h-2 iff claimed;
#   (c) the CONSTRUCTION: for n=8..100 and EVERY t in [0, W(n)], take the canonical
#       (greedy) Lagrange four-square decomposition t = a1^2+a2^2+a3^2+a4^2, perturb a
#       near-regular score vector (lower 4 scores by a_j, raise 4 others by a_j; upper
#       half only when n even), then verify the DEFINITION of a Landau sequence
#       (every prefix sum >= C(k,2), total = C(n,2), 0 <= s <= n-1) and that
#       f = m_nearreg(n) + t exactly (hence c3 = M(n) - t covers [M(n-1), M(n)]).
#   (d) stronger: for n=16..40 check that EVERY four-square decomposition works
#       (the proof text invokes an arbitrary Lagrange decomposition).
# All code fresh; nothing reused from the worker.
from math import comb, isqrt

def M(n):
    return (n**3 - n) // 24 if n % 2 == 1 else (n**3 - 4 * n) // 24

def m_nearreg(n):
    h = n // 2
    if n % 2 == 1:
        return n * comb(h, 2)
    return h * comb(h - 1, 2) + h * comb(h, 2)

print("=== (a) exact algebra n<=5000 ===")
alg_ok = True
for n in range(4, 5001):
    h = n // 2
    W = M(n) - M(n - 1)
    expect = h * (h + 1) // 2 if n % 2 == 1 else h * (h - 1) // 2
    if W != expect:
        print(f"W algebra FAIL at n={n}: {W} != {expect}")
        alg_ok = False
    if comb(n, 3) - m_nearreg(n) != M(n):
        print(f"M = C(n,3)-m_nearreg FAIL at n={n}")
        alg_ok = False
    if comb(n, 3) - comb(n - 1, 2) != comb(n - 1, 3):
        print(f"L6 identity FAIL at n={n}")
        alg_ok = False
print("algebra:", "PASS" if alg_ok else "FAIL")

print()
print("=== (b) threshold sqrt(W(n)) <= h-2 ===")
thr_ok = True
for n in range(4, 2001):
    h = n // 2
    W = M(n) - M(n - 1)
    holds = (W <= (h - 2) ** 2) and h >= 2   # a_max <= isqrt(t) <= isqrt(W) <= h-2
    claimed = (n >= 16) if n % 2 == 1 else (n >= 12)
    # the proof only NEEDS n>=16; record where it first holds for each parity
    if n % 2 == 1 and n >= 17 and not holds:
        print(f"odd threshold FAIL at n={n}")
        thr_ok = False
    if n % 2 == 0 and n >= 16 and not holds:
        print(f"even threshold FAIL at n={n}")
        thr_ok = False
first_odd = next(n for n in range(5, 100, 2) if M(n) - M(n - 1) <= (n // 2 - 2) ** 2)
first_even = next(n for n in range(4, 100, 2) if M(n) - M(n - 1) <= (n // 2 - 2) ** 2)
print(f"sqrt(W)<=h-2 first holds: odd n={first_odd}, even n={first_even}; "
      f"so all n>=16 covered: {'PASS' if thr_ok else 'FAIL'}")

def all_four_squares(t):
    res = []
    for x in range(isqrt(t), -1, -1):
        if 4 * x * x < t:
            break
        r1 = t - x * x
        for y in range(min(x, isqrt(r1)), -1, -1):
            if 3 * y * y < r1:
                break
            r2 = r1 - y * y
            for z in range(min(y, isqrt(r2)), -1, -1):
                if 2 * z * z < r2:
                    break
                r3 = r2 - z * z
                w = isqrt(r3)
                if w * w == r3 and w <= z:
                    res.append((x, y, z, w))
    return res

QUADCACHE = {}
def quads(t):
    if t not in QUADCACHE:
        QUADCACHE[t] = all_four_squares(t)
    return QUADCACHE[t]

def construction_ok(n, t, quad):
    """Build perturbed near-regular score vector; verify Landau DEFINITION + f value."""
    h = n // 2
    if n % 2 == 1:
        if n < 9:
            return False
        scores = [h] * n
        low_idx = [0, 1, 2, 3]
        high_idx = [4, 5, 6, 7]
    else:
        if h < 8:
            return False
        scores = [h - 1] * h + [h] * h
        low_idx = [h, h + 1, h + 2, h + 3]
        high_idx = [h + 4, h + 5, h + 6, h + 7]
    for j, a in enumerate(quad):
        scores[low_idx[j]] -= a
        scores[high_idx[j]] += a
    s = sorted(scores)
    if s[0] < 0 or s[-1] > n - 1:
        return False
    ps = 0
    for k, v in enumerate(s, start=1):
        ps += v
        if ps < k * (k - 1) // 2:
            return False
    if ps != n * (n - 1) // 2:
        return False
    return sum(v * (v - 1) // 2 for v in s) == m_nearreg(n) + t

print()
print("=== (c) canonical construction, n=8..100, every t in [0, W(n)] ===")
fail_by_n = {}
n_checked = 0
for n in range(8, 101):
    W = M(n) - M(n - 1)
    fails = []
    for t in range(W + 1):
        qs = quads(t)
        assert qs, f"Lagrange failure?! t={t}"
        n_checked += 1
        if not construction_ok(n, t, qs[0]):
            fails.append(t)
    if fails:
        fail_by_n[n] = fails
print(f"total (n,t) pairs checked: {n_checked}")
for n in sorted(fail_by_n):
    fs = fail_by_n[n]
    print(f"  n={n}: {len(fs)} failing t values "
          f"({fs if len(fs) <= 12 else str(fs[:6]) + ' ... ' + str(fs[-3:])})")
bad_big = [n for n in fail_by_n if n >= 16]
print("no failures for n>=16:", "PASS" if not bad_big else f"FAIL {bad_big}")

print()
print("=== (d) ALL four-square decompositions, n=16..40, every t ===")
alldec_ok = True
pairs = 0
for n in range(16, 41):
    W = M(n) - M(n - 1)
    for t in range(W + 1):
        for q in quads(t):
            pairs += 1
            if not construction_ok(n, t, q):
                print(f"FAIL n={n} t={t} quad={q}")
                alldec_ok = False
print(f"(n,t,quad) triples checked: {pairs}; all pass: "
      f"{'PASS' if alldec_ok else 'FAIL'}")

print()
print("OVERALL window-step verification:",
      "PASS" if (alg_ok and thr_ok and not bad_big and alldec_ok) else "FAIL")
