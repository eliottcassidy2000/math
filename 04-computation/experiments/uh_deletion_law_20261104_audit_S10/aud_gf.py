"""aud_gf.py -- independent audit code (uh_deletion_law_20261104_audit_S10). Shares no code with the session.

Checks THM-4874 (3) against THM-4869's core counts:
  per block of order w: cores 2^(d-w-2) at every d >= w + 2, plus one core at d = w + 1 when A != W (THM-4869 B4);
  rooted prefixes R(2) = 1, R(m) = F(2m - 5) (THM-4868 (F));
  classes with that block at order n = sum_d cores(d) R(n - d + 1);
compares with F(2n - 2w - 4) (A = W) and F(2n - 2w - 3) (A != W) on the stated ranges, sums them by type and compares with
THM-4869's c_k(d) and (C4) (u_3 = F(2n-9), u_4 = 2F(2n-11), u_5 = 4F(2n-12) - 3F(2n-14)), and checks
L(2n - 12) = F(2n - 12) + 2F(2n - 13) = F(2n - 13) + F(2n - 11). Also prints the first n at which each per-block formula
holds, to locate the true lower ends of the ranges."""


def F(m):
    if m < 0:
        return (-1) ** (m + 1) * F(-m)
    a, b = 0, 1
    for _ in range(m):
        a, b = b, a + b
    return a


def L(m):
    return F(m - 1) + F(m + 1)


def R(m):
    return 1 if m == 2 else (F(2 * m - 5) if m >= 3 else 0)


# (name, type k, order w, A == W)
BLOCKS = [("C3 all", 3, 3, True), ("W4 {c,d}", 3, 4, False), ("W4 {a,c,d}", 4, 4, False), ("W4 {b,c,d}", 4, 4, False),
          ("W4 all", 5, 4, True), ("AT5 {4,5}", 5, 5, False), ("W5 A1", 5, 5, False), ("W5 A2", 5, 5, False)]


def cores(w, AW, d):
    c = 2 ** (d - w - 2) if d >= w + 2 else 0
    if not AW and d == w + 1:
        c += 1
    return c


def classes(w, AW, n):
    return sum(cores(w, AW, d) * R(n - d + 1) for d in range(3, n))


def ck(k, d):   # THM-4869 B5 / note section 2.2 core counts
    if k == 3:
        return 2 if d == 5 else (3 * 2 ** (d - 6) if d >= 6 else 0)
    if k == 4:
        return 2 if d == 5 else (2 ** (d - 5) if d >= 6 else 0)
    if k == 5:
        return 0 if d <= 5 else (4 if d == 6 else 5 * 2 ** (d - 7))


ok = True
NMAX = 40
print("per-block counts vs (3):")
for name, k, w, AW in BLOCKS:
    first = None
    for n in range(w + 1, NMAX + 1):
        pred = F(2 * n - 2 * w - 4) if AW else F(2 * n - 2 * w - 3)
        got = classes(w, AW, n)
        if got == pred and first is None:
            first = n
        if got != pred and first is not None:
            print("  FAIL after first agreement", name, n, got, pred)
            ok = False
    lo = w + 2
    print(f"  {name}: formula {'F(2n-2w-4)' if AW else 'F(2n-2w-3)'} with w = {w} holds for all n in [{first}, {NMAX}]"
          f" (theorem's range starts at n = {6 if k == 3 else 7}); counts n=6..11: {[classes(w, AW, n) for n in range(6, 12)]}")
    if first is None or first > (6 if k == 3 else 7):
        ok = False
print("core counts per type vs THM-4869 c_k(d):")
for k in (3, 4, 5):
    bad = [d for d in range(3, NMAX) if sum(cores(w, AW, d) for _, kk, w, AW in BLOCKS if kk == k) != ck(k, d)]
    print(f"  k = {k}: mismatching d: {bad}")
    ok &= not bad
print("sums by type vs THM-4869 (C4):")
for n in range(7, NMAX + 1):
    u3 = sum(classes(w, AW, n) for _, kk, w, AW in BLOCKS if kk == 3)
    u4 = sum(classes(w, AW, n) for _, kk, w, AW in BLOCKS if kk == 4)
    u5 = sum(classes(w, AW, n) for _, kk, w, AW in BLOCKS if kk == 5)
    if (u3, u4, u5) != (F(2 * n - 9), 2 * F(2 * n - 11), 4 * F(2 * n - 12) - 3 * F(2 * n - 14)):
        print("  FAIL", n, u3, u4, u5)
        ok = False
u3_6 = sum(classes(w, AW, 6) for _, kk, w, AW in BLOCKS if kk == 3)
print(f"  u_3(6) from blocks = {u3_6} (F(3) = {F(3)})")
ok &= u3_6 == F(3)
print("Lucas identity L(2n-12) = F(2n-12) + 2F(2n-13) = F(2n-13) + F(2n-11), n = 7..40:",
      all(L(2 * n - 12) == F(2 * n - 12) + 2 * F(2 * n - 13) == F(2 * n - 13) + F(2 * n - 11) for n in range(7, 41)))
print("k = 5 split: 3-deletion classes F(2n-13), 2-deletion classes L(2n-12), n = 7..14:",
      [(F(2 * n - 13), L(2 * n - 12)) for n in range(7, 15)])
print("HYP-9467 table rows n = 10..15 (C3, W4{c,d}, W4 type 4 each, W4 all, each 5-block):")
for n in range(10, 16):
    print(f"  n={n}: {classes(3, True, n)}, {classes(4, False, n)}, {classes(4, False, n)}, {classes(4, True, n)}, {classes(5, False, n)};"
          f" u3={F(2*n-9)} u4={2*F(2*n-11)} u5={4*F(2*n-12)-3*F(2*n-14)}")
print("ALL OK:", ok)
