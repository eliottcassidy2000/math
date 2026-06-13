# taxicab_meetmiddle_kpc1_verify.py — ADVERSARIAL VERIFIER (kind-pasteur-2026-06-10-S1)
# Independent check of worker claims A4 + A5 (and Ta(2) part of A6):
#   Enumerate ALL n = x^3 + y^3 <= 10^12, 1 <= y <= x (so x <= 10^4).
#   - count pairs, count n with >= 2 positive reps (worker: 44164746 / 40426)
#   - doubly-primitive n (>= 2 reps with gcd(x,y)=1): worker 5464, with 10931
#     primitive reps total
#   - among them: gcd of primitive d's > 1: worker 3406
#   - n containing an inert prime (p == 2 mod 3): worker 1952; first = 4342914
#     = 2 * 3^2 * 31 * 43 * 181; histogram of inert primes
#   - Cor C audit on every primitive rep: v_p(d) = v_p(n) for every inert p|n;
#     v3(d) >= v3(n) - 1; p^{v_p(n)} | gcd(d_1,...,d_k)
#   - first 25 doubly-primitive numbers + factorizations; first five completely split
#   - Ta(2) = 1729 (smallest n with two positive reps)
#
# FRESH METHOD (different from worker): numpy GLOBAL SORT of all 44M pair-sums
# (int64, exact: max value 2*10^12 << 2^63), runs-of-equal extraction; per-n
# y recovered by exact integer cube root; factorization by Miller-Rabin +
# Pollard rho (not trial division).
import math
import numpy as np
from collections import Counter

LIMIT = 10**12
XMAX = 10**4

def icbrt(n):
    if n < 0:
        return -1
    r = round(n ** (1.0 / 3.0))
    while r * r * r > n:
        r -= 1
    while (r + 1) ** 3 <= n:
        r += 1
    return r

# ---------- enumerate all pairs ----------
chunks_n = []
chunks_x = []
for x in range(1, XMAX + 1):
    x3 = x * x * x
    rem = LIMIT - x3
    if rem < 1:
        break
    ymax = min(x, icbrt(rem))
    if ymax < 1:
        continue
    ys = np.arange(1, ymax + 1, dtype=np.int64)
    chunks_n.append(ys * ys * ys + x3)
    chunks_x.append(np.full(ymax, x, dtype=np.int32))
Narr = np.concatenate(chunks_n)
Xarr = np.concatenate(chunks_x)
del chunks_n, chunks_x
npairs = len(Narr)
assert int(Narr.max()) <= LIMIT and int(Narr.max()) < 2**62  # exactness guard
print(f"[enum] pairs (1<=y<=x, n<=10^12): {npairs}")

order = np.argsort(Narr, kind="stable")
Ns = Narr[order]
Xs = Xarr[order]
del Narr, Xarr, order

bnd = np.flatnonzero(Ns[1:] != Ns[:-1])
starts = np.concatenate(([0], bnd + 1))
ends = np.concatenate((bnd + 1, [len(Ns)]))
lens = ends - starts
multi = np.flatnonzero(lens >= 2)
print(f"[enum] n with >= 2 positive reps: {len(multi)}")
print(f"[enum] Ta(2) candidate (smallest multi-rep n): {int(Ns[starts[multi[0]]])}")
print(f"[enum] rep-multiplicity histogram: {Counter(int(l) for l in lens[multi])}")

# ---------- collect reps for multi-rep n ----------
twos = []
for idx in multi:
    s, e = int(starts[idx]), int(ends[idx])
    n = int(Ns[s])
    reps = []
    for xv in (int(v) for v in Xs[s:e]):
        y3 = n - xv * xv * xv
        yv = icbrt(y3)
        assert yv * yv * yv == y3 and 1 <= yv <= xv
        reps.append((xv, yv))
    twos.append((n, reps))
del Ns, Xs

# ---------- Miller-Rabin + Pollard rho ----------
def is_prime(n):
    if n < 2:
        return False
    for p in (2, 3, 5, 7, 11, 13, 17, 19, 23, 29, 31, 37):
        if n % p == 0:
            return n == p
    d = n - 1
    r = 0
    while d % 2 == 0:
        d //= 2
        r += 1
    for a in (2, 3, 5, 7, 11, 13, 17, 19, 23, 29, 31, 37):
        x = pow(a, d, n)
        if x == 1 or x == n - 1:
            continue
        for _ in range(r - 1):
            x = x * x % n
            if x == n - 1:
                break
        else:
            return False
    return True

def pollard(n):
    if n % 2 == 0:
        return 2
    c = 1
    while True:
        x = y = 2
        d = 1
        while d == 1:
            x = (x * x + c) % n
            y = (y * y + c) % n
            y = (y * y + c) % n
            d = math.gcd(abs(x - y), n)
        if d != n:
            return d
        c += 1

def factor(n):
    if n == 1:
        return Counter()
    if is_prime(n):
        return Counter({n: 1})
    d = pollard(n)
    f = factor(d)
    f.update(factor(n // d))
    return f

def vP(n, p):
    v = 0
    while n % p == 0:
        n //= p
        v += 1
    return v

# ---------- doubly-primitive analysis ----------
doubly = []
for n, reps in twos:
    prim = [(x, y) for (x, y) in reps if math.gcd(x, y) == 1]
    if len(prim) >= 2:
        doubly.append((n, prim, reps))
nprim_total = sum(len(p) for _, p, _ in doubly)
print(f"[dp] doubly-primitive n (>=2 primitive reps): {len(doubly)}")
print(f"[dp] total primitive reps across them: {nprim_total}")

gcd_gt1 = 0
inert_carriers = []
audit_vp = audit_v3 = audit_gcd = audit_both = 0
inert_membership = Counter()
min_inert_hist = Counter()
for n, prim, reps in doubly:
    ds = [x + y for (x, y) in prim]
    g = math.gcd(*ds) if len(ds) > 1 else ds[0]
    for d2 in ds[2:]:
        g = math.gcd(g, d2)
    if g > 1:
        gcd_gt1 += 1
    fac = factor(n)
    inerts = sorted(p for p in fac if p % 3 == 2)
    v3n = fac.get(3, 0)
    # Cor C audits over every primitive rep
    for (x, y) in prim:
        d = x + y
        for p in inerts:
            if vP(d, p) != fac[p]:
                audit_vp += 1
        if vP(d, 3) < v3n - 1:
            audit_v3 += 1
    for p in inerts:
        if vP(g, p) < fac[p]:
            audit_gcd += 1
        if g % p != 0:
            audit_both += 1
    if inerts:
        inert_carriers.append((n, fac))
        for p in inerts:
            inert_membership[p] += 1
        min_inert_hist[inerts[0]] += 1

print(f"[dp] gcd of primitive d's > 1: {gcd_gt1}")
print(f"[dp] n with an inert prime (p==2 mod 3): {len(inert_carriers)}")
if inert_carriers:
    n0, f0 = min(inert_carriers)
    print(f"[dp] first inert carrier: {n0} = "
          + " * ".join(f"{p}^{e}" if e > 1 else f"{p}"
                       for p, e in sorted(f0.items())))
print(f"[audit CorC] v_p(d) != v_p(n) violations: {audit_vp}")
print(f"[audit CorC] v3(d) < v3(n)-1 violations: {audit_v3}")
print(f"[audit CorC] p^(v_p(n)) not dividing gcd(d's): {audit_gcd}")
print(f"[audit CorC] inert p not dividing BOTH d's: {audit_both}")
print(f"[hist] inert-prime membership counts (top 10): "
      f"{inert_membership.most_common(10)}")
print(f"[hist] minimal-inert-prime counts (top 10): "
      f"{min_inert_hist.most_common(10)}")

print("[dp] first 25 doubly-primitive n with factorizations:")
def split_status(fac):
    if all(p % 3 == 1 for p in fac):
        return "completely split"
    return "has 3/inert"
for n, prim, reps in doubly[:25]:
    fac = factor(n)
    fs = " * ".join(f"{p}^{e}" if e > 1 else f"{p}" for p, e in sorted(fac.items()))
    print(f"   {n} = {fs}   reps {prim}   [{split_status(fac)}]")

first5_split = all(all(p % 3 == 1 for p in factor(n)) for n, _, _ in doubly[:5])
print(f"[dp] first five doubly-primitive all completely split: {first5_split}")
print(f"[dp] first five: {[n for n, _, _ in doubly[:5]]}")

ok = (audit_vp == audit_v3 == audit_gcd == audit_both == 0)
print("VERDICT(CorC audits):", "PASS" if ok else "FAIL")
