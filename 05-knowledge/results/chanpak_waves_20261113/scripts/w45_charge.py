"""Charge-adapted coordinate at a root x: r = sum over even-order paths P from x of H(T-V(P)),
s = same over odd-order (incl. the singleton), g = H(T-x), a = 2r+1-h (charge), b = s-r-g (readiness defect).
Claims: rooted append (z beats only x): r' = 2r + b, b' = 0, a' = a + 2b, h' = 2h-1+a (= h+2r).
Top shield (C3 over everything): r' = 3r, s' = 3s, g' = 3g, h' = 3h, so a' = 3a-2, b' = 3b."""
import random
from w45_hk import ham_split, sub
def regs(n, out, x):
    H = ham_split(n, out)[0]
    rest = [v for v in range(n) if v != x]
    g = ham_split(*sub(out, rest))[0]
    r = s = 0
    def ext(p, used):
        nonlocal r, s
        h = ham_split(*sub(out, [v for v in range(n) if not used >> v & 1]))[0]
        if len(p) % 2 == 0: r += h
        else: s += h
        for w in range(n):
            if not used >> w & 1 and out[p[-1]] >> w & 1:
                ext(p + [w], used | 1 << w)
    ext([x], 1 << x)
    return H, g, r, s
def append(n, out, x):
    o = out[:] + [0]
    for v in range(n):
        if v == x: o[n] |= 1 << v
        else: o[v] |= 1 << n
    return n + 1, o
def shield(n, out):
    o = out[:] + [0, 0, 0]; X = [n, n + 1, n + 2]
    for t in range(3):
        o[X[t]] |= 1 << X[(t + 1) % 3]
        for v in range(n): o[X[t]] |= 1 << v
    return n + 3, o
rng = random.Random(9); bad = 0; tests = 0
for trial in range(120):
    n = rng.randint(2, 6); out = [0] * n
    for i in range(n):
        for j in range(i + 1, n):
            if rng.random() < .5: out[i] |= 1 << j
            else: out[j] |= 1 << i
    x = rng.randrange(n)
    h, g, r, s = regs(n, out, x); a = 2 * r + 1 - h; b = s - r - g
    n2, o2 = append(n, out, x); h2, g2, r2, s2 = regs(n2, o2, x)
    ok1 = (r2 == 2 * r + b and s2 - r2 - g2 == 0 and 2 * r2 + 1 - h2 == a + 2 * b and h2 == 2 * h - 1 + a)
    n3, o3 = shield(n, out); h3, g3, r3, s3 = regs(n3, o3, x)
    ok2 = (r3 == 3 * r and s3 == 3 * s and g3 == 3 * g and h3 == 3 * h)
    tests += 1; bad += (not ok1) + (not ok2)
print(f"charge identities (append and top shield): {tests} random rooted tournaments, failures {bad}")
