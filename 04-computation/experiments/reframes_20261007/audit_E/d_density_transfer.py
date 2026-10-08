#!/usr/bin/env python3
"""Audit E, task D: statements 3 (unequal times), 5 (integer density transfer) and 6 (conjugation x -> r x).

 (D1) For p in {5,7,9}, several e, and every residue r mod 2^K: chain absorption by K, driven by the ACTUAL integer
      parities of T^t(n) for a big lift n = r + 2^K L, versus the direct integer test 'T^t(n) = T^t(n+e) for some t <= K'.
      A second lift n' = r + 2^K L' must give the same chain verdict (residue-determinedness). Densities reported exactly.
 (D2) Small integers 1 <= n <= 2^16 (cycles possible): every direct equal-time merge that is NOT a chain absorption
      (i.e. equal time, unequal odd counts) is listed: these are the finitely many exceptions the written proof of 5
      does not mention.
 (D3) Statement 5's sentence 'For p = 5 this density is below 0.0087 for every K' tested for e != 1.
 (D4) Statement 6: T_{p,r}(r x) = r T_{p,1}(x) on random 2-adic x (mod 2^M), r odd incl. negative; and the residue
      count of merges by K for px + r with translation r e equals that of px + 1 with translation e.
 (D5) Unequal-time meetings among small integers (they exist, inside cycles), for the Reading's caveat."""
import random
rnd = random.Random(4606)

def T(p, x, r=1):
    return x >> 1 if x % 2 == 0 else (p * x + r) >> 1

def chain_step(p, k, E, b):
    s = E & 1
    if s == 0:
        if b == 0: return k, E >> 1
        if k >= 0: return k, (p * E + 1 - p ** k) >> 1
        return k, (p * E + p ** (-k) - 1) >> 1
    if b == 0:
        if k >= 0: return k + 1, (p * E + 1) >> 1
        return k + 1, (E + p ** (-k - 1)) >> 1
    if k >= 1: return k - 1, (E - p ** (k - 1)) >> 1
    return k - 1, (p * E - 1) >> 1

def chain_absorbed_by(p, e, n, K):
    k, E = 0, e
    v = n
    for t in range(1, K + 1):
        k, E = chain_step(p, k, E, v & 1)
        v = T(p, v)
        if k == 0 and E == 0: return t
    return None

def direct_merge_by(p, e, n, K):
    a, b = n, n + e
    for t in range(1, K + 1):
        a, b = T(p, a), T(p, b)
        if a == b: return t
    return None

def d1():
    tot_bad = 0
    for p, e, K in ((5, 1, 16), (7, 1, 16), (9, 1, 16), (5, -1, 14), (5, 2, 14), (5, 3, 14), (7, -3, 14), (11, 1, 16), (13, 1, 16)):
        L1 = rnd.getrandbits(200) | (1 << 199); L2 = rnd.getrandbits(200) | (1 << 199)
        cnt = 0; bad_dir = 0; bad_lift = 0
        for r in range(1 << K):
            n1 = r + (L1 << K); n2 = r + (L2 << K)
            c1 = chain_absorbed_by(p, e, n1, K)
            c2 = chain_absorbed_by(p, e, n2, K)
            d = direct_merge_by(p, e, n1, K)
            if (c1 is None) != (c2 is None) or c1 != c2: bad_lift += 1
            if c1 != d: bad_dir += 1
            cnt += c1 is not None
        tot_bad += bad_dir + bad_lift
        print(f"(D1) p={p:2d} e={e:2d} K={K}: absorbed residues {cnt}/2^{K} = {cnt / 2**K:.6f};"
              f" chain-vs-direct-integer disagreements {bad_dir}; lift-dependence {bad_lift}")
    return tot_bad == 0

def d2():
    for p, e in ((5, 1), (7, 1), (9, 1), (5, 3)):
        K = 40
        ex = []
        for n in range(1, (1 << 16) + 1):
            a, b = n, n + e
            k, E = 0, e
            for t in range(1, K + 1):
                k, E = chain_step(p, k, E, a & 1)
                a, b = T(p, a), T(p, b)
                if a == b:
                    if not (k == 0 and E == 0):
                        ex.append((n, t, k, a))
                    break
                if k == 0 and E == 0:
                    print("IMPOSSIBLE: absorbed without integer merge", p, e, n, t); break
        print(f"(D2) p={p} e={e}: n in [1, 2^16], t <= 40: equal-time integer merges with unequal odd counts (k_t != 0): "
              f"{len(ex)}; first few (n, t, k_t, value): {ex[:6]}")

def d3():
    # densities of {n : T^t(n) = T^t(n+e) for some t <= K} for p = 5 and several e, exactly by residues mod 2^K
    p = 5
    L = rnd.getrandbits(200) | (1 << 199)
    for e in (1, 2, 3, -3, 5, 6, 13):
        row = []
        for K in (5, 8, 12, 14):
            cnt = sum(1 for r in range(1 << K) if chain_absorbed_by(p, e, r + (L << K), K) is not None)
            row.append(f"K={K}: {cnt}/2^{K}={cnt / 2**K:.5f}")
        flag = "  <-- exceeds 0.0087" if cnt / 2**14 > 0.0087 else ""
        print(f"(D3) p=5, e={e:3d}: " + "; ".join(row) + flag)
    # the explicit witness class for e = 3
    w = [n for n in range(1, 200) if direct_merge_by(5, 3, n, 5) == 5]
    print(f"(D3) witness: for p = 5, e = 3 the class n = 16 mod 32 merges at t = 5: {w[:6]} ...; e.g. T^5(16) = T^5(19) = "
          f"{[T(5, T(5, T(5, T(5, T(5, 16)))))][0]}")

def d4():
    M = 512; mod = 1 << M
    bad = 0
    for _ in range(20000):
        p = rnd.choice([3, 5, 7, 9, 11])
        r = rnd.choice([1, -1, 3, -3, 5, 7, -7, 9, 15, 101])
        x = rnd.getrandbits(M)
        lhs = T(p, (r * x) % mod, r) % (mod >> 1)
        rhs = (r * T(p, x, 1)) % (mod >> 1)
        if lhs != rhs: bad += 1
    print(f"(D4) T_(p,r)(r x) = r T_(p,1)(x) mod 2^{M-1} on 20000 random 2-adic x, r odd incl. negative: failures {bad}")
    K = 12
    for p, r, e in ((5, 3, 1), (5, -1, 1), (7, 5, 1), (5, 3, 2)):
        L = rnd.getrandbits(100) | (1 << 99)
        cnt_r = 0; cnt_1 = 0
        for res in range(1 << K):
            n = res + (L << K)
            a, b = n, n + r * e
            for t in range(1, K + 1):
                a, b = T(p, a, r), T(p, b, r)
                if a == b: cnt_r += 1; break
            a, b = n, n + e
            for t in range(1, K + 1):
                a, b = T(p, a), T(p, b)
                if a == b: cnt_1 += 1; break
        print(f"(D4) residues mod 2^{K} merging by K: p={p}, px+{r} with translation {r*e}: {cnt_r}; px+1 with translation {e}: {cnt_1}"
              f" -> {'equal' if cnt_r == cnt_1 else 'DIFFERENT'}")
    return bad == 0

def d5():
    p = 5
    found = []
    for n in range(1, 400):
        orb_a = {}
        a = n
        for m in range(60):
            orb_a.setdefault(a, m); a = T(p, a)
        b = n + 1
        for m2 in range(60):
            if b in orb_a and orb_a[b] != m2:
                found.append((n, orb_a[b], m2, b)); break
            b = T(p, b)
        if len(found) >= 5: break
    print(f"(D5) p=5 integers: unequal-time meetings T^m(n) = T^m'(n+1), m != m' (inside cycles): first examples (n, m, m', value) {found}")

if __name__ == '__main__':
    ok = d1()
    d2()
    d3()
    ok &= d4()
    d5()
    print("D: chain/direct agreement and conjugation checks PASS" if ok else "SOME D CHECK FAILED")
