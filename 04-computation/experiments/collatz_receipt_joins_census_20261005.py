#!/usr/bin/env python3
"""
Receipt joins, multiplicative fusion, and the first deficit of multiples of 3 (opus, 2026-10-05).

Objects (universal weak receipts note, collatz_universal_weak_receipts_20261005.md): U(u) = oddpart(3u+1) on odd
u >= 1; a weak receipt for n is a nonnegative edge multiset c with prod (u/U(u))^c(u) = n; its boundary
d c = sum c(u)([u] - [U(u)]); defect D_n(c) = d c - [n] + [1]; grounded-deficit rule G1.

Part A  (Q1) product joins: for odd a <= b, does the orbit of ab contain a PRODUCT x*y of orbit points of a and b
        (x, y > 1)?  Ordinary joins (x > 1 = y) are the familiar common futures; product joins would transport a
        fusion obligation R(a,b) to R(x,y) at a smaller composite.  Also: synchronous joins U^i(ab) = U^i(a) U^i(b).
Part B  synchronous no-go search to larger a, b, i.
Part C  (Q2) multiples of 3: the exact law of the first deficit of the source-anchored receipt (n = 3m must start
        with the edge 3m -> U(3m)); the stopping time k(n) = least k with U^k(n) < n; the 3-parent variant
        k'(n) = least k with pi(U^k(n)) < n, pi(v) = (2^a(v) v - 1)/3, a(v) = least a >= 1 with 2^a v = 1 mod 9
        (pi(v) is the smallest multiple of 3 with U(pi(v)) = v), which closes an induction inside multiples of 3.
Part D  (Q1) depth-2 joins U^2(ab) = U(a): the families b = [2^A (3a+1) - 3 - 2^alpha]/(9a), beyond F2 (depth 1).
Part E  the Applegate--Lagarias multiplier identity: for x = -1 mod 2^k, j <= k odd, m_j = (2^j+1)/3,
        U(m_j x) = oddpart(x + (x+1)/2^j); m_j is the 3x-1 trunk ((3 m_j - 1) = 2^j); the descent is from m_j x, not x.
Part F  (Q4) shadow mismatch: the least admissible exponent e of the S1 shadow for small n, and r/n.
Usage: python collatz_receipt_joins_census_20261005.py [Amax=301] [Bmax=1001] [Nmult3=3000000]
"""
import sys, math, time
from functools import lru_cache

def U(n):
    n = 3 * n + 1
    while n % 2 == 0:
        n //= 2
    return n

def orbit(n, cap=100000):
    out = [n]
    while n != 1 and len(out) < cap:
        n = U(n); out.append(n)
    return out

def trunk(x):
    """x = (4^e - 1)/3 ?"""
    y = 3 * x + 1
    return y & (y - 1) == 0 and (y.bit_length() - 1) % 2 == 0

def main():
    Amax = int(sys.argv[1]) if len(sys.argv) > 1 else 301
    Bmax = int(sys.argv[2]) if len(sys.argv) > 2 else 1001
    N3 = int(sys.argv[3]) if len(sys.argv) > 3 else 3_000_000
    out = []
    def P(*a):
        s = " ".join(str(x) for x in a); print(s, flush=True); out.append(s)
    t0 = time.time()

    # ---------------- Part A ----------------
    P(f"# Part A: product joins, odd 3 <= a <= b <= {Amax}")
    odds = list(range(3, Amax + 1, 2))
    orb = {a: orbit(a) for a in odds}
    pairs = 0; ordinary_above_trunk = 0; ordinary_any = 0; product_pairs = 0; sync_pairs = 0; l0_pairs = 0
    examples = []; sync_examples = []; l0_examples = []
    for ia, a in enumerate(odds):
        Oa = orb[a]
        for b in odds[ia:]:
            Ob = orb[b]
            Oab = orbit(a * b)
            pos = {z: i for i, z in enumerate(Oab) if i >= 1 and z != 1}   # depth >= 1, above 1
            pairs += 1
            # ordinary joins: orbit point of a (depth >= 0, > 1) on ab's orbit at depth >= 1
            ord_hit = [z for z in Oa if z > 1 and z in pos] + [z for z in Ob if z > 1 and z in pos]
            if ord_hit:
                ordinary_any += 1
                if any(not trunk(z) for z in ord_hit):
                    ordinary_above_trunk += 1
            found = False; sync = False; l0 = False
            for j, x in enumerate(Oa):
                if x == 1: break
                for l, y in enumerate(Ob):
                    if y == 1: break
                    z = x * y
                    if z in pos:
                        i = pos[z]
                        found = True
                        if i == j == l: sync = True
                        if l == 0 or j == 0: l0 = True
                        if len(examples) < 40:
                            examples.append((a, b, i, j, l, x, y, z))
                        if i == j == l and len(sync_examples) < 20:
                            sync_examples.append((a, b, i, x, y, z))
                        if (l == 0 or j == 0) and len(l0_examples) < 20:
                            l0_examples.append((a, b, i, j, l, x, y, z))
            product_pairs += found; sync_pairs += sync; l0_pairs += l0
    P(f"pairs {pairs}; ordinary joins (an orbit point >1 of a or b lies on the orbit of ab at depth >= 1): {ordinary_any} "
      f"({100*ordinary_any/pairs:.1f}%), of which above the trunk ray (4^e-1)/3: {ordinary_above_trunk} ({100*ordinary_above_trunk/pairs:.1f}%)")
    P(f"product joins U^i(ab) = U^j(a) U^l(b) with all three factors > 1: pairs {product_pairs} ({100*product_pairs/pairs:.2f}%); "
      f"synchronous (i=j=l): {sync_pairs}; with one factor at depth 0 (y = b or x = a): {l0_pairs}")
    P("  examples (a, b, i, j, l, x, y, z=xy):")
    for e in examples[:40]:
        P("   ", e)
    P("  synchronous examples:", sync_examples if sync_examples else "none")
    P("  depth-0 factor examples (a,b,i,j,l,x,y,z):", l0_examples if l0_examples else "none")
    P(f"  [{time.time()-t0:.0f}s]")

    # ---------------- Part B ----------------
    P(f"# Part B: synchronous search U^i(ab) = U^i(a) U^i(b), i >= 1, all > 1, odd 3 <= a <= b <= {Bmax}")
    odds2 = list(range(3, Bmax + 1, 2))
    orb2 = {a: orbit(a) for a in odds2}
    hits = []
    for ia, a in enumerate(odds2):
        Oa = orb2[a]
        for b in odds2[ia:]:
            Ob = orb2[b]
            z = a * b
            i = 0
            L = min(len(Oa), len(Ob))
            while True:
                i += 1
                if i >= L: break
                z = U(z)
                if z == 1: break
                x, y = Oa[i], Ob[i]
                if x == 1 or y == 1: break
                if z == x * y:
                    hits.append((a, b, i, x, y, z))
    P(f"synchronous hits: {len(hits)}", hits[:20])
    P(f"  [{time.time()-t0:.0f}s]")

    # ---------------- Part C ----------------
    P(f"# Part C: multiples of 3, n = 3m <= {N3}, m odd")
    def a9(v):
        a = 1; p = 2 % 9
        while (p * v) % 9 != 1:
            a += 1; p = (p * 2) % 9
            if a > 7: raise RuntimeError
        return a
    hist_k = {}; hist_kp = {}; law_ok = 0; law_n = 0; maxk = (0, 0); maxkp = (0, 0)
    firstdef_below = 0; cnt = 0
    class_count = {1: [0, 0], 3: [0, 0]}
    for m in range(1, N3 // 3 + 1, 2):
        n = 3 * m
        # stopping time in odd steps
        v = U(n); k = 1
        while v >= n:
            v = U(v); k += 1
        hist_k[k] = hist_k.get(k, 0) + 1
        if k > maxk[0]: maxk = (k, n)
        # exact law: k = 1 iff m = 3 mod 4
        law_n += 1
        if (k == 1) == (m % 4 == 3): law_ok += 1
        cl = m % 4
        class_count[cl][0] += 1
        if k == 1: class_count[cl][1] += 1
        # 3-parent variant
        v = U(n); kp = 1
        while True:
            a = a9(v)
            par = ((1 << a) * v - 1) // 3
            if par < n or v == 1:
                break
            v = U(v); kp += 1
        hist_kp[kp] = hist_kp.get(kp, 0) + 1
        if kp > maxkp[0]: maxkp = (kp, n)
        cnt += 1
    P(f"count {cnt}; law 'first deficit U(3m) < 3m iff m = 3 mod 4' holds in {law_ok}/{law_n} cases; "
      f"per class m mod 4: {{1: {class_count[1]}, 3: {class_count[3]}}} (members, k=1)")
    ks = sorted(hist_k)
    P("stopping-time k (odd steps to a deficit below n) distribution: " + ", ".join(f"{k}:{hist_k[k]}" for k in ks[:25]) + (" ..." if len(ks) > 25 else ""))
    tail = sum(c for k, c in hist_k.items() if k > 1)
    P(f"  share with k > 1: {tail/cnt:.4f}; k > 10: {sum(c for k,c in hist_k.items() if k>10)/cnt:.5f}; k > 50: {sum(c for k,c in hist_k.items() if k>50)/cnt:.6f}; max k {maxk}")
    kps = sorted(hist_kp)
    P("3-parent rank k' (odd steps until the deficit's smallest 3-multiple preimage is below n): " + ", ".join(f"{k}:{hist_kp[k]}" for k in kps[:25]) + (" ..." if len(kps) > 25 else ""))
    P(f"  share with k' > 1: {sum(c for k,c in hist_kp.items() if k>1)/cnt:.4f}; k' > 10: {sum(c for k,c in hist_kp.items() if k>10)/cnt:.5f}; max k' {maxkp}")
    P(f"  [{time.time()-t0:.0f}s]")

    # ---------------- Part D ----------------
    P("# Part D: depth-2 joins U^2(ab) = U(a): b = [2^A (3a+1) - 3 - 2^alpha]/(9a), checked by direct iteration")
    for a in range(3, 52, 2):
        sols = []
        for alpha in range(1, 14):
            for A in range(alpha + 1, alpha + 61):
                num = (1 << A) * (3 * a + 1) - 3 - (1 << alpha)
                if num <= 0 or num % (9 * a):
                    continue
                b = num // (9 * a)
                if b % 2 == 0 or b <= 1:
                    continue
                ab = a * b
                v1 = 3 * ab + 1; al = (v1 & -v1).bit_length() - 1
                if al != alpha:
                    continue
                if U(U(ab)) == U(a):
                    sols.append((alpha, A, b))
        if sols:
            sols.sort(key=lambda t: t[2])
            P(f"  a={a}: {len(sols)} solutions with alpha<=13, A<=alpha+60; smallest b: {sols[0]} ; next: {sols[1:3]}")
        else:
            P(f"  a={a}: none in range")
    P(f"  [{time.time()-t0:.0f}s]")

    # ---------------- Part E ----------------
    P("# Part E: Applegate--Lagarias multiplier identity (x = -1 mod 2^k, odd j <= k, m_j = (2^j+1)/3)")
    ok = 0; tot = 0; ex = []
    for k in range(3, 16):
        for j in range(1, k + 1, 2):
            m = ((1 << j) + 1) // 3
            if (3 * m - 1) != (1 << j):
                raise RuntimeError
            for t in range(1, 6):
                x = (t << k) - 1
                lhs = U(m * x)
                rhs = x + (x + 1) // (1 << j)
                while rhs % 2 == 0:
                    rhs //= 2
                tot += 1; ok += (lhs == rhs)
                if len(ex) < 4:
                    ex.append((k, j, m, x, lhs, rhs))
    P(f"  U(m_j x) = oddpart(x + (x+1)/2^j) in {ok}/{tot} checks; m_j = (2^j+1)/3 satisfies 3 m_j - 1 = 2^j (the 3x-1 trunk); examples {ex}")
    P("  descent factor of the one-step image: (x + (x+1)/2^j)/x -> 1 + 2^-j: the multiplier resets the binary tail (x' = -1 mod 2^(k-j) only).")

    # ---------------- Part F ----------------
    P("# Part F: shadow mismatch (S1 with m = 1, h = 1): least admissible e, r, r/n")
    def shadow(n, m, h):
        # prefix of length m from n (positive odd)
        d = n; A = 0; B = 0
        for _ in range(m):
            d = U(d)
        # exact affine data: d = (3^m n + B)/2^A
        P3 = 3 ** m
        # recover A, B from the actual word
        x = n; A = 0; B = 0
        for _ in range(m):
            y = 3 * x + 1; a = (y & -y).bit_length() - 1
            A += a; x = y >> a
        B = (1 << A) * d - P3 * n
        mod = 3 ** (m + h + 1)
        target = (3 * d + 1) % mod
        e = 1; p = 4 % mod
        while p != target:
            e += 1; p = (p * 4) % mod
            if e > mod: return None
        while not (e >= 2 and 2 * e > ((3 * d + 1) & -(3 * d + 1)).bit_length() - 1):
            e += 3 ** (m + h)
        r = ((1 << A) * ((1 << (2 * e)) - 1) - 3 * B) // (3 * P3)
        return e, r
    for n in (27, 97, 255, 999, 7):
        res = shadow(n, 1, 1)
        if res:
            e, r = res
            P(f"  n={n}: e={e}, r has {r.bit_length()} bits, r/n ~ 2^{math.log2(r/n):.1f}; r = n mod 3: {r % 3 == n % 3}")
    P(f"  [{time.time()-t0:.0f}s total]")
    with open(__file__.replace(".py", ".out"), "w", encoding="utf-8") as fh:
        fh.write("\n".join(out) + "\n")

if __name__ == "__main__":
    main()
