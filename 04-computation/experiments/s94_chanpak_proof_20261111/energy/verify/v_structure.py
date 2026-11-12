# v_structure.py - independent checks of Lemmas 1-4, Cor. 3a/3b (verifier, lane energy)
# Independent implementation: q by its recursion, everything else from the definitions. Exact integers.
import sys, random
import numpy as np
from numba import njit

@njit(cache=True)
def q_table(n):
    N = 1 << n
    q = np.zeros(N, dtype=np.int64)
    for h in range(1, N):
        p = h >> 1
        if h & 1:
            q[h] = q[p] + (p + 1) // 2      # ceil(p/2)
        else:
            q[h] = q[p]
    return q

def kappa(d):
    return 2 if d == 1 else 1 << (d - 1)

def bits(h, n):
    return [(h >> i) & 1 for i in range(n)]

def rev(h, n):
    r = 0
    for i in range(n):
        if (h >> i) & 1:
            r |= 1 << (n - 1 - i)
    return r

def main():
    NMAX = 13
    q = q_table(18)
    eta = lambda h: 2 * h + 1 - 2 * int(q[h])
    # ---- Lemma 1
    ok1 = True
    for h in range(1 << NMAX):
        n = max(h.bit_length(), 1)
        b = bits(h, n)
        tot = 0
        for i in range(n):
            if b[i]:
                for j in range(i + 1, n):
                    if b[j]:
                        tot += kappa(j - i)
        if tot != 2 * int(q[h]):
            ok1 = False; print("L1 FAIL", h); break
    print("Lemma 1 (2q = sum kappa h_i h_j), all h < 2^%d: %s" % (NMAX, ok1))
    # ---- Lemma 3: 2^n (eta - 1 - rho) == h c(h),  rho exact as integer / 2^n ; bounds 0 < rho < 3b/4 <= 3(n+1)/8
    ok3 = True; okb = True; maxratio = 0.0; minrho = 1e9
    for n in range(1, 15):
        for h in range(1 << (n - 1), 1 << n):
            b = bits(h, n)
            num = 0  # rho * 2^n
            blocks = 0
            for i in range(n):
                if b[i]:
                    t = (1 << (n - 1))            # 1/2
                    if i >= 1 and b[i - 1]:
                        t -= 3 * (1 << (n - 2)) if n >= 2 else 0
                    else:
                        blocks += 1
                    for j in range(0, i - 1):
                        if b[j]:
                            t += 1 << (n + j - i - 1)
                    num += t
            c = (1 << (n + 1)) - rev(h, n)
            if (eta(h) - 1) * (1 << n) - num != h * c:
                ok3 = False; print("L3 FAIL", n, h); break
            # bounds: 0 < num, 4 num < 3 blocks 2^n, blocks <= (n+1)/2
            if not (num > 0 and 4 * num < 3 * blocks * (1 << n) and 2 * blocks <= n + 1):
                okb = False; print("L3 bound FAIL", n, h, num, blocks); break
            maxratio = max(maxratio, num / (1 << n) / (3 * (n + 1) / 8))
            minrho = min(minrho, num / (1 << n))
    print("Lemma 3 identity (exact), all n <= 14: %s ; bounds 0 < rho < 3b/4 <= 3(n+1)/8: %s ; max rho/(3(n+1)/8) = %.4f ; min rho = %.3g" % (ok3, okb, maxratio, minrho))
    # ---- Corollary 3a on all colliding pairs n = 12: |h c - h' c'| < 2^n 3(n+1)/8
    n = 12
    vals = {}
    for h in range(1 << (n - 1), 1 << n):
        vals.setdefault(eta(h), []).append(h)
    ok3a = True; cnt = 0; worst = 0.0
    for v, L in vals.items():
        for a in range(len(L)):
            for b_ in range(a + 1, len(L)):
                h, hp = L[a], L[b_]
                d = abs(h * ((1 << (n + 1)) - rev(h, n)) - hp * ((1 << (n + 1)) - rev(hp, n)))
                cnt += 1
                worst = max(worst, d / (1 << n))
                if not (8 * d < 3 * (n + 1) * (1 << n)):
                    ok3a = False
    print("Cor. 3a on all %d colliding pairs at n = 12: %s (max |Phi - Phi'| = %.3f, bound %.3f)" % (cnt, ok3a, worst, 3 * (n + 1) / 8))
    # ---- Lemma 4: eta(2^(n-1)+r) == eta(2^(n-1)+rc(r)) for all n <= 18
    ok4 = True
    for n in range(2, 19):
        m = n - 1
        for r in range(1 << m):
            rc = rev((~r) & ((1 << m) - 1), m)
            if eta((1 << m) + r) != eta((1 << m) + rc):
                ok4 = False; print("L4 FAIL", n, r); break
        if not ok4: break
    print("Lemma 4 involution, all n <= 18: %s" % ok4)
    # ---- Lemma 2 on random pairs and both lambda formulas; Remark (i) bound |lambda_j - X/2^(j+1) - Xi/2^(w-j)| < 3
    random.seed(20261111)
    ok2 = True; okl = True; okr = True
    for trial in range(100000):
        n = random.randint(2, 17)
        h = random.randrange(1 << (n - 1), 1 << n)
        hp = random.randrange(1 << (n - 1), 1 << n)
        if h == hp: continue
        bh = bits(h, n); bp = bits(hp, n)
        eps = [bh[i] - bp[i] for i in range(n)]
        S = [i for i in range(n) if eps[i] != 0]
        D = h - hp
        L0x2 = 4 * D
        for a in range(len(S)):
            for b_ in range(a + 1, len(S)):
                i, j = S[a], S[b_]
                L0x2 -= kappa(j - i) * (eps[i] + eps[j])
        tot2 = L0x2
        for j in range(n):
            if eps[j] == 0 and bh[j]:
                lam = sum(kappa(abs(i - j)) * eps[i] for i in S)
                Xgt = sum(eps[i] * (1 << (i - j - 1)) for i in range(j + 1, n))
                Xilt = sum(eps[i] * (1 << (j - 1 - i)) for i in range(0, j))
                lam2 = Xgt + Xilt + (eps[j + 1] if j + 1 < n else 0) + (eps[j - 1] if j >= 1 else 0)
                if lam != lam2: okl = False
                tot2 -= 2 * lam
        if tot2 != 2 * (eta(h) - eta(hp)):
            ok2 = False; print("L2 FAIL", n, h, hp); break
        # remark (i): window coordinates
        s = S[0]; top = S[-1]; w = top - s + 1
        X = sum(eps[s + i] * (1 << i) for i in range(w)); Xi = sum(eps[s + i] * (1 << (w - 1 - i)) for i in range(w))
        for jj in range(1, w - 1):
            j = s + jj
            if eps[j] == 0:
                lam = sum(kappa(abs(i - j)) * eps[i] for i in S)
                dev = abs(lam - X / 2 ** (jj + 1) - Xi / 2 ** (w - jj))
                if not dev < 3: okr = False
    print("Lemma 2 (100000 random pairs, n <= 17): identity %s ; lambda_j formula %s ; Remark (i) |.| < 3: %s" % (ok2, okl, okr))
    # ---- Cor. 3b: B(K,eps) = [D(K)D*(eps)+D*(K)D(eps)]/2^n + E, |E| < 3 min(supp)
    ok3b = True; worstE = 0.0
    for trial in range(100000):
        n = random.randint(2, 24)
        K = [0] * n; eps = [0] * n
        for i in range(n):
            r = random.random()
            if r < 0.3: K[i] = 1
            elif r < 0.6: eps[i] = random.choice([-1, 1])
        B = 0
        for j in range(n):
            if K[j]:
                for i in range(n):
                    if eps[i]:
                        B += kappa(abs(i - j)) * eps[i]
        DK = sum(K[i] << i for i in range(n)); DsK = sum(K[i] << (n - 1 - i) for i in range(n))
        De = sum(eps[i] * (1 << i) for i in range(n)); Dse = sum(eps[i] * (1 << (n - 1 - i)) for i in range(n))
        # E * 2^n exactly
        E2 = B * (1 << n) - (DK * Dse + DsK * De)
        mn = min(sum(K), sum(1 for e in eps if e))
        if mn == 0:
            if E2 != 0: ok3b = False
        else:
            if not abs(E2) < 3 * mn * (1 << n): ok3b = False
            worstE = max(worstE, abs(E2) / (1 << n) / mn)
    print("Cor. 3b (100000 random K, eps, n <= 24): %s (max |E|/min supp = %.3f)" % (ok3b, worstE))

if __name__ == "__main__":
    main()
