#!/usr/bin/env python3
"""hyp_w.py (Q-10, Q-11): for strong D of orders 3-6, is phi(D) = h(D) + P_D(2) = h(T) + 2^k B_t(T) (k >= 1) for some
T = D + w, where the new vertex w has out-neighbourhood X in D (all 2^n choices), some root t of T?  Also with |X| = 1."""
exec(open('nearjoin_check.py').read().split('# (1)')[0])
def rev(T): n = len(T); return [[T[j][i] for j in range(n)] for i in range(n)]
def Bvec(T):
    n = len(T); N = 1 << n; dp, H = ham_all(T); dpr, _ = ham_all(rev(T))
    return H[N-1], [sum(dpr[m][t] * H[(N-1) & ~m] for m in range(1, N) if (m >> t) & 1) for t in range(n)]
def P2(D):
    n = len(D); N = 1 << n; dp, H = ham_all(D); return sum(H[S] * H[(N-1) & ~S] for S in range(N)), H[N-1]
tot = ok = ok1 = 0; fails = []; wits = []
for n in range(3, 7):
    for D in classes(n):
        if not strong(D): continue
        P, h = P2(D); v = h + P; tot += 1; found = False; found1 = False; wit = None
        for X in range(1 << n):
            T = [r[:] + [0 if (X >> i) & 1 else 1] for i, r in enumerate(D)] + [[1 if (X >> j) & 1 else 0 for j in range(n)] + [0]]
            hT, BT = Bvec(T)
            for t, b in enumerate(BT):
                k = 1
                while hT + (b << k) <= v:
                    if hT + (b << k) == v:
                        found = True
                        if bin(X).count('1') == 1: found1 = True
                        if wit is None: wit = (bin(X).count('1'), hT, b, k, t == n)
                    k += 1
        ok += found; ok1 += found1
        if not found: fails.append((n, h, P, v))
        else: wits.append((n, v) + wit)
print(f"strong D of orders 3-6: {tot}; phi(D) rewritten as a near-join over some D + w: {ok}; with w a sole exit: {ok1}")
print("failures (n, h, P, phi):", fails)
print("witnesses (n, phi, |X|, h(T), B, k, root is w):", wits[:20])
