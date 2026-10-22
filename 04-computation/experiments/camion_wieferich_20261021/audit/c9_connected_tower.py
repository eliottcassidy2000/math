# c9_connected_tower.py -- second-pass audit of the remark "each Hurwitz-closed structure has its own lift and therefore its own
# list" (THM-4750 5b, connected labelled graphs).  For connected graphs the circulant target is 2^((n-1)/2) - 1 (the empty
# connection set is disconnected), and with E(n) = b(n) - 2^((n-1)/2) + 1 the Session25 law b(n + phi(p^r)) = eps b(n) (mod p^r),
# n >= 2r, gives the DRIFTED tower  E(p^r) = eps E(p^(r-1)) + (1 - eps)  (mod p^r)  (r >= 2, p^(r-1) >= 2r).
# For eps = (2/p) = -1 this forces E(p^2) = 2 - E(p) = 2 (mod p): p^2 is then never a "connected Camion pseudoprime",
# whatever p^2 | E(p) says.  Both listed "connected-CW" primes 5 and 197 have (2/p) = -1.  Exact integers, n <= 3^5, 5^3, 7^2, 11^2, 13^2.
# Second-pass auditor, 2026-10-21.
from math import comb
N = 243
g = [2 ** (n * (n - 1) // 2) for n in range(N + 1)]; b = [0] * (N + 1)
for n in range(1, N + 1): b[n] = g[n] - sum(comb(n - 1, k - 1) * b[k] * g[n - k] for k in range(1, n))
E = lambda n: b[n] - 2 ** ((n - 1) // 2) + 1
for p in (3, 5, 7, 11, 13):
    eps = 1 if p % 8 in (1, 7) else -1
    ks = [k for k in range(1, 7) if p ** k <= N]
    row = []
    for k in ks:
        n = p ** k; e = E(n); v = 0; x = e
        while x % p == 0 and v < 12: x //= p; v += 1
        row.append(f"E({p}^{k}) mod {p}^{k} = {e % n} (v_p = {v})")
    drift = []
    for k in ks[1:]:
        n = p ** k
        if p ** (k - 1) >= 2 * k:
            drift.append(f"r={k}: E(p^r) - eps E(p^(r-1)) - (1 - eps) = 0 mod p^r: {(E(n) - eps * E(n // p) - (1 - eps)) % n == 0}")
    print(f"p={p:2d} eps={eps:+d}: " + "; ".join(row) + ("  | drifted tower " + "; ".join(drift) if drift else ""))
