# AUDIT a17 (independent auditor, 2026-10-21): is s(m) quasi-periodic MODULO p^2 with period phi(p^2) = p(p-1) and
# multiplier u = 2^(p(p-1)/2) mod p^2 (= +-1 mod p^2), as Mani-Stones (SIAM J. Discrete Math. 2016) prove for connected
# labelled graphs?  For each p, report the set of m in [1, M] where s(m+p(p-1)) != u s(m) (mod p^2), the same for
# D = s - c, and whether D is periodic with period 2p(p-1) from some point on.  Also connected graphs (control).
import subprocess
from math import comb
def smod_seq(N, M):
    out = subprocess.run(["./a2_smod", "seq", str(N), str(M)], capture_output=True, text=True).stdout
    s = [0] * (N + 1)
    for line in out.split("\n"):
        if line.strip():
            m, v = map(int, line.split()); s[m] = v
    return s
def conn_mod(N, M):
    g = [pow(2, j * (j - 1) // 2, M) for j in range(N + 1)]; c = [0] * (N + 1); row = [1]
    for m in range(1, N + 1):
        if m > 1: row = [1] + [(row[i - 1] + row[i]) % M for i in range(1, m - 1)] + [1]
        c[m] = (g[m] - sum(row[k - 1] * c[k] * g[m - k] for k in range(1, m))) % M
    return c
for p in [3, 5, 7, 11, 13]:
    M = p * p; L = p * (p - 1); u = pow(2, L // 2, M)
    N = 6 * L + 20
    s = smod_seq(N + L, M)
    cc = lambda m: (1 if m == 1 else (pow(2, (m - 1) // 2, M) if m % 2 else 0))
    bad_s = [m for m in range(1, N + 1) if (s[m + L] - u * s[m]) % M]
    bad_D = [m for m in range(1, N + 1) if ((s[m + L] - cc(m + L)) - u * (s[m] - cc(m))) % M]
    bad_per = [m for m in range(1, N - L + 1) if ((s[m + 2 * L] - cc(m + 2 * L)) - (s[m] - cc(m))) % M]
    print(f"p={p}: u = 2^(p(p-1)/2) mod p^2 = {u}; m<= {N} with s(m+p(p-1)) != u s(m) mod p^2: {bad_s}; same for D: {bad_D}; "
          f"m with D(m+2p(p-1)) != D(m): {bad_per[:10]}")
    if p <= 7:
        cg = conn_mod(N + L, M)
        bad_c = [m for m in range(1, N + 1) if (cg[m + L] - u * cg[m]) % M]
        print(f"      connected graphs: m with cg(m+p(p-1)) != u cg(m) mod p^2: {bad_c}")
