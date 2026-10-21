# AUDIT a21 (independent auditor, 2026-10-21): Theorem 1 at the primes of 2^26 - 1 = 3 * 2731 * 8191 (and 43, 127 of
# 16383), m <= 3p, plus the residues r_p of n = 67108863 and 16383; and every 2p, 3p <= 2000 passes directly.
import subprocess
from sympy import isprime
def seq(N, M):
    out = subprocess.run(["./a2_smod", "seq", str(N), str(M)], capture_output=True, text=True).stdout
    s = [0] * (N + 1)
    for line in out.split("\n"):
        if line.strip():
            m, v = map(int, line.split()); s[m] = v
    return s
for p in [43, 127, 2731, 8191]:
    eps = 1 if pow(2, (p - 1) // 2, p) == 1 else -1
    s = seq(4 * p, p)
    bad = [m for m in range(1, 3 * p + 1) if (s[m + p - 1] - eps * s[m]) % p]
    print(f"p={p} eps={eps:+d}: violations of s(m+p-1) = eps s(m) for 1<=m<=3p: {len(bad)}; D(3) mod p = {(s[3] - 2) % p}, D(1) = {(s[1] - 1) % p}")
for n, ps in [(67108863, [3, 2731, 8191]), (16383, [3, 43, 127])]:
    print(n, "residues r_p:", {p: (n - 1) % (p - 1) + 1 for p in ps})
direct = {}
for line in open("a2_direct_4_2000.out"):
    a, ok, *_ = line.split(); direct[int(a)] = int(ok)
twop = [2 * p for p in range(3, 1001) if isprime(p)]
threep = [3 * p for p in range(5, 667) if isprime(p)]
print("2p <= 2000 failing directly:", [n for n in twop if not direct[n]], f"(of {len(twop)})")
print("3p <= 2000 failing directly:", [n for n in threep if not direct[n]], f"(of {len(threep)})")
