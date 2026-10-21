# AUDIT a17b (independent auditor, 2026-10-21): extended test of the auditor's lemma
#    s(m + phi(p^k)) = eps * s(m)  (mod p^k),  eps = (2/p),  for all m >= m0(k)
# (proved below for k = 2 with m0 = 3 in the audit report), for k = 2, p <= 31 and k = 3, p <= 7.
import subprocess
def smod_seq(N, M):
    out = subprocess.run(["./a2_smod", "seq", str(N), str(M)], capture_output=True, text=True).stdout
    s = [0] * (N + 1)
    for line in out.split("\n"):
        if line.strip():
            m, v = map(int, line.split()); s[m] = v
    return s
for k, plist in [(2, [3, 5, 7, 11, 13, 17, 19, 23, 29, 31]), (3, [3, 5, 7]), (4, [3])]:
    for p in plist:
        M = p ** k; L = p ** (k - 1) * (p - 1); eps = 1 if pow(2, (p - 1) // 2, p) == 1 else -1
        N = 3 * L + 10
        s = smod_seq(N + L, M)
        bad = [m for m in range(1, N + 1) if (s[m + L] - eps * s[m]) % M]
        print(f"k={k} p={p}: modulus {M}, shift phi(p^k)={L}, eps={eps:+d}: m in [1,{N}] violating s(m+L) = eps s(m): {bad}")
