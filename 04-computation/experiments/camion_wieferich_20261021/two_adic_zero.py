# 2-adic zero of the strong-tournament class sum on the class n = 2 (mod 3).
# f(n) = s(n)/n! = sum over strong classes of 1/|Aut| in Z_2.  Observed: v2(f(n)) = 2 + v2(n - n*) for n = 2 mod 3, n >= 5,
# with a 2-adic integer n*.  This script computes v2(f(n)) for n <= NMAX and reads off n* mod 2^k: the unique residue r mod 2^k
# such that every tested n = 2 (mod 3), n >= 5, has v2(f(n)) >= 2 + k iff n = r (mod 2^k); it reports the consistency of the law.
# mac-mini-2026-10-21-fanopsp.
import sys
NMAX = int(sys.argv[1]) if len(sys.argv) > 1 else 1500
N = NMAX + 80; MOD = 1 << N
def v2(x): return (x & -x).bit_length() - 1 if x else 10**9
s = [0] * (NMAX + 1); fv = [0] * (NMAX + 1); row = [1]
for n in range(1, NMAX + 1):
    row = [1] + [(row[i - 1] + row[i]) % MOD for i in range(1, n)] + [1]
    acc = (1 << (n * (n - 1) // 2)) if n * (n - 1) // 2 < N else 0
    for k in range(max(1, n - 60), n):                 # 2^C(n-k,2) >= 2^N beyond this
        e = (n - k) * (n - k - 1) // 2
        if e < N: acc -= row[k] * s[k] << e
    s[n] = acc % MOD; fv[n] = fv[n - 1] + v2(n)
E = {n: v2(s[n]) - fv[n] for n in range(5, NMAX + 1) if n % 3 == 2}
# law: E(n) = 2 + v2(n - n*).  Determine n* mod 2^k greedily from the data.
r, k = 0, 0
while True:
    cands = [c for c in (r, r + (1 << k)) if all((E[n] >= 3 + k) == ((n - c) % (1 << (k + 1)) == 0) for n in E if E[n] >= 2 + k)]
    if len(cands) != 1: break
    r = cands[0]; k += 1
    if k > 40: break
print(f"n <= {NMAX}: {len(E)} values n = 2 (mod 3); n* determined mod 2^{k}: n* = {r} (mod 2^{k})  [binary {bin(r)}]")
viol = [(n, E[n], 2 + v2(n - r)) for n in E if E[n] != 2 + min(v2(n - r), k)] if k else []
print("violations of E(n) = 2 + v2(n - n*) (with v2 capped at k):", viol[:10], "count", len(viol))
print("largest valuations:", sorted(((E[n], n) for n in E), reverse=True)[:8])
