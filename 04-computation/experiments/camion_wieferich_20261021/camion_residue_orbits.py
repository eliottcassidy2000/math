# Combinatorial reading of the Camion residue (THM-4750 §5b): rho_C(2) = -eps * ns(p)/p (mod p), where ns(p) = 2^C(p,2) - s(p)
# counts non-strong labelled tournaments and ns(p)/p is the number of their rotation orbits on Z_p.  Checked against
# camion_residue_2000.out.  mac-mini-2026-10-21-fanopsp.
from math import comb
import re
rho = {}
for line in open('camion_residue_2000.out'):
    m = re.match(r'p=(\d+) eps=([+-]1) rho_C=(\d+)', line)
    if m: rho[int(m.group(1))] = (int(m.group(2)), int(m.group(3)))
def s_mod(N, q):
    T = [pow(2, j * (j - 1) // 2, q) for j in range(N + 1)]; s = [0] * (N + 1); row = [1]
    for n in range(1, N + 1):
        row = [1] + [(row[i - 1] + row[i]) % q for i in range(1, n)] + [1]
        s[n] = (T[n] - sum(row[k] * s[k] * T[n - k] for k in range(1, n))) % q
    return s
bad = 0; checked = 0
for p in [q for q in sorted(rho) if q < 400]:
    q = p * p; s = s_mod(p, q); ns = (pow(2, comb(p, 2), q) - s[p]) % q
    assert ns % p == 0
    eps, r = rho[p]
    checked += 1; bad += ((-eps * (ns // p)) % p != r)
print(f"rho_C(2) = -eps * ns(p)/p (mod p): {checked} primes 5 <= p < 400, failures {bad}")
