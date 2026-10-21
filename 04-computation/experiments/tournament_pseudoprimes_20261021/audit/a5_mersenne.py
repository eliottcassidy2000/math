# AUDIT a5 (independent auditor, 2026-10-21): which 2^k - 1 (2 <= k <= 64) are Camion pseudoprimes?
# For each odd prime p | n: r_p = ((n-1) mod (p-1)) + 1 and the NECESSARY condition p | D(r_p) (valid for every n,
# squarefree or not, by THM-4746 Thm 1).  D(r) mod p from exact values (r <= 200) or from a2_smod (r <= RMAX).
# Squarefree n passing at every p => PASS (criterion).  Non-squarefree n passing at every p => needs the p^2 part:
# decided here by DIRECT computation of s(n) mod n when n <= 262143 (a2_smod single), else reported 'open'.
import subprocess
from sympy import factorint, isprime

D = {}
for line in open("a1_exact_D.txt"):
    r, d = line.split(); D[int(r)] = int(d)
RMAX = 200000
direct_known = {}
try:
    for line in open("a2_single_262143.out"):
        if line.startswith("n=262143"):
            direct_known[262143] = int(line.strip().split("pass=")[1][0])
except FileNotFoundError:
    pass

def Dmod(r, p):
    if r in D:
        return D[r] % p
    if r <= RMAX:
        out = subprocess.run(["./a2_smod", "single", str(r), str(p)], capture_output=True, text=True).stdout
        sv = int(out.split("s(n) mod M = ")[1].split()[0]); cv = int(out.split("c(n) mod M = ")[1].split()[0])
        return (sv - cv) % p
    return None

for k in range(2, 65):
    n = 2 ** k - 1
    if isprime(n):
        continue
    f = factorint(n)
    sqf = all(e == 1 for e in f.values())
    verdict = None; info = []
    undecided = False
    for p in sorted(f):
        r = (n - 1) % (p - 1) + 1
        d = Dmod(r, p)
        info.append((p, f[p], r, None if d is None else d == 0))
        if d is None:
            undecided = True
        elif d != 0:
            verdict = f"fail at {p}"; break
    if verdict is None:
        if undecided:
            verdict = "undecided (residue too large)"
        elif sqf:
            verdict = "PASS (criterion)"
        else:
            if n <= 4095:
                out = subprocess.run(["./a2_smod", "single", str(n)], capture_output=True, text=True).stdout
                verdict = ("PASS" if "pass=1" in out else "fail") + " (direct; non-squarefree)"
            elif n in direct_known:
                verdict = ("PASS" if direct_known[n] else "fail") + " (direct; non-squarefree)"
            else:
                verdict = "passes every prime; p^2 part OPEN (non-squarefree)"
    print(k, "squarefree" if sqf else "NON-squarefree", verdict, info)
