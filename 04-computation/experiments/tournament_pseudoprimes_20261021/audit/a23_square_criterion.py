# AUDIT a23 (independent auditor, 2026-10-21): full Korselt-type criterion for n whose odd part is cube-free,
# derived from THM-4746 Thm 1 (mod p) and the auditor's mod-p^2 lemma (m >= 3):
#   p || n : need p   | D(r),  r = ((n-1) mod (p-1)) + 1;
#   p^2 || n: need p^2 | D(r),  r = ((n-3) mod p(p-1)) + 3   (so 3 <= r <= p(p-1)+2);   2-part automatic.
# Checked against direct s(n) mod n for every composite n <= 2500 whose odd part is cube-free and not squarefree.
import subprocess
from sympy import factorint
direct = {}
for f in ["a2_direct_4_2000.out", "a2_direct_2001_2500.out"]:
    for line in open(f):
        n, ok, *_ = line.split(); direct[int(n)] = int(ok)
cache = {}
def Dmod(r, M):
    key = (r, M)
    if key not in cache:
        out = subprocess.run(["./a2_smod", "single", str(r), str(M)], capture_output=True, text=True).stdout
        sv = int(out.split("s(n) mod M = ")[1].split()[0]); cv = int(out.split("c(n) mod M = ")[1].split()[0])
        cache[key] = (sv - cv) % M
    return cache[key]
agree = dis = 0; examples = []
for n in range(4, 2501):
    if n not in direct: continue
    f = factorint(n)
    odd = {p: e for p, e in f.items() if p > 2}
    if any(e > 2 for e in odd.values()) or all(e == 1 for e in odd.values()): continue
    ok = True
    for p, e in odd.items():
        if e == 1:
            r = (n - 1) % (p - 1) + 1; ok &= (Dmod(r, p) == 0)
        else:
            L = p * (p - 1); r = (n - 3) % L + 3; ok &= (Dmod(r, p * p) == 0)
    if ok == bool(direct[n]): agree += 1
    else: dis += 1; print("DISAGREE", n, f)
    if ok and len(examples) < 40: examples.append(n)
print(f"composite n <= 2500, odd part cube-free but not squarefree: criterion agrees {agree}, disagrees {dis}")
print("passing examples:", examples)
