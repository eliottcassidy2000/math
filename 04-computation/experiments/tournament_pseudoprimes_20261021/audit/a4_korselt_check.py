# AUDIT a4 (independent auditor, 2026-10-21): THM-4746 statement 4 (Korselt-type criterion) against DIRECT
# computation of s(n) mod n (a2_smod direct, file a2_direct_4_2000.out) for every composite 4 <= n <= 2000.
# Criterion: odd part squarefree  =>  [s(n) = c(n) mod n  <=>  p | D(r_p) for all odd p | n], r_p = ((n-1) mod (p-1)) + 1.
# Uses the exact D(r) of a1_exact_D.txt (r <= 200 suffices here since r_p <= p-1 and r_p <= ...; checked below).
import sys
from sympy import factorint, isprime

D = {}
for line in open("a1_exact_D.txt"):
    r, d = line.split(); D[int(r)] = int(d)
direct = {}
for line in open("a2_direct_4_2000.out"):
    n, ok, sv, cv = line.split(); direct[int(n)] = int(ok)

def Dmod(r, p):
    if r in D:
        return D[r] % p
    raise ValueError(r)

def criterion(n):
    f = factorint(n); odd = [p for p in f if p > 2]
    sqf = all(f[p] == 1 for p in odd)
    conds = {}
    for p in odd:
        r = (n - 1) % (p - 1) + 1
        conds[p] = (r, Dmod(r, p) == 0)
    return sqf, all(v[1] for v in conds.values()), conds

agree = dis = 0; nsq_total = 0; nsq_pass = 0; nec_viol = 0
sqf_odd = sqf_even = 0
nonsqf_pass = []
for n in sorted(direct):
    sqf, pred, conds = criterion(n)
    if sqf:
        if n % 2: sqf_odd += 1
        else: sqf_even += 1
        if pred == bool(direct[n]): agree += 1
        else:
            dis += 1; print("DISAGREE", n, pred, direct[n], conds)
    else:
        nsq_total += 1
        if direct[n]:
            nsq_pass += 1; nonsqf_pass.append(n)
            if not pred:
                nec_viol += 1; print("NECESSITY VIOLATED", n, conds)
print(f"composite n in [4,2000] with squarefree odd part: {agree + dis} ({sqf_odd} odd, {sqf_even} even): "
      f"criterion agrees {agree}, disagrees {dis}")
print(f"composite n with non-squarefree odd part: {nsq_total}; direct passes {nsq_pass}; "
      f"necessary condition violated by a passing n: {nec_viol}")
print("passing n <= 2000 with NON-squarefree odd part:", nonsqf_pass)
odd_pass_1200 = [n for n in sorted(direct) if n % 2 and n <= 1200 and direct[n]]
print(f"odd composite Camion pseudoprimes <= 1200 (direct, {len(odd_pass_1200)}):", odd_pass_1200)
# compare with the session's camion_psp_1200.out (odd composites with camion=1)
sess = []
for line in open("../camion_psp_1200.out"):
    if line.startswith("n=") and "camion=1" in line:
        sess.append(int(line.split()[0][2:]))
print("session odd list (camion_psp_1200.out) identical:", sess == odd_pass_1200)
even_pass = [n for n in sorted(direct) if n % 2 == 0 and direct[n]]
print("even composite passes <= 400 :", [n for n in even_pass if n <= 400])
print("composites <= 600 passing (all):", [n for n in sorted(direct) if n <= 600 and direct[n]])
# specific cases requested
for n in [21, 1023, 561, 1105, 1729, 341, 85, 25, 27, 9, 169, 49, 81, 121, 125, 63, 186, 195, 399, 1443, 1239]:
    if n in direct:
        sqf, pred, conds = criterion(n)
        print(f"n={n}: direct pass={direct[n]}  odd-part squarefree={sqf}  residues {conds}")
