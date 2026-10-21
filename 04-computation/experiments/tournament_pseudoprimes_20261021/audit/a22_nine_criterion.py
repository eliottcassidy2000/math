# AUDIT a22 (independent auditor, 2026-10-21): consequence of the auditor's mod-p^2 lemma at p = 3
# (s(m+6) = -s(m) mod 9 for m >= 3): 9 | D(n) for every n >= 3 with n = 0 mod 3, so for n = 9 m', 3 !| m',
# m' with squarefree odd part, the 9-part never obstructs and
#     n is a Camion pseudoprime  <=>  p | D(r_p) for every prime p >= 5 dividing n.
# Checked against direct s(n) mod n for every such n <= 2500.
from sympy import factorint
D = {}
for line in open("a1_exact_D.txt"):
    r, d = line.split(); D[int(r)] = int(d)
direct = {}
for f in ["a2_direct_4_2000.out", "a2_direct_2001_2500.out"]:
    for line in open(f):
        n, ok, *_ = line.split(); direct[int(n)] = int(ok)
s9 = {}
for line in open("a2_seq_mod45.txt"):
    m, v = map(int, line.split()); s9[m] = v % 9
c9 = lambda m: (1 if m == 1 else (pow(2, (m - 1) // 2, 9) if m % 2 else 0))
print("9 | D(m) for all 3 <= m <= 1401 with 3 | m:", all((s9[m] - c9(m)) % 9 == 0 for m in range(3, 1402, 3)))
agree = dis = 0
for n in range(9, 2501, 9):
    f = factorint(n)
    if f.get(3) != 2 or any(e > 1 for p, e in f.items() if p > 3):
        continue
    pred = all(D[(n - 1) % (p - 1) + 1] % p == 0 for p in f if p >= 5)
    if pred == bool(direct[n]): agree += 1
    else: dis += 1; print("DISAGREE", n)
print(f"n <= 2500 with 9 || n and the rest squarefree (odd and even): criterion without the 9-part agrees {agree}, disagrees {dis}")
