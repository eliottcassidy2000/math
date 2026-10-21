# AUDIT a18 (independent auditor, 2026-10-21): THM-4747 checks on the auditor's circulant data a11_circ_all.out
# (all multiplier classes of Cay(Z_n,S), odd 3 <= n <= 23) + comparison with the session's circ_*.out.
import re
from collections import Counter, defaultdict
from sympy import isprime, factorint
rows = []
for line in open("a11_circ_all.out"):
    m = re.match(r"n=(\d+) S=\{([\d,]+)\} mu=(\d+) A=(\d+) H=(\d+)(?: Hfull=(\d+))? C=(\d+) rho0=(\d+)", line)
    if m:
        n, S, mu, A, H, Hf, C, R = m.groups()
        rows.append((int(n), tuple(map(int, S.split(","))), int(mu), int(A), int(H), int(C), int(R)))
per_n = Counter(r[0] for r in rows)
print("multiplier classes per odd n:", dict(sorted(per_n.items())), " total", len(rows))
viol = []
naive = defaultdict(lambda: [0, 0])
for n, S, mu, A, H, C, R in rows:
    if H % n: viol.append(("H mod n", n, S))
    if R % 2 == 0: viol.append(("rho0 even", n, S))
    if isprime(n):
        if H % (n * mu): viol.append(("H mod n mu", n, S))
        if (C - (n - 1) // 2) % (n * mu): viol.append(("C cong", n, S))
        if (R - mu) % (2 * mu): viol.append(("rho0 cong", n, S))
        if mu % 2 == 0: viol.append(("mu even", n, S))
    else:
        naive[n][1] += 1
        if (C - A) % n == 0: naive[n][0] += 1
print("violations of THM-4747 (2),(3), H = 0 mod n, rho0 odd:", viol)
print("composite n: classes with C = A (mod n):", {n: f"{a}/{b}" for n, (a, b) in sorted(naive.items())})
# compare with session data (multisets of (H, C, rho0, mu, A) per n)
sess = []
for f in ["circ_small.out", "circ_13_17.out", "circ_19.out", "circ_21.out", "circ_23.out"]:
    for line in open("../" + f):
        m = re.match(r"n=(\d+) S=\{([\d,]+)\} mult=(\d+) H=(\d+) H/n=(\d+) C=(\d+) rho0=(\d+) A=(\d+)", line)
        if m:
            n, S, mu, H, h0, C, R, A = m.groups()
            sess.append((int(n), int(H), int(C), int(R), int(mu), int(A)))
mine = Counter((n, H, C, R, mu, A) for n, S, mu, A, H, C, R in rows)
print("session circ_*.out (without paley19/23 duplicates): rows", len(sess), "; multiset equal to auditor's:", Counter(sess) == mine)
# circulant maxima per n and Paley data
for n in sorted(per_n):
    best = max((r for r in rows if r[0] == n), key=lambda r: r[4])
    rot = tuple(range(1, (n - 1) // 2 + 1))
    import math
    print(f"n={n}: max H over circulants {best[4]} at S={best[1]} (rotational: {best[1] == rot}); ratio H*2^(n-1)/n! = {best[4] * 2 ** (n - 1) / math.factorial(n):.4f}")
for q in [3, 7, 11, 19, 23]:
    QR = tuple(sorted({x * x % q for x in range(1, q)}))
    # find the class row isomorphic to P_q: mu = (q-1)/2
    r = [x for x in rows if x[0] == q and x[2] == (q - 1) // 2]
    assert len(r) == 1
    n, S, mu, A, H, C, R = r[0]
    aut = q * (q - 1) // 2; t = (C // ((q - 1) // 2) - 1) // q
    assert C == ((q - 1) // 2) * (1 + q * t)
    tc = H // aut
    print(f"P_{q}: H={H} C={C} tc=H/(q(q-1)/2)={tc} {factorint(tc)} tc mod q = {tc % q}; rho0={R}; t={t}; "
          f"rotational H={[x[4] for x in rows if x[0] == q and x[1] == tuple(range(1, (q - 1) // 2 + 1))]}")
