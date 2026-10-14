"""
collatz_irrationality_exponent_zeta_20261014.py -- S31 companion: the owner's hypothesis that the odd zeta values and the
logarithms have irrationality exponent 2, read through THM-4632 (the formal law R_x(s) = 1/(1 - e^(-s)/2 - e^(s(x-1))/2)
with parameter x in (1, 2) has boundary growth exponent 2(mu(x) - 1), and has a pole on its line iff x is rational).
zeta(3) = 1.2020569..., zeta(5) = 1.0369277..., zeta(7) = 1.0083492..., log 2 + 1 = 1.6931..., log 3 = 1.0986..., and
theta = log2 3 all lie in (1, 2).  For each: 2600 partial quotients validated at two precisions, the largest local
exponent e_n = log q_(n+1)/log q_n beyond n = 100 and n = 1000 (mu = 1 + limsup e_n), the Gauss-Kuzmin frequencies, and
the known upper bounds for mu (CITED from memory: Rhin-Viola 2001 mu(zeta(3)) <= 5.513891; Marcovecchio 2009
mu(log 2) <= 3.57455; Salikhov mu(log 3) <= 5.125; zeta(5), zeta(7): irrationality itself open).  NUMERICAL only: a
continued fraction to 2600 terms cannot distinguish mu = 2 from mu = 2 + epsilon; it can only show that no large
partial quotient appears in the computed range.
"""
import math
from mpmath import mp, log, zeta, floor, mpf

def cf(x, N):
    a = []
    for _ in range(N):
        fl = int(floor(x)); a.append(fl)
        frac = x - fl
        if frac == 0: break
        x = 1 / frac
    return a

def analyse(name, make, known):
    mp.dps = 4500; a1 = cf(make(), 2600)
    mp.dps = 5200; a2 = cf(make(), 2600)
    common = 0
    while common < min(len(a1), len(a2)) and a1[common] == a2[common]: common += 1
    a = a1[:common]
    q = [0, 1]
    for n in range(1, len(a)): q.append(a[n] * q[-1] + q[-2])
    qn = q[1:]                                   # q_n, n = 0, 1, ...
    e = {n: math.log(qn[n + 1]) / math.log(qn[n]) for n in range(2, len(qn) - 1)}
    out = []
    for N in (100, 1000):
        sub = {n: v for n, v in e.items() if n >= N}
        nmax = max(sub, key=sub.get); out.append((N, sub[nmax], nmax, a[nmax + 1]))
    gk = [sum(1 for t in a[1:] if t == k) / (len(a) - 1) for k in (1, 2, 3)]
    big = sorted(range(1, len(a)), key=lambda i: -a[i])[:4]
    print(f"   {name:12s}: {common} terms validated; max e_n beyond n=100: {out[0][1]:.5f} (n={out[0][2]}, a={out[0][3]}); beyond n=1000: {out[1][1]:.5f} (n={out[1][2]}, a={out[1][3]}); "
          f"largest partial quotients {[(i, a[i]) for i in sorted(big)]}; Gauss-Kuzmin a_n=1,2,3: {[round(g, 4) for g in gk]} (0.415, 0.170, 0.093); known: {known}")

print("# local exponents e_n = log q_(n+1)/log q_n (mu = 1 + limsup e_n; growth exponent of the formal law 2 e_n) for six parameters in (1, 2)")
analyse("log2 3", lambda: log(3) / log(2), "mu open; Rhin-type effective bound ~14.3 (S29)")
analyse("zeta(3)", lambda: zeta(3), "mu <= 5.513891 (Rhin-Viola 2001)")
analyse("zeta(5)", lambda: zeta(5), "irrationality open")
analyse("zeta(7)", lambda: zeta(7), "irrationality open")
analyse("1 + log 2", lambda: 1 + log(2), "mu(log 2) <= 3.57455 (Marcovecchio 2009); = 2 claimed by arXiv:2610.10192")
analyse("log 3", lambda: log(3), "mu(log 3) <= 5.125 (Salikhov); = 2 claimed by arXiv:2610.10192")
print("   (THM-4632 statement 3: the formal law with parameter zeta(5) has a pole on its line iff zeta(5) is rational; the computed partial quotients show no sign of termination.)")
print("done")
