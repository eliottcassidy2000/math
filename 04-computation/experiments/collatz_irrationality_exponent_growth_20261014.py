"""
collatz_irrationality_exponent_growth_20261014.py -- S31 (owner: positive integers as irrationality measures; 1 <-> 2 the
only abnormality; mu(log) = 2 hypothesis; Fisher information).  Computations behind THM-4632 and HYP-9262.

THM-4632: for the trivial character the Collatz transfer function R(s) = 1/(1 - e^(-s)/2 - e^(s(theta-1))/2) on the line
Re s = 0 has growth exponent limsup log|R(2 pi i f)|/log f = 2 (mu(theta) - 1), mu = irrationality exponent of theta = log2 3
(upper bound: the S29 line bound theta^2/(4 theta_f^2) and the definition of mu; lower bound: the S29 peak law at the
convergents, |R| ~ (2 - theta)^2/(2 pi^2 theta_f^2) with theta_f = Q_n theta - P_n ~ 1/Q_(n+1)).  The same holds for the
formal law with any irrational parameter x in (1, 2) in place of theta.

Part 1.  Convergents of log2 3 to n = 2600 (continued fraction validated at two precisions): local exponents
         e_n = log Q_(n+1)/log Q_n (so mu_n = 1 + e_n and the local growth exponent is 2 e_n); the largest e_n beyond
         n = 10, 100, 1000; the same along the tau-coloured convergents (Q odd, P even) only.
Part 2.  Exact peaks: for n = 1..30, the value |R| at the drifted frequency Q_n + Im(w)/(2 pi) (w the exact root of
         e^(-w) + e^(2 pi i theta_f) e^(w(theta-1)) = 2 continued from 0), against (2 - theta)^2/(2 pi^2 theta_f^2) and
         against Q_(n+1)^2; the local exponent log|R_n|/log Q_n against 2 e_n.
Part 3.  Jarnik laws: for mu = 3, 4, 5 a continued fraction with a_(n+1) = ceil(q_n^(mu-2)) (irrationality exponent
         exactly mu), shifted into (1, 2); the exact peaks of the formal law at its convergents and the local exponents
         against 2 (mu - 1) = 4, 6, 8.  mu = 2 control: the golden-type number [1; 2, 1, 1, 1, ...] = 3 - phi.
Part 4.  The characters modulo m up to sign against the cusps of the principal congruence curve X(m): J_2(m)/2 against
         the cusp count; the distribution of the 2599 convergents over the 12 cusp classes of X(5) and over the 6 points
         of P^1(F_5).
Part 5.  Artin primes for the base 2 below 10^5 (one doubling cycle per level in THM-4631): count against Artin's
         constant; Wieferich primes below 10^5.
Part 6.  Fisher information of a line: |R|^2 near a line is a Lorentzian of half-width gamma = |Re w|/(2 pi) in f; the
         Fisher information of the Cauchy location family is 1/(2 gamma^2); table for the first lines; the Dirac limit.
"""
import math
import numpy as np
from mpmath import mp, mpf, mpc, log, exp, floor, pi, findroot, nstr, sqrt
from sympy import primerange, n_order

def cf(x, N):
    a = []
    for _ in range(N):
        fl = int(floor(x)); a.append(fl)
        frac = x - fl
        if frac == 0: break
        x = 1 / frac
    return a

def convergents(a):
    P = [1, a[0]]; Q = [0, 1]
    for n in range(1, len(a)):
        P.append(a[n] * P[-1] + P[-2]); Q.append(a[n] * Q[-1] + Q[-2])
    return [(Q[i + 1], P[i + 1]) for i in range(len(a))]     # (Q_n, P_n), n = 0, 1, ...

print("# Part 1: local exponents e_n = log Q_(n+1)/log Q_n along the convergents of log2 3")
mp.dps = 4500; a1 = cf(log(3) / log(2), 2600)
mp.dps = 5200; a2 = cf(log(3) / log(2), 2600)
common = 0
while common < min(len(a1), len(a2)) and a1[common] == a2[common]: common += 1
a = a1[:common]; conv = convergents(a)
e = [math.log(conv[n + 1][0]) / math.log(conv[n][0]) for n in range(2, len(conv) - 1)]   # n >= 2 so Q_n >= 2
ns = list(range(2, len(conv) - 1))
for N in (10, 100, 1000, 2000):
    sub = [(ns[i], e[i]) for i in range(len(e)) if ns[i] >= N]
    best = max(sub, key=lambda t: t[1])
    print(f"   n >= {N:4d}: max e_n = {best[1]:.6f} at n = {best[0]} (a_(n+1) = {a[best[0]+1]}), so mu_n = {1+best[1]:.6f}, local growth exponent 2 e_n = {2*best[1]:.6f}; mean e_n = {np.mean([t[1] for t in sub]):.6f}")
tau = [(n, e[i]) for i, n in enumerate(ns) if conv[n][0] % 2 == 1 and conv[n][1] % 2 == 0]
print(f"   tau-coloured convergents (Q odd, P even): {len(tau)} of {len(ns)}; max e_n beyond n = 100 among them: {max(t[1] for t in tau if t[0] >= 100):.6f}")
print(f"   first 12 e_n (n = 2..13): {[round(x, 4) for x in e[:12]]}")

print("\n# Part 2: exact peaks of the trivial-character law at the convergents, against the peak law and Q_(n+1)^2")
mp.dps = 80
TH = log(3) / log(2); L = TH - 1; delta = 2 - TH
def peak_at(Q, P, theta, Lm, dl):
    tf = Q * theta - P
    if abs(tf) < mpf('0.05'):
        w0 = 2j * pi * tf / dl - 4 * pi**2 * tf**2 / dl**3
        w = findroot(lambda w: exp(-w) + exp(2j * pi * tf) * exp(w * Lm) - 2, mpc(w0))
        s = 2j * pi * Q + 1j * w.imag                       # on the line, at the drifted frequency
        R = 1 / (1 - exp(-s) / 2 - exp(s * Lm) / 2)
        return tf, w, abs(R)
    # large defect: the root expansion is not valid; maximise |R| on the line near f = Q by a two-stage grid
    best = (mpf(0), None)
    for width, pts in ((mpf('0.5'), 2001), (mpf('0.001'), 2001)):
        centre = Q if best[1] is None else best[1]
        for k in range(pts):
            f = centre - width + 2 * width * k / (pts - 1)
            s = 2j * pi * f
            val = abs(1 / (1 - exp(-s) / 2 - exp(s * Lm) / 2))
            if val > best[0]: best = (val, f)
    return tf, mpc(0, 2 * pi * (best[1] - Q)), best[0]
print("    n     Q_n          theta_f          |R| at drifted f     (2-theta)^2/(2 pi^2 theta_f^2)   ratio     log|R|/log Q_n   2 e_n")
for n in range(1, 31):
    Q, P = conv[n]; Qn1 = conv[n + 1][0]
    tf, w, Rabs = peak_at(mpf(Q), mpf(P), TH, L, delta)
    pred = delta**2 / (2 * pi**2 * tf**2)
    en = math.log(Qn1) / math.log(Q) if Q > 1 else float('nan')
    print(f"   {n:2d}  {Q:12d}  {nstr(tf, 6):>14s}  {nstr(Rabs, 8):>18s}  {nstr(pred, 8):>18s}  {nstr(Rabs/pred, 6):>8s}   {float(log(Rabs)/log(Q)) if Q > 1 else float('nan'):8.4f}   {2*en:8.4f}")

print("\n# Part 3: Jarnik laws -- formal transfer law with parameter x of prescribed irrationality exponent mu")
def jarnik(mu, levels, dps):
    mp.dps = dps
    a = [1, 2]                     # x = 1 + 1/(2 + 1/(a_2 + ...)) in (1, 2)
    q = [1, 2]; p = [1, 3]
    for _ in range(levels):
        an = int(math.ceil(q[-1] ** (mu - 2))) if mu > 2 else 1
        a.append(an); q.append(an * q[-1] + q[-2]); p.append(an * p[-1] + p[-2])
    # the number: evaluate the continued fraction
    x = mpf(0)
    for an in reversed(a[1:]): x = 1 / (an + x)
    x = a[0] + x
    return a, [(q[i], p[i]) for i in range(len(a))], x
for mu, levels, dps in ((2, 40, 60), (3, 7, 400), (4, 6, 1000), (5, 5, 1700)):
    a_j, conv_j, x = jarnik(mu, levels, dps)
    Lj = x - 1; dj = 2 - x
    print(f"   mu = {mu}: x = {nstr(x, 12)}, partial quotients a_1.. = {[str(t)[:24] for t in a_j[1:6]]}")
    rows = []
    rng = range(2, min(len(conv_j) - 1, 8)) if mu > 2 else range(5, 40, 6)
    for n in rng:
        Q, P = conv_j[n]; Qn1 = conv_j[n + 1][0]
        tf, w, Rabs = peak_at(mpf(Q), mpf(P), x, Lj, dj)
        en = math.log(Qn1) / math.log(Q)
        rows.append((n, len(str(Q)), float(log(Rabs) / log(Q)), 2 * en)); assert mp.dps > 2.2 * len(str(Qn1)) + 40, (mu, n, mp.dps, len(str(Qn1)))
    print("      n, digits(Q_n), log|R_n|/log Q_n, 2 e_n (target 2(mu-1) = %d): " % (2 * (mu - 1)) + "; ".join(f"({n}, {d}, {g:.4f}, {t:.4f})" for n, d, g, t in rows))

print("\n# Part 4: characters modulo m up to sign against the cusps of X(m); the convergents of log2 3 over the cusps of X(5)")
def J2(m):
    out = m * m
    for p in set(f for f in range(2, m + 1) if m % f == 0 and all(f % r for r in range(2, int(f ** 0.5) + 1))):
        out = out * (p * p - 1) // (p * p)
    return out
def cusps(m): return 3 if m == 2 else J2(m) // 2
print("   m, J_2(m), J_2(m)/2 (signed pairs up to sign), cusps of X(m):", [(m, J2(m), J2(m) // 2 if m > 2 else 3, cusps(m)) for m in range(2, 13)])
print("   Platonic counts: X(2) 3 (dihedron), X(3) 4 (tetrahedron), X(4) 6 (octahedron), X(5) 12 (icosahedron); X(7) 24 (Klein quartic)")
classes = {}
for (Q, P) in conv[1:]:
    key = min((Q % 5, P % 5), ((-Q) % 5, (-P) % 5)); classes[key] = classes.get(key, 0) + 1
print(f"   X(5): {len(classes)} cusp classes hit by {sum(classes.values())} convergents; counts {sorted(classes.values())} (uniform would be {sum(classes.values())/12:.1f})")
pts = {}
for (Q, P) in conv[1:]:
    key = "inf" if Q % 5 == 0 else (P * pow(Q, -1, 5)) % 5; pts[key] = pts.get(key, 0) + 1
print(f"   P^1(F_5) (6 points, the exotic S_5 = PGL(2,5) of the outer automorphism of S_6): counts {dict(sorted(pts.items(), key=lambda kv: str(kv[0])))}")

print("\n# Part 5: Artin primes for the base 2 below 10^5 (THM-4631 multiplicity one at every level) and Wieferich primes")
primes = list(primerange(3, 100000))
artin = [p for p in primes if n_order(2, p) == p - 1]
wief = [p for p in primes if pow(2, p - 1, p * p) == 1]
print(f"   primes 3 <= p < 10^5: {len(primes)}; 2 a primitive root: {len(artin)} = {len(artin)/len(primes):.5f} of them (Artin's constant 0.3739558); Wieferich below 10^5: {wief}")
print(f"   simple fast spectrum at every level (primitive root mod p and not Wieferich): {len([p for p in artin if p not in wief])}")

print("\n# Part 6: Fisher information of a line (trivial character): Lorentzian half-width gamma = |Re w|/(2 pi), I = 1/(2 gamma^2)")
mp.dps = 60
for n in range(4, 14):
    Q, P = conv[n]
    tf, w, Rabs = peak_at(mpf(Q), mpf(P), TH, L, delta)
    gamma = abs(w.real) / (2 * pi); I = 1 / (2 * gamma**2)
    print(f"   n={n:2d} (Q,P)=({Q},{P}): theta_f = {nstr(tf, 5)}, gamma = {nstr(gamma, 5)}, Fisher information 1/(2 gamma^2) = {nstr(I, 6)}, ratio I / theta_f^(-4) = {nstr(I * tf**4, 6)} (law: (2-theta)^6/(8 pi^2) = {nstr(delta**6/(8*pi**2), 6)})")
print("done")
