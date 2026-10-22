#!/usr/bin/env python3
"""(b) Paper 2: printed cumulant formulas, the rank-43 witness pair, Delta c5,
identification of c_m as log-coefficients of the orthogonal-root-frame generating function,
Euler-Boole kernel, Coxeter-defect identities."""
import sys; sys.path.insert(0, '.')
from fractions import Fraction as Fr
import numpy as np, sympy as sp
from math import comb, factorial
from rootsys import roots, gram, count_tuples

def cA(k):
    return [Fr(k*(k+1)), Fr(-k*(k+1)*(2*k-1)), Fr(4, 3)*k*(k+1)*(5*k**2-7*k+3),
            Fr(-2)*k*(k+1)*(14*k**3-37*k**2+39*k-15), Fr(4, 5)*k*(k+1)*(42*k**4-175*k**3+301*k**2-245*k+85)]
def cD(k):
    return [Fr(2*k*(k-1)), Fr(-2*k*(k-1)*(4*k-7)), Fr(8, 3)*k*(k-1)*(20*k**2-74*k+69),
            Fr(-4)*k*(k-1)*(112*k**3-672*k**2+1378*k-963), Fr(16, 5)*k*(k-1)*(336*k**4-2800*k**3+8876*k**2-12705*k+6920)]
cE = {6: [Fr(x) for x in (72, -1512, 50976, -2011824, 87400704)],
      7: [Fr(x) for x in (126, -4158, 223272, -14196924, 990424560)],
      8: [Fr(x) for x in (240, -13680, 1281600, -143445600, 17540208000)]}
def cprint(t, k): return cA(k) if t == 'A' else cD(k) if t == 'D' else cE[k]

def sig(ms, f=cprint):
    s = [Fr(0)]*5
    for t, k, m in ms:
        for i, x in enumerate(f(t, k)): s[i] += m*x
    return s

R1 = [('A',1,2),('A',4,2),('A',5,1),('A',7,2),('A',8,1),('D',6,1)]
R2 = [('A',2,2),('A',3,1),('A',6,2),('A',9,1),('D',5,3)]
print("rank R1 =", sum(k*m for _,k,m in R1), " rank R2 =", sum(k*m for _,k,m in R2))
s1, s2 = sig(R1), sig(R2)
print("printed-formula signature R1:", [str(x) for x in s1])
print("printed-formula signature R2:", [str(x) for x in s2])
print("c1..c4 equal:", s1[:4] == s2[:4], " (c1,c2,3c3,c4) =", (s1[0], s1[1], 3*s1[2], s1[3]))
print("c5(R1) =", s1[4], " c5(R2) =", s2[4], " Delta c5 =", s2[4]-s1[4], " 480*Delta c5 =", 480*(s2[4]-s1[4]))

# ---- orthogonal-frame EGF: Z_X(z) = sum_m O_m(X) z^m/m!, O_m = # ordered m-tuples of mutually orthogonal roots
M = 9
z = sp.symbols('z')
def series_mul(a, b): return [sum(a[i]*b[n-i] for i in range(n+1)) for n in range(M)]
from functools import lru_cache
@lru_cache(None)
def Z(t, k):
    """EGF coefficients O_m/m! as Fractions, via Z_X' = N_X * Z_{X cap alpha^perp} (W transitive on roots)."""
    if t == 'A':
        if k <= 0: return tuple([Fr(1)] + [Fr(0)]*(M-1))
        N, sub = k*(k+1), [Z('A', k-2)]
    elif t == 'D':
        if k == 1: return tuple([Fr(1)] + [Fr(0)]*(M-1))
        if k == 2: return tuple(series_mul(list(Z('A', 1)), list(Z('A', 1))))
        if k == 3: return Z('A', 3)
        N, sub = 2*k*(k-1), [Z('D', k-2), Z('A', 1)]
    else:
        N = {6: 72, 7: 126, 8: 240}[k]
        sub = [Z('A', 5)] if k == 6 else [Z('D', 6)] if k == 7 else [Z('E', 7)]
    prod = [Fr(1)] + [Fr(0)]*(M-1)
    for s in sub: prod = series_mul(prod, list(s))
    # integrate: coefficient of z^m in Z is N * [z^{m-1}]prod / m
    out = [Fr(1)] + [N*prod[m-1]/m for m in range(1, M)]
    return tuple(out)

def logcoeffs(a):
    # log of power series with a[0]=1: l' = a'/a
    l = [Fr(0)]*M
    for n in range(1, M):
        # n a_n = sum_{j=1}^{n} j l_j a_{n-j}
        l[n] = (n*a[n] - sum(j*l[j]*a[n-j] for j in range(1, n)))/n
    return l

# brute-force check of O_m for a few systems
print("\nbrute-force check of orthogonal-frame counts O_m (Z_X recursion vs explicit enumeration):")
for t, k in [('A',5),('A',7),('D',4),('D',6),('E',6),('E',7),('E',8)]:
    R, sc = roots(t, k); IP = gram(R, sc)
    for m in (2, 3, 4):
        G = 2*np.eye(m, dtype=np.int64)
        bf = R.shape[0]*count_tuples(IP, G, first_candidates=[0])
        rec = Z(t, k)[m]*factorial(m)
        assert bf == rec, (t, k, m, bf, rec)
    print(f"  {t}{k}: O_1..O_4 = {[int(Z(t,k)[m]*factorial(m)) for m in range(1,5)]}  (brute force agrees for m=2,3,4)")

print("\nlog-EGF coefficients l_m vs printed c_m:")
bad5 = []
allmatch14 = True
for t, ks in [('A', range(1, 16)), ('D', range(4, 16)), ('E', (6, 7, 8))]:
    for k in ks:
        l = logcoeffs(list(Z(t, k)))[1:6]
        c = cprint(t, k)
        m14 = all(l[i] == c[i] for i in range(4))
        allmatch14 &= m14
        if l[4] != c[4]: bad5.append((f"{t}{k}", str(l[4]), str(c[4])))
print("  l_1..l_4 == printed c_1..c_4 for A1..A15, D4..D15, E6, E7, E8:", allmatch14)
print("  l_5 vs printed c_5 mismatches (type, true l5, printed c5):")
for b in bad5[:8]: print("    ", b)
print("   ... total mismatching types:", len(bad5), "of", 15+12+3)

# fit the true l5 polynomials
kk = sp.symbols('k')
ptsA = [(k, logcoeffs(list(Z('A', k)))[5]) for k in range(1, 12)]
ptsD = [(k, logcoeffs(list(Z('D', k)))[5]) for k in range(4, 16)]
polyA = sp.factor(sp.interpolate([(a, sp.Rational(b.numerator, b.denominator)) for a, b in ptsA], kk))
polyD = sp.factor(sp.interpolate([(a, sp.Rational(b.numerator, b.denominator)) for a, b in ptsD], kk))
print("  true l5(A_k) =", polyA)
print("  true l5(D_k) =", polyD)
print("  true l5(E6,E7,E8) =", [str(logcoeffs(list(Z('E', k)))[5]) for k in (6, 7, 8)])
# also check printed c1..c4 polynomials are the interpolants (degree m+1)
for m in range(1, 5):
    pA = sp.expand(sp.interpolate([(k, sp.Rational(*(lambda f: (f.numerator, f.denominator))(logcoeffs(list(Z('A', k)))[m]))) for k in range(1, 12)], kk))
    print(f"  l{m}(A_k) interpolant = {sp.factor(pA)}")

def sig_true(ms):
    s = [Fr(0)]*6
    for t, k, m in ms:
        l = logcoeffs(list(Z(t, k)))
        for i in range(1, 6): s[i] += m*l[i]
    return s[1:]
t1, t2 = sig_true(R1), sig_true(R2)
print("\nTRUE log-EGF signature R1:", [str(x) for x in t1])
print("TRUE log-EGF signature R2:", [str(x) for x in t2])
print("true Delta l5 =", t2[4]-t1[4])
# orthogonal frame counts O_m of R1, R2
def Zsum(ms):
    a = [Fr(1)] + [Fr(0)]*(M-1)
    for t, k, m in ms:
        for _ in range(m): a = series_mul(a, list(Z(t, k)))
    return [a[i]*factorial(i) for i in range(M)]
O1, O2 = Zsum(R1), Zsum(R2)
print("ordered orthogonal m-frame counts O_m, m=1..6:")
print("  R1:", [int(x) for x in O1[1:7]])
print("  R2:", [int(x) for x in O2[1:7]])
print("  Delta O_5 = O_5(R2)-O_5(R1) =", int(O2[5]-O1[5]), " = 120 * true Delta l5:", O2[5]-O1[5] == 120*(t2[4]-t1[4]))

# ---- Euler-Boole: kernel of the 5x6 matrix (rank, c1..c4) on A1..A6
V = sp.Matrix([[k] + [sp.Rational(x.numerator, x.denominator) for x in cA(k)[:4]] for k in range(1, 7)]).T
ker = V.nullspace()
print("\nEuler-Boole: kernel of (rank,c1..c4) on A1..A6 has dim", len(ker), "spanned by", list((ker[0]/ker[0][5]).T))
d = [(-1)**(6-j)*comb(6, j) for j in range(1, 7)]
print("  6th finite-difference vector:", d, " pos-part rank:", sum(j*x for j, x in zip(range(1,7), d) if x > 0),
      " neg-part rank:", sum(-j*x for j, x in zip(range(1,7), d) if x < 0))
s1 = sig([('A',1,6),('A',3,20),('A',5,6)]); s2 = sig([('A',2,15),('A',4,15),('A',6,1)])
print("  EB pair A1^6A3^20A5^6 vs A2^15A4^15A6: c1..c4 equal:", s1[:4] == s2[:4], " Delta c5 =", s2[4]-s1[4],
      " (6!*leading coeff of c5(A_k) = 720*(4/5)*42 =", Fr(720*4*42, 5), ")")

# ---- Coxeter-defect identities (3.3)-(3.6)
h = sp.symbols('h')
def cAs(k): return [sp.Rational(x.numerator, x.denominator) if isinstance(x, Fr) else x for x in
                    [k*(k+1), -k*(k+1)*(2*k-1), sp.Rational(4,3)*k*(k+1)*(5*k**2-7*k+3), -2*k*(k+1)*(14*k**3-37*k**2+39*k-15)]]
def cDs(k): return [2*k*(k-1), -2*k*(k-1)*(4*k-7), sp.Rational(8,3)*k*(k-1)*(20*k**2-74*k+69), -4*k*(k-1)*(112*k**3-672*k**2+1378*k-963)]
a = [sp.simplify(x/(h-1)) for x in cAs(h-1)]
dd = [sp.simplify(x/(h/2+1)) for x in cDs(h/2+1)]
print("\nper-rank A(h): ", [sp.factor(x) for x in a])
print("per-rank D(h): ", [sp.factor(x) for x in dd])
print("delta4(h) =", sp.factor(a[3]-dd[3]))

# ---- the small-M witness pairs of Table 3.1 (printed formulas)
def chk(Ra, Rb, Mm):
    sa, sb = sig(Ra), sig(Rb)
    ra, rb = sum(k*m for _,k,m in Ra), sum(k*m for _,k,m in Rb)
    return ra, rb, sa[:Mm] == sb[:Mm], sa[Mm] != sb[Mm]
T = {
 'EB M1': ([('A',1,3),('A',3,1)], [('A',2,3)], 1),
 'EB M2': ([('A',1,4),('A',3,4)], [('A',2,6),('A',4,1)], 2),
 'EB M3': ([('A',1,5),('A',3,10),('A',5,1)], [('A',2,10),('A',4,5)], 3),
 'EB M4': ([('A',1,6),('A',3,20),('A',5,6)], [('A',2,15),('A',4,15),('A',6,1)], 4),
 'A M2': ([('A',1,2),('A',6,2)], [('A',3,1),('A',4,1),('A',7,1)], 2),
 'A M3': ([('A',1,3),('A',5,1),('A',6,2),('A',8,1)], [('A',2,2),('A',3,1),('A',7,3)], 3),
 'A M4': ([('A',1,2),('A',8,2),('A',15,2)], [('A',3,1),('A',5,1),('A',11,1),('A',13,1),('A',16,1)], 4),
 'ADE M2': ([('A',2,4),('D',4,1)], [('A',1,4),('A',4,2)], 2),
 'ADE M3': ([('A',2,2),('A',4,1),('D',4,1)], [('A',1,1),('A',3,2),('A',5,1)], 3),
 'ADE M4': (R1, R2, 4)}
print("\nTable 3.1 witness pairs (rank1, rank2, first M moments equal, moment M+1 differs):")
for k, (Ra, Rb, Mm) in T.items(): print("  ", k, chk(Ra, Rb, Mm))
