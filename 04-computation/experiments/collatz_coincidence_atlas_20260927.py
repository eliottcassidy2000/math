#!/usr/bin/env python3
"""collatz_coincidence_atlas_20260927.py -- exact checks behind the coincidence atlas (session collatz-posets-zeta5-20260927,
opus, 2026-09-27, third note).

 (1) The owner's triangle {1},{2,1},{3,3,1},{4,6,4,1},{5,10,9,5,1},... is T(n,j) = P_(j+2)(n-j), the (n-j)-th (j+2)-gonal number;
     edges: P_2(n) = n, P_k(1) = 1, P_k(2) = k; hidden region = negative index, P_k(-m) = P_k(m) + (k-4) m; apex P_5(-1) = 2;
     pentagonal column: P_5(n) = sum_(j<n) (3j+1), P_5(-n) = sum_(j<=n) (3j-1) (the two Collatz sheets' affine forms),
     24 P_5(n) + 1 = (6n-1)^2, 24 P_5(-n) + 1 = (6n+1)^2; Syracuse predecessor valuations: v odd for m = 6k-1, v even for m = 6k+1.
 (2) Euler's pentagonal theorem prod (1 - q^n) = sum_k (-1)^k q^(k(3k-1)/2) checked to q^400; the eta indexing by 6k -+ 1 with the
     character chi_12.
 (3) Golden facts: phi^n = (L_n + F_n sqrt 5)/2, round(phi^n) = L_n for n >= 2; zeta_5 + conj = 1/phi; cube-in-dodecahedron
     volume ratio 2/(2 + phi) (exact in Q(sqrt 5)); spectra of the dodecahedron and icosahedron graphs; their antipodal quotients
     are the Petersen graph and K_6 (spectral identification; Petersen is determined by its spectrum).
 (4) The number 139: 3^7 - 2^11; the -17 cycle (7 odd steps, 11 halvings, S = 17 * 139); 139 = 3 mod 4 (a Paley tournament
     P_139 exists, |Aut| = 139 * 69); genus of X_0(p) for p = 7, 11, 139; Ihara discriminants 9 - 4n of K_n (Heegner at n = 3,4,5,7,13).
Usage: python3 collatz_coincidence_atlas_20260927.py
"""
import math, cmath, itertools
from fractions import Fraction


def P(k, n):
    return ((k - 2) * n * n - (k - 4) * n) // 2


def part1():
    print("== (1) the polygonal triangle and its hidden region ==")
    rows = [[1], [2, 1], [3, 3, 1], [4, 6, 4, 1], [5, 10, 9, 5, 1], [6, 15, 16, 12, 6, 1], [7, 21, 25, 22, 15, 7, 1]]
    print(" T(n,j) = P_(j+2)(n-j):", all(rows[n - 1][j] == P(j + 2, n - j) for n in range(1, 8) for j in range(n)))
    print(" edges: P_2(n) = n, P_k(1) = 1, P_k(2) = k:", all(P(2, n) == n and P(k, 1) == 1 and P(k, 2) == k for n in range(1, 40) for k in range(2, 40)))
    print(" hidden region P_k(-m) = P_k(m) + (k-4) m:", all(P(k, -m) == P(k, m) + (k - 4) * m for k in range(2, 30) for m in range(1, 30)))
    print(" apex P_5(-1) = %d; hidden pentagonal 2, 7, 15, 26 = visible 1, 5, 12, 22 plus m" % P(5, -1))
    print(" P_5(n) = sum_(j<n)(3j+1):", all(P(5, n) == sum(3 * j + 1 for j in range(n)) for n in range(1, 60)),
          "; P_5(-n) = sum_(1<=j<=n)(3j-1):", all(P(5, -n) == sum(3 * j - 1 for j in range(1, n + 1)) for n in range(1, 60)))
    print(" 24 P_5(n) + 1 = (6n-1)^2, 24 P_5(-n) + 1 = (6n+1)^2:", all(24 * P(5, n) + 1 == (6 * n - 1) ** 2 and 24 * P(5, -n) + 1 == (6 * n + 1) ** 2 for n in range(1, 60)))
    def pred_vals(m, vmax=13):
        return [v for v in range(1, vmax) if (2 ** v * m - 1) % 3 == 0]
    ok = all(all(v % 2 == 1 for v in pred_vals(m)) for m in range(5, 400, 6)) and all(all(v % 2 == 0 for v in pred_vals(m)) for m in range(7, 400, 6)) and all(pred_vals(m) == [] for m in range(3, 400, 6))
    print(" Syracuse predecessors (2^v m - 1)/3: v odd iff m = 6k-1, v even iff m = 6k+1, none for 3 | m:", ok)


def part2():
    print("== (2) Euler's pentagonal theorem and the eta indexing ==")
    N = 400
    c = [0] * (N + 1); c[0] = 1
    for n in range(1, N + 1):
        for e in range(N, n - 1, -1):
            c[e] -= c[e - n]
    expo = {}
    for k in range(-30, 31):
        e = k * (3 * k - 1) // 2
        if 0 <= e <= N:
            expo[e] = (-1) ** k
    ok = all(c[e] == expo.get(e, 0) for e in range(N + 1))
    print(" prod_(n<=%d) (1 - q^n) = sum_k (-1)^k q^(k(3k-1)/2) to order %d: %s" % (N, N, ok))
    # eta indexing: 24 e + 1 = m^2 with m = 6k - 1 (k > 0) or 6|k| + 1 (k < 0); sign = chi_12(m) = +1 for m = +-1 mod 12, -1 for m = +-5 mod 12
    def chi12(m):
        r = m % 12
        return 1 if r in (1, 11) else (-1 if r in (5, 7) else 0)
    ok2 = True
    for k in range(1, 30):
        e1, e2 = k * (3 * k - 1) // 2, k * (3 * k + 1) // 2
        m1, m2 = 6 * k - 1, 6 * k + 1
        ok2 &= (24 * e1 + 1 == m1 * m1 and 24 * e2 + 1 == m2 * m2 and (-1) ** k == chi12(m1) == chi12(m2))
    print(" exponents e with 24e + 1 = m^2, m = 6k -+ 1, and (-1)^k = chi_12(m) (the classical eta(tau) = sum chi_12(m) q^(m^2/24)): %s" % ok2)
    print(" so eta's exponents are indexed by the odd numbers coprime to 3 (the Syracuse core, the nodes with odd predecessors),")
    print(" the two signs of k being the two predecessor classes 6k-1 (v odd) and 6k+1 (v even)")


def part3():
    print("== (3) golden facts, exact ==")
    F = [0, 1]; L = [2, 1]
    for _ in range(40):
        F.append(F[-1] + F[-2]); L.append(L[-1] + L[-2])
    phi = (1 + 5 ** 0.5) / 2
    ok = all(abs(phi ** n - (L[n] + F[n] * 5 ** 0.5) / 2) < 1e-6 for n in range(0, 30))
    okr = all(round(phi ** n) == L[n] for n in range(2, 30))
    print(" phi^n = (L_n + F_n sqrt5)/2 (n < 30): %s; round(phi^n) = L_n for 2 <= n < 30: %s (phi^0 = 1 vs L_0 = 2, phi vs L_1 = 1: the two reversed)" % (ok, okr))
    z = cmath.exp(2j * math.pi / 5)
    print(" zeta_5 + conj(zeta_5) = %.12f = 1/phi = phi - 1 = %.12f (a golden integer; [Q(zeta_5):Q(sqrt5)] = 2)" % ((z + z.conjugate()).real, phi - 1))
    # cube in dodecahedron: edge ratio phi, volume ratio phi^3 / ((15 + 7 sqrt5)/4) = 2/(2 + phi)  -- exact in Q(sqrt5): represent a + b sqrt5
    class Q5:
        def __init__(s, a, b): s.a = Fraction(a); s.b = Fraction(b)
        def __mul__(s, o): return Q5(s.a * o.a + 5 * s.b * o.b, s.a * o.b + s.b * o.a)
        def __add__(s, o): return Q5(s.a + o.a, s.b + o.b)
        def inv(s):
            n = s.a * s.a - 5 * s.b * s.b
            return Q5(s.a / n, -s.b / n)
        def __eq__(s, o): return s.a == o.a and s.b == o.b
        def __repr__(s): return "%s + %s sqrt5" % (s.a, s.b)
    ph = Q5(Fraction(1, 2), Fraction(1, 2))
    vol_dodec = Q5(Fraction(15, 4), Fraction(7, 4))
    ratio = (ph * ph * ph) * vol_dodec.inv()
    target = Q5(2, 0) * (Q5(2, 0) + ph).inv()
    print(" cube edge = phi * dodecahedron edge (pentagon diagonal); volume ratio phi^3 / ((15 + 7 sqrt5)/4) = %s; 2/(2 + phi) = %s; equal: %s" % (ratio, target, ratio == target))
    # dodecahedron and icosahedron graphs from coordinates; spectra; antipodal quotients
    import numpy as np
    def spectrum(A):
        w = np.linalg.eigvalsh(A)
        return sorted([round(float(x), 6) for x in w], reverse=True)
    def graph_from_points(pts, edge_len):
        n = len(pts); A = np.zeros((n, n))
        for i in range(n):
            for j in range(i + 1, n):
                if abs(np.linalg.norm(np.array(pts[i]) - np.array(pts[j])) - edge_len) < 1e-6:
                    A[i, j] = A[j, i] = 1
        return A
    def antipodal_quotient(pts, A):
        n = len(pts); rep = {}
        cls = []
        for i in range(n):
            j = next(k for k in range(n) if np.allclose(np.array(pts[k]), -np.array(pts[i])))
            key = min(i, j)
            if key not in rep:
                rep[key] = len(cls); cls.append(key)
        m = len(cls); Q = np.zeros((m, m))
        for i in range(n):
            for j in range(n):
                if A[i, j]:
                    a = rep[min(i, next(k for k in range(n) if np.allclose(np.array(pts[k]), -np.array(pts[i]))))]
                    b = rep[min(j, next(k for k in range(n) if np.allclose(np.array(pts[k]), -np.array(pts[j]))))]
                    if a != b:
                        Q[a, b] = 1
        return Q
    ph_ = phi
    ico = [p for s1 in (1, -1) for s2 in (1, -1) for p in ((0, s1, s2 * ph_), (s1, s2 * ph_, 0), (s2 * ph_, 0, s1))]
    A_ico = graph_from_points(ico, 2.0)
    dod = [(s1, s2, s3) for s1 in (1, -1) for s2 in (1, -1) for s3 in (1, -1)]
    dod += [p for s1 in (1, -1) for s2 in (1, -1) for p in ((0, s1 / ph_, s2 * ph_), (s1 / ph_, s2 * ph_, 0), (s2 * ph_, 0, s1 / ph_))]
    A_dod = graph_from_points(dod, 2.0 / ph_)
    print(" icosahedron graph: %d vertices, degree %d, spectrum %s" % (len(ico), int(A_ico.sum(1)[0]), spectrum(A_ico)))
    print(" dodecahedron graph: %d vertices, degree %d, spectrum %s" % (len(dod), int(A_dod.sum(1)[0]), spectrum(A_dod)))
    Qd = antipodal_quotient(dod, A_dod); Qi = antipodal_quotient(ico, A_ico)
    # Petersen = Kneser K(5,2)
    pairs = list(itertools.combinations(range(5), 2)); A_pet = np.zeros((10, 10))
    for i, a in enumerate(pairs):
        for j, b in enumerate(pairs):
            if not set(a) & set(b):
                A_pet[i, j] = 1
    print(" dodecahedron / antipodal: %d vertices, spectrum %s; Petersen spectrum %s; equal (Petersen is determined by its spectrum): %s" % (len(Qd), spectrum(Qd), spectrum(A_pet), spectrum(Qd) == spectrum(A_pet)))
    print(" icosahedron / antipodal: %d vertices, spectrum %s = K_6's [5, -1 x5]: %s" % (len(Qi), spectrum(Qi), spectrum(Qi) == [5.0] + [-1.0] * 5))
    print(" reading: the golden eigenvalues +-sqrt5 (multiplicity 3, the icosahedral 3-dimensional representations) are exactly the part killed by the antipodal quotient")


def genus_X0(p):
    # prime level: g = (p+1)/12 - 1/4 (1 + (-1/p)) - 1/3 (1 + (-3/p)) + 1 ... use the standard formula
    def leg(a, p):
        return 0 if a % p == 0 else (1 if pow(a % p, (p - 1) // 2, p) == 1 else -1)
    mu = p + 1
    nu2 = 1 + leg(-1, p) if p != 2 else 1
    nu3 = 1 + leg(-3, p) if p != 3 else 1
    cusps = 2
    return 1 + Fraction(mu, 12) - Fraction(nu2, 4) - Fraction(nu3, 3) - Fraction(cusps, 2)


def part4():
    print("== (4) the number 139 across the threads ==")
    print(" 3^7 - 2^11 = %d; the -17 cycle: 7 odd steps, 11 halvings, -17 = S/(2^11 - 3^7), S = %d = 17 * 139" % (3 ** 7 - 2 ** 11, 17 * 139))
    print(" 139 mod 4 = %d: a Paley tournament on 139 vertices exists, |Aut| = 139 * 69 = %d; 7 mod 4 = 3: the Paley heptagon, |Aut| = 21" % (139 % 4, 139 * 69))
    for p in (7, 11, 139):
        print(" genus of X_0(%d) = %s" % (p, genus_X0(p)))
    print(" Ihara quadratic factor of K_n has discriminant 9 - 4n:", {n: 9 - 4 * n for n in (3, 4, 5, 7, 13)}, "(Heegner discriminants -3, -7, -11, -19, -43; K_5 gives -11)")
    print(" tournament thread (cited, not recomputed): flip-rank k(7) = 12; the best 11-arc configuration's 2^11 completions reach 454 of the 456 classes and miss the Paley heptagon")
    print(" verdict: 7 and 11 are the cycle's (p, A); 2^11 < 3^7 is a near-miss of log_2 3 (11/7 a semiconvergent); the modular-curve and Paley facts are properties of 7, 11, 139 as integers, not of the Collatz map")


def main():
    part1(); part2(); part3(); part4()


if __name__ == '__main__':
    main()
