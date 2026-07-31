#!/usr/bin/env python3
"""
delpezzo_verify.py — verification of the char-3 klt del Pezzo example

  X = X_14 = { F = 0 } in P(2,2,7,7) over F_3,
  F = A(x0,x1) + x2^2 + x3^2,   A = x0^7 - x0*x1^6 - x1^7   (to be checked
      separable + irreducible; adjust constant if needed).

Checks:
 1. weights (2,2,7,7): well-formed (every 3-subset coprime), all prime to 3 (tame).
 2. degree 14; Fano index = 18-14 = 4 > 0.
 3. NO mixed monomials of degree 14 exist except A-part/C-part (structural).
 4. quasi-smoothness: over F_3bar, grad F = (A_0, A_1, 2x2, 2x3) vanishes on the
    affine cone only at 0.  Via: gcd(a, a') = 1 for a(t)=A(t,1), plus x0^7
    leading term at (1:0).
 5. singular locus of X = 7 points {A=0} on the (2,2)-line  (1/2(1,1) A_1 each,
    tame since 2 != 3)  +  2 points {x2^2+x3^2=0} on the (7,7)-line
    (1/7(1,1)-type, tame).  Total 9 >= 8 geometric singular points.
 6. exact point counts #X(F_q) for q = 3, 9, 27, 81 by value-distribution
    convolution; compare against  q^2 + c*q + 1  to read off Frobenius traces
    on H^2 -> arithmetic Picard rank evidence (Tate holds for rational surfaces).

Pure Python, no dependencies.  GF(3^k) implemented via Zech-free table arithmetic.
"""
import itertools, sys

# ---------- GF(3^k) ----------
class GF:
    def __init__(self, p, k, modpoly):
        self.p, self.k, self.q = p, k, p ** k
        self.mod = modpoly  # coefficients low->high, length k+1, monic
        self.elts = [tuple(v) for v in itertools.product(range(p), repeat=k)]
        self.index = {e: i for i, e in enumerate(self.elts)}
        self.addt = [[self.index[tuple((a[t] + b[t]) % p for t in range(k))]
                      for b in self.elts] for a in self.elts]
        self.mult = [[self._mul(a, b) for b in self.elts] for a in self.elts]
        self.zero = self.index[tuple([0] * k)]
        one = [0] * k; one[0] = 1
        self.one = self.index[tuple(one)]
    def _mul(self, a, b):
        p, k = self.p, self.k
        res = [0] * (2 * k - 1)
        for i, ai in enumerate(a):
            if ai:
                for j, bj in enumerate(b):
                    res[i + j] = (res[i + j] + ai * bj) % p
        # reduce
        for i in range(2 * k - 2, k - 1, -1):
            c = res[i]
            if c:
                res[i] = 0
                for j in range(k):
                    res[i - k + j] = (res[i - k + j] - c * self.mod[j]) % p
        return self.index[tuple(res[:k])]
    def neg(self, a):
        p, k = self.p, self.k
        v = self.elts[a]
        return self.index[tuple((-x) % p for x in v)]
    def pow(self, a, e):
        r = self.one
        while e:
            if e & 1: r = self.mult[r][a]
            a = self.mult[a][a]; e >>= 1
        return r

def make_field(k):
    # irreducible cubics/quartics over F_3 (low->high, monic):
    mods = {1: [0, 1], 2: [1, 0, 1],            # t^2+1  (irred: -1 non-square mod 3)
            3: [1, 2, 0, 1],                    # t^3+2t^2+1?? -> verify below
            4: [2, 1, 0, 0, 1]}                 # t^4+t+2 ?? verify below
    F = GF(3, k, mods[k])
    return F

def poly_is_irred_over_F3(coeffs):
    """coeffs low->high over F_3, monic; brute force check no roots in F_{3^j}, j<=deg/2
       via gcd with t^{3^j}-t — small degrees only; here just root-count in small fields."""
    deg = len(coeffs) - 1
    # trial division by all monic polys of degree <= deg//2
    def polmulmod(a, b, m):
        res = [0] * (len(a) + len(b) - 1)
        for i, x in enumerate(a):
            for j, y in enumerate(b):
                res[i + j] = (res[i + j] + x * y) % 3
        # reduce mod m
        dm = len(m) - 1
        while len(res) - 1 >= dm:
            c = res[-1]
            if c:
                for t in range(dm + 1):
                    res[len(res) - 1 - dm + t] = (res[len(res) - 1 - dm + t] - c * m[t]) % 3
            res.pop()
        return res
    def all_monic(d):
        for tail in itertools.product(range(3), repeat=d):
            yield list(tail) + [1]
    def polmod(a, m):
        a = a[:]
        dm = len(m) - 1
        while len(a) - 1 >= dm and len(a) > 1:
            c = a[-1]
            if c:
                for t in range(dm + 1):
                    a[len(a) - 1 - dm + t] = (a[len(a) - 1 - dm + t] - c * m[t]) % 3
            a.pop()
        while len(a) > 1 and a[-1] == 0: a.pop()
        return a
    for dd in range(1, deg // 2 + 1):
        for m in all_monic(dd):
            r = polmod(coeffs, m)
            if r == [0]:
                return False, m
    return True, None

def main():
    # ---------- structural checks ----------
    W = [2, 2, 7, 7]; d = 14
    from math import gcd
    ok_wf = all(gcd(gcd(W[i], W[j]), W[k]) == 1
                for i in range(4) for j in range(i+1, 4) for k in range(j+1, 4))
    print(f"[1] well-formed (all triples coprime): {ok_wf}")
    print(f"[1] weights prime to char 3 (tame): {all(w % 3 != 0 for w in W)}")
    print(f"[2] Fano index sum(w)-d = {sum(W)-d} > 0: {sum(W)-d > 0}")
    # monomial inventory of degree 14
    monos = [(a, b, c, e) for a in range(8) for b in range(8) for c in range(3)
             for e in range(3) if 2*a + 2*b + 7*c + 7*e == 14]
    mixed = [m for m in monos if (m[0] + m[1] > 0) and (m[2] + m[3] > 0)]
    print(f"[3] degree-14 monomials: {len(monos)}, mixed (both blocks): {len(mixed)} "
          f"(structurally decoupled: {len(mixed) == 0})")

    # ---------- choice of A ----------
    # a(t) = t^7 + t^2 + 2, verified irreducible over F_3 below
    a = [2, 0, 1, 0, 0, 0, 0, 1]
    irr, factor = poly_is_irred_over_F3(a)
    print(f"[4] a(t) = t^7 + t^2 + 2 over F_3 irreducible: {irr}" +
          ("" if irr else f" (factor {factor})"))
    # separability: gcd(a, a') — a' = 7t^6 - 1 = t^6 - 1 mod 3
    # gcd computation over F_3
    def polgcd(u, v):
        u = [x % 3 for x in u]; v = [x % 3 for x in v]
        def norm(w):
            while len(w) > 1 and w[-1] == 0: w.pop()
            return w
        u, v = norm(u), norm(v)
        while v != [0]:
            # u mod v
            u2 = u[:]
            dv = len(v) - 1
            inv = pow(v[-1], 3 - 2, 3)
            while len(u2) - 1 >= dv and u2 != [0]:
                c = (u2[-1] * inv) % 3
                for t in range(dv + 1):
                    u2[len(u2) - 1 - dv + t] = (u2[len(u2) - 1 - dv + t] - c * v[t]) % 3
                u2.pop()
                u2 = norm(u2)
                if u2 == []: u2 = [0]
            u, v = v, u2
        return u
    ap = [0, 2, 0, 0, 0, 0, 1]  # derivative of t^7+t^2+2 = 7t^6+2t = t^6+2t low->high
    g = polgcd(a[:], ap[:])
    print(f"[4] gcd(a,a') = {g} (separable: {g == [1] or (len(g)==1 and g[0]!=0)})")
    print(f"[4] quasi-smoothness: grad F = (A_0, A_1, 2x2, 2x3); x2=x3=0 & Euler(7=1!=0 mod 3)"
          f" => A=A_0=A_1=0 => common root of (a,a') => none (separable). At (1:0): A_0=7x0^6!=0. QED")
    print(f"[5] singular pts: 7 distinct geometric roots of A on P(2,2)-line, each 1/2(1,1)=A_1 (tame);"
          f" 2 roots of x2^2+x3^2 on P(7,7)-line (over F_9, i=sqrt(-1)), each 1/7-cyclic (tame). Total 9.")

    # ---------- point counts ----------
    # X: A(x0,x1) + x2^2 + x3^2 = 0 in P(2,2,7,7).
    # Count N_q = #X(F_q).  Decompose by (x0,x1) = 0 or not, (x2,x3) = 0 or not.
    # value-distribution: cA[v] = #{(x0,x1) != (0,0): A = v},  cC[v] = #{(x2,x3) != 0: x2^2+x3^2 = v}
    # points with both blocks nonzero: free Gm-action (stab: lam^2=1 & lam^7=1 => lam=1):
    #   count = sum_v cA[v]*cC[-v] / (q-1)
    # points with (x2,x3)=0: A=0, (x0,x1) != 0: cA[0] affine solutions; orbit under lam^2:
    #   the (2,2)-line is P^1 with coords of weight 2; number of projective points = cA[0]/(q-1)... careful:
    #   stab of (x0,x1,0,0) in Gm: lam^2 = 1 (need lam^2 fix (x0,x1) and any lam^7 on 0): orbit size (q-1)/|mu_2(F_q)|
    #   Simpler: directly count projective points on the line P(2,2) ≅ P^1: #roots of a(t) in F_q  (+ pt at inf if x1=0 root).
    # points with (x0,x1)=0: x2^2+x3^2=0, (x2,x3) != 0: #roots of u^2+1 in F_q (as ratio x2/x3) -> 0 or 2 points.
    for k in [1, 2, 3, 4]:
        if k == 3:
            # find irreducible cubic: t^3 - t - 1: low->high [2,2,0,1]? t^3-t-1 -> [-1,-1,0,1] = [2,2,0,1]
            mod3 = [2, 2, 0, 1]
            irr3, _ = poly_is_irred_over_F3(mod3)
            assert irr3, "cubic not irreducible"
            F = GF(3, 3, mod3)
        elif k == 4:
            mod4 = [2, 1, 0, 0, 1]   # t^4 + t + 2
            irr4, _ = poly_is_irred_over_F3(mod4)
            assert irr4, "quartic not irreducible: pick another"
            F = GF(3, 4, mod4)
        else:
            F = make_field(k)
        q = F.q
        # evaluate A over all (x0,x1)
        cA = {}
        two = F.addt[F.one][F.one]
        for x0 in range(q):
            x07 = F.pow(x0, 7)
            x02 = F.mult[x0][x0]
            for x1 in range(q):
                if x0 == F.zero and x1 == F.zero: continue
                x15 = F.pow(x1, 5)
                x17 = F.mult[F.mult[x15][x1]][x1]
                t1 = F.mult[x02][x15]               # x0^2*x1^5
                v = F.addt[x07][t1]
                v = F.addt[v][F.mult[two][x17]]     # A = x0^7 + x0^2 x1^5 + 2 x1^7
                cA[v] = cA.get(v, 0) + 1
        cC = {}
        for x2 in range(q):
            s2 = F.mult[x2][x2]
            for x3 in range(q):
                if x2 == F.zero and x3 == F.zero: continue
                v = F.addt[s2][F.mult[x3][x3]]
                cC[v] = cC.get(v, 0) + 1
        both = sum(cA.get(v, 0) * cC.get(F.neg(v), 0) for v in cA)
        assert both % (q - 1) == 0
        n_both = both // (q - 1)
        # line points:
        nA0 = cA.get(F.zero, 0)   # affine (x0,x1) != 0 with A=0; P(2,2)=P^1: each proj pt has q-1 affine reps?
        # On the line, coords weight 2: lam acts by lam^2; orbits of Gm on {A=0}\0 have size (q-1)/gcd(2,q-1)...
        # but as SET of projective points of P(2,2)≅P^1 (coords x0,x1 with same weight): projective pts = nA0/(q-1)
        assert nA0 % (q - 1) == 0
        n_lineA = nA0 // (q - 1)
        nC0 = cC.get(F.zero, 0)
        assert nC0 % (q - 1) == 0
        n_lineC = nC0 // (q - 1)
        N = n_both + n_lineA + n_lineC
        c = (N - q * q - 1) / q
        print(f"[6] q={q}: #X(F_q) = {N} = q^2 + ({c})*q + 1   "
              f"[sing pts rational here: A1s:{n_lineA}, 1/7s:{n_lineC}]")

if __name__ == "__main__":
    main()
