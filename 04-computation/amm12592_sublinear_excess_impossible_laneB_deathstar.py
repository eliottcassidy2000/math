#!/usr/bin/env python3
"""
laneB_referee.py -- Lane B (C = 1 impossibility) referee for AMM 12592 frontier.

Verifies EVERY finite computation used in laneB_draft.md with exact arithmetic
(fractions.Fraction / sympy).  The classical theorems (Polya-Carlson, Szego,
Fatou's lemma, Kronecker) are cited, not machine-checked; everything else in
the draft is either a self-contained proof or is corroborated here.

Checks:
  C1  Decided-tree polynomial fact, BOTH directions, for depth d = 0,1,2,3:
      { G_rule : rule a stopping rule of depth <= d } equals
      { sum_k w_k p^(d-k) q^k : 0 <= w_k <= C(d,k) integer }.
  C2  Vandermonde refinement bound: sum_j C(d',j) C(d-d',k-j) = C(d,k), d<=8.
  C3  Coefficient bound: |coeff_j( (1-p) * sum_k w_k p^(d-k)(1-p)^k )| <= 2*3^d
      for all 0 <= w_k <= C(d,k), d <= 6 (via vertex enumeration; coeffs are
      linear in w so the box max is attained at a vertex).
  C4  {z in C : |z| = 1 and |1-z| = 1} = {(1 +- i sqrt3)/2} = primitive 6th
      roots of unity; 1 - zeta = conj(zeta) = zeta^5; both are roots of
      Phi_6(p) = p^2 - p + 1 = cyclotomic_poly(6).
  C5  Phi_6(1-p) = Phi_6(p); power series of 1/Phi_6 has period-6 coefficients
      (1,1,0,-1,-1,0); coefficients of (b1*p+b0)/Phi_6 are b0*a_n + b1*a_(n-1)
      (symbolic, to order 24); coefficients of 1/Phi_6^j (j = 1,2,3) restricted
      to each residue class mod 6 agree with a polynomial of degree <= j-1
      (exact Lagrange fit, checked to n = 72).
  C6  Positive control for the spine telescoping identity (S): the classical
      von Neumann pairs rule has W_m = 1 (m odd), q/2 (m even), V_m = 0
      (m odd), q + p/2 (m even), restart value pq/(1-p^2-q^2) = 1/2, and
      sum_m p^m q W_m + sum_m q^m p V_m == 1/2 exactly (symbolic).
  C7  d = 0 exhaustive: for every eventually periodic S (pre-period L <= 2,
      period P <= 8) the companion T forced by q*S(p) + p*T(q) = 1/2 fails to
      have {0,1} power-series coefficients (checked exactly to order 60).
  C8  d = 1 exhaustive: for every purely periodic (W_m) with period P <= 5,
      W_m in W_1 = {0, p, q, 1}, the forced companion series
      G(u) = 1/2 - F(1-u) is NOT representable as sum_m u^m (1-u) V_m(u) with
      V_m in {0, u, 1-u, 1}: exact sliding-window DP dies before depth 60.
  C9  Newton backward-propagation: g(t) = sum_i Delta^i g(T) * C(t-T, i)
      (symbolic identity, deg <= 5) and C(x,i) is an integer for ALL integers
      x (spot check incl. negatives) -- the lemma that integer values on a
      tail of an arithmetic progression propagate to the whole progression.

Exit code 0 and final line "LANE-B REFEREE: ALL CHECKS PASS" on success.
"""

import sys
import itertools
from fractions import Fraction
from math import comb

import sympy as sp

p, q, u, t_sym = sp.symbols('p q u t')

FAIL = []


def report(name, ok, detail=""):
    print(f"[{'PASS' if ok else 'FAIL'}] {name}" + (f"  {detail}" if detail else ""))
    if not ok:
        FAIL.append(name)


# ---------------------------------------------------------------- Fraction poly helpers
def padd(a, b):
    n = max(len(a), len(b))
    return [ (a[i] if i < len(a) else Fraction(0)) + (b[i] if i < len(b) else Fraction(0)) for i in range(n) ]

def pmul(a, b):
    r = [Fraction(0)] * (len(a) + len(b) - 1) if a and b else [Fraction(0)]
    for i, ai in enumerate(a):
        if ai:
            for j, bj in enumerate(b):
                if bj:
                    r[i + j] += ai * bj
    return r

def pscal(c, a):
    return [Fraction(c) * x for x in a]

def ptrim(a):
    while len(a) > 1 and a[-1] == 0:
        a = a[:-1]
    return a

def pcompose_1mu(a):
    """a(1-u) as polynomial in u, exact (Horner)."""
    res = [Fraction(0)]
    one_minus = [Fraction(1), Fraction(-1)]
    for c in reversed(a):
        res = padd(pmul(res, one_minus), [Fraction(c)])
    return res

def series_div(num, den, N):
    """First N+1 coefficients of num/den as power series (den[0] != 0)."""
    assert den[0] != 0
    out = []
    for n in range(N + 1):
        s = num[n] if n < len(num) else Fraction(0)
        for j in range(1, min(n, len(den) - 1) + 1):
            s -= den[j] * out[n - j]
        out.append(s / den[0])
    return out


# ---------------------------------------------------------------- C1 decided-tree fact
def check_C1():
    for d in range(4):
        # target set: expand sum_k w_k p^(d-k) (1-p)^k as coeff tuple (deg d)
        basis = []
        for k in range(d + 1):
            poly = [Fraction(0)] * (d - k) + [Fraction(1)]     # p^(d-k)
            for _ in range(k):
                poly = pmul(poly, [Fraction(1), Fraction(-1)])  # * (1-p)
            poly = poly + [Fraction(0)] * (d + 1 - len(poly))
            basis.append(poly[:d + 1])
        target = set()
        for w in itertools.product(*[range(comb(d, k) + 1) for k in range(d + 1)]):
            acc = [Fraction(0)] * (d + 1)
            for k, wk in enumerate(w):
                acc = padd(acc, pscal(wk, basis[k]))
            target.add(tuple(acc[:d + 1] + [Fraction(0)] * (d + 1 - len(acc))))

        # all stopping rules of depth <= d:  Leaf(0), Leaf(1), Node(rule0, rule1)
        # G(Leaf(b)) = b ;  G(Node(r0, r1)) = p*G(r0) + (1-p)*G(r1)
        # represent G as coeff tuple padded to degree d
        def all_G(depth):
            if depth == 0:
                return {(Fraction(0),), (Fraction(1),)}
            sub = all_G(depth - 1)
            out = {(Fraction(0),), (Fraction(1),)}
            for g0 in sub:
                pg0 = pmul([Fraction(0), Fraction(1)], list(g0))
                for g1 in sub:
                    qg1 = pmul([Fraction(1), Fraction(-1)], list(g1))
                    out.add(tuple(padd(pg0, qg1)))
            return out

        got = set()
        for g in all_G(d):
            gg = list(g) + [Fraction(0)] * (d + 1 - len(g))
            got.add(tuple(gg[:d + 1]))
        report(f"C1 decided-tree fact d={d}", got == target,
               f"|rules G-set|={len(got)} |target|={len(target)}")


# ---------------------------------------------------------------- C2 Vandermonde
def check_C2():
    ok = True
    for d in range(9):
        for dp_ in range(d + 1):
            for k in range(d + 1):
                s = sum(comb(dp_, j) * comb(d - dp_, k - j)
                        for j in range(0, dp_ + 1) if 0 <= k - j <= d - dp_)
                if s != comb(d, k):
                    ok = False
    report("C2 Vandermonde refinement bound d<=8", ok)


# ---------------------------------------------------------------- C3 coefficient bound
def check_C3():
    ok = True
    worst = Fraction(0)
    for d in range(7):
        basis = []
        for k in range(d + 1):
            poly = [Fraction(0)] * (d - k) + [Fraction(1)]
            for _ in range(k):
                poly = pmul(poly, [Fraction(1), Fraction(-1)])
            basis.append(poly)
        for verts in itertools.product(*[(0, comb(d, k)) for k in range(d + 1)]):
            acc = [Fraction(0)]
            for k, wk in enumerate(verts):
                acc = padd(acc, pscal(wk, basis[k]))
            acc = pmul([Fraction(1), Fraction(-1)], acc)   # * (1-p)
            m = max(abs(c) for c in acc)
            worst = max(worst, Fraction(m, 3 ** d))
            if m > 2 * 3 ** d:
                ok = False
    report("C3 |coeff((1-p)W)| <= 2*3^d for d<=6", ok,
           f"max ratio |coeff|/3^d over vertices = {worst}")


# ---------------------------------------------------------------- C4 zeta6 lemma
def check_C4():
    x, y = sp.symbols('x y', real=True)
    sols = sp.solve([x**2 + y**2 - 1, (1 - x)**2 + y**2 - 1], [x, y])
    expect = {(sp.Rational(1, 2), -sp.sqrt(3) / 2), (sp.Rational(1, 2), sp.sqrt(3) / 2)}
    ok1 = set((sp.nsimplify(a), sp.nsimplify(b)) for a, b in sols) == expect
    zeta = sp.Rational(1, 2) + sp.I * sp.sqrt(3) / 2
    ok2 = sp.simplify(zeta**6 - 1) == 0
    ok3 = all(sp.simplify(zeta**k - 1) != 0 for k in range(1, 6))
    ok4 = sp.simplify((1 - zeta) - sp.conjugate(zeta)) == 0
    ok5 = sp.simplify((1 - zeta) - zeta**5) == 0
    phi6 = sp.cyclotomic_poly(6, p)
    ok6 = sp.expand(phi6 - (p**2 - p + 1)) == 0
    ok7 = sp.simplify(phi6.subs(p, zeta)) == 0
    report("C4 {|z|=1,|1-z|=1} = primitive 6th roots; 1-z = conj(z) = z^5; Phi6 roots",
           all([ok1, ok2, ok3, ok4, ok5, ok6, ok7]))


# ---------------------------------------------------------------- C5 Phi6 structure
def check_C5():
    ok1 = sp.expand(((1 - p)**2 - (1 - p) + 1) - (p**2 - p + 1)) == 0

    # series of 1/Phi6: period-6 pattern 1,1,0,-1,-1,0
    N = 72
    inv = series_div([Fraction(1)], [Fraction(1), Fraction(-1), Fraction(1)], N)
    pattern = [1, 1, 0, -1, -1, 0]
    ok2 = all(inv[n] == pattern[n % 6] for n in range(N + 1))
    ok2b = all(inv[n] == inv[n - 1] - inv[n - 2] for n in range(2, N + 1))

    # (b1*p + b0)/Phi6 coefficients = b0*a_n + b1*a_(n-1) symbolically
    b0, b1 = sp.symbols('b0 b1')
    ser = sp.series((b1 * p + b0) / (p**2 - p + 1), p, 0, 25).removeO()
    ok3 = True
    alpha = lambda n: pattern[n % 6] if n >= 0 else 0
    for n in range(25):
        cn = ser.coeff(p, n)
        if sp.simplify(cn - (b0 * alpha(n) + b1 * alpha(n - 1))) != 0:
            ok3 = False

    # 1/Phi6^j classwise polynomial of degree <= j-1, j = 1,2,3 (exact fit)
    ok4 = True
    for j in (1, 2, 3):
        den = [Fraction(1)]
        for _ in range(j):
            den = pmul(den, [Fraction(1), Fraction(-1), Fraction(1)])
        ser_j = series_div([Fraction(1)], den, N)
        for r in range(6):
            idx = [n for n in range(N + 1) if n % 6 == r]
            pts = idx[:j]                       # fit through first j points
            def lag(x):
                s = Fraction(0)
                for a_i in pts:
                    li = Fraction(1)
                    for b_i in pts:
                        if b_i != a_i:
                            li *= Fraction(x - b_i, a_i - b_i)
                    s += ser_j[a_i] * li
                return s
            if not all(ser_j[n] == lag(n) for n in idx):
                ok4 = False
    report("C5 Phi6(1-p)=Phi6(p); 1/Phi6 period-6; linear/Phi6 coeffs; quasi-poly classes",
           all([ok1, ok2, ok2b, ok3, ok4]))


# ---------------------------------------------------------------- C6 von Neumann control
def check_C6():
    qq = 1 - p
    restart = sp.cancel(p * qq / (1 - p**2 - qq**2) - sp.Rational(1, 2))
    ok0 = restart == 0
    # sum_{m odd} p^m q * 1  +  sum_{m even>=2} p^m q * (q/2)
    # + sum_{m odd} q^m p * 0 +  sum_{m even>=2} q^m p * (q + p/2)  == 1/2
    S_odd_p = p / (1 - p**2)
    S_even_p = p**2 / (1 - p**2)
    S_even_q = qq**2 / (1 - qq**2)
    total = qq * S_odd_p + qq * (qq / 2) * S_even_p + p * (qq + p / 2) * S_even_q
    ok1 = sp.cancel(sp.together(total - sp.Rational(1, 2))) == 0

    # independent numeric spot checks with exact Fractions
    ok2 = True
    for pv in (Fraction(1, 2), Fraction(1, 3), Fraction(7, 11), Fraction(97, 101)):
        qv = 1 - pv
        tot = (qv * pv / (1 - pv**2) + qv * qv / 2 * pv**2 / (1 - pv**2)
               + pv * (qv + pv / 2) * qv**2 / (1 - qv**2))
        if tot != Fraction(1, 2):
            ok2 = False
    report("C6 von Neumann positive control for identity (S)", all([ok0, ok1, ok2]))


# ---------------------------------------------------------------- C7 d=0 exhaustive
def check_C7():
    N = 60
    bad = []
    tested = 0
    for L in range(3):
        for P in range(1, 9):
            for bits in itertools.product((0, 1), repeat=L + P):
                s_pre, s_per = bits[:L], bits[L:]
                tested += 1
                # S(p) = sum_{i<=L} s_i p^i + p^L * sum_j sigma_j p^j / (1-p^P)
                den = [Fraction(1)] + [Fraction(0)] * (P - 1) + [Fraction(-1)]
                numS = [Fraction(0)] * (L + P + 1)
                for i, si in enumerate(s_pre, start=1):
                    if si:
                        numS[i] += 1
                numS = pmul(numS, den) if False else None
                # build num_S directly: (sum s_i p^i)*(1-p^P) + p^L * sum sigma_j p^j
                pre_poly = [Fraction(0)] * (L + 1)
                for i, si in enumerate(s_pre, start=1):
                    pre_poly[i] = Fraction(si)
                per_poly = [Fraction(0)] * (L + P + 1)
                for j, sj in enumerate(s_per, start=1):
                    per_poly[L + j] = Fraction(sj)
                numS = padd(pmul(pre_poly, den), per_poly)
                # R(p) = 1/2 - (1-p) S(p)  =  (den/2 - (1-p) numS)/den
                numR = padd(pscal(Fraction(1, 2), den),
                            pscal(-1, pmul([Fraction(1), Fraction(-1)], numS)))
                # T(u) = R(1-u) / (1-u) ; den_T(u) = den(1-u)*(1-u) both sides /u
                numT = pcompose_1mu(numR)
                denT = pmul(pcompose_1mu(den), [Fraction(1), Fraction(-1)])
                if ptrim(numT) == [Fraction(0)]:
                    tser = [Fraction(0)] * (N + 1)
                else:
                    # both must vanish at u=0 (cancel factor u)
                    if numT[0] != 0 or denT[0] != 0:
                        bad.append((L, P, bits, "no common factor u"))
                        continue
                    tser = series_div(numT[1:], denT[1:], N)
                ok_coeffs = tser[0] == 0 and all(c in (0, 1) for c in tser[:N + 1])
                if ok_coeffs:
                    bad.append((L, P, bits, "companion T valid to order 60"))
    report("C7 d=0: no eventually periodic S (L<=2, P<=8) admits a valid companion T",
           len(bad) == 0, f"tested {tested} configs, flags = {bad[:3]}")


# ---------------------------------------------------------------- C8 d=1 DP search
def check_C8():
    N = 60
    # vectors (1-u)*V for V in {0, u, 1-u, 1}: coeffs (c0,c1,c2)
    OPTIONS = [(Fraction(0), Fraction(0), Fraction(0)),
               (Fraction(0), Fraction(1), Fraction(-1)),
               (Fraction(1), Fraction(-2), Fraction(1)),
               (Fraction(1), Fraction(-1), Fraction(0))]
    flags = []
    tested = 0
    Wopts = [(0, 0), (0, 1), (1, 0), (1, 1)]      # (w0, w1): W = w0*p + w1*q
    for P in range(1, 6):
        for cfg in itertools.product(Wopts, repeat=P):
            tested += 1
            den = [Fraction(1)] + [Fraction(0)] * (P - 1) + [Fraction(-1)]
            numF = [Fraction(0)]
            for r, (w0, w1) in enumerate(cfg, start=1):
                # (1-p) * (w0 p + w1 (1-p)) * p^r
                Wp = padd(pscal(w0, [Fraction(0), Fraction(1)]),
                          pscal(w1, [Fraction(1), Fraction(-1)]))
                term = pmul([Fraction(1), Fraction(-1)], Wp)
                term = pmul(term, [Fraction(0)] * r + [Fraction(1)])
                numF = padd(numF, term)
            # R = 1/2 - F ; G(u) = R(1-u) with den(1-u); cancel factor u
            numR = padd(pscal(Fraction(1, 2), den), pscal(-1, numF))
            numG = pcompose_1mu(numR)
            denG = pcompose_1mu(den)
            if numG[0] != 0 or denG[0] != 0:
                flags.append((P, cfg, "no factor u"))
                continue
            g = series_div(numG[1:], denG[1:], N)
            if g[0] != 0 or any(c != int(c) for c in g):
                continue   # instantly infeasible (needs g_0 = 0, integer coeffs)
            g = [int(c) for c in g]
            # DP: states (pend1, pend2); position m finalized at step m
            states = {(0, 0)}
            depth_reached = 0
            for m in range(1, N + 1):
                new = set()
                for (p1, p2) in states:
                    for (v0, v1, v2) in OPTIONS:
                        if p1 + v0 == g[m]:
                            new.add((p2 + v1, v2))
                states = new
                if not states:
                    break
                depth_reached = m
            if states:
                flags.append((P, cfg, f"DP alive at depth {N}"))
    report("C8 d=1: no purely periodic (W_m), P<=5, admits ANY companion (V_m) to depth 60",
           len(flags) == 0, f"tested {tested} configs, flags = {flags[:3]}")


# ---------------------------------------------------------------- C9 Newton backward
def check_C9():
    ok1 = True
    T0 = sp.Symbol('T0', integer=True)
    for deg in range(6):
        cs = sp.symbols(f'c0:{deg + 1}')
        g = sum(cs[i] * t_sym**i for i in range(deg + 1))
        # forward differences at T0
        def delta(f, i):
            for _ in range(i):
                f = f.subs(t_sym, t_sym + 1) - f
            return f.subs(t_sym, T0)
        newton = sum(delta(g, i) * sp.binomial(t_sym - T0, i) for i in range(deg + 1))
        if sp.simplify(sp.expand(newton - g)) != 0:
            ok1 = False
    ok2 = all(float(sp.binomial(x, i)).is_integer()
              for x in range(-8, 9) for i in range(6))
    report("C9 Newton backward-propagation identity (deg<=5) and integer C(x,i), x in Z",
           ok1 and ok2)


def main():
    check_C1()
    check_C2()
    check_C3()
    check_C4()
    check_C5()
    check_C6()
    check_C7()
    check_C8()
    check_C9()
    if FAIL:
        print(f"LANE-B REFEREE: {len(FAIL)} CHECK(S) FAILED: {FAIL}")
        sys.exit(1)
    print("LANE-B REFEREE: ALL CHECKS PASS")


if __name__ == "__main__":
    main()
