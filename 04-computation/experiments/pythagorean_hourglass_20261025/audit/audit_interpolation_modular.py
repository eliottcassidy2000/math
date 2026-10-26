# Audit checks D (re-solve the interpolation systems from the audit's own data) and G (HYP-9390 and THM-4780 (6) with
# modular data from PARI/GP's mf package, modular_coeffs_gp.out, not from eta/theta products).  Audit 2026-10-25.
import re, sys
from fractions import Fraction
import numpy as np
from sympy import primerange, Matrix
sys.path.insert(0, '.')
from audit_graph_sums import S_norm, W4, SIX

def load_mine():
    tt5, tt6, pp = {}, {}, {}
    for l in open('audit_tt_values.out'):
        m = re.match(r'p=(\d+) TT5=(\d+)(?: TT6=(\d+))?', l)
        if m:
            p = int(m.group(1)); tt5[p] = int(m.group(2))
            if m.group(3) is not None: tt6[p] = int(m.group(3))
        m = re.match(r'q=(\d+) TT5=(\d+) TT6=(\d+)', l)
        if m: pp[int(m.group(1))] = (int(m.group(2)), int(m.group(3)))
    return tt5, tt6, pp
TT5, TT6, PP = load_mine()
gp = {}
for l in open('modular_coeffs_gp.out'):
    m = re.match(r'p=(\d+) a8=(-?\d+) bre=(-?\d+) bim=(-?\d+)', l)
    if m: gp[int(m.group(1))] = (int(m.group(2)), int(m.group(3)), int(m.group(4)))

out = []
# ---------------------------------------------------------------- D: re-solve the interpolation systems
P5 = [p for p in primerange(7, 200) if p % 4 == 3]
SW4 = {p: S_norm(p, 5, W4) for p in P5}; SW4[3] = 0            # T_W4(3) = 0: no injective 5-tuples in F_3
P6 = [p for p in primerange(7, 104) if p % 4 == 3]
SKP = {p: S_norm(p, 6, SIX['K6-P6']) for p in P6}; SKP[3] = 0
def solve_exact(rows, rhs):
    A = Matrix(rows); b = Matrix(rhs)
    rk = A.rank()
    aug = A.row_join(b)
    if aug.rank() != rk: return rk, None
    sol = (A.T * A).LUsolve(A.T * b)
    assert A * sol == b
    return rk, sol
rows = [[q ** i for i in range(6)] + [q * (q - 1) * SW4[q]] for q in [3] + P5]
rk, sol = solve_exact(rows, [1024 * TT5[q] for q in [3] + P5])
out.append(f'D. TT5: {len(rows)} equations, 7 unknowns, rank {rk}; solution {list(sol) if sol is not None else None}')
rows = [[q ** i for i in range(7)] + [q * (q - 1) * SW4[q], q * q * (q - 1) * SW4[q], q * (q - 1) * SKP[q]] for q in [3] + P6]
rk, sol = solve_exact(rows, [32768 * TT6[q] for q in [3] + P6])
out.append(f'D. TT6: {len(rows)} equations, 10 unknowns, rank {rk}; solution {list(sol) if sol is not None else None}')
# the q = 3 row is needed? drop it and check rank
rows_no3 = rows[1:]
out.append(f'   rank without the q = 3 row: {Matrix(rows_no3).rank()}')

# ---------------------------------------------------------------- the symbolic formulas, evaluated
def F5(q, s): return q * (q - 1) * ((q - 3) * (q * q - 21 * q + 105) + s)
def F6(q, sw, skp): return q * (q - 1) * ((q - 3) * (q ** 3 - 46 * q * q + 682 * q - 3239) + 7 * (q - 15) * sw - 8 * skp)
ok5 = all(F5(q, SW4[q]) == 1024 * TT5[q] for q in [3] + P5)
ok6 = all(F6(q, SW4[q], SKP[q]) == 32768 * TT6[q] for q in [3] + P6)
out.append(f'THM-4780 (5) with the audit S values: TT5 at {len(P5) + 1} q: {ok5}; TT6 at {len(P6) + 1} q: {ok6}')

# ---------------------------------------------------------------- G: modular forms (PARI)
def a8(p): return gp[p][0]
def absb2(p): return gp[p][1] ** 2 + gp[p][2] ** 2

def HYP9390(q, aq, b2): return q * (q - 1) * (q ** 4 - 49 * q ** 3 + 828 * q * q - 5313 * q + 10193 - 7 * (q - 15) * aq - 8 * b2)
def TT5mod(q, aq): return q * (q - 1) * (q ** 3 - 24 * q * q + 168 * q - 319 - aq)
bad5 = [p for p in TT5 if p in gp and TT5mod(p, a8(p)) != 1024 * TT5[p]]
bad6 = [p for p in TT6 if p in gp and HYP9390(p, a8(p), absb2(p)) != 32768 * TT6[p]]
out.append(f'G. THM-4780 (6) TT5(P_p) = p(p-1)(p^3-24p^2+168p-319-a_p)/2^10 with PARI a_p: tested at {sum(1 for p in TT5 if p in gp)} primes '
           f'(3..1499, p = 3 mod 4), failures {bad5}')
out.append(f'G. HYP-9390 (a) with PARI a_p and |b(p)|^2: tested at {sum(1 for p in TT6 if p in gp)} primes (3..751), failures {bad6}')
# algebra of the substitution, symbolically
import sympy as sp
q, sw, skp, aq, B2 = sp.symbols('q S_W4 S_KP a_q B2')
expr = (q - 3) * (q ** 3 - 46 * q ** 2 + 682 * q - 3239) + 7 * (q - 15) * sw - 8 * skp
sub = sp.expand(expr.subs({sw: -4 - aq, skp: B2 - q ** 2 - 7}))
tgt = sp.expand(q ** 4 - 49 * q ** 3 + 828 * q ** 2 - 5313 * q + 10193 - 7 * (q - 15) * aq - 8 * B2)
out.append(f'G. algebra THM-4780 (5) -> HYP-9390 (a): {sp.simplify(sub - tgt) == 0}')
expr5 = sp.expand(((q - 3) * (q ** 2 - 21 * q + 105) + sw).subs(sw, -4 - aq))
out.append(f'   algebra THM-4780 (5) -> (6): {sp.expand(expr5 - (q ** 3 - 24 * q ** 2 + 168 * q - 319 - aq)) == 0}')
Mq = sp.symbols('M_q')
exprM = sp.expand(expr.subs(skp, Mq + 2 * q - 6))
out.append(f'   M_q form (q-3)(q^3-46q^2+682q-3255) + 7(q-15)S_W4 - 8M_q: '
           f'{sp.expand(exprM - ((q - 3) * (q ** 3 - 46 * q ** 2 + 682 * q - 3255) + 7 * (q - 15) * sw - 8 * Mq)) == 0}')

# prime powers: Frobenius traces from PARI's a_p and b(p)
def power_sum(s1, prod, r):
    s = [2, s1]
    for _ in range(2, r + 1): s.append(s1 * s[-1] - prod * s[-2])
    return s[r]
for (p, r) in [(3, 3), (3, 5), (7, 3), (11, 3), (3, 7)]:
    Q = p ** r
    Aq = power_sum(a8(p), p ** 3, r)
    # b(p) = i h(p) at p = 3 mod 4; Satake roots of b: alpha + beta = b(p), alpha beta = chi_-4(p) p^2 = -p^2.
    # Trace of Frob^r: alpha^r + beta^r = i^r (u^r + v^r), u + v = h, u v = p^2.
    h = gp[p][2]; assert gp[p][1] == 0
    hr = power_sum(h, p * p, r)
    # cross-check with complex Satake roots
    al, be = np.roots([1, -1j * h, -p * p])
    tr = al ** r + be ** r
    assert abs(abs(tr) ** 2 - hr * hr) < 1e-6 * hr * hr + 1e-6
    t5, t6 = PP[Q]
    ok5 = TT5mod(Q, Aq) == 1024 * t5
    ok6 = HYP9390(Q, Aq, hr * hr) == 32768 * t6
    # Sym^2 trace S_q of HYP-9276 (2) and the relation |b(q)|^2 = q^2 - S_q
    Sq = (al ** (2 * r) + (al * be) ** r + be ** (2 * r)).real
    out.append(f'G. q = {p}^{r} = {Q}: A_q = {Aq}, |b(q)|^2 := h_r^2 = {hr * hr}; TT5 formula {ok5}, TT6 formula {ok6}; '
               f'Sym^2 trace S_q = {round(Sq)}, q^2 - S_q = {round(Q * Q - Sq)}, S_q - q^2 = {round(Sq - Q * Q)}')
print('\n'.join(out))
with open('audit_interpolation_modular.out', 'w') as fh: fh.write('\n'.join(out) + '\n')
