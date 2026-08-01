"""
FC(2)-homogeneous <=> the Lebesgue polynomial moment problem (Pakovich PMP, weight 1).

L(f^m)=int_[0,inf)^2 f^m e^{-x-y} for homogeneous f of degree d equals (dm+1)! int_0^1 phi^m da,
phi(a)=f(a,1-a). So all moments vanish <=> int_0^1 phi^m da = 0 for all m; and phi ranges over all
polynomials of degree <= d (Bernstein basis). No non-constant solution => FC(2)-homogeneous holds.
Exact through degree 4. WARNING: double-precision int phi^m via coefficient sums catastrophically
cancels and fakes solutions; use exact / high precision.
"""
import mpmath as mp
import sympy as sp

# --- 1. verify the reduction on a sample homogeneous f (deg 2) ---
mp.mp.dps = 30
def Lfm(m):
    return mp.quad(lambda x: mp.quad(lambda y: (x*x-3*x*y+y*y)**m*mp.e**(-x-y), [0, mp.inf]), [0, mp.inf])
def red(m, d=2):
    return mp.factorial(d*m+1)*mp.quad(lambda a: (a*a-3*a*(1-a)+(1-a)**2)**m, [0, 1])
print("reduction L(f^m) = (dm+1)! int_0^1 phi^m da  (f = x^2-3xy+y^2):")
for m in [1, 2, 3]:
    print("  m=%d: L=%s  reduction=%s  match=%s" % (m, mp.nstr(Lfm(m), 10), mp.nstr(red(m), 10), mp.nstr(Lfm(m)-red(m), 3)))

# --- 2. EXACT: no non-constant deg-d solution of int_0^1 phi^m da = 0 (Lebesgue PMP) ---
a = sp.symbols('a')
print("\nLebesgue PMP int_0^1 phi^m da = 0 (leading coeff = 1), exact solve:")
for d in [2, 3, 4]:
    cs = sp.symbols('c0:%d' % d)
    phi = sum(cs[k]*a**k for k in range(d)) + a**d
    eqs = [sp.expand(sp.integrate(sp.expand(phi**m), (a, 0, 1))) for m in range(1, d+3)]
    sol = sp.solve(eqs, list(cs), dict=True)
    print("  deg %d: %s" % (d, sol if sol else "NONE (no non-constant solution => FC(2)-homogeneous holds)"))

# --- 3. the catastrophic-cancellation guard, demonstrated ---
print("\nGUARD: a phi whose double-precision moments fake ~1e-12 but are truly ~1e-6:")
c = [mp.mpc('-0.00995', '0.00147'), mp.mpc('0.08791', '0.11067'), mp.mpc('-0.03275', '-0.53493'),
     mp.mpc('-0.28012', '0.69883'), mp.mpc('0.23471', '-0.26604')]
mp.mp.dps = 50
def phi_val(x): return sum(c[k]*x**k for k in range(5))
print("  high-precision |int_0^1 phi^1 da| =", mp.nstr(abs(mp.quad(lambda x: phi_val(x), [0, 1])), 6),
      " (double precision had reported ~6.8e-12 for the full residual -- WRONG)")
