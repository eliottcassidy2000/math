"""
The RADIAL bridge: [int_0^1 e^Q transcendental, nonconstant algebraic Q] => FC(2)-homogeneous.
(Distinct from HYP-9078's failed log-bridge: here the exponent Q=s*phi is a POLYNOMIAL in the integration var.)

FC(2)-homog = [0,1] Lebesgue PMP (Liu-Sun 2020). Moment EGF: sum_m (int_0^1 phi^m) s^m/m! = int_0^1 e^{s phi}.
If FC-homog fails (phi!=0, all moments 0) => int_0^1 e^{s phi}=1 for all s => at algebraic s0, int e^{s0 phi}=1
(rational) with nonconstant algebraic exponent => contradicts transcendence. So the claim => FC(2)-homog.
Full FC(2) (saddle-descent) needs int_0^1 e^{s phi_D + phi_{D-1}/(D phi_D)} = int e^{poly+RATIONAL} (weighted).
"""
import mpmath as mp
mp.mp.dps = 30

# EGF identity (the bridge is a polynomial-exponent integral)
phi = lambda a: a*a - a + mp.mpf('0.3')
s = mp.mpf('0.7')
lhs = mp.quad(lambda a: mp.e**(s*phi(a)), [0, 1])
rhs = mp.nsum(lambda m: s**m*mp.quad(lambda a: phi(a)**int(m), [0, 1])/mp.factorial(int(m)), [0, mp.inf])
print("radial bridge EGF: int_0^1 e^{s phi} da =", mp.nstr(lhs, 12),
      " = sum s^m<phi^m>/m! =", mp.nstr(rhs, 12), " (Q=s*phi is a POLYNOMIAL in a)")

# degree-1 transcendence (Lindemann-Weierstrass): int_0^1 e^{at+b} = e^b(e^a-1)/a, algebraic iff e^a=1 (a transc)
a, b = mp.mpc('0.5', '1.3'), mp.mpc('0.2', '-0.4')
print("deg-1: int_0^1 e^{at+b} =", mp.nstr(mp.e**b*(mp.e**a-1)/a, 10), " (transcendental unless e^a=1, i.e. a in 2pi i Z)")

# inhomogeneous weight: full FC(2) EGF = int_0^1 e^{s phi_D + phi_{D-1}/(D phi_D)} da (poly + rational exponent)
print("full FC(2): leak EGF = int_0^1 e^{s phi_D + phi_{D-1}/(D phi_D)} da = int e^{poly+RATIONAL} (weighted generalization)")
print("golden AMM constant C = 1 + 2 log(phi)/log 5 =", mp.nstr(1+2*mp.log((1+mp.sqrt(5))/2)/mp.log(5), 12), "(klein THM-3027 = opus two-ray)")
