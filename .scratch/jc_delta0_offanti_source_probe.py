#!/usr/bin/env python3
"""Scratch-only exact probe of the Delta=0 off-anti exact-M=9 source."""

import sympy as sp


s, p = sp.symbols("s p")
Phi, Theta, eta, zeta = sp.symbols("Phi Theta eta zeta")
u = sp.symbols("u")
K = sp.Rational(2848, 45)
epsilon = -sp.Rational(1376, 135)
t = p - s**2
H = (
    -3*p + sp.Rational(8, 3)*p**2 + epsilon*p**3
    + K*s**2*p**2 + Phi*s*p**3 + Theta*s**2*p**3
    + eta*s*p**4 + zeta*s**3*p**3
)
G = -s**2/(2*t) + H
A = sp.cancel((-s*p + t**2*sp.diff(H, s))/p)
C0 = sp.expand(s**2 + 2*t**2*sp.diff(H, p))
B = sp.cancel((C0 + s*A)/t**2)
print("den", sp.denom(A), sp.denom(B), flush=True)
print("degrees", sp.degree(A, s), sp.degree(B, s), flush=True)
print("LC", sp.factor(sp.Poly(A, s).LC()), sp.factor(sp.Poly(B, s).LC()), flush=True)
res = sp.factor(sp.resultant(A, B, s))
v = min(k[0] for k in sp.Poly(res, p).as_dict())
R = sp.Poly(sp.cancel(res/p**v), p)
print("valuation_degree", v, R.degree(), flush=True)
print("TC", sp.factor(R.TC()), flush=True)
print("LC", sp.factor(R.LC()), flush=True)
# Full coefficients are intentionally suppressed; only endpoint and wall
# cascades are useful for this boundary audit.

# Exact chart of the first bottom-endpoint wall
# 4*Theta*K^2-27*zeta^2=0, with zeta!=0:
# Theta=3*u^2, zeta=2*K*u/3.
wall_sub = {Theta: 3*u**2, zeta: 2*K*u/3}
wall_res = sp.factor(res.subs(wall_sub))
wall_v = min(k[0] for k in sp.Poly(wall_res, p).as_dict())
wall_R = sp.Poly(sp.cancel(wall_res/p**wall_v), p)
print("wall_valuation_degree", wall_v, wall_R.degree(), flush=True)
print("wall_TC", sp.factor(wall_R.TC()), flush=True)
print("wall_LC", sp.factor(wall_R.LC()), flush=True)
for j in range(1, 5):
    print(f"wall_low{j}", sp.factor(wall_R.nth(j)), flush=True)

phi_J = sp.factor((1215*u**3 + 22784*u)/8544)
wall2_res = sp.factor(wall_res.subs(Phi, phi_J))
wall2_v = min(k[0] for k in sp.Poly(wall2_res, p).as_dict())
wall2_R = sp.Poly(sp.cancel(wall2_res/p**wall2_v), p)
print("wall2_valuation_degree", wall2_v, wall2_R.degree(), flush=True)
print("wall2_TC", sp.factor(wall2_R.TC()), flush=True)
print("wall2_LC", sp.factor(wall2_R.LC()), flush=True)
for j in range(1, 4):
    print(f"wall2_low{j}", sp.factor(wall2_R.nth(j)), flush=True)

S = 547499520*eta + u*(2460375*u**4 - 204543360*u**2 + 5580439552)
eta_S = sp.factor(-u*(2460375*u**4 - 204543360*u**2 + 5580439552)
                  / 547499520)
wall3_res = sp.factor(wall2_res.subs(eta, eta_S))
wall3_v = min(k[0] for k in sp.Poly(wall3_res, p).as_dict())
wall3_R = sp.Poly(sp.cancel(wall3_res/p**wall3_v), p)
print("wall3_valuation_degree", wall3_v, wall3_R.degree(), flush=True)
print("wall3_TC", sp.factor(wall3_R.TC()), flush=True)
print("wall3_LC", sp.factor(wall3_R.LC()), flush=True)
for j in range(1, 3):
    print(f"wall3_low{j}", sp.factor(wall3_R.nth(j)), flush=True)
T0 = 27064125*u**4 - 5739517440*u**2 + 47239069696
eta_zero = 2460375*u**4 - 204543360*u**2 + 5580439552
anti_zero = 18225*u**4 - 1515136*u**2 - 129777664
print("terminal_gcd_T0_low1", sp.gcd(sp.Poly(T0, u), sp.Poly(wall3_R.nth(1), u)), flush=True)
print("hostile_gcd_T0_eta0", sp.gcd(sp.Poly(T0, u), sp.Poly(eta_zero, u)), flush=True)
print("hostile_gcd_T0_anti", sp.gcd(sp.Poly(T0, u), sp.Poly(anti_zero, u)), flush=True)
