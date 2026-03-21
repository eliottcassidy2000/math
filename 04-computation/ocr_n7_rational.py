#!/usr/bin/env python3
"""Find exact rational OCR(7) from the C output sums."""
from math import gcd
from fractions import Fraction

N = 2097152
sum_H = 165150720
sum_S2 = 22020096
sum_H2 = 16386048000
sum_S22 = 286261248
sum_HS2 = 1321205760

cov_num = N * sum_HS2 - sum_H * sum_S2
var_S2 = N * sum_S22 - sum_S2 * sum_S2
var_H = N * sum_H2 - sum_H * sum_H

r2_num = cov_num ** 2
r2_den = var_S2 * var_H

g = gcd(r2_num, r2_den)
p = r2_num // g
q = r2_den // g

print(f"R^2 numerator (cov^2): {r2_num}")
print(f"R^2 denominator (var_S2 * var_H): {r2_den}")
print(f"GCD: {g}")
print(f"R^2 = {p} / {q}")
print(f"1 - R^2 = {q-p} / {q}")
print(f"OCR(7) = {p/q:.12f}")
print()

# Factor the denominator
def factor(n):
    factors = {}
    d = 2
    while d * d <= abs(n):
        while n % d == 0:
            factors[d] = factors.get(d, 0) + 1
            n //= d
        d += 1
    if abs(n) > 1:
        factors[abs(n)] = factors.get(abs(n), 0) + 1
    return factors

print(f"Numerator factors: {factor(p)}")
print(f"Denominator factors: {factor(q)}")
print(f"(q-p) factors: {factor(q-p)}")
print()

# The OCR sequence so far:
print("THE EXACT OCR SEQUENCE:")
print("  n=3: OCR = 1/1, 1-OCR = 0")
print("  n=4: OCR = 1/1, 1-OCR = 0")
print("  n=5: OCR = 18/19, 1-OCR = 1/19")
print("  n=6: OCR = 12/13, 1-OCR = 1/13")
print(f"  n=7: OCR = {p}/{q}, 1-OCR = {q-p}/{q}")
print()

# Also compute exact Var(H), Var(c3), Cov(H, c3)
E_H = Fraction(sum_H, N)
E_S2 = Fraction(sum_S2, N)
Var_H = Fraction(sum_H2, N) - E_H**2
Var_S2 = Fraction(sum_S22, N) - E_S2**2
Cov_HS2 = Fraction(sum_HS2, N) - E_H * E_S2

print(f"E[H] = {E_H} = {float(E_H)}")
print(f"E[S2] = {E_S2} = {float(E_S2)}")
print(f"Var(H) = {Var_H} = {float(Var_H)}")
print(f"Var(S2) = {Var_S2} = {float(Var_S2)}")
print(f"Cov(H, S2) = {Cov_HS2} = {float(Cov_HS2)}")
print()

R2_exact = Cov_HS2**2 / (Var_H * Var_S2)
print(f"R^2 = {R2_exact} = {float(R2_exact):.12f}")

# Check: Cov(H, S2) = -2*Cov(H, c3) = -2*(Cov(2c3+eps, c3))
# S2 = -2*c3 + const, so Cov(H, S2) = -2*Cov(H, c3)
# R^2(S2, H) = Cov(S2,H)^2 / (Var(S2)*Var(H))
# = 4*Cov(c3,H)^2 / (4*Var(c3)*Var(H))
# = Cov(c3,H)^2 / (Var(c3)*Var(H))
# = R^2(c3, H)  [consistent]

# The key: what is Var(H)/Var(c3)?
ratio = Var_H / (4 * Var_S2)  # Since Var(c3) = Var(S2)/4 (S2 = -2c3+const)
# Wait: S2 = -2*c3 + const. Var(S2) = 4*Var(c3). So Var(c3) = Var(S2)/4.
# R^2 = Cov^2/(Var_S2*Var_H) = Cov_S2_H^2/(Var_S2*Var_H)
# Alternative: R^2 = [Cov(c3,H)]^2 / (Var(c3)*Var(H))
#  where Cov(c3,H) = -Cov(S2,H)/2

Var_c3 = Var_S2 / 4
Cov_c3_H = -Cov_HS2 / 2
print(f"\nDerived: Var(c3) = {Var_c3}, Cov(c3,H) = {Cov_c3_H}")
print(f"Check R^2 = Cov(c3,H)^2/(Var(c3)*Var(H)) = {Cov_c3_H**2 / (Var_c3 * Var_H)}")
print()

# What fraction of Var(H) is explained by the linear term 2*c3?
# If H = 1 + 2*c3 + eps: Var(H) = 4*Var(c3) + 4*Cov(c3,eps) + Var(eps)
# And Cov(H, c3) = 2*Var(c3) + Cov(eps, c3)
Cov_eps_c3 = Cov_c3_H - 2*Var_c3
print(f"Cov(eps, c3) = Cov(c3,H) - 2*Var(c3) = {Cov_eps_c3}")
print(f"  = {float(Cov_eps_c3):.6f}")
print(f"  This is {'ZERO' if Cov_eps_c3 == 0 else 'NONZERO'} => eps {'is' if Cov_eps_c3 == 0 else 'is NOT'} uncorrelated with c3")
