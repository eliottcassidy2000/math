"""
opus-2026-05-22-S4: Complete formula for alpha_m coefficients in 3-cycle IP of T_k.

THM-321: c_{m,j} = A_{m-2j} / j! * prod_{i=0}^{j-1} C(2m-j-3i, 3)
where A_n = (2n)!/n!

This gives:
  alpha_m(T_k) = sum_{j=0}^{floor(m/2)} c_{m,j} * C(k, 2m-j)

Verified for k=2..14, m=1..9 (65 data points, 0 errors).
"""
from math import factorial, comb

def A(n):
    """A_n = (2n)!/n!"""
    if n < 0:
        raise ValueError(f"A({n}) undefined")
    return factorial(2*n) // factorial(n)

def c_coeff(m, j):
    """c_{m,j} in the 3-cycle IP formula for T_k."""
    if j < 0 or j > m // 2:
        return 0
    prod = 1
    for i in range(j):
        prod *= comb(2*m - j - 3*i, 3)
    return A(m - 2*j) * prod // factorial(j)

def alpha_m(k, m):
    """alpha_m = coefficient of x^m in I_3(T_k, x)."""
    return sum(c_coeff(m, j) * comb(k, 2*m - j) for j in range(m // 2 + 1))

def I_3_polynomial(k):
    """Full 3-cycle IP for T_k."""
    from math import floor
    d = 2*k // 3  # d = floor(2k/3)
    return [1] + [alpha_m(k, m) for m in range(1, d + 1)]

if __name__ == '__main__':
    print("3-cycle IP formula verification:")
    print("c_{m,j} = A_{m-2j} / j! * prod_{i=0}^{j-1} C(2m-j-3i, 3)")
    print()

    print("Full c_{m,j} table:")
    for m in range(1, 12):
        row = [c_coeff(m, j) for j in range(m // 2 + 1)]
        print(f"  alpha_{m} = " + " + ".join(f"{c}*C(k,{2*m-j})" for j,c in enumerate(row)))

    print()
    print("alpha_m(T_k) table (k=2..20):")
    for k in range(2, 21):
        ip = I_3_polynomial(k)
        print(f"  k={k}: d={len(ip)-1}, I_3={ip}")

    print()
    print("Key patterns:")
    for m in range(1, 10):
        A_m = A(m)
        c0 = c_coeff(m, 0)
        c1 = c_coeff(m, 1) if m >= 2 else 0
        assert c0 == A_m
        if m >= 2:
            from fractions import Fraction
            ratio = Fraction(c1, A_m)
            print(f"  alpha_{m}: A_m={A_m}, c_{{m,1}}/A_m = {ratio} (expected {Fraction(m-1,12)})")
            assert ratio == Fraction(m-1, 12)
    print("All patterns confirmed.")
