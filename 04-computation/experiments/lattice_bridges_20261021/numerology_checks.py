# 1729 and the cube pattern H(P_q)/|Aut(P_q)| - 1 = 0, 8, 1728 at q = 3, 7, 11 (wide-net spark 28); the cannonball identity
# for the transitive tournament (sum of squared scores of TT_n is a square only for n = 1, 2, 25); E_8 roots via A_8.
# mac-mini-2026-10-21-fanopsp.
from math import comb, isqrt
def icbrt(n):
    x = round(abs(n) ** (1 / 3))
    for y in (x - 1, x, x + 1):
        if y ** 3 == abs(n): return y if n >= 0 else -y
    return None
data = {3: (3, 3), 7: (189, 21), 11: (95095, 55), 19: (1172695746915, 171), 23: (15760206976379349, 253),
        27: (418425945688463432031, 1053), 31: (19876729650849384218015865, 465)}
for q, (H, aut) in data.items():
    v = H // aut; assert v * aut == H
    print(f"q={q}: H(P_q)/|Aut| = {v}; minus 1 = {v - 1}; integer cube root = {icbrt(v - 1)}")
print("1729 = 7*13*19 =", 7 * 13 * 19, "; 1729 = 1^3+12^3 = 9^3+10^3:", 1 + 12 ** 3 == 9 ** 3 + 10 ** 3 == 1729)
sq = [n for n in range(1, 10 ** 6) if isqrt((n - 1) * n * (2 * n - 1) // 6) ** 2 == (n - 1) * n * (2 * n - 1) // 6]
print("TT_n with sum of squared scores a perfect square, n < 10^6:", sq, " (n = 25: sum = 4900 = 70^2)")
print("E_8 roots via A_8: 9*8 =", 9 * 8, "+ 2*C(9,3) =", 2 * comb(9, 3), "=", 72 + 2 * comb(9, 3), "; dim e_8 = 80 + 84 + 84 =", 80 + 84 + 84)
