"""Wave 59.  Ordinal/C3 class D: species equation A = S/(1-S), S = x + (A^3 + 2A(x^3))/3 (session284_structural_births, credited)."""
N = 18
from fractions import Fraction as F
A = [F(0)] * (N + 1); S = [F(0)] * (N + 1)
def mul(P, Q):
    R = [F(0)] * (N + 1)
    for i, a in enumerate(P):
        if a:
            for j, b in enumerate(Q[:N + 1 - i]):
                if b: R[i + j] += a * b
    return R
for it in range(N + 2):
    A3 = mul(mul(A, A), A); Ax3 = [F(0)] * (N + 1)
    for i in range(N + 1):
        if 3 * i <= N: Ax3[3 * i] = A[i]
    S = [F(0)] * (N + 1); S[1] = F(1)
    for i in range(N + 1): S[i] += (A3[i] + 2 * Ax3[i]) / 3
    newA = [F(0)] * (N + 1); P = S[:]
    for k in range(1, N + 1):
        for i in range(N + 1): newA[i] += P[i]
        P = mul(P, S)
    A = newA
print("ordinal/C3 classes by order 1..18:", [int(a) for a in A[1:]])
print("strong (C3-rooted) ones:          ", [int(s) for s in S[1:]])
