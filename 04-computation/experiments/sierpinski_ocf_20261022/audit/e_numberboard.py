# Audit of the note's section 7 (number board) and section 2 arithmetic.
from itertools import permutations, combinations
from math import comb, isqrt
def collatz_stats(n):
    steps = odd = 0
    while n != 1:
        if n % 2: n = 3 * n + 1; odd += 1
        else: n //= 2
        steps += 1
    return steps, odd, steps - odd
def U(n):  # Syracuse (odd-to-odd)
    m = 3 * n + 1
    while m % 2 == 0: m //= 2
    return m
print("27: (steps, odd, halvings) =", collatz_stats(27))
print("U(9) =", U(9), " U(21) =", U(21), " C(9) = 3*9+1 =", 3 * 9 + 1)
print("3^k = k^3 for integer 1<=k<=200:", [k for k in range(1, 201) if 3 ** k == k ** 3])
phi6 = lambda k: k * k - k + 1
print("Phi_6(3), (3^3+1)/4, Phi_6(5), (5^3+1)/6:", phi6(3), (27 + 1) // 4, phi6(5), (125 + 1) // 6,
      "; k^3+1 = (k+1)Phi_6(k) and Phi_6 = 1+2C(k,2) for k<100:",
      all(k ** 3 + 1 == (k + 1) * phi6(k) and phi6(k) == 1 + 2 * comb(k, 2) for k in range(100)))
print("111_2, 111_3, 111_4 =", int('111', 2), int('111', 3), int('111', 4), "; (3x+1)^3(0) =", 3 * (3 * (3 * 0 + 1) + 1) + 1)
# c_3 of a regular tournament on 2N+1 vertices = sum_{i<=N} i^2 ; Watson
ok = all(comb(2 * N + 1, 3) - (2 * N + 1) * comb(N, 2) == N * (N + 1) * (2 * N + 1) // 6 for N in range(0, 500))
sq = [N for N in range(0, 2 * 10 ** 6) if isqrt(N * (N + 1) * (2 * N + 1) // 6) ** 2 == N * (N + 1) * (2 * N + 1) // 6]
print("c_3(regular, 2N+1) = sum i^2 (N<500):", ok, "; sum_{i<=N} i^2 square for N < 2e6:", sq, "; sum_{i<=24} i^2 =", sum(i * i for i in range(25)), "= 70^2:", 70 ** 2)
# H values attained at n <= 6 (is 13 an H value?)
def H(n, A):
    return sum(all(A[p[i]][p[i + 1]] for i in range(n - 1)) for p in permutations(range(n)))
vals = {}
for n in (5, 6):
    E = [(i, j) for i in range(n) for j in range(i + 1, n)]
    s = set()
    for code in range(1 << len(E)):
        A = [[0] * n for _ in range(n)]
        for b, (i, j) in enumerate(E):
            if code >> b & 1: A[i][j] = 1
            else: A[j][i] = 1
        s.add(H(n, A))
    vals[n] = sorted(s)
print("H values n=5:", vals[5]); print("H values n=6:", vals[6])
# Section 2: profiles alpha for H = 7, 21 (alpha_k = # independent k-sets; alpha_2 <= C(alpha_1,2), alpha_3>=1 forces alpha_2>=3, alpha_1>=3)
def profiles(Hv):
    out = []
    for a1 in range(0, Hv):
        for a2 in range(0, comb(a1, 2) + 1):
            for a3 in range(0, comb(a1, 3) + 1):
                if a3 >= 1 and (a2 < 3 or a1 < 3): continue
                if 1 + 2 * a1 + 4 * a2 + 8 * a3 == Hv: out.append((a1, a2, a3))
    return out
print("profiles H=7:", profiles(7), "  H=21:", profiles(21))
print("3^3+1 = 28 = 7+21:", 28 == 7 + 21, "; 3^3-1 = 26 = 25+1:", 26 == 25 + 1)
