"""Wave 16: the terrace corridor kernel is the backward tribonacci sequence.

Checks (exact integer arithmetic; sympy only for Bernoulli numbers and factoring):
 1. K_q entries A=T(q+3)-T(q+1), B=T(q+3)-T(q), D=T(q+3) and det K_q = -T(-q-3).
 2. Zero set of the two-sided tribonacci sequence on [-6000, 6000] and a Skolem certificate at p=83
    (T != 0 mod p outside the zero classes; v_p(T(z+period)-T(z)) = 1 in each zero class).
 3. Dual/doubling identities with S_n=T(-n), the forced squares, and all GPs T_a T_c = T_b^2.
 4. Two-sided k-bonacci zeros: structured family versus sporadic zeros, k=3..15.
 5. Periods mod 3 and 5 (Singer cycles of PG(2,3), PG(2,5)).
 6. Eisenstein constants 4k/|B_2k| with numerator 2^a (2^2k - 1).
Run: python3 corridor_tribonacci.py
"""
from math import isqrt
from sympy import bernoulli, Rational, factorint

R = 6000
T = {0: 0, 1: 0, 2: 1}
for n in range(3, R + 40): T[n] = T[n-1] + T[n-2] + T[n-3]
for n in range(-1, -R - 40, -1): T[n] = T[n+3] - T[n+2] - T[n+1]
S = lambda n: T[-n]

def series(num, N):
    a = [0] * N
    for q in range(N): a[q] = (num[q] if q < len(num) else 0) + sum(a[q-k] for k in (1, 2, 3) if q - k >= 0)
    return a
A, B, D = series([1, 0, 1], 200), series([1, 1], 200), series([1, 1, 1], 200)
ok1 = all(A[q] == T[q+3]-T[q+1] and B[q] == T[q+3]-T[q] and D[q] == T[q+3] and A[q]*D[q]-B[q]**2 == -T[-q-3] for q in range(200))
print("1. entries and det K_q = -T(-q-3) for 0<=q<200:", ok1, "| singular q:", [q for q in range(200) if A[q]*D[q] == B[q]**2])
print("   K_13 =", (A[13], B[13], D[13]), "= (47^2, 47*56, 56^2):", (A[13], B[13], D[13]) == (47**2, 47*56, 56**2))

zeros = [n for n in range(-R, R) if T[n] == 0]
p = 83; st = (0, 0, 1); s = st; per = 0
while True:
    s = (s[1], s[2], (s[0]+s[1]+s[2]) % p); per += 1
    if s == st: break
res = {z % per for z in zeros}
other = [r for r in range(per) if T[r] % p == 0 and r not in res]
def vp(x, p):
    v = 0
    while x % p == 0: x //= p; v += 1
    return v
lin = {z: vp(T[z+per]-T[z], p) for z in zeros}
print(f"2. zeros of T on [-{R},{R}):", zeros, f"| p={p} period {per}, extra zero residues {other}, v_p(T(z+per)-T(z)) {lin}")
print("   certificate valid:", len(res) == len(zeros) and not other and all(v == 1 for v in lin.values()))

idents = {
 "T(q-1)T(q+1)-T(q)^2 = -T(1-q)": lambda q: T[q-1]*T[q+1]-T[q]**2 == -T[1-q],
 "T(q-3)T(q+3)-T(q)^2 = -T(-q-3)": lambda q: T[q-3]*T[q+3]-T[q]**2 == -T[-q-3],
 "T(q-4)T(q+4)-T(q)^2 = -4T(-q-2)": lambda q: T[q-4]*T[q+4]-T[q]**2 == -4*T[-q-2],
 "S(q-1)S(q+1)-S(q)^2 = -T(q+1)": lambda q: S(q-1)*S(q+1)-S(q)**2 == -T[q+1],
 "S(q-3)S(q+3)-S(q)^2 = -T(q-3)": lambda q: S(q-3)*S(q+3)-S(q)**2 == -T[q-3],
 "S(q-4)S(q+4)-S(q)^2 = -4T(q-2)": lambda q: S(q-4)*S(q+4)-S(q)**2 == -4*T[q-2],
 "S(2n) = S(n-1)^2+S(n)^2+2S(n)S(n+1)": lambda n: S(2*n) == S(n-1)**2+S(n)**2+2*S(n)*S(n+1),
}
for name, f in idents.items(): print("3.", name, all(f(q) for q in range(-400, 400)))
print("   forced squares: T10=S13^2", T[10] == S(13)**2, " T16=S15^2", T[16] == S(15)**2, " T18=S17^2", T[18] == S(17)**2,
      " S30=47^2+56^2", S(30) == 47**2+56**2)
print("   square T_n, 0<=n<6000:", [(n, T[n]) for n in range(R) if isqrt(T[n])**2 == T[n]])
print("   GPs T_a T_c = T_b^2, 3<=a<b<c<300:", [(a, b, c) for b in range(4, 300) for a in range(3, b) for c in range(b+1, 300) if T[a]*T[c] == T[b]**2])
sm = []
for n in range(2, 3000):
    x = T[n]
    for q in (2, 3, 7):
        while x % q == 0: x //= q
    if x == 1: sm.append(T[n])
print("   {2,3,7}-smooth T_n, n<3000:", sm)

def kzeros(k, lo):
    a = {i: 0 for i in range(k-1)}; a[k-1] = 1
    for n in range(k, k+3): a[n] = sum(a[n-j] for j in range(1, k+1))
    for n in range(-1, lo-1, -1): a[n] = a[n+k] - sum(a[n+j] for j in range(1, k))
    z = {n for n in a if a[n] == 0 and n < 0}
    struct = {-3-j*(k+1)-i for j in range(0, k-2) for i in range(0, k-2-j)}
    return sorted(z - struct), sorted(struct - z), len(z)
print("4. k-bonacci negative zeros (sporadic, missing structured, count), index >= -4000:")
for k in range(3, 16): print("   k =", k, kzeros(k, -4000))

def period(p):
    st = (0, 0, 1); s = st; n = 0
    while True:
        s = (s[1], s[2], (s[0]+s[1]+s[2]) % p); n += 1
        if s == st: return n
print("5. tribonacci periods mod 2,3,5,7,13:", {q: period(q) for q in (2, 3, 5, 7, 13)})

hits = []
for k2 in range(2, 61, 2):
    c = abs(Rational(2*k2) / bernoulli(k2)); M = 2**k2 - 1; y = c.p // M if c.p % M == 0 else 0
    if y and y & (y-1) == 0: hits.append((k2, c))
print("6. 2k<=60 with numerator of 4k/|B_2k| equal to 2^a(2^2k-1):", hits)
