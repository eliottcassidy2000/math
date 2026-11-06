import numpy as np, cmath, math
# Reading v's bits after the top 1, MSB first. State (P mod q, last bit). Bit b=0 adds P (the prefix value before the bit,
# = floor(v/2^(i+1))) and adds 1 if last bit was 1 (#10); g(v) = 1 + total.  Twisted operator T_a: weight e(a*add/q).
def T(q, a):
    n = 2 * q; M = np.zeros((n, n), dtype=complex)
    for P in range(q):
        for l in (0, 1):
            i = 2 * P + l
            for b in (0, 1):
                add = (P + (1 if l == 1 else 0)) if b == 0 else 0
                j = 2 * ((2 * P + b) % q) + b
                M[j, i] += cmath.exp(2j * math.pi * a * add / q)
    return M
def S_exact(q, a, k):
    M = T(q, a); vec = np.zeros(2 * q, dtype=complex); vec[2 * (1 % q) + 1] = 1      # top bit: P = 1, last = 1
    for _ in range(k): vec = M @ vec
    return abs(cmath.exp(2j * math.pi * a / q) * vec.sum()) / 2 ** k
for q, a in [(3, 1), (5, 1), (5, 2), (7, 1), (7, 2), (7, 3), (9, 1), (21, 1), (49, 1), (11, 1), (13, 1), (31, 1), (73, 1)]:
    ev = np.linalg.eigvals(T(q, a)); top = max(ev, key=abs)
    print(f"q={q:2d} a={a}: |lambda|/2 = {abs(top)/2:.4f}, arg/(2pi) = {cmath.phase(top)/(2*math.pi):+.4f}; exact |S|/V at k=12,20: {S_exact(q,a,12):.4f}, {S_exact(q,a,20):.4f}")
