from fractions import Fraction as Fr
from itertools import product, accumulate
from collections import defaultdict
import sympy as sp
def f(word, x=Fr(-1)):
    for c in word: x = (3*x + 1) / 2**c
    return x
# (1) the sporadic collision, symbolic in the last letter a
a = sp.symbols('a', integer=True, positive=True)
def fsym(word, x):
    for c in word: x = (3*x + 1) / 2**c
    return sp.simplify(x)
lhs = fsym((2, 2, 10, a), sp.Integer(-1)); rhs = fsym((6, 3, 2, 1, a + 2), sp.Integer(-1))
print("f_(2,2,10,a)(-1) =", sp.simplify(lhs), "; equal to f_(6,3,2,1,a+2)(-1):", sp.simplify(lhs - rhs) == 0, "; equals 8207/2^(13+a):", sp.simplify(lhs - sp.Integer(8207)/2**(13 + a)) == 0)
# (2) the real switch n = 53803, m = 26901
def U(n):
    n = 3*n + 1; e = 0
    while n % 2 == 0: n //= 2; e += 1
    return n, e
def orbit(n, k):
    ws = []; xs = [n]
    for _ in range(k):
        n, e = U(n); ws.append(e); xs.append(n)
    return xs, ws
for n in (53803, 26901):
    xs, ws = orbit(n, 6); print(n, "orbit", xs, "exponents", ws)
n = 53803; m = (n + 1)//2 - 1
print("m = (n+1)/2 - 1 =", m, "; U^5 equal:", orbit(n, 5)[0][-1] == orbit(m, 5)[0][-1], "; first meeting index:",
      next(j for j in range(1, 50) if orbit(n, j)[0][-1] == orbit(m, j)[0][-1]))
# (3) census: N-values of prefixes (reduced: first letter >= 2), prefix length <= 4 (word length <= 5)
def N_of(pre):
    p = len(pre) + 1
    S = list(accumulate(pre))
    return -3**(p-1) + sum(3**(p-1-i) * 2**(S[i-1] - 1) for i in range(1, p))
# sanity: f_w(-1) = N/2^(A-1)
for w in [(3,), (2, 1), (8, 5), (4, 1, 1, 7), (2, 2, 10, 4), (6, 3, 2, 1, 6)]:
    assert f(w) == Fr(N_of(w[:-1]), 2**(sum(w) - 1)), w
def root_partner(pre):
    if len(pre) >= 1 and pre[0] >= 3: return (2, pre[0] - 2) + pre[1:]
    if len(pre) >= 2 and pre[0] == 2: return (pre[1] + 2,) + pre[2:]
    return None
for L in (8, 22):
    groups = defaultdict(list)
    for k in range(0, 5):
        for pre in product(range(1, L + 1), repeat=k):
            if k and pre[0] < 2: continue
            groups[N_of(pre)].append(pre)
    spor = []
    for N, P in groups.items():
        S = set(P); seen = set(); ncls = 0
        for p in P:
            if p in seen: continue
            q = root_partner(p); seen.add(p)
            if q in S: seen.add(q)
            ncls += 1
        if ncls >= 2: spor.append((N, sorted(P)))
    print(f"letters <= {L}: sporadic N-values (>= 2 root classes): {len(spor)}")
    if L == 22:
        new = [s for s in spor if any(max(p, default=0) > 8 for p in s[1])]
        fam = [s for s in new if any(len(p) == 3 and p[:2] == (2, 6) for p in s[1])]
        print("   with a letter > 8:", len(new), "; of which in the (2,6,c)~(4,1,1,c+2) family:", len(fam))
        print("   the others:", [s for s in new if s not in fam])
