"""Near-join with an arbitrary tail D: T = (A => D) + t, t beats only t0 in A.  Claims:
 (1) sum over subsets S of V(D) of (-1)^|S| H(D[S]) H(D[V-S]) = 0 for nonempty D  (so F_e(D) = F_o(D));
 (2) H(T) = h(A) h(D) + B(A,t0) * P_D(2),  P_D(2) := sum over S of H(D[S]) H(D[V-S])  (= 2^k for transitive D)."""
import random, sys, itertools
sys.path.insert(0, '.')
from w45_hk import ham_split, sub
def H(n, out, keep):
    return ham_split(*sub(out, keep))[0] if keep else 1
def rt(n, rng):
    out = [0] * n
    for i in range(n):
        for j in range(i + 1, n):
            if rng.random() < .5: out[i] |= 1 << j
            else: out[j] |= 1 << i
    return out
def Bsum(n, out, t):
    tot = 0
    def ext(p, used):
        nonlocal tot
        tot += H(n, out, [v for v in range(n) if not used >> v & 1])
        for w in range(n):
            if not used >> w & 1 and out[p[-1]] >> w & 1: ext(p + [w], used | 1 << w)
    ext([t], 1 << t); return tot
rng = random.Random(56); bad1 = bad2 = 0
for trial in range(120):
    m = rng.randint(1, 7); D = rt(m, rng); V = list(range(m))
    alt = 0; P = 0
    for r in range(m + 1):
        for S in itertools.combinations(V, r):
            Sc = [v for v in V if v not in S]
            x = H(m, D, list(S)) * H(m, D, Sc); P += x; alt += (-1) ** r * x
    bad1 += alt != 0
    a = rng.randint(1, 5); A = rt(a, rng); t0 = rng.randrange(a)
    n = a + m + 1; out = [0] * n
    for i in range(a): out[i] = A[i]
    for i in range(m): out[a + i] = D[i] << a
    for i in range(a):
        for j in range(m): out[i] |= 1 << (a + j)
    t = n - 1
    for v in range(n - 1):
        if v == t0: out[t] |= 1 << v
        else: out[v] |= 1 << t
    lhs = ham_split(n, out)[0]
    rhs = ham_split(a, A)[0] * ham_split(m, D)[0] + Bsum(a, A, t0) * P
    bad2 += lhs != rhs
print("alternating identity failures:", bad1, "of 120 | tail formula h(A)h(D) + B*P_D(2) failures:", bad2, "of 120")
