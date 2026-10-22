# Check the OCF form of the alternating forward-step count:
#   F(T,-1) = sum_sigma (-1)^(#backward steps) = (-1)^((n-1)/2) sum_S prod_(c in S) [(-1)^((|c|-1)/2) 2^|c|] E_(n - phi(S)),
# S over collections of vertex-disjoint directed odd cycles of T (independent sets of the OCF conflict graph), phi(S) =
# sum (|c| - 1), E_k = tangent numbers (Euler zigzag, k odd).  It follows from F(T,t) = sum_S 2^|S| (1-t)^phi(S) A~_(n-phi(S))(t)
# (Redei-Berge / deformed Eulerian form, cf. THM-062/063) at t = -1, with A~_k(-1) = (-1)^((k-1)/2) E_k for odd k.
# Also checks the forest form F(T,x) = sum_(L even) N_L (n-L)! ((x-1)/2)^L ((x+1)/2)^(n-1-L).  mac-mini-2026-10-22-sierp.
import itertools, random
from math import factorial
from fractions import Fraction

def zigzag(N):
    # Euler zigzag numbers via boustrophedon: E_0..E_N
    E = [1]; row = [1]
    for k in range(1, N + 1):
        new = [0]
        for x in reversed(row): new.append(new[-1] + x)
        row = new; E.append(row[-1])
    return E

def directed_odd_cycles(A, n):
    cyc = []
    for k in range(3, n + 1, 2):
        for sub in itertools.combinations(range(n), k):
            first = sub[0]
            for perm in itertools.permutations(sub[1:]):
                c = (first,) + perm
                if all(A[c[i]][c[(i + 1) % k]] for i in range(k)): cyc.append(frozenset(sub) if False else c)
    return cyc

def check(A, n, E):
    F = 0
    for w in itertools.permutations(range(n)):
        F += (-1) ** sum(1 for i in range(n - 1) if not A[w[i]][w[i + 1]])
    cyc = directed_odd_cycles(A, n)
    sets = [frozenset(c) for c in cyc]
    tot = 0
    # enumerate collections of pairwise vertex-disjoint cycles
    def rec(start, used, wt, phi):
        nonlocal tot
        tot += wt * (-1) ** ((n - phi - 1) // 2) * E[n - phi]
        for i in range(start, len(cyc)):
            if not (sets[i] & used):
                k = len(cyc[i]); rec(i + 1, used | sets[i], wt * 2 ** k * 1, phi + k - 1)
    rec(0, frozenset(), 1, 0)
    # forest form at x = -1/2 (a generic point): compare F(T,x) by brute force
    x = Fraction(-1, 2)
    Fx = sum(x ** sum(1 for i in range(n - 1) if A[w[i]][w[i + 1]]) for w in itertools.permutations(range(n)))
    S = [[0 if i == j else (1 if A[i][j] else -1) for j in range(n)] for i in range(n)]
    # N_L by brute force over orderings: 2^(n-1)... use the identity N_L (n-L)! = sum_w e_L(s_1..s_(n-1)) (elementary symm.)
    NL = [0] * n
    for w in itertools.permutations(range(n)):
        s = [S[w[i]][w[i + 1]] for i in range(n - 1)]
        e = [1] + [0] * (n - 1)
        for v in s:
            for L in range(n - 1, 0, -1): e[L] += e[L - 1] * v
        for L in range(n): NL[L] += e[L]
    # NL[L] here = sum_w e_L = N_L (n-L)! * ... ; forest form: F(T,x) = sum_L (sum_w e_L) ((x-1)/2)^L ((x+1)/2)^(n-1-L)
    Ff = sum(NL[L] * ((x - 1) / 2) ** L * ((x + 1) / 2) ** (n - 1 - L) for L in range(n))
    odd_zero = all(NL[L] == 0 for L in range(1, n, 2))
    return F, tot, Fx == Ff, odd_zero

E = zigzag(15)
random.seed(1)
ok = True
for n in (3, 5, 7):
    for trial in range(6 if n < 7 else 3):
        A = [[0] * n for _ in range(n)]
        for i in range(n):
            for j in range(i + 1, n):
                if random.random() < 0.5: A[i][j] = 1
                else: A[j][i] = 1
        F, tot, ff, oz = check(A, n, E)
        ok &= (F == tot) and ff and oz
        print(f"n={n}: F(T,-1) = {F}, OCF-tangent sum = {tot}, forest form ok: {ff}, odd N_L vanish: {oz}")
print('ALL OK' if ok else 'MISMATCH')
