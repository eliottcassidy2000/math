"""Independent check of the containment census for targets N < 45 (orders <= 7 suffice,
since a strong tournament with H < 45 = f(8) has at most 7 vertices).  H of every induced
subtournament is computed from the odd-cycle formula H = I(Omega, 2) (THM-002), not by a
path DP: for |S| <= 7 a pack has at most two cycles, so H(S) = 1 + 2 c(S) + 4 d(S)."""
import itertools, subprocess, collections, sys
def isprime(x): return x > 1 and all(x % d for d in range(2, int(x**.5) + 1))
def tours(n):
    out = subprocess.run(['gentourng', '-q', str(n)], capture_output=True, text=True).stdout.split()
    for s in out:
        A = [[False]*n for _ in range(n)]; k = 0
        for i in range(n):
            for j in range(i+1, n):
                if s[k] == '1': A[i][j] = True
                else: A[j][i] = True
                k += 1
        yield A
def odd_cycles(A, n):
    cyc = []
    for L in (3, 5, 7):
        for S in itertools.combinations(range(n), L):
            m = 0
            first = S[0]
            for perm in itertools.permutations(S[1:]):
                seq = (first,) + perm
                if all(A[seq[i]][seq[(i+1) % L]] for i in range(L)): m += 1
            for _ in range(m): cyc.append(frozenset(S))
    return cyc
real = collections.defaultdict(set)
for n in range(3, 8):
    for A in tours(n):
        cyc = odd_cycles(A, n)
        masks = [sum(1 << v for v in c) for c in cyc]
        def Hs(S):
            inside = [m for m in masks if m & ~S == 0]
            d = sum(1 for i in range(len(inside)) for j in range(i+1, len(inside)) if inside[i] & inside[j] == 0)
            return 1 + 2*len(inside) + 4*d
        V = (1 << n) - 1
        N = Hs(V)
        if N < 23 or N >= 45: continue
        pres = {Hs(S) for S in range(1, V)}
        for c in (1, -1, -3):
            s = N + c
            for p in range(3, s//2 + 1, 2):
                q = s - p
                if 7 in (p, q) or not (isprime(p) and isprime(q)): continue
                if p in pres and q in pres: real[N].add((p, q, c))
for N in range(23, 45, 2):
    S = sorted(real.get(N, set()), key=lambda x: (x[2], x[0]))
    print(N, ' '.join(f"{p}+{q}{'-1' if c==1 else ('+1' if c==-1 else '+3')}" for p, q, c in S) or 'NONE')
