# opus-2026-07-15-S312 -- debug: take the 8 real qf triples at n=5, examine
# witness pi: cycle type, pi p0 = p0?, p0 in im(pi+1)?, solve S, verify Y anti.
import itertools
from collections import defaultdict

n = 5
pairs = [(u, v) for u in range(n) for v in range(u+1, n)]
P = len(pairs)
pidx = {p: i for i, p in enumerate(pairs)}
perms = list(itertools.permutations(range(n)))
ONE = (1 << P) - 1

def apply_perm(pi, X):
    Y = 0
    for i, (u, v) in enumerate(pairs):
        pu, pv = pi[u], pi[v]
        a, b = (pu, pv) if pu < pv else (pv, pu)
        bit = (X >> i) & 1
        if pu > pv: bit ^= 1
        Y |= bit << pidx[(a, b)]
    return Y

def apply_perm_pairs(pi, S):
    R = 0
    for i, (u, v) in enumerate(pairs):
        if (S >> i) & 1:
            pu, pv = pi[u], pi[v]
            a, b = (pu, pv) if pu < pv else (pv, pu)
            R |= 1 << pidx[(a, b)]
    return R

def cycle_type(pi):
    seen = [False]*n; ct = []
    for v in range(n):
        if not seen[v]:
            l = 0; w = v
            while not seen[w]:
                seen[w] = True; w = pi[w]; l += 1
            ct.append(l)
    return tuple(sorted(ct))

p0_pairs = 0
for k in range(1, n):
    p0_pairs |= 1 << pidx[(k-1, k)]

def contains_std_path(X):
    return all(not (X >> pidx[(k-1, k)]) & 1 for k in range(1, n))

rho = tuple(n-1-v for v in range(n))
found = 0
for X in range(1 << P):
    if not contains_std_path(X): continue
    target = X ^ ONE ^ p0_pairs
    for pi in perms:
        if apply_perm(pi, X) == target:
            found += 1
            blue = apply_perm(rho, X) == X ^ ONE
            pip0 = apply_perm_pairs(pi, p0_pairs)
            # solve (pi+1)S = p0: brute force
            sols = [S for S in range(1 << P) if S ^ apply_perm_pairs(pi, S) == p0_pairs]
            print(f"X={X:04x} blue={blue} pi={pi} ct={cycle_type(pi)} "
                  f"pi(p0)==p0: {pip0 == p0_pairs}  #S-solutions: {len(sols)}")
            if sols:
                S = sols[0]; Y = X ^ S
                print(f"    Y anti-check: {apply_perm(pi, Y) == Y ^ ONE}")
print(f"total qf triples found: {found}")
