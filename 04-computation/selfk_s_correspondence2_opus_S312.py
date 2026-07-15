# opus-2026-07-15-S312 -- HYP-6895: the S-correspondence, correctly scoped.
# GLOBAL IDENTITY (data: n=5: 960=960; n=6,7 via weighted-black=SC):
#   #{(X,p,pi): p directed Ham path of X, pi X = X^p, (X,p) black}
#     = #{(Y,pi): pi Y = rev Y}  = n! * SC(n)     [n >= 5]
# Correspondence: (Y,pi,S) -> (X = Y+S, p = (pi+1)S, pi); fibers = ker(pi+1) torsor.
# TEST at n=5: per anti-pair (Y,pi), histogram of
#   #{S : (pi+1)S is a Ham-path pair-set, X=Y+S orients it end-to-end, black}
# plus: per qf triple, is p in im(pi+1)?  orbit counts, cycle types.
import itertools, time
from collections import defaultdict

n = 5
t0 = time.time()
pairs = [(u, v) for u in range(n) for v in range(u+1, n)]
P = len(pairs)
pidx = {p: i for i, p in enumerate(pairs)}
perms = list(itertools.permutations(range(n)))
ONE = (1 << P) - 1
rho = tuple(n-1-v for v in range(n))

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

def path_of_pairset(S):
    # return vertex sequence if S is a Ham path pair-set, else None
    if bin(S).count('1') != n - 1: return None
    deg = [0]*n; adjp = defaultdict(list)
    for i, (u, v) in enumerate(pairs):
        if (S >> i) & 1:
            deg[u] += 1; deg[v] += 1
            adjp[u].append(v); adjp[v].append(u)
    ends = [v for v in range(n) if deg[v] == 1]
    if len(ends) != 2 or any(d > 2 for d in deg): return None
    seq = [ends[0]]; prev = -1
    while len(seq) < n:
        nxts = [w for w in adjp[seq[-1]] if w != prev]
        if len(nxts) != 1: return None
        prev = seq[-1]; seq.append(nxts[0])
    return seq

def directed_along(X, seq):
    # X orients seq[0] <- seq[1] <- ... ?  accept either end-to-end direction
    fwd = bwd = True
    for a, b in zip(seq, seq[1:]):
        u, v = (a, b) if a < b else (b, a)
        bit = (X >> pidx[(u, v)]) & 1      # 1 = u -> v
        arc_ab = (bit == 1) if a == u else (bit == 0)   # True if a -> b
        if not arc_ab: bwd = False
        else: pass
        if arc_ab: pass
        else: fwd = False
    return fwd or bwd

def tiling_of(X, seq):
    # relabel so seq (directed hi->lo) becomes n-1 -> ... -> 0; return relabeled X
    # find direction: make seq run source->sink along arcs
    a, b = seq[0], seq[1]
    u, v = (a, b) if a < b else (b, a)
    bit = (X >> pidx[(u, v)]) & 1
    if not ((bit == 1) if a == u else (bit == 0)):
        seq = seq[::-1]
    lam = [0]*n
    for j, vtx in enumerate(seq): lam[vtx] = n - 1 - j
    return apply_perm(tuple(lam), X)

def is_blue(X, seq):
    t = tiling_of(X, seq)
    return apply_perm(rho, t) == t ^ ONE

# anti-pairs
anti = []
for Y in range(1 << P):
    for pi in perms:
        if apply_perm(pi, Y) == Y ^ ONE:
            anti.append((Y, pi))
print(f"anti-pairs: {len(anti)}")

def orbits_on_pairs(pi):
    seen = [False]*P; cnt = 0
    for i in range(P):
        if not seen[i]:
            cnt += 1; j = i
            while not seen[j]:
                seen[j] = True
                (u, v) = pairs[j]
                pu, pv = pi[u], pi[v]
                a, b = (pu, pv) if pu < pv else (pv, pu)
                j = pidx[(a, b)]
    return cnt

def cycle_type(pi):
    seen = [False]*n; ct = []
    for v in range(n):
        if not seen[v]:
            l = 0; w = v
            while not seen[w]:
                seen[w] = True; w = pi[w]; l += 1
            ct.append(l)
    return tuple(sorted(ct))

hist = defaultdict(int); tot = 0
by_ct = defaultdict(lambda: defaultdict(int))
for (Y, pi) in anti:
    cnt = 0
    for S in range(1 << P):
        p = S ^ apply_perm_pairs(pi, S)
        seq = path_of_pairset(p)
        if seq is None: continue
        X = Y ^ S
        if not directed_along(X, seq): continue
        if is_blue(X, seq): continue
        cnt += 1
    hist[cnt] += 1; tot += cnt
    by_ct[cycle_type(pi)][cnt] += 1
print(f"S-correspondence per-anti-pair histogram: {dict(sorted(hist.items()))}")
for ct, h in sorted(by_ct.items()):
    print(f"   cycle type {ct}: {dict(sorted(h.items()))}, orbits example: "
          f"{orbits_on_pairs(next(pi for (Y,pi) in anti if cycle_type(pi)==ct))}")
print(f"total black (Y,pi,S) = {tot};  n!*SC = {len(anti)}")
print(f"({time.time()-t0:.0f}s)")
