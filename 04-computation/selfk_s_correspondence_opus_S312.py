# opus-2026-07-15-S312 -- HYP-6895: the S-correspondence test.
# DERIVATION: flips = GF(2) vectors on pair-positions. rev = +ONE (all pairs);
# X^p = X + ONE + p. qf witness: pi X = X + ONE + p0.
# If S solves (pi+1)S = p0 then Y := X + S has pi Y = rev Y  (Y self-converse!).
# Conversely (Y, pi, S) with (pi+1)S = p0 gives X = Y + S qf.
# TEST: for each anti-pair (Y, pi), how many S give a valid qf triple
# (X = Y+S must CONTAIN the directed standard path p0)?  If exactly-1 on the
# black side: M_black = n! * SC  -- the proof skeleton.
import itertools, time
from collections import defaultdict

def run(n):
    t0 = time.time()
    pairs = [(u, v) for u in range(n) for v in range(u+1, n) ]
    P = len(pairs)
    pidx = {p: i for i, p in enumerate(pairs)}
    perms = list(itertools.permutations(range(n)))
    # X encoded as bitmask over pairs: bit i = 1 means arc u->v (u<v), else v->u
    def apply_perm(pi, X):
        # (pi X)(u,v) = X(pi^-1 u, pi^-1 v); build permuted mask
        Y = 0
        for i, (u, v) in enumerate(pairs):
            pu, pv = pi[u], pi[v]
            a, b = (pu, pv) if pu < pv else (pv, pu)
            j = pidx[(a, b)]
            bit = (X >> i) & 1
            if pu > pv: bit ^= 1     # orientation flips when order swaps
            Y |= bit << j
        return Y
    ONE = (1 << P) - 1
    # standard path p0: arcs k -> k-1 for k = n-1..1 (0-indexed v: n-1 -> ... -> 0)
    # pair {k, k-1}, orientation: higher -> lower = arc (u=k-1 < v=k): v->u: bit 0
    p0_pairs = 0
    for k in range(1, n):
        p0_pairs |= 1 << pidx[(k-1, k)]
    def contains_std_path(X):
        for k in range(1, n):
            if (X >> pidx[(k-1, k)]) & 1: return False   # bit1 = u->v = k-1->k WRONG dir
        return True
    def is_ham_path_of(X, S_pairs):
        # S_pairs: pair-mask; check it forms a Ham path and orientations in X consistent?
        # for the correspondence we only need p = (pi+1)S to BE p0 (std path pairs)
        return S_pairs == p0_pairs
    # sigma-blue test: rho = reversal relabel v -> n-1-v; t blue iff rho is anti-automorphism
    rho = tuple(n-1-v for v in range(n))

    # 1. enumerate qf triples directly:  pi X = X ^ ONE ^ p0_pairs  (X contains p0)
    M_black = M_blue = 0
    qf_black = set()
    witnesses = defaultdict(list)
    for X in range(1 << P):
        if not contains_std_path(X): continue
        target = X ^ ONE ^ p0_pairs
        for pi in perms:
            if apply_perm(pi, X) == target:
                blue = apply_perm(rho, X) == X ^ ONE   # rho anti-automorphism
                if blue: M_blue += 1
                else:
                    M_black += 1
                    qf_black.add(X)
                    witnesses[X].append(pi)
    # 2. anti-pairs
    anti = []
    for Y in range(1 << P):
        for pi in perms:
            if apply_perm(pi, Y) == Y ^ ONE:
                anti.append((Y, pi))
    fact = 1
    for i in range(2, n+1): fact *= i
    SC = len(anti) // fact
    print(f"n={n}: qf triples: black={M_black}, blue={M_blue}; "
          f"anti-pairs={len(anti)} = {fact}*{SC} (SC={SC})")
    print(f"   M_black == n!*SC ? {M_black} vs {fact*SC}: {M_black == fact*SC}")
    print(f"   M_total == n!*(weighted qf) sanity: M_tot={M_black+M_blue}")

    # 3. the S-correspondence: per anti-pair count S with (pi+1)S = p0 and X=Y^S contains p0
    #    (pi+1)S = S ^ pi(S) as PAIR-SETS (positions only, orientation-free on S)
    def apply_perm_pairs(pi, Smask):
        R = 0
        for i, (u, v) in enumerate(pairs):
            if (Smask >> i) & 1:
                pu, pv = pi[u], pi[v]
                a, b = (pu, pv) if pu < pv else (pv, pu)
                R |= 1 << pidx[(a, b)]
        return R
    hist = defaultdict(int)
    tot_valid = 0
    for (Y, pi) in anti:
        cnt = 0
        for S in range(1 << P):
            if S ^ apply_perm_pairs(pi, S) == p0_pairs:
                X = Y ^ S
                if contains_std_path(X):
                    blue = apply_perm(rho, X) == X ^ ONE
                    if not blue: cnt += 1
        hist[cnt] += 1
        tot_valid += cnt
    print(f"   S-correspondence multiplicity histogram over anti-pairs: {dict(sorted(hist.items()))}")
    print(f"   total valid (Y,pi,S) black triples = {tot_valid} (M_black = {M_black})")
    print(f"   ({time.time()-t0:.0f}s)")

run(4)
run(5)
