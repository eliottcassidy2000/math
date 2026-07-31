#!/usr/bin/env python3
"""
Independent second path for the sensitivity-vs-degree exhaustive claims.

Different algorithm & representation from sens_search.c:
  - work over the *Moebius coefficient* space directly, level by level
    (c_S for |S|=2 in {-2,-1} forced by f(e_i+e_j)=2+c_ij in {0,1}, etc.)
  - meet constraints via dictionary of partial level values, DFS by level,
    no mask-numeric ordering, pure Python big-int arithmetic.

Claim checked:  EXISTS f:{0,1}^n->{0,1} multilinear with deg<=d, f(0)=0,
                f(e_i)=1 for all i?    (n,d) in {(4,2),(7,3), ...}

f(x) = sum_{S subseteq supp(x), |S|<=d} c_S   with c_emptyset=0, c_{i}=1.
Booleanness on level k means: for every |X|=k: sum_{S subset X,|S|<=d} c_S in {0,1}.
DFS over c_S level by level: level-2 coefficients, then level-3, ..., level-d;
after level d all higher-level values are determined linear forms: check all.
Within level k, we choose c_S for each S; the value f(X)=k + partial sums must
land in {0,1} for |X|=k, which pins c_S = f(S) - (sum over proper subsets),
f(S) in {0,1}: so the branching is exactly over Boolean values f(S) at levels
2..d — mirrors the C search but organized by level, implemented independently.
"""
from itertools import combinations
import sys

def decide(n, d, find_all=False):
    subsets = {}  # frozenset -> f value
    subsets[frozenset()] = 0
    for i in range(n):
        subsets[frozenset([i])] = 1
    levels = [[frozenset(c) for c in combinations(range(n), k)] for k in range(n + 1)]
    sols = []

    # precompute proper-subset lists per set (all sizes)
    def fval(X, assigned):
        # determined value at X (|X|>d): f(X) = -sum_{T proper subset X} (-1)^{|X|-|T|} f(T)
        s = 0
        lx = len(X)
        Xl = sorted(X)
        for r in range(lx):
            for c in combinations(Xl, r):
                T = frozenset(c)
                s += (-1) ** (lx - r) * assigned[T]
        return -s

    order = []
    for k in range(2, d + 1):
        order.extend(levels[k])

    def check_forced(assigned):
        for k in range(d + 1, n + 1):
            for X in levels[k]:
                v = fval(X, assigned)
                if v not in (0, 1):
                    return False
                assigned[X] = v
        return True

    def undo_forced(assigned):
        for k in range(d + 1, n + 1):
            for X in levels[k]:
                assigned.pop(X, None)

    count = [0]
    def dfs(idx, assigned):
        if idx == len(order):
            ok = check_forced(assigned)
            if ok:
                sols.append(dict(assigned))
                count[0] += 1
            undo_forced(assigned)
            return ok and not find_all
        X = order[idx]
        # partial pruning: none (keep the two paths maximally different);
        # rely on final check only for small n, but that is 2^56 for (3,7)...
        # so DO add the *level-local* pruning that is mathematically forced:
        # when |X| = k <= d is assigned we require nothing (free), but any
        # completed superset at level k+1..d is also free... so the only
        # pruning available before the end is the forced-level check on
        # complete lower sets: implement partial forced check on sets whose
        # subsets are all assigned:
        for v in (0, 1):
            assigned[X] = v
            if partial_ok(X, assigned):
                if dfs(idx + 1, assigned):
                    del assigned[X]
                    return True
            del assigned[X]
        return False

    # partial forced check: after assigning X (level k), any forced set Y
    # (level > d) ALL of whose free subsets are assigned can be evaluated.
    # For independence keep it simple: check all forced Y contained in the
    # union of assigned sets support and with all free subsets assigned.
    idx_of = {X: i for i, X in enumerate(order)}
    def partial_ok(X, assigned):
        i = idx_of[X]
        # a forced Y is fully determined iff every free subset of Y has index <= i
        # equivalently the max index among free subsets of Y equals i and X ⊆ Y.
        # Enumerate candidate Y = X ∪ extra with |Y| in d+1..n, extra from earlier vars…
        # For simplicity (still exact): check all Y ⊇ X with |Y| <= n such that all
        # free subsets of Y are assigned.
        base = sorted(X)
        rest = [v for v in range(n) if v not in X]
        from itertools import combinations as comb
        for add in range(max(0, d + 1 - len(X)), n - len(X) + 1):
            for extra in comb(rest, add):
                Y = frozenset(base + list(extra))
                if len(Y) <= d:
                    continue
                ready = True
                Yl = sorted(Y)
                for r in range(2, d + 1):
                    for c in comb(Yl, r):
                        if frozenset(c) not in assigned:
                            ready = False
                            break
                    if not ready:
                        break
                if not ready:
                    continue
                v = fval_partial(Y, assigned)
                if v not in (0, 1):
                    return False
        return True

    def fval_partial(Y, assigned):
        # evaluate forced value at Y given all free subsets assigned;
        # forced intermediate levels between d+1 and |Y|-1 are themselves
        # determined; compute recursively with memo
        memo = {}
        def g(Z):
            if Z in assigned:
                return assigned[Z]
            if Z in memo:
                return memo[Z]
            s = 0
            Zl = sorted(Z)
            lz = len(Z)
            for r in range(lz):
                for c in combinations(Zl, r):
                    s += (-1) ** (lz - r) * g(frozenset(c))
            memo[Z] = -s
            return -s
        return g(Y)

    found = dfs(0, subsets)
    return (len(sols) > 0, count[0])

if __name__ == "__main__":
    n, d = int(sys.argv[1]), int(sys.argv[2])
    find_all = len(sys.argv) > 3 and sys.argv[3] == "all"
    sat, cnt = decide(n, d, find_all=find_all)
    print(f"INDEPENDENT n={n} d={d}: {'SAT' if sat else 'UNSAT'}"
          + (f" (solutions={cnt})" if find_all else ""))
