#!/usr/bin/env python3
"""spine_numeration.py -- the spine family (A |> v1 ... vk) + x of session124 (every vertex of A beats every spine
vertex, earlier spine vertices beat later ones, x -> v_i iff w_i = 1), re-derived and explored.

State of a parent (A, x): h = H(A), F = H(A + x), t = #Hamiltonian paths of A ending at a vertex that beats x,
q = #ordered pairs (P, Q) of vertex-disjoint directed paths with V(P) u V(Q) = V(A), Q nonempty, x -> first(Q), P may be
empty.  Transfer (session124): bit 1: (F, q, t) -> (F, 2q + h, 0); bit 0: (F, q, t) -> (F + q + h - t, 2q, h).
Part 1 checks the transfer by brute force; part 2 computes the achievable H-sets and their coverage of the odd numbers.
"""
import sys, itertools
from functools import lru_cache

def ham_paths(n, adj):
    """number of Hamiltonian paths of a tournament on 0..n-1 (adj[u] = bitmask of out-neighbours), subset DP."""
    full = (1 << n) - 1
    dp = [[0] * n for _ in range(1 << n)]
    for v in range(n): dp[1 << v][v] = 1
    for S in range(1, 1 << n):
        for v in range(n):
            c = dp[S][v]
            if not c: continue
            nxt = adj[v] & ~S
            while nxt:
                u = (nxt & -nxt).bit_length() - 1; nxt &= nxt - 1
                dp[S | (1 << u)][u] += c
    return sum(dp[full])

def build(nA, adjA, xA, w):
    """A on 0..nA-1 (adjA out-masks), x's out-mask into A = xA, spine of length k = len(w): vertices nA..nA+k-1, x = nA+k."""
    k = len(w); n = nA + k + 1; x = nA + k
    adj = [0] * n
    for u in range(nA):
        adj[u] = adjA[u] | sum(1 << (nA + i) for i in range(k))        # A beats the spine
        if not (xA >> u) & 1: adj[u] |= 1 << x                          # u -> x unless x -> u
    for i in range(k):
        v = nA + i
        adj[v] = sum(1 << (nA + j) for j in range(i + 1, k))             # earlier beats later
        if w[i] == 1: adj[x] |= 1 << v
        else: adj[v] |= 1 << x
    adj[x] |= xA
    return n, adj

def parent_state(nA, adjA, xA):
    """(h, F, q, t) by brute force."""
    h = ham_paths(nA, adjA)
    n1, adj1 = build(nA, adjA, xA, [])
    F = ham_paths(n1, adj1)
    # t: Hamiltonian paths of A ending at a vertex u with u -> x (i.e. x does not beat u)
    full = (1 << nA) - 1
    dp = [[0] * nA for _ in range(1 << nA)]
    for v in range(nA): dp[1 << v][v] = 1
    for S in range(1, 1 << nA):
        for v in range(nA):
            c = dp[S][v]
            if not c: continue
            nxt = adjA[v] & ~S
            while nxt:
                u = (nxt & -nxt).bit_length() - 1; nxt &= nxt - 1
                dp[S | (1 << u)][u] += c
    t = sum(dp[full][v] for v in range(nA) if not (xA >> v) & 1)
    # q: ordered pairs (P, Q), P possibly empty, Q nonempty, x -> first(Q); paths = Hamiltonian paths of the induced subtournaments
    paths_of = {}
    for S in range(1, 1 << nA):
        sub = [v for v in range(nA) if (S >> v) & 1]
        # count Hamiltonian paths of A[S] by start vertex
        idx = {v: i for i, v in enumerate(sub)}; m = len(sub)
        sadj = [sum(1 << idx[u] for u in sub if (adjA[v] >> u) & 1) for v in sub]
        d = [[0] * m for _ in range(1 << m)]
        for i in range(m): d[1 << i][i] = 1
        for T in range(1, 1 << m):
            for i in range(m):
                c = d[T][i]
                if not c: continue
                nx = sadj[i] & ~T
                while nx:
                    j = (nx & -nx).bit_length() - 1; nx &= nx - 1
                    d[T | (1 << j)][j] += c
        tot = sum(d[(1 << m) - 1])
        # paths starting at a given vertex: count via reversed DP; simpler: paths of reversed tournament ending there
        radj = [sum(1 << idx[u] for u in sub if (adjA[u] >> v) & 1) for v in sub]
        e = [[0] * m for _ in range(1 << m)]
        for i in range(m): e[1 << i][i] = 1
        for T in range(1, 1 << m):
            for i in range(m):
                c = e[T][i]
                if not c: continue
                nx = radj[i] & ~T
                while nx:
                    j = (nx & -nx).bit_length() - 1; nx &= nx - 1
                    e[T | (1 << j)][j] += c
        starts = {sub[i]: e[(1 << m) - 1][i] for i in range(m)}
        paths_of[S] = (tot, starts)
    q = 0
    for S in range(1, 1 << nA):            # S = V(Q)
        comp = full & ~S
        hP = 1 if comp == 0 else paths_of[comp][0]
        tot, starts = paths_of[S]
        q += hP * sum(c for v, c in starts.items() if (xA >> v) & 1)
    return h, F, q, t

def transfer(state, w, h):
    F, q, t = state
    for b in w:
        if b == 1: F, q, t = F, 2 * q + h, 0
        else: F, q, t = F + q + h - t, 2 * q, h
    return F

def reachable(h, F, q, t, kmax, limit):
    """all H values reachable with words of length <= kmax (values <= limit), by DFS over states (F, q, t)."""
    out = set()
    stack = [(F, q, t, 0)]
    while stack:
        F, q, t, k = stack.pop()
        if F <= limit: out.add(F)
        if k == kmax or F > limit: continue
        stack.append((F, 2 * q + h, 0, k + 1))
        stack.append((F + q + h - t, 2 * q, h, k + 1))
    return out

if __name__ == "__main__":
    # ---- part 1: brute-force check of the transfer on all tournaments A of order <= 4, all cuts x, all words of length <= 5
    bad = 0; checked = 0
    for nA in range(1, 5):
        pairs = [(i, j) for i in range(nA) for j in range(i + 1, nA)]
        for bits in itertools.product([0, 1], repeat=len(pairs)):
            adjA = [0] * nA
            for (i, j), b in zip(pairs, bits):
                if b: adjA[i] |= 1 << j
                else: adjA[j] |= 1 << i
            for xA in range(1 << nA):
                h, F, q, t = parent_state(nA, adjA, xA)
                for k in range(0, 6):
                    for w in itertools.product([0, 1], repeat=k):
                        n, adj = build(nA, adjA, xA, list(w))
                        H = ham_paths(n, adj); Ht = transfer((F, q, t), w, h); checked += 1
                        if H != Ht: bad += 1
    print("part 1: transfer checked on %d (A, x, w) triples, orders of A <= 4, |w| <= 5: mismatches = %d" % (checked, bad), flush=True)
    # ---- part 2: coverage by the transitive parent and by all parents of order <= 3
    LIMIT = int(sys.argv[1]) if len(sys.argv) > 1 else 200000
    def oddmiss(S, lim):
        return [v for v in range(1, lim + 1, 2) if v not in S]
    # transitive parent: A = single vertex a, x -> a: (h, F, q, t) = (1, 1, 1, 0)
    R = reachable(1, 1, 1, 0, 40, LIMIT)
    miss = oddmiss(R, LIMIT)
    print("transitive parent (1,1,1,0): odd values reachable <= %d: %d of %d; missed (first 40): %s" % (LIMIT, len([v for v in R if v % 2]), (LIMIT + 1) // 2, miss[:40]), flush=True)
    print("   last missed below the limit:", miss[-5:] if miss else None, " number missed:", len(miss))
    # union over all parents A of order <= 3 and all cuts
    U = set(); states = {}
    for nA in range(1, 4):
        pairs = [(i, j) for i in range(nA) for j in range(i + 1, nA)]
        for bits in itertools.product([0, 1], repeat=len(pairs)):
            adjA = [0] * nA
            for (i, j), b in zip(pairs, bits):
                if b: adjA[i] |= 1 << j
                else: adjA[j] |= 1 << i
            for xA in range(1 << nA):
                st = parent_state(nA, adjA, xA)
                states[st] = (nA, bits, xA)
    print("distinct parent states (h,F,q,t) from orders <= 3:", sorted(states))
    for st in states:
        h, F, q, t = st
        U |= reachable(h, F, q, t, 40, LIMIT)
    miss = oddmiss(U, LIMIT)
    print("union over parents of order <= 3: odd missed <= %d: %d; first 40: %s" % (LIMIT, len(miss), miss[:40]), flush=True)
    print("   last missed:", miss[-8:])
