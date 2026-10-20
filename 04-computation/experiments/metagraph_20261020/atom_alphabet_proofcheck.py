#!/usr/bin/env python3
r"""Proof-check for THM-4730 (strong pure-blue atoms are K1, C3, T5).
(A) odd n >= 7: for EVERY circulant tournament on Z_n with 1 in S, the proof's explicit Hamiltonian path
    (three-interval path (X,Z,Y) with x = a-1, y = b-1 for a != b in S\{1}, a+b <= n+1; or the descending-pairs path
    for S = {1} u [(n+1)/2, n-2]) is verified to be a Hamiltonian path that is NOT palindromic
    by testing positional reversal on every arc. This is finite evidence only;
    the original all-order inference from nonconstant partner sums has a proof gap.
(B) n = 2 mod 4, n >= 6: for EVERY (S', D) bicirculant (S' a circulant set on Z_k, k = n/2 odd; D a set of odd residues
    with 1 in D, -1 not in D) the proof's cycle C' = (0, 2s, ..., -2s, d-2s, d, ..., d-4s) with d in D, d-4s not in D,
    s = +-1 in S', is verified to be a Hamiltonian cycle containing the even-even arc 0 -> 2s.
(C) direct brute force: no strong pure-blue tournament in the n = 6 and n = 10 bicirculant families (all Ham paths tested)."""
import itertools, sys
def check_path(n, path, arc):
    assert len(path) == n and len(set(path)) == n, path
    for u, v in zip(path, path[1:]): assert arc(u, v), (path, u, v)
def palindromic(n, path, arc):
    # Positional reversal can be nonaffine in a composite-order circulant.
    return all(arc(path[i], path[j]) == arc(path[n - 1 - j], path[n - 1 - i])
               for i in range(n) for j in range(i + 1, n))
print("(A) circulants of odd order with 1 in S")
for n in range(7, 24, 2):
    total = 0; caseA = 0; caseB = 0
    for choice in itertools.product([0, 1], repeat=(n - 3) // 2):
        S = {1} | {(d if ch == 0 else n - d) for d, ch in zip(range(2, (n + 1) // 2), choice)}
        arc = lambda u, v: (v - u) % n in S
        small = sorted(d for d in range(2, (n + 1) // 2) if d in S)
        pair = None
        if len(small) >= 2: pair = (small[0], small[1])
        elif len(small) == 1:
            a = small[0]; d = a - 1 if a >= 3 else 3
            assert 2 <= d <= (n - 1) // 2 and d != a and (n - d) in S
            pair = (a, n - d)
        if pair:
            a, b = pair; assert a != b and a + b <= n + 1 and a in S and b in S and 1 not in (a, b)
            x, y = a - 1, b - 1
            path = list(range(0, x)) + list(range(x + y, n)) + list(range(x, x + y)); caseA += 1
        else:
            assert S == {1} | set(range((n + 1) // 2, n - 1)), S
            path = [0, 1]
            for j in range(1, (n - 3) // 2 + 1): path += [n - 2 * j, n - 2 * j + 1]
            path += [2]; caseB += 1
        check_path(n, path, arc); assert not palindromic(n, path, arc), (n, S, path)
        total += 1
    print(f"  n = {n}: {total} connection sets, three-interval {caseA}, descending-pairs {caseB}: all have a non-palindromic Hamiltonian path", flush=True)
print("(B) bicirculants of order n = 2 mod 4")
for n in range(6, 31, 4):
    k = n // 2; total = 0
    for chS in itertools.product([0, 1], repeat=(k - 1) // 2):
        Sp = {(y if ch == 0 else k - y) for y, ch in zip(range(1, (k + 1) // 2), chS)}
        s = 1 if 1 in Sp else -1; assert (s % k) in Sp
        odd = [r for r in range(1, n, 2) if r not in (1, n - 1)]
        for chD in itertools.product([0, 1], repeat=len(odd)):
            D = {1} | {r for r, ch in zip(odd, chD) if ch}
            def arc(u, v):
                diff = (v - u) % n
                if diff % 2 == 0: return (diff // 2) % k in Sp
                if u % 2 == 0: return diff in D
                return ((u - v) % n) not in D   # odd -> even iff (o - e) not in D
            ds = [d for d in D if (d - 4 * s) % n not in D]; assert ds, (n, Sp, D)
            d = ds[0]
            cyc = [(2 * s * i) % n for i in range(k)] + [(d - 2 * s + 2 * s * i) % n for i in range(k)]
            check_path(n, cyc, arc); assert arc(cyc[-1], cyc[0]); assert cyc[0] % 2 == 0 and cyc[1] % 2 == 0
            total += 1
    print(f"  n = {n}: {total} bicirculants (S', D): the mixing Hamiltonian cycle C' exists in every one", flush=True)
print("(C) brute force: all Hamiltonian paths palindromic? in the n = 6, 10 bicirculant families")
import networkx as nx
from networkx.algorithms.isomorphism import DiGraphMatcher
def all_ham_paths(n, adj):
    full = (1 << n) - 1; out = []
    def rec(v, mask, path):
        if mask == full: out.append(tuple(path)); return
        for w in adj[v]:
            if not (mask >> w) & 1:
                path.append(w); rec(w, mask | (1 << w), path); path.pop()
    for v in range(n): rec(v, 1 << v, [v])
    return out
for n in (6, 10):
    k = n // 2; pb = 0; total = 0
    for chS in itertools.product([0, 1], repeat=(k - 1) // 2):
        Sp = {(y if ch == 0 else k - y) for y, ch in zip(range(1, (k + 1) // 2), chS)}
        odd = [r for r in range(1, n, 2) if r not in (1, n - 1)]
        for chD in itertools.product([0, 1], repeat=len(odd)):
            D = {1} | {r for r, ch in zip(odd, chD) if ch}
            def arc(u, v):
                diff = (v - u) % n
                if diff % 2 == 0: return (diff // 2) % k in Sp
                if u % 2 == 0: return diff in D
                return ((u - v) % n) not in D
            adjm = [[1 if u != v and arc(u, v) else 0 for v in range(n)] for u in range(n)]
            assert all(adjm[u][v] + adjm[v][u] == 1 for u in range(n) for v in range(n) if u != v)
            adj = [[v for v in range(n) if adjm[u][v]] for u in range(n)]
            D_ = nx.DiGraph(); D_.add_nodes_from(range(n)); D_.add_edges_from((u, v) for u in range(n) for v in range(n) if adjm[u][v])
            R = D_.reverse(copy=True)
            antis = [[mp[i] for i in range(n)] for mp in DiGraphMatcher(D_, R).isomorphisms_iter()]
            paths = all_ham_paths(n, adj)
            pal = sum(1 for p in paths if any(all(sig[p[x]] == p[n - 1 - x] for x in range(n)) for sig in antis))
            total += 1
            if pal == len(paths): pb += 1
    print(f"  n = {n}: {total} tournaments, pure-blue among them: {pb}", flush=True)
