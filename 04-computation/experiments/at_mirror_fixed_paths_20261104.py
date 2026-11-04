"""AT_n (transitive order 1 < ... < n with the extreme arc reversed, n -> 1) and its converse involution.

sigma(i) = n + 1 - i is an isomorphism AT_n -> converse(AT_n). It acts on Hamiltonian paths by
iota(P) = reverse(sigma(P)). This script enumerates every Hamiltonian path by DFS (no formula used) and checks:
  (1) H(AT_n) = 2^(n-2) + 1;
  (2) iota is a well-defined involution on the paths;
  (3) Fix(iota) = 1 + 2^((n-2)/2) for even n and 1 for odd n, so Fix(AT_(2k+2)) = H(AT_(k+2));
  (4) the 2^(n-2) back-arc paths A n 1 B are distributed by |A| as the binomial row C(n-2, |A|), and iota maps
      |A| = j to |A| = n - 2 - j (the mirror k <-> m - k of the row about its peak).
Also prints the Fermat rows n = 2^i + 2.
kind-pasteur-2026-11-04-S10.
"""
import sys
from math import comb

def at(n):
    """adjacency sets, vertices 1..n"""
    out = {v: set() for v in range(1, n + 1)}
    for i in range(1, n + 1):
        for j in range(i + 1, n + 1):
            if (i, j) == (1, n):
                out[n].add(1)
            else:
                out[i].add(j)
    return out

def ham_paths(n, out):
    paths = []
    def dfs(path, seen):
        if len(path) == n:
            paths.append(tuple(path)); return
        for w in sorted(out[path[-1]]):
            if w not in seen:
                seen.add(w); path.append(w); dfs(path, seen); path.pop(); seen.discard(w)
    for s in range(1, n + 1):
        dfs([s], {s})
    return paths

def main(nmax):
    for n in range(3, nmax + 1):
        out = at(n)
        # sigma is an anti-automorphism: x -> y iff sigma(y) -> sigma(x)
        sig = lambda v: n + 1 - v
        assert all((y in out[x]) == (sig(x) in out[sig(y)]) for x in out for y in out if x != y)
        P = ham_paths(n, out)
        S = set(P)
        iota = {p: tuple(reversed([sig(v) for v in p])) for p in P}
        assert all(q in S for q in iota.values()) and all(iota[iota[p]] == p for p in P)
        fix = sum(1 for p in P if iota[p] == p)
        m = n - 2
        exp_fix = 1 + 2 ** (m // 2) if m % 2 == 0 else 1
        # |A| distribution for paths through the back arc n -> 1
        dist = [0] * (m + 1); mirror_ok = True
        for p in P:
            if any(p[t] == n and p[t + 1] == 1 for t in range(n - 1)):
                a = p.index(n)
                dist[a] += 1
                q = iota[p]; b = q.index(n)
                mirror_ok &= (b == m - a)
        ok = (len(P) == 2 ** m + 1 and fix == exp_fix and dist == [comb(m, j) for j in range(m + 1)] and mirror_ok)
        fermat = " <- Fermat row: n = 2^i + 2" if m > 0 and (m & (m - 1)) == 0 else ""
        print(f"n={n:2d} H={len(P):5d} (2^{m}+1) Fix(iota)={fix:4d} expected={exp_fix:4d} "
              f"H(AT_(n/2+1))={'%d' % (2 ** (m // 2) + 1) if m % 2 == 0 else '-':>4s} "
              f"binomial row + mirror k<->m-k: {dist == [comb(m, j) for j in range(m + 1)] and mirror_ok}  ok={ok}{fermat}")
        assert ok

if __name__ == "__main__":
    main(int(sys.argv[1]) if len(sys.argv) > 1 else 14)
