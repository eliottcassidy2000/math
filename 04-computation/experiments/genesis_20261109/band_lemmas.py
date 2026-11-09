"""Finite checks of the band lemmas (genesis lane, mac-mini-2026-11-09).

band_s(n): vertices 0..n-1; for i < j the arc is j -> i when j - i <= s (a "back arc"), else i -> j.
Checks:
 (W) two-step walk lemma: the digraph on {0..n-1} with arcs i -> i+a and i -> i-b (when in range), gcd(a,b) = 1,
     is strongly connected for every n >= a + b   (all a, b <= 25, n <= 80);
 (L) lambda(band_s(n)) = s for 2s+1 <= n <= NMAX (arc-strong connectivity by max-flow, all pairs from vertex 0);
 (P) band_s(n) is prime for 2s+3 <= n <= NMAX and not prime for n = 2s+2;
 (D) the s arc-disjoint spanning strong subdigraphs D_r = {back span r} + {forward span r*s+1}, r = 1..s, exist
     (each strongly connected) for n >= s*s + s + 1.
"""
import sys
from math import gcd
from collections import deque


def strongly_connected(n, arcs_out):
    def reach(adj):
        seen = [False] * n
        seen[0] = True
        dq = deque([0])
        while dq:
            v = dq.popleft()
            for w in adj[v]:
                if not seen[w]:
                    seen[w] = True
                    dq.append(w)
        return all(seen)
    rev = [[] for _ in range(n)]
    for v in range(n):
        for w in arcs_out[v]:
            rev[w].append(v)
    return reach(arcs_out) and reach(rev)


def check_W(amax=25, nmax=80):
    bad = 0
    for a in range(1, amax + 1):
        for b in range(1, amax + 1):
            if gcd(a, b) != 1:
                continue
            for n in range(a + b, nmax + 1):
                adj = [[w for w in (v + a, v - b) if 0 <= w < n] for v in range(n)]
                if not strongly_connected(n, adj):
                    bad += 1
    return bad


def band_out(n, s):
    out = [set() for _ in range(n)]
    for i in range(n):
        for j in range(i + 1, n):
            if j - i <= s:
                out[j].add(i)
            else:
                out[i].add(j)
    return out


def maxflow(n, out, src, dst):
    cap = {}
    for v in range(n):
        for w in out[v]:
            cap[(v, w)] = cap.get((v, w), 0) + 1
            cap.setdefault((w, v), 0)
    nbr = [[] for _ in range(n)]
    for (v, w) in cap:
        nbr[v].append(w)
    flow = 0
    while True:
        par = {src: None}
        dq = deque([src])
        while dq and dst not in par:
            v = dq.popleft()
            for w in nbr[v]:
                if w not in par and cap[(v, w)] > 0:
                    par[w] = v
                    dq.append(w)
        if dst not in par:
            return flow
        w = dst
        while par[w] is not None:
            v = par[w]
            cap[(v, w)] -= 1
            cap[(w, v)] += 1
            w = v
        flow += 1


def lam(n, out):
    best = None
    for v in range(1, n):
        for a, b in ((0, v), (v, 0)):
            f = maxflow(n, out, a, b)
            best = f if best is None else min(best, f)
    return best


def is_prime(n, out):
    full = set(range(n))
    for u in range(n):
        for w in range(u + 1, n):
            M = {u, w}
            ch = True
            while ch:
                ch = False
                for x in range(n):
                    if x in M:
                        continue
                    k = len(out[x] & M)
                    if 0 < k < len(M):
                        M.add(x)
                        ch = True
            if M != full:
                return False
    return True


def check_D(s, n):
    out = band_out(n, s)
    for r in range(1, s + 1):
        f = r * s + 1
        adj = [[w for w in (v - r, v + f) if 0 <= w < n] for v in range(n)]
        for v in range(n):
            for w in adj[v]:
                assert w in out[v], (s, n, r, v, w)
        if not strongly_connected(n, adj):
            return False
    return True


if __name__ == '__main__':
    NMAX = int(sys.argv[1]) if len(sys.argv) > 1 else 26
    print('(W) two-step walk failures (a,b<=25, a+b<=n<=80, gcd=1):', check_W())
    for s in range(1, 7):
        lams = {}
        pr = {}
        for n in range(2 * s + 1, NMAX + 1):
            out = band_out(n, s)
            lams[n] = lam(n, out)
            pr[n] = is_prime(n, out)
        badL = [n for n in lams if lams[n] != s]
        badP = [n for n in pr if (n >= 2 * s + 3) != pr[n]]
        dok = all(check_D(s, n) for n in range(s * s + s + 1, max(s * s + s + 1, NMAX) + 1))
        print('s=%d: lambda==s for n in [%d,%d]: %s; prime exactly for n>=2s+3: %s (n=2s+2 prime? %s); '
              'D_r strong for n in [%d,%d]: %s' % (s, 2 * s + 1, NMAX, 'yes' if not badL else 'NO %s' % badL,
                                                     'yes' if not badP else 'NO %s' % badP, pr.get(2 * s + 2),
                                                     s * s + s + 1, max(s * s + s + 1, NMAX), dok))
