"""
LRC coloring / flow-duality core.

Setup: distinct positive integer speeds v_1..v_m, observer 0, n=m+1, threshold 1/n.
At time t, points P(t) = {0} U {v_i t mod 1}. Circular distance d(a,b) = min(|a-b|, 1-|a-b|).
Danger graph G(t): n vertices (observer + m runners), edge iff d < 1/n.

Exact arithmetic via Fraction. Time circle [0,1) decomposed into cells: the dynamics are
piecewise constant in combinatorial type between "events". Events occur when:
 - two points collide: positions equal  -> distance 0
 - a pairwise circular distance crosses exactly 1/n (edge appears/disappears)
We collect all critical t in [0,1), sort, and sample the OPEN midpoint of each cell to get
the generic combinatorial type on that cell. (Endpoints are measure-zero; cell = open interval.)
"""

from fractions import Fraction as F
from math import gcd, lcm
from itertools import combinations


def circ_dist(a, b):
    """Circular distance on [0,1), a,b are Fractions in [0,1)."""
    d = abs(a - b)
    if d > F(1, 2):
        d = 1 - d
    return d


def positions(speeds, t):
    """Positions of observer(0) + runners at time t. Returns list of Fractions in [0,1)."""
    pts = [F(0)]
    for v in speeds:
        x = (v * t) % 1
        pts.append(x)
    return pts


def critical_times(speeds, n):
    """
    All t in [0,1) where some pairwise circular distance equals 0 or 1/n.
    Points indexed 0..m, point i moves at speed s_i with s_0 = 0.
    Pairwise relative position: (s_i - s_j) t mod 1. We need
      (s_i - s_j) t = k   (collision, dist 0), or
      (s_i - s_j) t = k +/- 1/n  (dist exactly 1/n)
    for integer k. Solve t in [0,1).
    """
    s = [0] + list(speeds)
    thr = F(1, n)
    crit = set()
    crit.add(F(0))
    for i, j in combinations(range(len(s)), 2):
        w = s[i] - s[j]          # integer relative speed
        if w == 0:
            continue
        aw = abs(w)
        # relative position r(t) = w*t mod 1; events when w*t = integer + offset, offset in {0, thr, -thr -> 1-thr}
        for offset in (F(0), thr, 1 - thr):
            # w*t = k + offset  => t = (k+offset)/w ; want t in [0,1)
            for k in range(-aw - 1, aw + 1):
                t = F(k) + offset
                t = t / w
                if 0 <= t < 1:
                    crit.add(t % 1)
    return sorted(crit)


def cell_midpoints(crit):
    """Open-cell representative times: midpoints between consecutive critical times (wrap)."""
    mids = []
    cl = crit + [F(1)]
    for a, b in zip(cl, cl[1:]):
        mids.append((a + b) / 2)
    return mids


def danger_graph(speeds, n, t):
    """Return adjacency (set of frozenset edges) and vertex count for danger graph at t."""
    pts = positions(speeds, t)
    N = len(pts)
    thr = F(1, n)
    edges = set()
    for i, j in combinations(range(N), 2):
        if circ_dist(pts[i], pts[j]) < thr:
            edges.add(frozenset((i, j)))
    return N, edges, pts


def clique_number_arc(speeds, n, t):
    """
    Clique number of the danger graph. The danger graph is a unit-circular-arc graph:
    vertex i = arc of half-width 1/(2n)? No -- edge iff point-distance < 1/n. Equivalent:
    place each point on circle; two points adjacent iff within 1/n. A clique = set of points
    pairwise within 1/n. Max clique = max number of points contained in some arc of length < 1/n
    ... but pairwise-within is NOT same as all-in-one-small-arc in general on a circle.
    We compute clique number directly (small n) via max set pairwise-adjacent.
    """
    N, edges, pts = danger_graph(speeds, n, t)
    adj = {i: set() for i in range(N)}
    for e in edges:
        a, b = tuple(e)
        adj[a].add(b); adj[b].add(a)
    best = 0
    # Bron-Kerbosch-ish for tiny graphs
    def expand(R, P):
        nonlocal best
        if not P:
            best = max(best, len(R))
            return
        if len(R) + len(P) <= best:
            return
        for v in list(P):
            expand(R | {v}, P & adj[v])
            P = P - {v}
    expand(set(), set(range(N)))
    return best


def chromatic_number(speeds, n, t):
    """
    Danger graph is a circular-arc graph (each point -> arc [x-1/(2n), x+1/(2n)]?).
    Actually edge iff distance<1/n. Model each point i as arc A_i = (x_i - 1/(2n), x_i + 1/(2n))
    on the circle; then A_i ∩ A_j != empty  <=>  |x_i - x_j|_circ < 1/n. Yes! exact.
    So G(t) is the intersection graph of n equal arcs of length 1/n on the circle => a
    *proper circular-arc graph* (unit circular-arc). Coloring such graphs: chi can exceed omega
    in general circular-arc graphs! (circular-arc graphs are NOT perfect.) So we compute chi exactly
    by greedy-with-backtracking / brute force for small n.
    """
    N, edges, pts = danger_graph(speeds, n, t)
    adj = {i: set() for i in range(N)}
    for e in edges:
        a, b = tuple(e)
        adj[a].add(b); adj[b].add(a)
    # exact chromatic number via increasing k, backtracking assignment
    order = sorted(range(N), key=lambda v: -len(adj[v]))
    for k in range(1, N + 1):
        color = {}
        def assign(idx):
            if idx == len(order):
                return True
            v = order[idx]
            used = {color[u] for u in adj[v] if u in color}
            for c in range(k):
                if c not in used:
                    color[v] = c
                    if assign(idx + 1):
                        return True
                    del color[v]
            return False
        if assign(0):
            return k, color
    return N, {}


def observer_danger_count(speeds, n, t):
    """N(t) = number of runners within 1/n of observer (the time-coloring level)."""
    thr = F(1, n)
    cnt = 0
    for v in speeds:
        x = (v * t) % 1
        if circ_dist(F(0), x) < thr:
            cnt += 1
    return cnt
