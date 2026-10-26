"""Implement the reverse valuation potential on the production ordinary A000568 state graph (port_residual_integral counting/README,
formulas I5-I7 with e=0, one channel). Build the structurally reachable terminal-reaching graph, compute suffix minima V,
run the transformed recurrence Z modulo 2^K from Z_initial=1, and compare a(N)/2^B mod 2^K with the exact count."""
import sys, json, time
from math import gcd, factorial, comb
from fractions import Fraction
from functools import lru_cache
from collections import Counter
from itertools import combinations
from math import prod

def vfact(k): return k - k.bit_count()
@lru_cache(None)
def oddfact(k): return factorial(k) >> vfact(k)
def bound(defect):
    return max(((value + 2*i - 1)//(2*i) for i, value in enumerate(defect) if i), default=0)
def intrinsic(ell, k): return k*(ell-1)//2 + ell*k*(k-1)//2
@lru_cache(None)
def transfer(u, ell, k):
    up = u + ell*k
    q, r = divmod(oddfact(up), oddfact(u) * ell**k * oddfact(k))
    assert r == 0
    return q   # positive odd (beta(k,0)=1)

def build_graph(n):
    """Vertices (ell, s, defect); edges (src, dst, T, J). Returns stages list, edges, initial, terminal."""
    top = n if n % 2 else n - 1
    initial = (top, n, (0,)*((top+1)//2))
    edges = []
    frontier = {initial}
    vertices = {initial}
    for ell in range(top, 0, -2):
        nxt = set()
        for (l, s, defect) in frontier:
            selected = defect[(ell-1)//2] if ell <= s else 0
            b = bound(defect)
            for k in range(s//ell + 1):
                rest = s - ell*k
                if ell == 1 and rest: continue
                dn = tuple(defect[(d-1)//2] + k*(d - gcd(ell, d)) for d in range(1, min(ell-2, rest)+1, 2))
                bp = bound(dn)
                J = intrinsic(ell, k) - vfact(k) + k*(ell*b - selected) + rest*(k + b - bp)
                assert J >= 0
                T = transfer(n - s, ell, k)
                dst = (ell-2, rest, dn)
                edges.append(((l, s, defect), dst, T, J))
                nxt.add(dst); vertices.add(dst)
        frontier = nxt
    terminal = (-1, 0, ())
    assert frontier == {terminal}, frontier
    return vertices, edges, initial, terminal

def ordinary_direct(n):
    def partitions(m, low=1):
        if m == 0: yield (); return
        for d in range(low, m+1, 2):
            for tail in partitions(m-d, d): yield (d,)+tail
    total = Fraction(0)
    for part in partitions(n):
        mult = Counter(part)
        expo = sum((d-1)//2 for d in part) + sum(gcd(a,b) for a,b in combinations(part,2))
        den = prod(d**k * factorial(k) for d,k in mult.items())
        total += Fraction(1 << expo, den)
    assert total.denominator == 1
    return int(total)

def run(n, K):
    t0 = time.time()
    vertices, edges, initial, terminal = build_graph(n)
    # reverse pass: V(terminal)=0; V(v)=min(J+V(w)); vertices are stage-decreasing so process by increasing stage order from terminal
    out = {}
    for (src, dst, T, J) in edges: out.setdefault(src, []).append((dst, T, J))
    V = {terminal: 0}
    order = sorted(vertices, key=lambda v: v[0])   # stage -1 first, then 1, 3, ...
    for v in order:
        if v == terminal: continue
        V[v] = min(J + V[w] for (w, T, J) in out[v])
    B = V[initial]
    reduced = [J + V[w] - V[v] for v in vertices if v != terminal for (w, T, J) in out[v]]
    assert all(j >= 0 for j in reduced)
    # transformed recurrence modulo 2^K, forward by decreasing stage
    mod = 1 << K
    Z = {initial: 1 % mod}
    for v in sorted(vertices, key=lambda v: -v[0]):
        if v == terminal: continue
        zv = Z.get(v, 0)
        if zv == 0: continue
        for (w, T, J) in out[v]:
            Jp = J + V[w] - V[v]
            Z[w] = (Z.get(w, 0) + T * pow(2, Jp, mod) * zv) % mod
    Zt = Z.get(terminal, 0)
    # exact comparison: Y_terminal = O(N) a(N); Z_terminal = Y_terminal / 2^B ; so a(N)/2^B = Z_t * O(N)^{-1} mod 2^K
    a = ordinary_direct(n)
    ON = oddfact(n)
    assert (ON * a) % (1 << B) == 0, 'baseline divides the terminal value'
    normalized_exact = (ON * a) >> B
    ok = (normalized_exact % mod) == Zt
    a_norm_mod = (Zt * pow(ON, -1, mod)) % mod
    ok2 = ((a >> 0) * pow(2, 0, mod)) % mod == ((a_norm_mod * pow(2, B, 2*mod*(1<<B))) ) % mod if False else True
    v2a = (a & -a).bit_length() - 1
    return {'n': n, 'vertices': len(vertices), 'edges': len(edges), 'B': B, 'floor_n_over_2': n//2, 'B_equals_floor': B == n//2 if n >= 3 else B == 0,
            'v2_a': v2a, 'Z_terminal_mod_2^K': Zt, 'matches_exact_normalized_value': ok,
            'a_over_2^B_mod_2^K': a_norm_mod, 'exact_a_over_2^B_mod_2^K': ((a >> B) % mod) if a % (1<<B) == 0 else None,
            'reduced_cost_histogram': dict(sorted(Counter(reduced).items())[:6])}

if __name__ == '__main__':
    K = 20
    for n in [3,4,5,6,7,8,9,10,12,15,20,25,30,35,40,45,50]:
        r = run(n, K)
        print(json.dumps(r))
