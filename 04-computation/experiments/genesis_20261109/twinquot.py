"""For each order-12 escaper: does it arise from a first creator Q of order 11 by doubling one vertex into a 2-vertex
module (substitution Q[x := TT_2]) and reversing ONE arc?  Equivalently: a pair {u,w} with exactly one distinguisher y;
reversing the arc between y and u (or w) makes {u,w} a module; the quotient Q = T' - w is tested for H(Q) not in S(10).
Also reports the 3-cut census class from cuts.py logic for cross-tabulation.  (genesis lane; exact per class)"""
import sys
from collections import Counter
from verify_escapers import adj_from, ham
from spectra import in_S

def flip(out, a, b):
    o = out[:]
    o[a] ^= 1 << b; o[b] ^= 1 << a
    return o

def quotient_fc(out, n):
    full = (1 << n) - 1
    found = False; qvals = set()
    for u in range(n):
        for w in range(u + 1, n):
            dist = [y for y in range(n) if y not in (u, w) and ((out[y] >> u) ^ (out[y] >> w)) & 1]
            if len(dist) != 1:
                continue
            y = dist[0]
            for z in (u, w):                      # reverse the arc between y and z
                o = flip(out, y, z)
                # {u,w} is now a module of o; quotient Q = o - w (keep u as representative)
                Q = full & ~(1 << w)
                h = ham(o, Q)
                if in_S(n - 2, h) is False:
                    found = True; qvals.add(h)
    return found, qvals

def fc_cut(out, n):
    full = (1 << n) - 1
    for S in range(1, full):
        k = bin(S).count('1')
        if k < 2 or k > n - 2: continue
        c = sum(bin(out[v] & ~S & full).count('1') for v in range(n) if S >> v & 1)
        if c != 3: continue
        A, B = full & ~S, S
        if in_S(n - k - 1, ham(out, A)) is False and in_S(k - 1, ham(out, B)) is False:
            return True
    return False

tab = Counter(); qv = Counter()
for line in open(sys.argv[1]):
    n, code = line.split()[:2]; n = int(n)
    out = adj_from(code, n)
    a = fc_cut(out, n)
    b, vals = quotient_fc(out, n)
    tab[(a, b)] += 1
    for v in vals: qv[v] += 1
print('(two first creators re-connected by 3 arcs, first creator with a doubled vertex and one arc reversed) -> classes')
for k, v in sorted(tab.items()): print('  ', k, v)
print('distinct order-11 first-creator quotients H:', len(qv), ' smallest', sorted(qv)[:8])
