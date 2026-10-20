#!/usr/bin/env python3
"""THM-4722 certificate. (1) Jordan lemma check: a 3-cycle and a 5-cycle sharing >= 1 point generate a group of even order
(all overlap patterns). (2) Hall violator for n = 5..9 from the census: S = tournaments with Aut containing a Sylow-3 or a
Sylow-5 subgroup of S_n; |S| = tau_T(P3) + tau_T(P5) (no tournament has both, by (1)); |N(S)| <= tau_E(P3) + tau_E(P5) - 1
(the empty graph has both). Values tau from subgroups_fix*.out."""
import itertools, re
from math import gcd
def mul(a, b): return tuple(a[b[i]] for i in range(len(a)))
def closure(gens, n):
    idt = tuple(range(n)); S = {idt}; fr = [idt]
    while fr:
        x = fr.pop()
        for g in gens:
            y = mul(g, x)
            if y not in S: S.add(y); fr.append(y)
    return S
def cyc(cycles, n):
    g = list(range(n))
    for c in cycles:
        for i in range(len(c)): g[c[i]] = c[(i + 1) % len(c)]
    return tuple(g)
# (1): 3-cycle on {0,1,2}; 5-cycle on 5 points overlapping in k = 1,2,3 points, all arrangements up to relabelling
for k in (1, 2, 3):
    n = 3 + 5 - k
    c = cyc([[0, 1, 2]], n)
    pts = list(range(3 - k, 3 - k + 5))   # overlap = {3-k..2}
    orders = set()
    for perm in itertools.permutations(pts):
        d = cyc([list(perm)], n); G = closure([c, d], n); orders.add(len(G))
    print(f"overlap {k}: |<c,d>| in {sorted(orders)} -> all even: {all(o % 2 == 0 for o in orders)}")
# (2): read tau values
vals = {}
for f in ['subgroups_fix.out', 'subgroups_fix_n8.out', 'subgroups_fix_n9.out']:
    cur = None
    for line in open(f):
        m = re.match(r"=== n = (\d+)", line)
        if m: cur = int(m.group(1)); continue
        m = re.match(r"\s+H = (.*?) \(\|H\|=(\d+).*tau = (\d+), eps = (\d+)", line)
        if m and cur: vals[(cur, m.group(1))] = (int(m.group(2)), int(m.group(3)), int(m.group(4)))
sylow = {5: ("Z3 type (3, 1, 1)", "Z5 type (5,)"), 6: ("Z3xZ3 = <(012),(345)>", "Z5 type (5, 1)"),
         7: ("Z3xZ3 = <(012),(345)>", "Z5 type (5, 1, 1)"), 8: ("Z3xZ3 = <(012),(345)>", "Z5 type (5, 1, 1, 1)"),
         9: ("Z3 wr Z3 (Sylow-3 of S_9, order 81)", "Z5 type (5, 1, 1, 1, 1)")}
print("n | tau_T(P3) tau_T(P5) | tau_E(P3) tau_E(P5) | |S| | |N(S)| <= | deficit >=")
for n, (p3, p5) in sylow.items():
    o3, t3, e3 = vals[(n, p3)]; o5, t5, e5 = vals[(n, p5)]
    assert t3 == e3 and t5 == e5
    S = t3 + t5; N = e3 + e5 - 1
    print(f"{n} | {t3} {t5} | {e3} {e5} | {S} | {N} | {S - N}")
