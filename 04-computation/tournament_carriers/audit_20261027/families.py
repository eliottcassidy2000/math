"""Referee: THM-4835 parts 6 and 8, plus counterexamples found during the audit.
Run: python families.py"""
import itertools
from math import comb

from core import automorphisms, PermGroup, chi_induced, fixed_set
from decide import decide


def bigon(k):
    n = 2 * k + 1
    s, t = 0, 1
    p1 = list(range(2, k + 1))
    p2 = list(range(k + 1, 2 * k))
    arcs = []
    for path in (p1, p2):
        seq = [s] + path + [t]
        arcs += [(seq[i], seq[i + 1]) for i in range(len(seq) - 1)]
    return n, arcs


def alt_cycle(k):
    n = 2 * k + 1
    arcs = []
    for i in range(0, 2 * k, 2):
        arcs.append((i, (i + 1) % (2 * k)))
        arcs.append((i, (i - 1) % (2 * k)))
    return n, arcs


def theta(p, k):
    u, w = 0, 1
    arcs = []
    nxt = 2
    for _ in range(p):
        seq = [u] + list(range(nxt, nxt + k - 1)) + [w]
        nxt += k - 1
        arcs += [(seq[i], seq[i + 1]) for i in range(len(seq) - 1)]
    n = nxt + (p - 1)
    return n, arcs


def tt4():
    return [(i, j) for i in range(4) for j in range(i + 1, 4)]


def show(name, n, arcs, elems=None, j=1, counts=True):
    if elems is None:
        elems = automorphisms(n, arcs)
    r = decide(n, arcs, elems, jdrop=j, want_counts=counts)
    print(f"{name}: n={n} |H|={r['order']} crit(C)={r['crit']} defects(|K|,delta)={r['defects']} "
          f"linear={r['linear']} weak={r['weak']} count={r.get('count')} "
          f"count_matches_some_X={r.get('count_matches')} linX={r['linear_X'][:2]} weakX={r['weak_X'][:2]}",
          flush=True)
    return r


if __name__ == "__main__":
    print("== 6(a) bigons B_k ==")
    for k in range(2, 6):
        n, arcs = bigon(k)
        r = show(f"B_{k}", n, arcs)
        print(f"   predicted count 2^(C(2k,2)-1) = {2 ** (comb(2 * k, 2) - 1)}; equal: {r['count'] == 2 ** (comb(2 * k, 2) - 1)}")
    print("== 6(b) alternating C_2k + K1 ==")
    for k in range(2, 6):
        n, arcs = alt_cycle(k)
        show(f"altC_{2 * k}+K1", n, arcs)
    print("== 6(c) dense core TT4 + K2 + 2K1 (n=8) ==")
    show("TT4+K2+2K1", 8, tt4() + [(4, 5)])
    print("== 6(d) odd thetas (p=3) ==")
    for k in (2, 3, 4):
        n, arcs = theta(3, k)
        show(f"theta_3^{k}+2K1", n, arcs)
    print("== 8 group dependence: TT4 + 3K1 ==")
    n, arcs = 7, tt4()
    full = automorphisms(n, arcs)
    print("   |Aut| =", len(full))
    Z3 = [tuple(range(4)) + c for c in [(4, 5, 6), (5, 6, 4), (6, 4, 5)]]
    Z2 = [tuple(range(7)), (0, 1, 2, 3, 5, 4, 6)]
    show("TT4+3K1 under Z3", n, arcs, Z3)
    show("TT4+3K1 under Z2", n, arcs, Z2)
    show("TT4+3K1 under S3 (full)", n, arcs, full)
    print("== audit: strong T4 + 3K1 (weak with odd-order defect) ==")
    strong4 = [(0, 1), (0, 2), (3, 0), (1, 2), (1, 3), (2, 3)]
    show("strongT4+3K1", 7, strong4)
    print("== audit: HYP-9435(3) counterexample, K4(C3+source) + K1 + K2 + dirP3 (n=10) ==")
    arcs = [(0, 1), (0, 2), (0, 3), (1, 2), (2, 3), (3, 1), (5, 6), (7, 8), (8, 9)]
    r = show("K4(C3+src)+K1+K2+P3", 10, arcs)
    elems = automorphisms(10, arcs)
    G = PermGroup(elems)
    print("   element orders:", sorted(G.order_of), " involutions:",
          sum(1 for o in G.order_of if o == 2))
    print("== audit: rigid diamond + 2K1 (n=6) ==")
    diamond = [(0, 1), (0, 2), (0, 3), (1, 2), (3, 1)]
    show("diamond+2K1", 6, diamond)
