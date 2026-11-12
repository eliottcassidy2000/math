#!/usr/bin/env python3
"""Verifier, small checks written from the definitions (pure python, no numba, no lane code):
 A. roots of rho_G(s) = 1 for the three gadget sets, from the constants printed by vk.exe neutral 8 / 9 (my recount);
 B. Theorem 2' : kappa_X({e, f}) = 2 Pi_2(X - ef) - Pi_2(X'), literal enumeration, all labelled tournaments of order <= 4
    and random ones of order 5, 6;
 C. Theorem 2 (a): W and kappa of (TT_a, P) for every P, a <= 6, literal enumeration of vertex orders;
 D. Theorem 2 (c): no unit transport except pairs at distance 1, 2, from the closed formulas, a <= 12, M up to 10^6
    (exact divisibility argument instead of a loop over M);
 E. twin equation counts (finisher's F5)."""
import itertools
import random
from math import comb


# ---------------------------------------------------------------- A
def read_units(fn):
    g = []
    for line in open(fn):
        if line.startswith("UNIT"):
            t = line.split()
            n, code, P, s, H, Pi, W, rX, kap, u, sig, u3, eps, nu = map(int, t[1:15])
            prim = int(t[15].split("=")[1])
            g.append(dict(n=n, H=H, rho=rX, u=u, sig=sig, u3=u3, c=2 * (2 * u - sig + u3), prim=prim))
    return g


def rho_of(gad, s, terms=120):
    t = 0.0
    for c, sig in gad:
        t += sum((c * 3.0 ** m - sig) ** (-s) for m in range(terms))
        t += c ** (-s) * 2.0 ** (-s) / ((1 - 2.0 ** (-s)) * (1 - 3.0 ** (-s)))
    return t


def root(gad):
    lo, hi = 0.5, 1.0
    for _ in range(80):
        mid = (lo + hi) / 2
        if rho_of(gad, mid) > 1:
            lo = mid
        else:
            hi = mid
    return lo


tp = [(10, 1), (26, 1)]                     # transitive pairs: c = 2 gamma, sigma = 1
g8 = [(g["c"], g["sig"]) for g in read_units("v_neutral_8.out") if g["prim"]]
u9 = read_units("v_neutral_9.out")
g9 = [(g["c"], g["sig"]) for g in u9 if g["prim"]]
print(f"A   transitive pairs only: root {root(tp):.6f}, rho(1) = {rho_of(tp, 1.0):.6f}   (PROOFS: 0.788507, 0.429511)")
print(f"A   + {len(g8)} gadgets of order 8: root {root(tp + g8):.6f}, rho(1) = {rho_of(tp + g8, 1.0):.6f}   (PROOFS: 0.792831; lane output 0.431174)")
print(f"A   + {len(g9)} primitive gadgets of order 9: root {root(tp + g8 + g9):.6f}, rho(1) = {rho_of(tp + g8 + g9, 1.0):.6f}   (PROOFS: 0.797044, 0.432524)")
print(f"A   contributions to rho(0.8): order 8 {rho_of(g8, 0.8):.4f}, order 9 {rho_of(g9, 0.8):.4f}; to rho(1): {rho_of(g8, 1.0):.4f}, {rho_of(g9, 1.0):.4f}   (PROOFS: 0.0169, 0.0177; 0.0017, 0.0013)")
p9 = [g for g in u9 if g["prim"]]
print("A   order-9 primitive unit neutral gadgets: r_X = -1:", sum(1 for g in p9 if g["rho"] == -1), "H in", (min(g["H"] for g in p9 if g["rho"] == -1), max(g["H"] for g in p9 if g["rho"] == -1)),
      "; r_X = +1:", sum(1 for g in p9 if g["rho"] == 1), "H in", (min(g["H"] for g in p9 if g["rho"] == 1), max(g["H"] for g in p9 if g["rho"] == 1)),
      "; min over all gadgets of c - sigma - u:", min(g["c"] - g["sig"] - g["u"] for g in u9))


# ---------------------------------------------------------------- literal counts
def literal(n, arc):
    """arc[i][j] = True iff i -> j.  Returns Pi_2, Q[e][f], q[e][f] (D3 possibly empty), Pi_2 with e, f not consecutive (dict by pair)."""
    def path(seq):
        return all(arc[seq[k]][seq[k + 1]] for k in range(len(seq) - 1))
    Pi = 0
    Q = [[0] * n for _ in range(n)]
    q = [[0] * n for _ in range(n)]
    noncons = {}
    for p in itertools.permutations(range(n)):
        pre = [True] * (n + 1)
        for k in range(2, n + 1):
            pre[k] = pre[k - 1] and arc[p[k - 2]][p[k - 1]]
        suf = [True] * (n + 1)
        for k in range(n - 2, -1, -1):
            suf[k] = suf[k + 1] and arc[p[k]][p[k + 1]]
        for k in range(n + 1):
            if pre[k] and suf[k]:
                Pi += 1
                if 1 <= k < n:
                    Q[p[k - 1]][p[k]] += 1
                # pairs consecutive inside one of the two paths
                cons = set()
                for j in range(n - 1):
                    if j + 1 != k:
                        cons.add(frozenset((p[j], p[j + 1])))
                for e in range(n):
                    for f in range(e + 1, n):
                        if frozenset((e, f)) not in cons:
                            noncons[(e, f)] = noncons.get((e, f), 0) + 1
            if pre[k] and k >= 1:
                for k2 in range(k + 1, n + 1):
                    if path(p[k:k2]) and suf[k2]:
                        q[p[k - 1]][p[k]] += 1
    return Pi, Q, q, noncons


def tour_from_bits(n, bits):
    arc = [[False] * n for _ in range(n)]
    b = 0
    for i in range(n):
        for j in range(i + 1, n):
            if (bits >> b) & 1:
                arc[i][j] = True
            else:
                arc[j][i] = True
            b += 1
    return arc


# ---------------------------------------------------------------- B
random.seed(20261111)
cases = bad = 0
todo = [(n, bits) for n in (2, 3, 4) for bits in range(1 << comb(n, 2))]
todo += [(5, random.getrandbits(10)) for _ in range(40)] + [(6, random.getrandbits(15)) for _ in range(8)]
for n, bits in todo:
    arc = tour_from_bits(n, bits)
    Pi, Q, q, noncons = literal(n, arc)
    for e in range(n):
        for f in range(n):
            if e != f and arc[e][f]:
                arc2 = [row[:] for row in arc]
                arc2[e][f] = False
                arc2[f][e] = True
                Pi2 = literal(n, arc2)[0]
                kap = Pi - 2 * (q[e][f] + q[f][e])
                key = (min(e, f), max(e, f))
                if kap != 2 * noncons.get(key, 0) - Pi2:
                    bad += 1
                cases += 1
print(f"B   Theorem 2': kappa_X({{e, f}}) = 2 Pi_2(X - ef) - Pi_2(X'), literal enumeration: {cases} arcs, mismatches {bad}")

# ---------------------------------------------------------------- C
cases = bad = 0
for a in range(1, 7):
    arc = [[i < j for j in range(a)] for i in range(a)]
    Pi, Q, q, _ = literal(a, arc)
    for P in range(1, 1 << a):
        vs = [v for v in range(a) if (P >> v) & 1]
        W = sum(Q[e][f] for e in vs for f in vs if e != f)
        kap = Pi - 2 * sum(q[e][f] for e in vs for f in vs if e != f)
        Wf = sum(2 ** (j - i - 1) + (1 if j - i == 1 else 0) for i in vs for j in vs if i < j)
        # kappa = 2^a (1 - sum phi(d)), phi(d) = (1 + 3^(d-1)) / 2^d, in integers
        kf = 2 ** a - sum(2 ** (a - (j - i)) * (1 + 3 ** (j - i - 1)) for i in vs for j in vs if i < j)
        if W != Wf or kap != kf or Pi != 2 ** a:
            bad += 1
        cases += 1
print(f"C   Theorem 2 (a): W_X(P) and kappa_X(P) for X = TT_a, a <= 6, every P, literal enumeration: {cases} gadgets, mismatches {bad}")

# ---------------------------------------------------------------- D
# r+ = (1 - W) r + kappa M.  |r+| = 1 with r = +-1 and M >= 1 (M >= 3 if r = -1): solve kappa M = t - (1 - W) r exactly.
viol = []
ng = 0
for a in range(1, 13):
    for P in range(0, 1 << a):
        vs = [v for v in range(a) if (P >> v) & 1]
        pairs = [(i, j) for i in vs for j in vs if i < j]
        W = sum(2 ** (j - i - 1) + (1 if j - i == 1 else 0) for i, j in pairs)
        kap = 2 ** a - sum(2 ** (a - (j - i)) * (1 + 3 ** (j - i - 1)) for i, j in pairs)
        ng += 1
        good = len(pairs) == 1 and pairs[0][1] - pairs[0][0] <= 2
        for r in (1, -1):
            for t in (1, -1):
                rhs = t - (1 - W) * r
                if kap == 0:
                    sol = rhs == 0          # every M
                else:
                    sol = rhs % kap == 0 and rhs // kap >= (1 if r == 1 else 3)
                if sol and not good:
                    viol.append((a, P, r, t))
                if good and not (kap == 0 and (1 - W) == -1):
                    viol.append(("good-fail", a, P))
print(f"D   Theorem 2 (c): TT_a, a <= 12, every P ({ng} gadgets), every M >= 1 (seed) / M >= 3 (anti-seed), exact: unit transport only by pairs at distance 1, 2: {not viol}")


# ---------------------------------------------------------------- E
def I_of(b):
    return sum(2 ** (k - i - 1) for i in range(len(b)) if b[i] == 1 for k in range(i + 1, len(b)) if b[k] == 0)


cnt = {}
for L in range(0, 9):
    for LY in range(0, L + 1):
        LZ = L - LY
        for Yv in range(1 << LY):
            Y = [(Yv >> (LY - 1 - i)) & 1 for i in range(LY)]
            if Y and Y[-1] == 0:
                continue
            for Zv in range(1 << LZ):
                Z = [(Zv >> (LZ - 1 - i)) & 1 for i in range(LZ)]
                if Z and Z[0] == 0:
                    continue
                y = Y + [0, 0] + Z
                J = sum(1 for k in range(len(y) - 1) if y[k] == 0 and y[k + 1] == 1)
                val = sum(Y[i] << (LY - 1 - i) for i in range(LY))
                C = (1 if y[0] == 1 else 0) + (1 if y[-1] == 0 else 0) + J + I_of(Y) + I_of(Z)
                if val == C + 1:
                    cnt[L] = cnt.get(L, 0) + 1
print("E   twin equation val(Y) = 1 + C, maximal 0-blocks, solutions by |Y| + |Z| = 2..8:", [cnt.get(L, 0) for L in range(2, 9)], "(finisher: 1, 3, 5, 7, 12, 17, 22); |Y| + |Z| < 2:", [cnt.get(L, 0) for L in (0, 1)])
