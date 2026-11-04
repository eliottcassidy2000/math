#!/usr/bin/env python3
"""audit_A (a): THM-4871 statement 1 bookkeeping, recomputed from scratch.
  * b_k(4) by Burnside over odd-cycle-type permutations (4^cycles * 2^(pair orbits) / z_lambda), k <= 8,
    and by explicit orbit enumeration (min over S_k) for k <= 4;
  * t_m (dilation classes of skew sets) and e_m (multiplier-even dilation classes of symmetric sets) by brute force
    for m = 3, 5, 7, 9, 11, 13, 25, and the prime formula t_p = (1/(p-1)) sum_{d | p-1, d odd} phi(d) 2^((p-1)/(2d));
  * c = (p+q)/2 + 2k + C(k,2) and |N| = p phi(p) q phi(q) k! against a direct closure of the normaliser for small n;
  * Z_9: 4 dilation classes of circulant tournaments -> 3 isomorphism types (pynauty); Z_25: 208 -> ? types;
  * the nauty keying subtlety (certificate alone merges coloured types: 82 vs 88 at k = 3).
"""
import itertools, math, sys
from fractions import Fraction
from collections import Counter, defaultdict
import pynauty

def partitions(n, maxpart=None):
    if maxpart is None: maxpart = n
    if n == 0:
        yield (); return
    for p in range(min(n, maxpart), 0, -1):
        for rest in partitions(n - p, p):
            yield (p,) + rest

def zlam(lam):
    r = 1
    for part, mult in Counter(lam).items():
        r *= part ** mult * math.factorial(mult)
    return r

def pair_orbits(lam):
    s = 0
    for l in lam:
        s += (l - 1) // 2 if l % 2 else l // 2
    for i in range(len(lam)):
        for j in range(i + 1, len(lam)):
            s += math.gcd(lam[i], lam[j])
    return s

def b_k4_burnside(k):
    tot = Fraction(0)
    for lam in partitions(k):
        if all(l % 2 for l in lam):
            tot += Fraction(4 ** len(lam) * 2 ** pair_orbits(lam), zlam(lam))
    assert tot.denominator == 1
    return int(tot)

def b_k4_orbits(k):
    pairs = [(i, j) for i in range(k) for j in range(i + 1, k)]
    perms = list(itertools.permutations(range(k)))
    idx = {pr: t for t, pr in enumerate(pairs)}
    seen = set()
    for bits in itertools.product([0, 1], repeat=len(pairs)):
        for col in itertools.product(range(4), repeat=k):
            best = None
            for g in perms:
                nb = [0] * len(pairs)
                for (i, j), b in zip(pairs, bits):
                    gi, gj = g[i], g[j]
                    # arc i->j if b == 1 else j->i ; image arc g(i)->g(j)
                    src, dst = (gi, gj) if b else (gj, gi)
                    nb[idx[(min(src, dst), max(src, dst))]] = 1 if src < dst else 0
                ncol = [0] * k
                for v in range(k): ncol[g[v]] = col[v]
                key = (tuple(ncol), tuple(nb))
                if best is None or key < best: best = key
            seen.add(best)
    return len(seen)

def units(m): return [a for a in range(1, m) if math.gcd(a, m) == 1]

def skew_sets(m):
    half = list(range(1, (m + 1) // 2))
    for bits in itertools.product([0, 1], repeat=len(half)):
        yield frozenset((d if b == 0 else m - d) for d, b in zip(half, bits))

def symmetric_sets(m):
    half = list(range(1, (m + 1) // 2))
    for bits in itertools.product([0, 1], repeat=len(half)):
        yield frozenset(x for d, b in zip(half, bits) if b for x in (d, m - d))

def dilation_classes(m, sets):
    seen, reps = set(), []
    for S in sets:
        if S in seen: continue
        orb = {frozenset((a * s) % m for s in S) for a in units(m)}
        seen |= orb; reps.append(S)
    return reps

def reversal_parity(g, edges):
    return sum(1 for (u, v) in edges if g[u] > g[v]) & 1

def multiplier_even(m, C):
    edges = sorted({(min(i, (i + c) % m), max(i, (i + c) % m)) for i in range(m) for c in C})
    for a in units(m):
        if frozenset((a * c) % m for c in C) == C:
            g = [(a * i) % m for i in range(m)]
            if reversal_parity(g, edges): return False
    return True

def t_formula(p):
    tot = Fraction(0)
    for d in range(1, p):
        if (p - 1) % d == 0 and d % 2 == 1:
            phi = sum(1 for x in range(1, d + 1) if math.gcd(x, d) == 1)
            tot += phi * 2 ** ((p - 1) // (2 * d))
    tot /= (p - 1)
    assert tot.denominator == 1
    return int(tot)

def circ_tournament_cert(m, S):
    adj = {i: [(i + s) % m for s in S] for i in range(m)}
    return pynauty.certificate(pynauty.Graph(m, directed=True, adjacency_dict=adj))

def normaliser_closure(p, q, k):
    n = p + q + k
    h = list(range(n))
    for i in range(p): h[i] = (i + 1) % p
    for j in range(q): h[p + j] = p + (j + 1) % q
    h = tuple(h)
    # generators: h, dilations, Sym(F)
    gens = [h]
    for (off, m) in ((0, p), (p, q)):
        for a in units(m):
            g = list(range(n))
            for i in range(m): g[off + i] = off + (a * i) % m
            gens.append(tuple(g))
    if k >= 2:
        g = list(range(n)); g[p + q], g[p + q + 1] = p + q + 1, p + q; gens.append(tuple(g))
        g = list(range(n))
        for i in range(k): g[p + q + i] = p + q + (i + 1) % k
        gens.append(tuple(g))
    seen = {tuple(range(n))}; frontier = [tuple(range(n))]
    while frontier:
        nxt = []
        for g in frontier:
            for s in gens:
                x = tuple(s[g[i]] for i in range(n))
                if x not in seen: seen.add(x); nxt.append(x)
        frontier = nxt
    # every element must normalise <h>: g h g^-1 is a power of h
    powers = set()
    x = tuple(range(n))
    for _ in range(p * q):
        powers.add(x); x = tuple(h[x[i]] for i in range(n))
    for g in seen:
        ginv = [0] * n
        for i in range(n): ginv[g[i]] = i
        conj = tuple(g[h[ginv[i]]] for i in range(n))
        assert conj in powers
    # and conversely the full normaliser has order pq k! phi(pq) (centraliser pq k! times phi(pq) conjugate generators)
    return len(seen)

def phi(m): return sum(1 for x in range(1, m + 1) if math.gcd(x, m) == 1)

if __name__ == "__main__":
    print("== b_k(4) by Burnside (odd cycle types) ==")
    bb = [b_k4_burnside(k) for k in range(0, 9)]
    print("k=0..8:", bb)
    print("claimed: [1, 4, 16, 88, 704, 8912, 188288, 6847200, 437069312]",
          "MATCH" if bb == [1, 4, 16, 88, 704, 8912, 188288, 6847200, 437069312] else "MISMATCH")
    print("== b_k(4) by explicit orbit enumeration (min over S_k), k <= 4 ==")
    bo = [b_k4_orbits(k) for k in range(0, 5)]
    print("k=0..4:", bo, "MATCH" if bo == bb[:5] else "MISMATCH")

    print("\n== t_m, e_m by brute force; formula for primes ==")
    for m in (3, 5, 7, 9, 11, 13, 25):
        tcl = dilation_classes(m, skew_sets(m))
        gcl = dilation_classes(m, symmetric_sets(m))
        ecl = [C for C in gcl if multiplier_even(m, C)]
        line = "m=%2d  t_m=%d  (sym classes %d)  e_m=%d" % (m, len(tcl), len(gcl), len(ecl))
        if all(m % d for d in range(2, m)):
            line += "  formula t_p=%d" % t_formula(m)
        # isomorphism types among the dilation classes of circulant tournaments
        certs = {circ_tournament_cert(m, S) for S in tcl}
        line += "  iso types of circulant tournaments=%d" % len(certs)
        print(line)

    print("\n== c and |N| ==")
    for (p, q, k) in [(3, 5, 0), (3, 5, 1), (3, 5, 2), (3, 5, 3), (5, 7, 2), (9, 5, 0), (9, 5, 1), (3, 7, 3)]:
        c = (p + q) // 2 + 2 * k + k * (k - 1) // 2
        Nord = p * phi(p) * q * phi(q) * math.factorial(k)
        Ncl = normaliser_closure(p, q, k)
        print("(p,q,k)=(%d,%d,%d) c=%d |N| formula=%d closure=%d %s" % (p, q, k, c, Nord, Ncl, "OK" if Nord == Ncl else "MISMATCH"))

    print("\n== nauty keying subtlety: coloured 3-tournaments keyed by certificate alone vs (sizes, certificate) ==")
    k = 3
    pairs = [(i, j) for i in range(k) for j in range(i + 1, k)]
    keys_cert, keys_full = set(), set()
    for bits in itertools.product([0, 1], repeat=3):
        adj = {v: [] for v in range(k)}
        for (i, j), b in zip(pairs, bits):
            if b: adj[i].append(j)
            else: adj[j].append(i)
        for col in itertools.product(range(4), repeat=k):
            cells = [set(v for v in range(k) if col[v] == cc) for cc in range(4)]
            present = tuple(cc for cc in range(4) if cells[cc])
            g = pynauty.Graph(k, directed=True, adjacency_dict=adj, vertex_coloring=[cells[cc] for cc in present])
            cert = pynauty.certificate(g)
            keys_cert.add(cert)
            keys_full.add((tuple((cc, len(cells[cc])) for cc in present), cert))
    print("certificate-only keys:", len(keys_cert), " (sizes, certificate) keys:", len(keys_full), " b_3(4) =", bb[3])
