#!/usr/bin/env python3
"""audit_A (b): adversarial search for a counterexample to THM-4871 statement 2 (disjoint-orbit lemma).
Part 1: for prime lengths r, s in {3, 5, 7} and every overlap 1 <= o < min(r, s) (and o < r when r = s), take the
  r-cycle g = (0 1 ... r-1) and EVERY s-cycle g' with support {r-o, ..., r-o+s-1}; the group <g, g'> must have even
  order (the lemma says an odd-order group cannot contain two prime cycles with overlapping unequal supports).
  Since every pair of overlapping cycles is Sym-conjugate to one of these, this is exhaustive for r, s <= 7.
  For (3,11), (5,11), (7,11), (11,11) a random sample of 300 s-cycles each.
Part 2: the sanity check quoted in the theorem: over all tournaments of orders 3-8, list every element of Aut(T)
  that is a single cycle of prime length (all other points fixed), count them, and check all supports pairwise.
"""
import itertools, random, subprocess, sys, time
from sympy.combinatorics import Permutation, PermutationGroup
import pynauty

def cycle_perm(seq, n):
    g = list(range(n))
    for i in range(len(seq)):
        g[seq[i]] = seq[(i + 1) % len(seq)]
    return Permutation(g)

def part1():
    random.seed(1)
    total = 0; odd_found = []
    for r in (3, 5, 7, 11):
        for s in (3, 5, 7, 11):
            if r > s: continue
            for o in range(1, min(r, s) + (0 if r == s else 1)):
                if r == s and o >= r: continue
                n = r + s - o
                support2 = list(range(r - o, r - o + s))
                g = cycle_perm(list(range(r)), n)
                # all s-cycles on support2 = cyclic orders = (s-1)! sequences starting at support2[0]
                rest = support2[1:]
                exhaustive = (s <= 7)
                if exhaustive:
                    seqs = [[support2[0]] + list(t) for t in itertools.permutations(rest)]
                else:
                    seqs = []
                    for _ in range(300):
                        t = rest[:]; random.shuffle(t); seqs.append([support2[0]] + t)
                cnt = 0; t0 = time.time()
                for seq in seqs:
                    gp = cycle_perm(seq, n)
                    G = PermutationGroup([g, gp])
                    order = G.order()
                    cnt += 1; total += 1
                    if order % 2 == 1:
                        odd_found.append((r, s, o, seq, order))
                print("r=%2d s=%2d overlap=%d degree=%2d: %d pairs tested (%s), all even order: %s  [%.1fs]" % (
                    r, s, o, n, cnt, "exhaustive" if exhaustive else "sampled",
                    "yes" if not any(x[0] == r and x[1] == s and x[2] == o for x in odd_found) else "NO",
                    time.time() - t0), flush=True)
    print("total pairs tested:", total, " odd-order groups found:", len(odd_found))
    for x in odd_found: print("  COUNTEREXAMPLE:", x)

def digraph6_to_adj(line):
    assert line[0] == '&'
    data = [ord(ch) - 63 for ch in line[1:]]
    n = data[0]; bits = []
    for d in data[1:]:
        for b in range(5, -1, -1): bits.append((d >> b) & 1)
    adj = {v: [] for v in range(n)}
    idx = 0
    for u in range(n):
        for v in range(n):
            if bits[idx]: adj[u].append(v)
            idx += 1
    return n, adj

def closure(gens, n):
    ident = tuple(range(n)); seen = {ident}; frontier = [ident]
    while frontier:
        nxt = []
        for g in frontier:
            for s in gens:
                x = tuple(s[g[i]] for i in range(n))
                if x not in seen: seen.add(x); nxt.append(x)
        frontier = nxt
    return seen

def is_prime(m): return m > 1 and all(m % d for d in range(2, int(m ** 0.5) + 1))

def cycle_type(g):
    n = len(g); seen = [False] * n; cyc = []
    for i in range(n):
        if not seen[i]:
            l = 0; j = i
            while not seen[j]: seen[j] = True; j = g[j]; l += 1
            cyc.append(l)
    return sorted(cyc, reverse=True)

def part2():
    grand = 0; overlaps = 0
    for n in range(3, 9):
        out = subprocess.run(["gentourng", "-z", "-q", str(n)], capture_output=True, text=True).stdout.split()
        cnt = 0; types = 0
        for line in out:
            nn, adj = digraph6_to_adj(line)
            types += 1
            gens = pynauty.autgrp(pynauty.Graph(nn, directed=True, adjacency_dict=adj))[0]
            G = closure([tuple(g) for g in gens], nn) if gens else {tuple(range(nn))}
            cycles = []
            for g in G:
                ct = cycle_type(g)
                if is_prime(ct[0]) and (len(ct) == 1 or ct[1] == 1):
                    cycles.append(frozenset(i for i in range(nn) if g[i] != i))
            cnt += len(cycles)
            for a, b in itertools.combinations(cycles, 2):
                if a != b and a & b: overlaps += 1
        print("n=%d: %d tournament types, %d prime-cycle elements in their automorphism groups" % (n, types, cnt))
        grand += cnt
    print("orders 3-8: %d prime-length cycles in automorphism groups, %d overlapping unequal supports" % (grand, overlaps))

if __name__ == "__main__":
    print("== Part 1: overlapping prime cycles generate even-order groups ==")
    part1()
    print("\n== Part 2: prime cycles in Aut(T), tournaments n <= 8 ==")
    part2()
