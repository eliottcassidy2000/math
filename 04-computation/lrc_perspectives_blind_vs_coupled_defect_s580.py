#!/usr/bin/env python3
r"""
lrc_perspectives_blind_vs_coupled_defect_s580.py    oracle-2026-06-03-S580o

THE "PERSPECTIVES" CURIOSITY (user, repeated rough wording), nailed down.

A "perspective" of a tournament = the tournament VIEWED FROM A DISTINGUISHED VERTEX.
The user's counts:
  n=3: 2 structures (transitive, cyclic); perspectives = 3 (transitive) + 1 (cyclic) = 4.
  n=4: 4 structures; perspectives = 4+4+2+2 = 12.
  "12 is the number of unique structures on 5 vertices."
i.e. the conjecture  #perspectives(n) = #structures(n+1) = A000568(n+1).

A "perspective" read as a VERTEX-ORBIT (how the tournament looks from a vertex, up to
automorphism) gives the ROOTED tournament count:
  rooted(n) = sum over iso classes T of (number of vertex-orbits of T).
This is the OBSERVER-BLIND count: it remembers the shape but FORGETS how the observer
connects to the others.

THE RESOLUTION (and the user's "you may have been misinterpreting it"):
  rooted(n)        = 2, 4, 12, 48, 296, 3040, ...    (blind / vertex-orbit)
  A000568(n+1)     = 2, 4, 12, 56, 456, 6880, ...    (coupled / full n+1 tournament)
  defect           = 0, 0, 0,  8, 160, 3840, ...
The clean identity "perspectives(n)=structures(n+1)" is a SMALL-n ACCIDENT: it holds
EXACTLY for n<=4 and BREAKS at n=5 with defect 8, then explodes. The defect is precisely
the OBSERVER-COUPLING information that the blind (vertex-orbit) reading drops -- the same
thing as the augmentation-index split (observer-blind j=0 vs observer-coupled j!=0) and
the repo's marked-observer thread (THM-381 LRC=observer-is-source; THM-385 observer score
= blocker count; HYP-1977 LRC = a projection-DEFECT over A000568; S369/S370 the 48-vs-56
"perspective gap"; S510/S511 "LRC is NOT a function of the unmarked half-turn class --
anchor on the observer").

LRC reading (where/why it works and doesn't): the runner tournament has m=n-1 vertices.
The defect is ZERO for m<=4 (LRC structurally clean) and turns on at m=5 then EXPLODES --
matching the structural era's difficulty ramp and its death after n=7. Pure unmarked
(blind) structure can settle LRC only while the observer-coupling defect is 0.
"""
from math import gcd, factorial
from itertools import combinations, product, permutations

# ---------- direct enumeration (small n): vertex-orbit (blind) perspectives ----------
def canon(adj, n):
    best = None
    for p in permutations(range(n)):
        key = tuple(adj[p[i]][p[j]] for i in range(n) for j in range(n) if i != j)
        if best is None or key < best:
            best = key
    return best

def vertex_orbits(adj, n):
    autos = [p for p in permutations(range(n))
             if all(adj[i][j] == adj[p[i]][p[j]] for i in range(n) for j in range(n) if i != j)]
    seen = set(); orbits = 0
    for v in range(n):
        if v in seen: continue
        orbits += 1
        for p in autos: seen.add(p[v])
    return orbits

def enumerate_perspectives(n):
    pairs = list(combinations(range(n), 2))
    classes = {}
    for bits in product((0, 1), repeat=len(pairs)):
        adj = [[0] * n for _ in range(n)]
        for (i, j), b in zip(pairs, bits):
            if b: adj[i][j] = 1
            else: adj[j][i] = 1
        c = canon(adj, n)
        if c not in classes: classes[c] = adj
    orbit_list = sorted((vertex_orbits(a, n) for a in classes.values()), reverse=True)
    return len(classes), sum(orbit_list), orbit_list

# ---------- Burnside (any n): rooted(n) and A000568(n+1) ----------
def odd_partitions(n, mx=None):
    if mx is None: mx = n
    if n == 0: yield []; return
    for p in range(min(n, mx), 0, -1):
        if p % 2 == 1:
            for rest in odd_partitions(n - p, p):
                yield [p] + rest

def class_size(n, parts):
    from collections import Counter
    c = Counter(parts); d = factorial(n)
    for l, m in c.items(): d //= (l ** m) * factorial(m)
    return d

def p_orbits(parts):
    s = sum((l - 1) // 2 for l in parts)
    s += sum(gcd(parts[a], parts[b]) for a, b in combinations(range(len(parts)), 2))
    return s

def tournaments(n):
    return sum(class_size(n, P) * 2 ** p_orbits(P) for P in odd_partitions(n)) // factorial(n)

def rooted(n):
    return sum(class_size(n, P) * sum(1 for l in P if l == 1) * 2 ** p_orbits(P)
               for P in odd_partitions(n)) // factorial(n)

def main():
    print("=" * 80)
    print("THE PERSPECTIVES CURIOSITY: observer-BLIND (vertex-orbit) vs COUPLED (n+1)")
    print("=" * 80)

    print("\n(1) Direct enumeration (vertex-orbit perspectives), small n:")
    for n in range(2, 6):
        ncl, tot, orbits = enumerate_perspectives(n)
        print(f"   n={n}: {ncl} structures, perspectives={tot}, per-class orbits={orbits},"
              f"  A000568(n+1)={tournaments(n+1)}  {'== (accident)' if tot==tournaments(n+1) else f'DEFECT {tournaments(n+1)-tot}'}")

    print("\n(2) Burnside defect sequence (all n):")
    print("    n   rooted(n)=blind     A000568(n+1)=coupled    defect       ratio")
    for n in range(2, 11):
        r = rooted(n); c = tournaments(n + 1)
        print(f"   {n:2d}   {r:11d}        {c:13d}     {c-r:11d}    {c/r:.4f}")

    print("\n" + "=" * 80)
    print("READING")
    print("=" * 80)
    print("""  A "perspective" = a tournament viewed from a distinguished vertex.  Read as a
  VERTEX-ORBIT it is OBSERVER-BLIND (forgets how the observer connects); the count is
  rooted(n) = 2,4,12,48,296,... .  The full count A000568(n+1) = 2,4,12,56,456,... is
  OBSERVER-COUPLED.  They coincide EXACTLY for n<=4 and split at n=5 (defect 8), then
  explode.  The clean recursion "perspectives(n)=structures(n+1)" is a small-n ACCIDENT.

  THE MISINTERPRETATION: the LRC runner tournament has been studied UNMARKED (the round
  body A000016; the 64 self-converse classes, S576/S577; the blind vertex-orbit reading).
  But LRC safety is NOT a function of the unmarked class (S511) -- it is observer-COUPLED:
  loneliness = the observer is a SOURCE (THM-381), observer score = blocker count
  (THM-385), LRC = a projection-DEFECT over A000568 (HYP-1977).  The defect computed here
  IS that coupling, = the augmentation-index split (blind j=0 vs coupled j!=0).

  WHERE/WHY LRC WORKS AND DOESN'T: the observer-coupling defect is 0 for vertex count <=4
  (LRC structurally trivial/clean) and turns on at 5, then explodes -- the exact place the
  structural era ramps up (n=5 first computer-assisted) and dies after 7.  Unmarked
  (blind) structure can decide LRC only while the defect is 0.""")

if __name__ == "__main__":
    main()
