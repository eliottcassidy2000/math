#!/usr/bin/env python3
"""
THE CARTESIAN-PRODUCT TRICHOTOMY  (monad-explorer-2026-06-07-S2)

Context.  opus-S699g unified Hadwiger-Nelson / unit-distance / Lonely-Runner as
the THREE classical invariants of one forbidden-distance Cayley graph G_d:
    HN  = chromatic number        chi(G_d)
    UD  = edge density / avgdeg    (max |E| on n vertices)
    LRC = independence density     alpha/|V|  (= the lonely/covering density)

My own S1 (unit_distance_product_cap_s1.py) found that AVERAGE DEGREE IS ADDITIVE
under the Cartesian product []:  avgdeg(G[]H) = avgdeg(G) + avgdeg(H).  That is
exactly the Erdos/Minkowski-sum product lower bound for u(n) (e.g. the n=21
extremal K3[]W7), and it is why a product family can build density toward 3N.

THIS SCRIPT establishes the sharp structural distinction that explains why the
product machinery is UNIQUE to the unit-distance invariant among the three:

   THE TRICHOTOMY (verified below; legs 1-2 are exact theorems, leg 3 a sandwich):
     1. avgdeg(G[]H)  = avgdeg(G) + avgdeg(H)            SUM     (UD)   AMPLIFY
        Proof: |E(G[]H)| = |E(G)||V(H)| + |V(G)||E(H)|, |V|=|V(G)||V(H)|.
     2. chi(G[]H)     = max(chi(G), chi(H))              MAX     (HN)   NEUTRAL  [Sabidussi 1957]
     3. i(G)*i(H) <= i(G[]H) <= min(i(G), i(H))          sub-MIN (LRC)  DEGRADE
        where i(G)=alpha(G)/|V(G)| is the independence density.
        Upper: for each h, {g:(g,h) in I} is independent in G => |I|<=alpha(G)|V(H)|.
        Lower: S_G x S_H is independent in G[]H => alpha(G[]H) >= alpha(G)alpha(H).
        For vertex-transitive G, chi_f(G)=1/i(G); so chi_f(G[]H) >= max(chi_f) and
        the upper bound CAN BE STRICT (C5[]Petersen: i=17/50 < 2/5 = min).

   THE THREE INVARIANTS RESPOND TO PRODUCTS IN THREE QUALITATIVELY DISTINCT WAYS:
     UD  edge density  -> AMPLIFIED (additive): products BUILD density.
     HN  chromatic     -> NEUTRAL  (exact max): products neither help nor hurt.
     LRC independence  -> DEGRADED (sub-min):  a product can be STRICTLY LESS
                          lonely than either factor (loneliness is fragile).

   Consequence:  only the unit-distance invariant is amplified, so the
   Erdos/Minkowski-sum product (= tie-induction) lower-bound machinery is an
   intrinsically UNIT-DISTANCE device.  LRC lower bounds (worry-sets) are NEVER
   products -- products would hurt -- matching practice.  This refines opus-S699g
   and explains S1's "N* in [25,28] is NON-product" (the true 3N-crossover is an
   irreducible Moser-lattice graph: a UD-only 'irreducibility premium').

We verify all three legs on the SAME vertex-transitive graphs, then cross-check
against the unit-distance product table.
"""

from itertools import product as iproduct
from fractions import Fraction

# ---------------------------------------------------------------------------
# graph utilities (adjacency as dict: vertex -> set of neighbours)
# ---------------------------------------------------------------------------

def cycle(n):
    g = {i: set() for i in range(n)}
    for i in range(n):
        g[i].add((i + 1) % n)
        g[i].add((i - 1) % n)
    return g

def complete(n):
    return {i: set(j for j in range(n) if j != i) for i in range(n)}

def petersen():
    # standard Kneser K(5,2)
    from itertools import combinations
    verts = list(combinations(range(5), 2))
    idx = {v: k for k, v in enumerate(verts)}
    g = {k: set() for k in range(len(verts))}
    for a in verts:
        for b in verts:
            if a < b and not (set(a) & set(b)):
                g[idx[a]].add(idx[b]); g[idx[b]].add(idx[a])
    return g

def circulant(n, conns):
    """Cayley graph on Z/n with connection set +-conns (undirected)."""
    g = {i: set() for i in range(n)}
    for i in range(n):
        for c in conns:
            g[i].add((i + c) % n); g[i].add((i - c) % n)
    for i in range(n):
        g[i].discard(i)
    return g

def cartesian(G, H):
    """G [] H : (g,h)~(g',h') iff (g=g' and h~h') or (h=h' and g~g')."""
    VG, VH = list(G), list(H)
    P = {(g, h): set() for g in VG for h in VH}
    for g in VG:
        for h in VH:
            for h2 in H[h]:
                P[(g, h)].add((g, h2))
            for g2 in G[g]:
                P[(g, h)].add((g2, h))
    return P

def nv(G):   return len(G)
def ne(G):   return sum(len(v) for v in G.values()) // 2
def avgdeg(G): return Fraction(2 * ne(G), nv(G))

# ---------------------------------------------------------------------------
# exact independence number (branch & bound) and chromatic number (backtrack)
# ---------------------------------------------------------------------------

def independence_number(G):
    verts = list(G)
    n = len(verts)
    idx = {v: i for i, v in enumerate(verts)}
    adj = [0] * n
    for v in verts:
        for w in G[v]:
            adj[idx[v]] |= (1 << idx[w])
    best = 0
    full = (1 << n) - 1

    def bb(cand, size):
        nonlocal best
        if size + bin(cand).count("1") <= best:
            return
        if cand == 0:
            best = max(best, size); return
        # pick lowest-index candidate vertex (pivot on it)
        v = (cand & -cand).bit_length() - 1
        # branch 1: include v
        bb(cand & ~(1 << v) & ~adj[v], size + 1)
        # branch 2: exclude v
        bb(cand & ~(1 << v), size)

    bb(full, 0)
    return best

def chromatic_number(G):
    verts = list(G)
    n = len(verts)
    idx = {v: i for i, v in enumerate(verts)}
    nbr = [[idx[w] for w in G[v]] for v in verts]
    # order by degree desc for speed
    order = sorted(range(n), key=lambda i: -len(nbr[i]))

    def colorable(k):
        color = [-1] * n
        def rec(pos):
            if pos == n:
                return True
            v = order[pos]
            used = set(color[w] for w in nbr[v] if color[w] != -1)
            maxc = min(k, max([color[order[p]] for p in range(pos)] + [-1]) + 2)
            for c in range(maxc):
                if c not in used:
                    color[v] = c
                    if rec(pos + 1):
                        return True
                    color[v] = -1
            return False
        return rec(0)

    k = 1
    while not colorable(k):
        k += 1
    return k

# ---------------------------------------------------------------------------
# verification
# ---------------------------------------------------------------------------

def report(name, G):
    a = avgdeg(G)
    al = independence_number(G)
    ch = chromatic_number(G)
    idns = Fraction(al, nv(G))
    return dict(name=name, V=nv(G), E=ne(G), avgdeg=a, alpha=al, idens=idns, chi=ch)

def show(r):
    print(f"  {r['name']:14s} V={r['V']:3d} E={r['E']:3d}  avgdeg={str(r['avgdeg']):>6s}"
          f"  alpha={r['alpha']:3d}  i={str(r['idens']):>6s}  chi={r['chi']}")

ATOMS = {
    "C5":        cycle(5),
    "C7":        cycle(7),
    "K3":        complete(3),
    "K4":        complete(4),
    "Petersen":  petersen(),
    "Circ7{1,2}":circulant(7, [1, 2]),     # a small Cayley/distance-type graph
}

print("=" * 78)
print("ATOMIC vertex-transitive graphs")
print("=" * 78)
atom_rep = {}
for nm, G in ATOMS.items():
    r = report(nm, G); atom_rep[nm] = r; show(r)

print()
print("=" * 78)
print("CARTESIAN PRODUCTS  G[]H  -- verifying the SUM / MAX / MIN trichotomy")
print("=" * 78)
print("  avgdeg(G[]H) =? avgdeg(G)+avgdeg(H)   [SUM, unit-distance -- AMPLIFY]")
print("  chi(G[]H)    =? max(chi(G),chi(H))     [MAX, Hadwiger-Nelson -- NEUTRAL]")
print("  i(G)i(H) <= i(G[]H) <= min(i(G),i(H))  [sub-MIN, lonely-runner -- DEGRADE]")
print()

PAIRS = [
    ("C5", "C5"),
    ("C5", "K3"),
    ("C5", "C7"),
    ("K3", "K3"),
    ("K3", "K4"),
    ("Petersen", "K3"),
    ("C7", "K3"),
    ("Circ7{1,2}", "K3"),
    ("C5", "Petersen"),
    ("C7", "C7"),
    ("Circ7{1,2}", "C5"),
]

ok_sum = ok_max = 0          # exact-equality legs
sandwich_ok = min_strict = 0 # sandwich leg
total = 0
for an, bn in PAIRS:
    GA, GB = ATOMS[an], ATOMS[bn]
    P = cartesian(GA, GB)
    rp = report(f"{an}[]{bn}", P)
    ra, rb = atom_rep[an], atom_rep[bn]

    pred_sum = ra['avgdeg'] + rb['avgdeg']
    pred_max = max(ra['chi'], rb['chi'])
    lo = ra['idens'] * rb['idens']        # i(G)i(H)
    hi = min(ra['idens'], rb['idens'])    # min(i(G),i(H))

    s_ok = (rp['avgdeg'] == pred_sum)
    m_ok = (rp['chi'] == pred_max)
    sand = (lo <= rp['idens'] <= hi)
    strict = (rp['idens'] < hi)
    ok_sum += s_ok; ok_max += m_ok; sandwich_ok += sand; min_strict += strict
    total += 1

    show(rp)
    tag = "  (= min)" if rp['idens'] == hi else "  (< min: DEGRADED)"
    print(f"      SUM avgdeg {rp['avgdeg']}=?{pred_sum} {'OK' if s_ok else 'FAIL'}"
          f"  | MAX chi {rp['chi']}=?{pred_max} {'OK' if m_ok else 'FAIL'}"
          f"  | sub-MIN  {lo} <= {rp['idens']} <= {hi} {'OK' if sand else 'FAIL'}{tag}")

print()
print(f"SUM  avgdeg additive (exact)        : {ok_sum}/{total}")
print(f"MAX  chi = max       (exact)        : {ok_max}/{total}")
print(f"sub-MIN sandwich i(G)i(H)<=i<=min(i): {sandwich_ok}/{total}"
      f"   (of which STRICTLY < min, i.e. DEGRADED: {min_strict})")

print()
print("=" * 78)
print("INTERPRETATION  (the unification, opus-S699g, made quantitative)")
print("=" * 78)
print("""
The three classical invariants of the forbidden-distance Cayley graph respond to
the Cartesian product in three QUALITATIVELY DISTINCT ways:
   UD  edge density / avgdeg  -> AMPLIFIED (exactly additive): products build density.
   HN  chromatic number       -> NEUTRAL  (exactly max): products neither help nor hurt.
   LRC independence density   -> DEGRADED (sub-min): a product can be STRICTLY LESS
                                 lonely than either factor (witness C5[]Petersen:
                                 i = 17/50 < 2/5 = min). Loneliness is fragile.

ONLY the unit-distance invariant is amplified.  So the Erdos/Minkowski-sum product
construction -- the tie-induction 'build density one factor at a time' machinery,
the n=21 extremal K3[]W7, the 3^3=27 cube tie (avgdeg 2+2+2=6) of S1 -- is an
intrinsically UNIT-DISTANCE device.  LRC lower bounds (worry-sets) are NEVER built
from products, because products would strictly hurt the lonely density -- which
matches practice and now has a structural reason.  HN chromatic forcing is product-
neutral, so it too needs irreducible graphs (Moser spindle / de Grey).

This is the SAME fault line as opus-S699g's 'integrality gap chi > chi_f = Vitali
wall':  the additive (product) lower bound is the easy/algebraic side; the true
optimum is an irreducible geometric object the additive machine cannot see.  S1's
'N* in [25,28] is NON-product, while products first beat 3N only at N=32' is the
unit-distance face of exactly this gap (the 'irreducibility premium' = 32 - N*).
""")
