#!/usr/bin/env python3
"""
WHY are regular maximizers hereditary?

At n=5: 24 regular maximizers (H=15, score (2,2,2,2,2)) all have
every vertex deletion giving H=5 = max H(4).

At n=7: 240 regular maximizers (H=189, score (3,3,3,3,3,3,3)) all have
every vertex deletion giving H=45 = max H(6).

Key hypothesis: For regular tournaments, vertex-transitivity means all
deletions are isomorphic. So if ONE deletion is optimal, all are.

But WHY is even one deletion optimal?

The deletion score at odd n: from regular (d,...,d) with d=(n-1)/2,
deleting v gives score depending on v's out-neighborhood.
For regular T, v has d out-neighbors and d in-neighbors.
The deletion score is sorted(s_u - indicator(u->v) for u != v).

For vertex-transitive T, the deletion score is the same for all v.
At n=5 (d=2): deletion gives score (1,1,2,2) = SC score = max H(4) score.
At n=7 (d=3): deletion gives score (2,2,2,3,3,3) = SC score = max H(6) score.

The question: among all tournaments with SC deletion score, why does
the Paley/regular deletion achieve MAX H?

Let's examine this by checking all tournaments with the right score
and comparing H values.

kind-pasteur-2026-03-06-S18g
"""
import sys
import os
sys.path.insert(0, os.path.join(os.path.dirname(os.path.abspath(__file__)), '..', '03-artifacts', 'code'))
from tournament_lib import tournament_from_bits, hamiltonian_path_count, find_odd_cycles, conflict_graph

MAX_H = {1: 1, 2: 1, 3: 3, 4: 5, 5: 15, 6: 45, 7: 189}

def score_seq(T):
    return tuple(sorted(sum(T[i]) for i in range(len(T))))

def delete_vertex(T, v):
    n = len(T)
    verts = [i for i in range(n) if i != v]
    return [[T[verts[i]][verts[j]] for j in range(n-1)] for i in range(n-1)]

def is_sc(T):
    """Check if tournament is self-complementary."""
    n = len(T)
    s = score_seq(T)
    return all(s[i] + s[n-1-i] == n-1 for i in range(n))

# ============================================================
# n=5: What makes the regular maximizer's deletion optimal?
# ============================================================
print("=" * 70)
print("n=5: REGULAR MAXIMIZER DELETION ANALYSIS")
print("=" * 70)

n = 5
m = n * (n - 1) // 2

# Find ALL n=5 maximizers and check their deletions
regular_maximizers = []
non_regular_maximizers = []
for bits in range(1 << m):
    T = tournament_from_bits(n, bits)
    h = hamiltonian_path_count(T)
    if h != MAX_H[n]:
        continue
    s = score_seq(T)
    if s == (2,2,2,2,2):
        regular_maximizers.append((bits, T))
    else:
        non_regular_maximizers.append((bits, T))

print(f"Regular maximizers: {len(regular_maximizers)}")
print(f"Non-regular maximizers: {len(non_regular_maximizers)}")

# For regular: check deletion details
T_reg = regular_maximizers[0][1]
print(f"\nRegular maximizer (bits={regular_maximizers[0][0]}):")
for v in range(n):
    sub = delete_vertex(T_reg, v)
    h_sub = hamiltonian_path_count(sub)
    s_sub = score_seq(sub)
    sc = is_sc(sub)
    print(f"  v={v}: H={h_sub}, score={s_sub}, SC={sc}")

# For non-regular: check deletion details
T_nreg = non_regular_maximizers[0][1]
print(f"\nNon-regular maximizer (bits={non_regular_maximizers[0][0]}):")
for v in range(n):
    sub = delete_vertex(T_nreg, v)
    h_sub = hamiltonian_path_count(sub)
    s_sub = score_seq(sub)
    sc = is_sc(sub)
    print(f"  v={v}: H={h_sub}, score={s_sub}, SC={sc}")

# ============================================================
# n=4: Among all SC tournaments with score (1,1,2,2),
# which achieve H=5 (max)?
# ============================================================
print(f"\n{'='*70}")
print("n=4: SC TOURNAMENT H DISTRIBUTION (score (1,1,2,2))")
print("=" * 70)

n = 4
m = n * (n - 1) // 2
h_dist = {}
for bits in range(1 << m):
    T = tournament_from_bits(n, bits)
    s = score_seq(T)
    if s != (1, 1, 2, 2):
        continue
    h = hamiltonian_path_count(T)
    sc = is_sc(T)
    key = (h, sc)
    h_dist[key] = h_dist.get(key, 0) + 1

print("H | SC? | count")
for (h, sc) in sorted(h_dist.keys(), reverse=True):
    print(f"  {h} | {sc} | {h_dist[(h, sc)]}")

# ============================================================
# n=6: Among all SC tournaments with score (2,2,2,3,3,3),
# which achieve H=45 (max)?
# ============================================================
print(f"\n{'='*70}")
print("n=6: SC TOURNAMENT H DISTRIBUTION (score (2,2,2,3,3,3))")
print("=" * 70)

n = 6
m = n * (n - 1) // 2
h_dist = {}
for bits in range(1 << m):
    T = tournament_from_bits(n, bits)
    s = score_seq(T)
    if s != (2, 2, 2, 3, 3, 3):
        continue
    h = hamiltonian_path_count(T)
    sc = is_sc(T)
    key = (h, sc)
    h_dist[key] = h_dist.get(key, 0) + 1

print("H | SC? | count")
for (h, sc) in sorted(h_dist.keys(), reverse=True):
    print(f"  {h} | {sc} | {h_dist[(h, sc)]}")

# ============================================================
# KEY: The Paley deletion at n=7 gives a subtournament that is
# NOT just any SC tournament — it's a SPECIFIC one.
# What makes this specific SC tournament achieve H=45?
# ============================================================
print(f"\n{'='*70}")
print("PALEY T_7 DELETION: STRUCTURAL ANALYSIS")
print("=" * 70)

def build_paley(p):
    qr = set()
    for x in range(1, p):
        qr.add((x * x) % p)
    T = [[0]*p for _ in range(p)]
    for i in range(p):
        for j in range(p):
            if i != j and (j - i) % p in qr:
                T[i][j] = 1
    return T

T7 = build_paley(7)
h7 = hamiltonian_path_count(T7)
print(f"T_7: H={h7}")

# Delete v=0
sub = delete_vertex(T7, 0)
h_sub = hamiltonian_path_count(sub)
s_sub = score_seq(sub)
print(f"T_7 - v=0: H={h_sub}, score={s_sub}")

# What are the cycles?
cycles_sub = find_odd_cycles(sub)
cg_sub = conflict_graph(cycles_sub)
c3 = len([c for c in cycles_sub if len(c) == 3])
c5 = len([c for c in cycles_sub if len(c) == 5])
print(f"  Cycles: c3={c3}, c5={c5}, total={len(cycles_sub)}")

# Independence polynomial
nbr = [0] * len(cg_sub)
for i in range(len(cg_sub)):
    for j in range(len(cg_sub)):
        if cg_sub[i][j]:
            nbr[i] |= 1 << j
coeffs = [0] * (len(cg_sub) + 1)
for mask in range(1 << len(cg_sub)):
    ok = True
    seen = 0
    temp = mask
    while temp:
        v = (temp & -temp).bit_length() - 1
        if nbr[v] & seen:
            ok = False
            break
        seen |= 1 << v
        temp &= temp - 1
    if ok:
        coeffs[bin(mask).count('1')] += 1
while len(coeffs) > 1 and coeffs[-1] == 0:
    coeffs.pop()
print(f"  IP = {coeffs}")
print(f"  Type: {'A' if coeffs == [1,14,4] else 'B' if coeffs == [1,20,1] else 'other'}")

# Are ALL H=45 n=6 tournaments with this IP value isomorphic to T_7-v?
# At least they should have the same cycle structure.
print(f"\n  Comparison: H=45 tournaments at n=6 by IP:")
n = 6
m = n * (n - 1) // 2
ip_dist = {}
for bits in range(1 << m):
    T = tournament_from_bits(n, bits)
    h = hamiltonian_path_count(T)
    if h != 45:
        continue
    cycles = find_odd_cycles(T)
    cg = conflict_graph(cycles)
    # Quick IP
    m_c = len(cg)
    nbr2 = [0] * m_c
    for i in range(m_c):
        for j in range(m_c):
            if cg[i][j]:
                nbr2[i] |= 1 << j
    co = [0] * (m_c + 1)
    for mask in range(1 << m_c):
        ok = True
        seen = 0
        temp = mask
        while temp:
            v = (temp & -temp).bit_length() - 1
            if nbr2[v] & seen:
                ok = False
                break
            seen |= 1 << v
            temp &= temp - 1
        if ok:
            co[bin(mask).count('1')] += 1
    while len(co) > 1 and co[-1] == 0:
        co.pop()
    ip_key = tuple(co)
    ip_dist[ip_key] = ip_dist.get(ip_key, 0) + 1

for ip, count in sorted(ip_dist.items()):
    print(f"    IP={list(ip)}: {count} tournaments")

print("\nDone.")
