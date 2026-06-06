#!/usr/bin/env python3
"""
S638 / HYP-2316 — The Rado graph as a tournament.

The countable RANDOM tournament (Fraisse limit of finite tournaments) is the
tournament analogue of the Rado graph: unique homogeneous universal countable
tournament, almost surely the result of orienting every pair by a fair coin.
The user's clue "a loop of the input causes a swap of the two outputs" = its
SELF-CONVERSENESS: reversing every arc (swap each vertex's out-set/in-set, the
two 'outputs') gives an isomorphic tournament. The involution is sigma: v -> -v.

This script:
  (A) verifies the finite avatar on 7 vertices = the Paley/QR tournament:
      connection set {1,2,4} = QR mod 7 = cube roots of unity mu_3; antipode -1
      is a non-residue (7 = 3 mod 4) => tournament; x->-x reverses all arcs.
  (B) checks the universal-tournament EXTENSION PROPERTY on small random
      tournaments (how the Rado/random tournament is characterized).
  (C) the perspective tie: # Aut-orbits on k-subsets of the homogeneous
      tournament = # iso types of k-vertex tournaments (A000568) -- the
      'perspective' counting (the user's recurring object) made exact.
No external libs.
"""

# ---------- (A) Paley tournament on F_7 ----------
print("=" * 68)
print("(A) Paley/QR tournament on F_7  (the finite Rado-tournament avatar)")
print("=" * 68)
p = 7
QR = sorted({(x * x) % p for x in range(1, p)})          # nonzero squares
cubes_of_1 = sorted({x for x in range(1, p) if pow(x, 3, p) == 1})
print(f"  nonzero squares mod 7  QR        = {QR}")
print(f"  cube roots of unity   {{x: x^3=1}} = {cubes_of_1}")
print(f"  QR == mu_3 (cube roots of unity)? {QR == cubes_of_1}")
print(f"  -1 mod 7 = {(-1) % p}; is it a residue? {((-1) % p) in QR}  (7=3 mod 4 -> NO)")
negQR = sorted({(-x) % p for x in QR})
print(f"  -QR = {negQR};  QR and -QR disjoint? {set(QR).isdisjoint(negQR)}"
      f";  QR u -QR = F_7*? {sorted(set(QR) | set(negQR)) == list(range(1, p))}")

def adj(x, y):                      # x beats y iff x-y is a QR
    return ((x - y) % p) in QR

# tournament property
ok_irr = all(not adj(x, x) for x in range(p))
ok_tot = all((adj(x, y) ^ adj(y, x)) for x in range(p) for y in range(p) if x != y)
print(f"  irreflexive? {ok_irr};  exactly-one-of(x->y, y->x)? {ok_tot}  => TOURNAMENT")

# the antipodal SWAP: sigma(x) = -x reverses every arc
swap_ok = all((adj((-x) % p, (-y) % p) == adj(y, x)) for x in range(p) for y in range(p))
print(f"  antipode sigma:x->-x reverses every arc  (adj(-x,-y)==adj(y,x))? {swap_ok}")
print("    => sigma is an ANTI-automorphism: it swaps every vertex's two outputs")
print("       (out-set<->in-set). 'A loop of the input swaps the two outputs.'")

# vertex-transitive (translations) and out-degree regular
outdeg = [sum(adj(x, y) for y in range(p) if y != x) for x in range(p)]
print(f"  out-degrees = {outdeg}  (regular = (7-1)/2 = 3: doubly-regular Paley)")

# ---------- (B) the extension property (what makes it the Rado tournament) ----------
print("\n" + "=" * 68)
print("(B) Extension property: for disjoint U,V exists w with U->w->V")
print("=" * 68)
# Check on the Paley-7 tournament for all small disjoint U,V (|U|+|V|<=2),
# the n=2 extension axiom (Paley graphs/tournaments satisfy bounded extension).
from itertools import combinations
def has_witness(U, V):
    return any(all(adj(w, u) for u in U) and all(adj(v, w) for v in V)
               for w in range(p) if w not in U and w not in V)
V_all = list(range(p)); cnt = 0; good = 0
for k in range(0, 3):
    for U in combinations(V_all, k):
        for j in range(0, 3 - k):
            for Vv in combinations([z for z in V_all if z not in U], j):
                if not U and not Vv:
                    continue
                cnt += 1; good += has_witness(set(U), set(Vv))
print(f"  Paley-7 satisfies the (|U|+|V|<=2) extension axiom on "
      f"{good}/{cnt} disjoint pairs  ({'ALL' if good==cnt else 'NOT all'})")
print("  (the countable random tournament satisfies it for ALL finite U,V =")
print("   the defining axiom of the universal homogeneous 'Rado' tournament.)")

# ---------- (C) perspective = orbits = iso-types (A000568) ----------
print("\n" + "=" * 68)
print("(C) Perspectives = Aut-orbits on k-subsets = iso-types of k-tournaments")
print("=" * 68)
# Number of isomorphism classes of tournaments on k vertices = A000568.
A000568 = [1, 1, 1, 2, 4, 12, 56, 456, 6880]   # k = 0,1,2,3,4,5,6,7,8
for k in range(1, 6):
    print(f"  k={k}: #iso-types of k-vertex tournaments = A000568({k}) = {A000568[k]}"
          f"  = #orbits of Aut(T_rado) on k-subsets")
print("  By homogeneity + Ryll-Nardzewski, the universal tournament's orbit")
print("  counts ARE the tournament iso-type counts -- the 'perspective' numbers")
print("  (lrc-perspective-key) made EXACT, with no finite coincidence-breaking.")
print("  n=3: the 2 types = transitive vs 3-CYCLE (the cube-root resonance).")
