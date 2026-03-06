#!/usr/bin/env python3
"""
SC Maximizer at n=8: Generate SC tournaments efficiently using fpf involution.

For sigma = (0,1)(2,3)(4,5)(6,7), the anti-aut constraint is:
  T[i][j] + T[sigma(i)][sigma(j)] = 1  for all i != j

Within-orbit pairs {i, sigma(i)}: constraint is T[i][sigma(i)] + T[sigma(i)][i] = 1,
which is automatically true. These arcs are FREE (4 bits).

Cross-orbit pairs: {i,j} and {sigma(i),sigma(j)} are linked.
T[i][j] determines T[sigma(i)][sigma(j)] = 1 - T[i][j].
24 cross-orbit pairs in 12 orbits (12 bits).

Total: 16 free bits, 2^16 = 65536 SC tournaments per sigma.

kind-pasteur-2026-03-06-S18f
"""
import sys
import os
import random
sys.path.insert(0, os.path.join(os.path.dirname(os.path.abspath(__file__)), '..', '03-artifacts', 'code'))
from tournament_lib import hamiltonian_path_count

random.seed(42)

def score_seq(T):
    return tuple(sorted(sum(T[i]) for i in range(len(T))))

def is_sc_score(s, n):
    return all(s[i] + s[n-1-i] == n-1 for i in range(n))

def generate_sc_tournaments_simple(n, sigma):
    """Generate all SC tournaments for given fpf involution sigma."""
    pairs = [(i, j) for i in range(n) for j in range(i+1, n)]

    # Classify pairs into free (within-orbit) and paired (cross-orbit)
    free = []
    paired_reps = []
    seen = set()

    for i, j in pairs:
        if (i, j) in seen:
            continue
        si, sj = sigma[i], sigma[j]
        img = (min(si, sj), max(si, sj))
        if img == (i, j):
            # Within-orbit pair (sigma swaps i,j) -- FREE
            free.append((i, j))
            seen.add((i, j))
        else:
            # Cross-orbit pair -- linked with img
            paired_reps.append(((i, j), img))
            seen.add((i, j))
            seen.add(img)

    total_bits = len(free) + len(paired_reps)

    tournaments = []
    for bits in range(1 << total_bits):
        T = [[0]*n for _ in range(n)]

        bit_idx = 0
        # Free arcs
        for i, j in free:
            val = (bits >> bit_idx) & 1
            T[i][j] = val
            T[j][i] = 1 - val
            bit_idx += 1

        # Paired arcs
        for (i, j), (si, sj) in paired_reps:
            val = (bits >> bit_idx) & 1
            T[i][j] = val
            T[j][i] = 1 - val
            # Anti-aut constraint: T[sigma(i)][sigma(j)] = 1 - T[i][j]
            T[sigma[i]][sigma[j]] = 1 - val
            T[sigma[j]][sigma[i]] = val
            bit_idx += 1

        tournaments.append(T)

    return tournaments

# ============================================================
n = 8
print("=" * 70)
print(f"SC MAXIMIZER TEST AT n={n}")
print("=" * 70)

sigma = [1, 0, 3, 2, 5, 4, 7, 6]  # (01)(23)(45)(67)
print(f"Sigma = {sigma}")

sc_tours = generate_sc_tournaments_simple(n, sigma)
print(f"Generated {len(sc_tours)} SC tournaments")

# Group by score, find max H in each SC score class
sc_by_score = {}
global_max_h = 0
count = 0

for T in sc_tours:
    s = score_seq(T)
    h = hamiltonian_path_count(T)
    count += 1

    if s not in sc_by_score:
        sc_by_score[s] = (h, 1)
    else:
        old_max, old_cnt = sc_by_score[s]
        sc_by_score[s] = (max(old_max, h), old_cnt + 1)

    if h > global_max_h:
        global_max_h = h

    if count % 10000 == 0:
        print(f"  Processed {count}/{len(sc_tours)}, current max H = {global_max_h}")

print(f"\nTotal SC tournaments: {count}")
print(f"Global max H among SC (sigma1): {global_max_h}")
print(f"OEIS A038375 max H at n=8: 661")
print(f"SC achieves global max: {global_max_h >= 661}")

# Show all SC score classes
print(f"\nSC score classes ({len(sc_by_score)} total):")
for s in sorted(sc_by_score.keys()):
    max_h, cnt = sc_by_score[s]
    is_sc = is_sc_score(s, n)
    if is_sc:
        print(f"  {s}: {cnt} tours, max H = {max_h} {'***' if max_h >= 661 else ''}")

# ============================================================
# Try another sigma to find more SC tournaments
# ============================================================
print(f"\n{'='*70}")
print("TRYING ADDITIONAL SIGMA")
print("=" * 70)

sigma2 = [2, 3, 0, 1, 6, 7, 4, 5]  # (02)(13)(46)(57)
sc_tours2 = generate_sc_tournaments_simple(n, sigma2)

for T in sc_tours2:
    s = score_seq(T)
    h = hamiltonian_path_count(T)
    if s not in sc_by_score:
        sc_by_score[s] = (h, 1)
    else:
        old_max, old_cnt = sc_by_score[s]
        sc_by_score[s] = (max(old_max, h), old_cnt + 1)
    if h > global_max_h:
        global_max_h = h

print(f"After sigma2: global max H = {global_max_h}")

sigma3 = [7, 6, 5, 4, 3, 2, 1, 0]  # (07)(16)(25)(34)
sc_tours3 = generate_sc_tournaments_simple(n, sigma3)

for T in sc_tours3:
    s = score_seq(T)
    h = hamiltonian_path_count(T)
    if s not in sc_by_score:
        sc_by_score[s] = (h, 1)
    else:
        old_max, old_cnt = sc_by_score[s]
        sc_by_score[s] = (max(old_max, h), old_cnt + 1)
    if h > global_max_h:
        global_max_h = h

print(f"After sigma3: global max H = {global_max_h}")

# Final summary
print(f"\n{'='*70}")
print("FINAL SUMMARY")
print("=" * 70)

print(f"Global max H among all SC tournaments found: {global_max_h}")
print(f"OEIS A038375 max at n=8: 661")
print(f"SC achieves or exceeds global max: {global_max_h >= 661}")

print(f"\nAll SC score classes with max H:")
for s in sorted(sc_by_score.keys()):
    max_h, cnt = sc_by_score[s]
    is_sc = is_sc_score(s, n)
    marker = " *** GLOBAL MAX ***" if max_h >= 661 else ""
    if is_sc:
        print(f"  {s}: max H = {max_h}{marker}")

print("\nDone.")
