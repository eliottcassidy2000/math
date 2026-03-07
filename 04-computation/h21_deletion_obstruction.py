"""
h21_deletion_obstruction.py — Analyze exactly WHY vertex deletions fail
to preserve cycle-richness. For each vertex v in a cycle-rich tournament:

A deletion T-v fails to be cycle-rich if:
  (a) T-v has a source or sink [score obstruction], or
  (b) Some vertex u in T-v is not in any 3-cycle [cycle obstruction]

For (a): u becomes source in T-v iff score(u) = 1 and u -> v.
         u becomes sink in T-v iff score(u) = n-2 and v -> u.
For (b): u loses all 3-cycles iff ALL 3-cycles through u also go through v.

Key question: can ALL n deletions fail simultaneously when max matching <= 2?

Author: opus-2026-03-07-S43
"""
import random
from itertools import combinations

def random_tournament(n):
    A = [[0]*n for _ in range(n)]
    for i in range(n):
        for j in range(i+1, n):
            if random.random() < 0.5:
                A[i][j] = 1
            else:
                A[j][i] = 1
    return A

def scores(A, n):
    return [sum(A[i]) for i in range(n)]

def is_cycle_rich(A, n):
    sc = scores(A, n)
    for s in sc:
        if s == 0 or s == n-1:
            return False
    for v in range(n):
        found = False
        for a in range(n):
            if a == v: continue
            for b in range(n):
                if b == v or b == a: continue
                if A[v][a] and A[a][b] and A[b][v]:
                    found = True
                    break
            if found: break
        if not found:
            return False
    return True

def get_3cycle_sets(A, n):
    """Return list of 3-cycle vertex sets."""
    cycles = []
    for a in range(n):
        for b in range(a+1, n):
            for c in range(b+1, n):
                if (A[a][b] and A[b][c] and A[c][a]) or (A[a][c] and A[c][b] and A[b][a]):
                    cycles.append((a,b,c))
    return cycles

def max_matching(cycles):
    """Greedy max matching of disjoint 3-cycle sets."""
    used = set()
    count = 0
    for c in cycles:
        if not any(v in used for v in c):
            used.update(c)
            count += 1
    return count

def deletion_analysis(A, n):
    """For each vertex v, determine why T-v fails or succeeds at being cycle-rich."""
    sc = scores(A, n)
    results = []

    for v in range(n):
        # Build T-v
        verts = [u for u in range(n) if u != v]
        n1 = n - 1

        # Check source/sink in T-v
        source_sink = False
        ss_vertices = []
        for u in verts:
            s_u_in_Tv = sum(A[u][w] for w in verts if w != u)
            if s_u_in_Tv == 0 or s_u_in_Tv == n1 - 1:
                source_sink = True
                ss_vertices.append((u, s_u_in_Tv))

        # Check cycle obstruction in T-v
        cycle_obs = False
        co_vertices = []
        for u in verts:
            found_3cyc = False
            for a in verts:
                if a == u: continue
                for b in verts:
                    if b == u or b == a: continue
                    if A[u][a] and A[a][b] and A[b][u]:
                        found_3cyc = True
                        break
                if found_3cyc: break
            if not found_3cyc:
                cycle_obs = True
                co_vertices.append(u)

        good = not source_sink and not cycle_obs
        results.append({
            'v': v, 'good': good, 'score_v': sc[v],
            'source_sink': source_sink, 'ss_verts': ss_vertices,
            'cycle_obs': cycle_obs, 'co_verts': co_vertices
        })

    return results

# Main analysis
random.seed(42)
n = 9
trials = 5000000
no_good_del_count = 0
no_good_del_and_low_mm = 0

# Track obstruction types
both_fail_reasons = {}

print(f"=== Deletion Obstruction Analysis at n={n} ===")
print(f"Sampling {trials} random tournaments...")

for trial in range(trials):
    A = random_tournament(n)
    if not is_cycle_rich(A, n):
        continue

    cycles = get_3cycle_sets(A, n)
    mm = max_matching(cycles)

    analysis = deletion_analysis(A, n)
    num_good = sum(1 for r in analysis if r['good'])

    if num_good == 0:
        no_good_del_count += 1
        t3 = len(cycles)

        # Classify obstruction pattern
        n_ss_only = sum(1 for r in analysis if r['source_sink'] and not r['cycle_obs'])
        n_co_only = sum(1 for r in analysis if r['cycle_obs'] and not r['source_sink'])
        n_both = sum(1 for r in analysis if r['source_sink'] and r['cycle_obs'])
        pattern = (n_ss_only, n_co_only, n_both)
        both_fail_reasons[pattern] = both_fail_reasons.get(pattern, 0) + 1

        if no_good_del_count <= 5:
            print(f"\n  NO GOOD DELETION #{no_good_del_count}: t3={t3}, mm={mm}")
            for r in analysis:
                flags = []
                if r['source_sink']: flags.append(f"SS:{r['ss_verts']}")
                if r['cycle_obs']: flags.append(f"CO:{r['co_verts']}")
                print(f"    v={r['v']} score={r['score_v']}: {' | '.join(flags)}")

        if mm <= 2:
            no_good_del_and_low_mm += 1

    if (trial + 1) % 1000000 == 0:
        print(f"  Progress: {trial+1}/{trials}, no_good={no_good_del_count}, mm<=2_and_no_good={no_good_del_and_low_mm}")

print(f"\n=== RESULTS ===")
print(f"Total cycle-rich found: (many)")
print(f"No good deletion: {no_good_del_count}")
print(f"No good deletion AND mm<=2: {no_good_del_and_low_mm}")
print(f"\nObstruction patterns (ss_only, co_only, both):")
for pat, cnt in sorted(both_fail_reasons.items(), key=lambda x: -x[1]):
    print(f"  {pat}: {cnt}")
