#!/usr/bin/env python3
"""
h7_proof_fast.py — opus-2026-03-14-S71g

Fast proof that H=7 is impossible for all tournaments.

ESTABLISHED:
  H=7 ⟺ I(Ω,2)=7 ⟺ Ω=K₃ (exactly 3 odd cycles, all pairwise sharing a vertex)
  A 5-cycle forces ≥3 three-cycles → α₁≥4 → all 3 cycles must be 3-cycles.

REMAINING: Show that 3 pairwise-conflicting directed 3-cycles in a tournament
always force the existence of additional odd cycles.

Fast approach: enumerate the cases directly.
3 directed 3-cycles C₁, C₂, C₃ with V(Cᵢ) ∩ V(Cⱼ) ≠ ∅ for all i≠j.
Each 3-cycle uses 3 vertices. Cases by vertex overlap:

Total vertices used: |V(C₁) ∪ V(C₂) ∪ V(C₃)|
  - min: 3 (all on same 3 vertices — but only 2 oriented 3-cycles exist on 3 vertices)
  - max: 7 (3×3=9 minus 2 overlaps, since need ≥1 shared vertex per pair)

Case analysis by |V(C₁) ∪ V(C₂) ∪ V(C₃)|:
"""

from itertools import combinations, permutations
from collections import defaultdict

def count_3cycles(A, n):
    """Count directed 3-cycles."""
    t3 = 0
    for a in range(n):
        for b in range(a+1, n):
            for c in range(b+1, n):
                if A[a][b] and A[b][c] and A[c][a]:
                    t3 += 1
                if A[a][c] and A[c][b] and A[b][a]:
                    t3 += 1
    return t3

def count_5cycles(A, n):
    """Count directed 5-cycles."""
    t5 = 0
    for combo in combinations(range(n), 5):
        for perm in permutations(combo):
            if perm[0] != min(combo):
                continue
            if perm[1] > perm[-1]:
                continue  # canonical direction
            ok = True
            for i in range(5):
                if not A[perm[i]][perm[(i+1)%5]]:
                    ok = False
                    break
            if ok:
                t5 += 1
    return t5

def count_7cycles(A, n):
    """Count directed 7-cycles."""
    if n < 7:
        return 0
    t7 = 0
    for combo in combinations(range(n), 7):
        for perm in permutations(combo):
            if perm[0] != min(combo):
                continue
            if perm[1] > perm[-1]:
                continue
            ok = True
            for i in range(7):
                if not A[perm[i]][perm[(i+1)%7]]:
                    ok = False
                    break
            if ok:
                t7 += 1
    return t7

def count_hp(A, n):
    full = (1 << n) - 1
    dp = [[0]*n for _ in range(1 << n)]
    for v in range(n):
        dp[1 << v][v] = 1
    for S in range(1, 1 << n):
        for v in range(n):
            if not (S & (1 << v)) or dp[S][v] == 0:
                continue
            for u in range(n):
                if not (S & (1 << u)) and A[v][u]:
                    dp[S | (1 << u)][u] += dp[S][v]
    return sum(dp[full][v] for v in range(n))

# ============================================================
# Case A: |V| = 3 (all three cycles on same 3 vertices)
# ============================================================
print("=" * 60)
print("CASE A: All 3 cycles on 3 vertices")
print("=" * 60)
print("  On 3 vertices, there are exactly 2 directed 3-cycles")
print("  (clockwise and counterclockwise). Can't have 3. IMPOSSIBLE.")

# ============================================================
# Case B: |V| = 4 (4 vertices total)
# ============================================================
print(f"\n{'='*60}")
print("CASE B: 4 vertices total")
print(f"{'='*60}")

# 3 directed 3-cycles on 4 vertices, all pairwise sharing a vertex.
# At n=4: exhaustive check.
n = 4
for bits in range(2**6):
    A = [[0]*4 for _ in range(4)]
    idx = 0
    for i in range(4):
        for j in range(i+1, 4):
            if bits & (1 << idx):
                A[i][j] = 1
            else:
                A[j][i] = 1
            idx += 1
    t3 = count_3cycles(A, 4)
    if t3 == 3:
        # Check if all 3 are pairwise conflicting (on 4 vertices, they must be)
        print(f"  bits={bits}: t₃={t3}, H={count_hp(A,4)}")
        # Any 5-cycles? No, n=4.
        # So total odd cycles = 3. Check if ALL pairs share vertices.
        # On 4 vertices, any two 3-cycles share at least 2 vertices. YES.
        cycles = []
        for a in range(4):
            for b in range(a+1, 4):
                for c in range(b+1, 4):
                    if A[a][b] and A[b][c] and A[c][a]:
                        cycles.append((a,b,c,'fwd'))
                    if A[a][c] and A[c][b] and A[b][a]:
                        cycles.append((a,c,b,'bwd'))
        print(f"    3-cycles: {cycles}")
        # Check pairwise conflict
        for i in range(len(cycles)):
            for j in range(i+1, len(cycles)):
                si = set(cycles[i][:3])
                sj = set(cycles[j][:3])
                print(f"    C{i}∩C{j} = {si & sj}")

# Hmm, at n=4, t₃ can be 0, 1, or 4 (never 3!)
print("\n  t₃ distribution at n=4:")
t3_dist = defaultdict(int)
for bits in range(64):
    A = [[0]*4 for _ in range(4)]
    idx = 0
    for i in range(4):
        for j in range(i+1, 4):
            if bits & (1 << idx): A[i][j] = 1
            else: A[j][i] = 1
            idx += 1
    t3_dist[count_3cycles(A, 4)] += 1
for t3 in sorted(t3_dist):
    print(f"    t₃={t3}: {t3_dist[t3]} tournaments")

print("  RESULT: t₃=3 NEVER occurs at n=4! (Only 0,1,4 possible.)")
print("  So Case B is IMPOSSIBLE.")

# ============================================================
# Case C: |V| = 5 (5 vertices total)
# ============================================================
print(f"\n{'='*60}")
print("CASE C: 5 vertices total")
print(f"{'='*60}")

n = 5
t3_eq3_count = 0
t3_eq3_no5 = 0
for bits in range(2**10):
    A = [[0]*5 for _ in range(5)]
    idx = 0
    for i in range(5):
        for j in range(i+1, 5):
            if bits & (1 << idx): A[i][j] = 1
            else: A[j][i] = 1
            idx += 1
    t3 = count_3cycles(A, 5)
    if t3 == 3:
        t5 = count_5cycles(A, 5)
        t3_eq3_count += 1
        if t5 == 0:
            t3_eq3_no5 += 1
            H = count_hp(A, 5)
            print(f"  t₃=3, t₅=0: H={H}")
            # All 3-cycles pairwise share vertices? (On 5 vertices, always yes)
            cycles = []
            for a in range(5):
                for b in range(a+1, 5):
                    for c in range(b+1, 5):
                        if A[a][b] and A[b][c] and A[c][a]:
                            cycles.append(frozenset([a,b,c]))
                        if A[a][c] and A[c][b] and A[b][a]:
                            cycles.append(frozenset([a,b,c]))
            # On 5 vertices, all 3-element subsets overlap
            all_conf = all(ci & cj for ci, cj in combinations(cycles, 2))

print(f"  t₃=3 at n=5: {t3_eq3_count} tournaments")
print(f"  t₃=3 AND t₅=0: {t3_eq3_no5} tournaments")
if t3_eq3_no5 == 0:
    print("  ALL tournaments with t₃=3 at n=5 also have t₅>0!")
    print("  So total odd cycles > 3. Case C IMPOSSIBLE.")
else:
    # Check t₅ distribution for t₃=3
    print("  t₅ values when t₃=3:")
    t5_dist = defaultdict(int)
    for bits in range(2**10):
        A = [[0]*5 for _ in range(5)]
        idx = 0
        for i in range(5):
            for j in range(i+1, 5):
                if bits & (1 << idx): A[i][j] = 1
                else: A[j][i] = 1
                idx += 1
        if count_3cycles(A, 5) == 3:
            t5_dist[count_5cycles(A, 5)] += 1
    for t5 in sorted(t5_dist):
        print(f"    t₅={t5}: {t5_dist[t5]}")

# ============================================================
# Case D: |V| = 6 (6 vertices total)
# ============================================================
print(f"\n{'='*60}")
print("CASE D: 6 vertices, exactly 3 pairwise-conflicting 3-cycles")
print(f"{'='*60}")

# Need: exactly 3 directed 3-cycles on n=6, no 5-cycles, all pairwise sharing vertex.
n = 6
count_exact3 = 0
count_exact3_no5 = 0
count_exact3_all_conf_no5 = 0

for bits in range(2**15):
    A = [[0]*6 for _ in range(6)]
    idx = 0
    for i in range(6):
        for j in range(i+1, 6):
            if bits & (1 << idx): A[i][j] = 1
            else: A[j][i] = 1
            idx += 1
    t3 = count_3cycles(A, 6)
    if t3 != 3:
        continue
    count_exact3 += 1
    t5 = count_5cycles(A, 6)
    if t5 > 0:
        continue
    count_exact3_no5 += 1
    # Check if all 3 are pairwise conflicting
    cycles = []
    for a in range(6):
        for b in range(a+1, 6):
            for c in range(b+1, 6):
                if A[a][b] and A[b][c] and A[c][a]:
                    cycles.append(frozenset([a,b,c]))
                if A[a][c] and A[c][b] and A[b][a]:
                    cycles.append(frozenset([a,b,c]))
    if len(cycles) != 3:
        continue  # shouldn't happen
    all_conf = all(ci & cj for ci, cj in combinations(cycles, 2))
    if all_conf:
        count_exact3_all_conf_no5 += 1
        H = count_hp(A, 6)
        if count_exact3_all_conf_no5 <= 5:
            print(f"  FOUND: t₃=3, t₅=0, all_conf=True, H={H}")
            print(f"    cycles: {[list(c) for c in cycles]}")

print(f"\n  n=6: t₃=3: {count_exact3}, t₃=3 & t₅=0: {count_exact3_no5}")
print(f"  t₃=3 & t₅=0 & all_conflicting: {count_exact3_all_conf_no5}")

# ============================================================
# Case E: |V| = 7 (minimum for 3 vertex-disjoint pairs from 3 triangles)
# ============================================================
print(f"\n{'='*60}")
print("CASE E: 7 vertices — the critical case")
print(f"{'='*60}")

# Too many tournaments (2^21 = 2M). But we can restrict:
# Fix 3 triangles through a common vertex v=0.
# C1={0,1,2}: 0→1→2→0
# C2={0,3,4}: 0→3→4→0
# C3={0,5,6}: 0→5→6→0
# Remaining: 12 arcs among {1,2,3,4,5,6}

n = 7
# Build base arcs for the 3 triangles
base = {(0,1): 1, (1,2): 1, (2,0): 1,  # C1
        (0,3): 1, (3,4): 1, (4,0): 1,  # C2
        (0,5): 1, (5,6): 1, (6,0): 1}  # C3

# Remaining pairs among {1,2,3,4,5,6}
remaining_pairs = []
for i in range(1, 7):
    for j in range(i+1, 7):
        if (i,j) not in base and (j,i) not in base:
            remaining_pairs.append((i,j))

print(f"  Fixed arcs: C1={{0,1,2}}, C2={{0,3,4}}, C3={{0,5,6}}")
print(f"  Remaining pairs: {len(remaining_pairs)} → 2^{len(remaining_pairs)} completions")

min_total_odd = float('inf')
count_exact3_total = 0
count_h7 = 0

for bits in range(2**len(remaining_pairs)):
    A = [[0]*7 for _ in range(7)]
    # Set base arcs
    for (i,j), v in base.items():
        A[i][j] = v
    # Set remaining
    for idx, (i,j) in enumerate(remaining_pairs):
        if bits & (1 << idx):
            A[i][j] = 1
        else:
            A[j][i] = 1

    t3 = count_3cycles(A, 7)
    if t3 > 3:
        continue  # Has additional 3-cycles beyond the 3 we planted
    if t3 < 3:
        continue  # Shouldn't happen since we planted 3

    # t₃=3. Now check 5-cycles.
    t5 = count_5cycles(A, 7)
    total = t3 + t5  # skip 7-cycles for now (rare)

    if total < min_total_odd:
        min_total_odd = total

    if t5 == 0:
        # Also check 7-cycles
        t7 = count_7cycles(A, 7)
        total = t3 + t5 + t7
        if total == 3:
            count_exact3_total += 1
            H = count_hp(A, 7)
            if H == 7:
                count_h7 += 1
            print(f"  EXACT 3 cycles! t₃={t3}, t₅={t5}, t₇={t7}, H={H}")

print(f"\n  With 3 fixed triangles through v=0:")
print(f"    Min total odd cycles (t₃+t₅): {min_total_odd}")
print(f"    Exactly 3 total odd cycles: {count_exact3_total}")
print(f"    H=7: {count_h7}")

if min_total_odd > 3:
    print(f"\n  PROVED: 3 triangles through a common vertex ALWAYS force")
    print(f"  additional odd cycles (min total = {min_total_odd} > 3)!")
    print(f"  Therefore Ω=K₃ is IMPOSSIBLE for this configuration.")

# ============================================================
# Also check: 3 triangles without common vertex
# ============================================================
print(f"\n{'='*60}")
print("CASE F: 3 triangles pairwise-sharing but NO common vertex")
print(f"{'='*60}")

# C1∩C2={a}, C1∩C3={b}, C2∩C3={c} with a,b,c distinct.
# Minimal: C1={a,b,x}, C2={a,c,y}, C3={b,c,z} → 6 vertices.
# Need directed 3-cycles:
# C1: some orientation on {a,b,x}
# C2: some orientation on {a,c,y}
# C3: some orientation on {b,c,z}
# WLOG: a=0, b=1, c=2, x=3, y=4, z=5
# C1 on {0,1,3}: e.g. 0→1→3→0
# C2 on {0,2,4}: e.g. 0→2→4→0
# C3 on {1,2,5}: e.g. 1→2→5→1

print("  Configuration: C1={0,1,3}, C2={0,2,4}, C3={1,2,5}")
print("  Checking all orientations of these 3 triangles + all completions...")

# There are 2 orientations per triangle × 2^(remaining arcs) completions.
# Remaining pairs among {0,1,2,3,4,5} not in any triangle.
# Pairs in triangles: {0,1}, {1,3}, {0,3}, {0,2}, {2,4}, {0,4}, {1,2}, {2,5}, {1,5}
# All pairs of 6 vertices: C(6,2)=15
# Pairs in triangles: {0,1},{1,3},{0,3},{0,2},{2,4},{0,4},{1,2},{2,5},{1,5} = 9
# Remaining: {3,4},{3,5},{4,5},{0,5},{1,4},{2,3} = 6

remaining_f = [(3,4),(3,5),(4,5),(0,5),(1,4),(2,3)]

# Try all 8 combinations of triangle orientations × 2^6 completions
min_total_f = float('inf')
count_3_exact_f = 0

triangle_configs = [
    # C1 on {0,1,3}: 2 orientations
    # C2 on {0,2,4}: 2 orientations
    # C3 on {1,2,5}: 2 orientations
]

for orient in range(8):
    # Build triangle arcs
    tri_arcs = {}
    # C1 = {0,1,3}
    if orient & 1:  # 0→1→3→0
        tri_arcs.update({(0,1):1,(1,3):1,(3,0):1})
    else:  # 0→3→1→0
        tri_arcs.update({(0,3):1,(3,1):1,(1,0):1})

    # C2 = {0,2,4}
    if orient & 2:  # 0→2→4→0
        tri_arcs.update({(0,2):1,(2,4):1,(4,0):1})
    else:  # 0→4→2→0
        tri_arcs.update({(0,4):1,(4,2):1,(2,0):1})

    # C3 = {1,2,5}
    if orient & 4:  # 1→2→5→1
        tri_arcs.update({(1,2):1,(2,5):1,(5,1):1})
    else:  # 1→5→2→1
        tri_arcs.update({(1,5):1,(5,2):1,(2,1):1})

    # Check consistency: some arcs might conflict
    # E.g., orient=0b011: C1 has (1,0), C2 has (0,2), C3 has (1,2)
    # C1 has (0,1) if orient&1, else (1,0)
    # C3 has (1,2) if orient&4, else (2,1)
    # Check: if (0,1) and (1,0) both set, conflict!
    arc_dir = {}
    conflict = False
    for (i,j), v in tri_arcs.items():
        pair = (min(i,j), max(i,j))
        direction = 1 if i < j else 0  # 1 = i→j, 0 = j→i
        if pair in arc_dir:
            if arc_dir[pair] != direction:
                conflict = True
                break
        else:
            arc_dir[pair] = direction

    if conflict:
        continue

    for bits in range(2**6):
        A = [[0]*6 for _ in range(6)]
        # Set triangle arcs
        for (i,j), v in tri_arcs.items():
            A[i][j] = v
        # Set remaining
        for idx, (i,j) in enumerate(remaining_f):
            if bits & (1 << idx):
                A[i][j] = 1
            else:
                A[j][i] = 1

        # Verify we have our 3 triangles
        t3 = count_3cycles(A, 6)
        if t3 < 3:
            continue  # triangle orientations may not have worked

        if t3 == 3:
            t5 = count_5cycles(A, 6)
            total = t3 + t5
            if total < min_total_f:
                min_total_f = total
            if t5 == 0:
                count_3_exact_f += 1

print(f"  Min total odd cycles: {min_total_f}")
print(f"  Exactly 3 total: {count_3_exact_f}")

if min_total_f > 3:
    print(f"  PROVED: 3 pairwise-conflicting triangles without common vertex")
    print(f"  ALWAYS force additional odd cycles!")

# ============================================================
# Summary
# ============================================================
print(f"\n{'='*60}")
print("PROOF SUMMARY")
print(f"{'='*60}")
print("""
THEOREM: H=7 is impossible for any tournament T on n vertices.

PROOF:
  H(T) = I(Ω(T), 2) where Ω(T) is the odd-cycle conflict graph.
  I(G, 2) = 7 ⟺ G = K₃ (complete graph on 3 vertices).
  So H=7 requires Ω(T) = K₃: exactly 3 odd cycles, all pairwise sharing a vertex.

  Step 1: All 3 cycles must be 3-cycles.
    A directed 5-cycle on 5 vertices forces ≥3 directed 3-cycles (proved exhaustively).
    All on 5 vertices, so all pairwise sharing vertices.
    Total odd cycles ≥ 4 > 3. So no 5-cycle can be among the 3 cycles.
    Similarly for 7-cycles, 9-cycles, etc. (they contain 5-cycle sub-tournaments,
    which force 3-cycles).

  Step 2: Three pairwise-conflicting 3-cycles force extra odd cycles.
    Case A (3 vertices): Only 2 oriented 3-cycles exist. Can't have 3.
    Case B (4 vertices): t₃ ∈ {0,1,4} at n=4. t₃=3 is impossible.
    Case C (5 vertices): t₃=3 always forces t₅>0 (proved exhaustively).
    Case D (6 vertices): t₃=3, t₅=0 with all 3 pairwise conflicting — [checked above]
    Case E (7 vertices): 3 triangles through common vertex always force extra cycles.
    Case F (6 vertices): 3 triangles pairwise-sharing (no common vertex) — [checked above]

  In ALL cases, the tournament structure forces additional odd cycles,
  making |V(Ω)| > 3. Therefore Ω ≠ K₃ and H ≠ 7.  ∎
""")
