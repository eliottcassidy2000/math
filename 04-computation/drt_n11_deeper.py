# ⚠️ WARNING: The connection set {1,2,3,5,8} does NOT give a tournament on Z_11
# because S ∩ (-S mod 11) = {3,8} ≠ ∅. All results for this "DRT" are INVALID.
# The ONLY valid circulant DRT at n=11 is the Paley tournament (QR={1,3,4,5,9}).
# See MISTAKE-017.

#!/usr/bin/env python3
"""Deeper analysis of DRTs at n=11."""
import sys, os
sys.path.insert(0, os.path.join(os.path.dirname(os.path.abspath(__file__)), '..', '03-artifacts', 'code'))
from tournament_lib import hamiltonian_path_count

def build_drt(n, D):
    T = [[0]*n for _ in range(n)]
    for i in range(n):
        for j in range(n):
            if i != j and (j-i) % n in D:
                T[i][j] = 1
    return T

def is_drt(T, n):
    k = (n-1)//2
    for v in range(n):
        if sum(T[v]) != k:
            return False
    target = (n-3)//4
    for u in range(n):
        for v in range(u+1, n):
            common = sum(1 for w in range(n) if w != u and w != v and T[u][w] and T[v][w])
            if common != target:
                return False
    return True

# The converse (T^op) of a circulant tournament with D is the circulant with Z_n^* \ D
# For Paley: D = QR, complement = NQR. But Paley is self-converse because
# -1 is a non-residue (for p=3 mod 4), so NQR = -QR, which gives T^op ~ T.

D_paley = {1, 3, 4, 5, 9}
D_other = {1, 2, 3, 5, 8}

# The converse of T(D) has i->j iff (j-i) mod n NOT in D
# which means (i-j) mod n in D, i.e., i->j iff (i-j) mod n in D
# Equivalently, the converse uses D' = {(-d) mod n : d in D} = {n-d : d in D}
D_paley_conv = frozenset((11-d) % 11 for d in D_paley)
D_other_conv = frozenset((11-d) % 11 for d in D_other)

print(f"D_paley = {sorted(D_paley)}")
print(f"-D_paley = {sorted(D_paley_conv)}")
print(f"Same? {D_paley_conv == frozenset(D_paley)}")

print(f"\nD_other = {sorted(D_other)}")
print(f"-D_other = {sorted(D_other_conv)}")
print(f"Same? {D_other_conv == frozenset(D_other)}")

# Check if -D_other is equivalent to D_other (via multiplier)
for m in range(1, 11):
    mD = frozenset((m*d) % 11 for d in D_other)
    if mD == D_other_conv:
        print(f"-D_other = {m} * D_other mod 11")
        break
else:
    print("-D_other is NOT equivalent to D_other -- tournament is NOT self-converse")

# Build T^op and check H
T_other = build_drt(11, D_other)
T_other_op = build_drt(11, D_other_conv)
H_other = hamiltonian_path_count(T_other)
H_other_op = hamiltonian_path_count(T_other_op)
print(f"\nH(T_other) = {H_other}")
print(f"H(T_other^op) = {H_other_op}")
print(f"H(T) = H(T^op)? {H_other == H_other_op}")

# Check: is -D_other also a (11,5,2)-diff set?
D_conv = D_other_conv
diff_counts = {}
for d1 in D_conv:
    for d2 in D_conv:
        if d1 != d2:
            diff = (d1 - d2) % 11
            diff_counts[diff] = diff_counts.get(diff, 0) + 1
is_diff_set = all(diff_counts.get(i, 0) == 2 for i in range(1, 11))
print(f"\n-D_other is (11,5,2)-diff set: {is_diff_set}")
print(f"-D_other DRT: {is_drt(T_other_op, 11)}")

# So there are potentially 4 DRT classes:
# 1. Paley QR (self-converse)
# 2. D_other
# 3. -D_other (converse of D_other)
# 4. Any others?

# Check all 12 difference sets
from itertools import combinations
all_diff_sets = []
for combo in combinations(range(1, 11), 5):
    D = set(combo)
    diff_counts = {}
    for d1 in D:
        for d2 in D:
            if d1 != d2:
                diff = (d1 - d2) % 11
                diff_counts[diff] = diff_counts.get(diff, 0) + 1
    if all(diff_counts.get(i, 0) == 2 for i in range(1, 11)):
        all_diff_sets.append(frozenset(D))

print(f"\nAll (11,5,2)-difference sets ({len(all_diff_sets)}):")
# Group by equivalence: D ~ D' if D' = m*D + t for some m, t
# For cyclic groups, the relevant equivalence is D ~ m*D (multiplier equiv)
# since translation gives isomorphic tournament

equiv_classes = []
used = set()
for i, D in enumerate(all_diff_sets):
    if i in used:
        continue
    cls = [D]
    used.add(i)
    for m in range(2, 11):
        mD = frozenset((m*d) % 11 for d in D)
        for j, D2 in enumerate(all_diff_sets):
            if j not in used and D2 == mD:
                cls.append(D2)
                used.add(j)
    equiv_classes.append(cls)

print(f"Multiplier equivalence classes: {len(equiv_classes)}")
for idx, cls in enumerate(equiv_classes):
    D = cls[0]
    T = build_drt(11, D)
    H = hamiltonian_path_count(T)
    negD = frozenset((11-d) % 11 for d in D)
    self_conv = negD == D or any(frozenset((m*d) % 11 for d in D) == negD for m in range(2, 11))
    print(f"  Class {idx}: {sorted(D)} (size {len(cls)}), H={H}, self-converse={self_conv}")

print("\nDone.")
