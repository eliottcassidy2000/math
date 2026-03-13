"""
General degree-d nullity pattern for Paley tournaments.
Compute nullity at ALL degrees for P_7, and at degrees 4-6 for P_11.

Looking for: nullity(d, m) formula.

Known: nullity(d<4) = 0, nullity(4) = m².

opus-2026-03-13-S71b
"""
import numpy as np
from itertools import product as iprod

def get_QR(p):
    QR = set()
    for x in range(1, p):
        QR.add(pow(x, 2, p))
    return QR

def build_diffseqs(p, d):
    """Build all d-diff-seqs for P_p."""
    QR = get_QR(p)
    QR_list = sorted(QR)
    m = (p - 1) // 2

    if d == 0:
        return [()]

    # Recursive construction
    results = []

    def backtrack(seq, ps_list, ps_set):
        if len(seq) == d:
            results.append(tuple(seq))
            return
        for s in QR_list:
            new_ps = (ps_list[-1] + s) % p
            if new_ps in ps_set:
                continue
            seq.append(s)
            ps_list.append(new_ps)
            ps_set.add(new_ps)
            backtrack(seq, ps_list, ps_set)
            seq.pop()
            ps_list.pop()
            ps_set.remove(new_ps)

    backtrack([], [0], {0})
    return results

def compute_nullity(p, d):
    """Compute constraint matrix nullity at degree d."""
    QR = get_QR(p)
    m = (p - 1) // 2
    Ad = build_diffseqs(p, d)

    if d <= 1:
        return 0, len(Ad), 0, len(Ad)

    # Interior faces i=1,...,d-1
    # face i merges positions i-1 and i: replaces (s_{i}, s_{i+1}) with (s_i + s_{i+1})
    # face i is junk iff s_i + s_{i+1} ∉ QR (the merged value is not in QR)
    # Sign of face i: (-1)^i

    junk_faces = {}  # dict mapping junk face tuple to row index
    junk_by_block = [set() for _ in range(d-1)]  # indexed 0..d-2 for faces 1..d-1

    for seq in Ad:
        for i in range(1, d):  # face i
            merged = (seq[i-1] + seq[i]) % p
            if merged not in QR and merged != 0:  # junk face
                face = list(seq)
                face[i-1] = merged
                del face[i]
                face_tuple = tuple(face)
                junk_by_block[i-1].add(face_tuple)
                if face_tuple not in junk_faces:
                    junk_faces[face_tuple] = len(junk_faces)
            elif merged == 0:
                # merged = 0 is also junk (0 ∉ QR)
                face = list(seq)
                face[i-1] = merged
                del face[i]
                face_tuple = tuple(face)
                junk_by_block[i-1].add(face_tuple)
                if face_tuple not in junk_faces:
                    junk_faces[face_tuple] = len(junk_faces)

    n_rows = len(junk_faces)
    n_cols = len(Ad)

    if n_rows == 0:
        return 0, n_cols, 0, n_cols

    # Build constraint matrix
    C = np.zeros((n_rows, n_cols), dtype=np.int8)
    for j, seq in enumerate(Ad):
        for i in range(1, d):
            merged = (seq[i-1] + seq[i]) % p
            if merged not in QR:  # junk face (including 0)
                face = list(seq)
                face[i-1] = merged
                del face[i]
                face_tuple = tuple(face)
                sign = (-1) ** i
                C[junk_faces[face_tuple], j] += sign

    rank = np.linalg.matrix_rank(C.astype(np.float64))
    nullity = n_rows - rank
    omega = n_cols - rank

    # Block analysis
    block_sizes = [len(b) for b in junk_by_block]

    return nullity, n_cols, rank, omega, n_rows, block_sizes

# Full analysis for P_7
print("P_7 (m=3):")
print(f"{'d':>3} {'|A_d|':>8} {'junk':>8} {'rank':>8} {'nullity':>8} {'Omega':>8}")
for d in range(7):
    result = compute_nullity(7, d)
    if len(result) == 4:
        nullity, ad, rank, omega = result
        print(f"{d:>3} {ad:>8} {0:>8} {rank:>8} {nullity:>8} {omega:>8}")
    else:
        nullity, ad, rank, omega, junk, block_sizes = result
        print(f"{d:>3} {ad:>8} {junk:>8} {rank:>8} {nullity:>8} {omega:>8}  blocks={block_sizes}")

print()

# P_11 degrees 0-6 (maybe 7)
print("P_11 (m=5):")
print(f"{'d':>3} {'|A_d|':>8} {'junk':>8} {'rank':>8} {'nullity':>8} {'Omega':>8}")
for d in range(8):
    try:
        result = compute_nullity(11, d)
        if len(result) == 4:
            nullity, ad, rank, omega = result
            print(f"{d:>3} {ad:>8} {0:>8} {rank:>8} {nullity:>8} {omega:>8}")
        else:
            nullity, ad, rank, omega, junk, block_sizes = result
            print(f"{d:>3} {ad:>8} {junk:>8} {rank:>8} {nullity:>8} {omega:>8}  blocks={block_sizes}")
    except MemoryError:
        print(f"{d:>3}  OOM")
        break

# Pattern analysis
print("\n\nNullity pattern analysis:")
print("P_7 nullities: ", end="")
nullities_7 = []
for d in range(7):
    result = compute_nullity(7, d)
    n = result[0]
    nullities_7.append(n)
    print(f"{n}", end=" ")
print()

# Check: are nullities related to binomial coefficients times m²?
m = 3
print(f"  m² = {m**2}")
for d in range(4, 7):
    n = nullities_7[d]
    ratio = n / m**2
    print(f"  d={d}: nullity/m² = {ratio:.4f}")

# Also check: nullity vs C(d-1,2)*m² or other patterns
for d in range(4, 7):
    n = nullities_7[d]
    from math import comb
    for k in range(1, d):
        if comb(d-1, k) * m**2 == n:
            print(f"  d={d}: nullity = C({d-1},{k}) * m² = {n} ✓")
        elif n % comb(d-1, k) == 0:
            q = n // comb(d-1, k)
            if q in [m, m**2, m*(m-1), m*(m+1)]:
                print(f"  d={d}: nullity = C({d-1},{k}) * {q}")
