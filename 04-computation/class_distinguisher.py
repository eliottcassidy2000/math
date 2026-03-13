#!/usr/bin/env python3
"""
class_distinguisher.py -- kind-pasteur-2026-03-13-S60

What algebraic invariant distinguishes H=93467 from H=92411?
Both have #QR ∈ {2,3} and QR_closure=13/25.

Try:
- Product of elements sum(s_i) mod p
- Subset sum properties
- Intersection with QR cosets
- Cayley table structure
"""

p = 11
m = (p - 1) // 2
S_qr_set = set(j for j in range(1, p) if pow(j, (p-1)//2, p) == 1)

class1 = [4,5,7,11,15,16,20,24,26,27]  # H=93467
class2 = [1,8,9,12,14,17,19,22,23,30]  # H=92411

def bits_to_S(bits):
    S = []
    for j in range(m):
        if bits & (1 << j):
            S.append(j+1)
        else:
            S.append(p-(j+1))
    return S

# Try: product of all elements of S (mod p)
print("Product of elements mod p:")
for label, bits_list in [("H=93467", class1), ("H=92411", class2)]:
    prods = set()
    for bits in bits_list:
        S = bits_to_S(bits)
        prod = 1
        for s in S:
            prod = (prod * s) % p
        prods.add(prod)
    print(f"  {label}: products = {sorted(prods)}")

# Try: sum of elements mod p
print("\nSum of elements mod p:")
for label, bits_list in [("H=93467", class1), ("H=92411", class2)]:
    sums = set()
    for bits in bits_list:
        S = bits_to_S(bits)
        sums.add(sum(S) % p)
    print(f"  {label}: sums = {sorted(sums)}")

# Try: S ∩ QR pattern structure
# Map each pair to 1 if QR chosen, 0 if NQR
# The pattern is a binary string of length 5
# What combinatorial property distinguishes the two classes?
print("\n\nBinary patterns (QR=1 at each pair position):")
patterns_1 = []
patterns_2 = []
for bits in class1:
    S = bits_to_S(bits)
    pat = tuple(1 if s in S_qr_set else 0 for s in S)
    # Actually let's use the pair indices
    pat = []
    for j in range(m):
        s = S[j]
        pat.append(1 if s in S_qr_set else 0)
    patterns_1.append(tuple(pat))

for bits in class2:
    S = bits_to_S(bits)
    pat = []
    for j in range(m):
        s = S[j]
        pat.append(1 if s in S_qr_set else 0)
    patterns_2.append(tuple(pat))

print(f"  H=93467: {sorted(patterns_1)}")
print(f"  H=92411: {sorted(patterns_2)}")

# Try looking at ADJACENCY between chosen QR positions
# The pairs correspond to j=1,...,5 i.e. differences {1,...,5} or {10,9,8,7,6}
# The QR element in pair j has value: 1,9,3,4,5 for j=1,2,3,4,5
qr_values = {1: 1, 2: 9, 3: 3, 4: 4, 5: 5}

print("\n\nQR pair indices chosen:")
print("  H=93467:")
for bits in class1:
    S = bits_to_S(bits)
    qr_idx = [j+1 for j in range(m) if S[j] in S_qr_set]
    print(f"    bits={bits:>2}: QR@{qr_idx}")

print("  H=92411:")
for bits in class2:
    S = bits_to_S(bits)
    qr_idx = [j+1 for j in range(m) if S[j] in S_qr_set]
    print(f"    bits={bits:>2}: QR@{qr_idx}")

# Hypothesis: maybe it's about whether the QR-chosen positions are
# "consecutive" or "separated" in some sense?
# Or maybe it's about the sum of pair indices?
print("\n\nSum and product of QR pair indices:")
for label, bits_list in [("H=93467", class1), ("H=92411", class2)]:
    idx_sums = []
    idx_prods = []
    for bits in bits_list:
        S = bits_to_S(bits)
        qr_idx = [j+1 for j in range(m) if S[j] in S_qr_set]
        idx_sums.append(sum(qr_idx))
        prod = 1
        for i in qr_idx:
            prod *= i
        idx_prods.append(prod)
    print(f"  {label}: sums={sorted(set(idx_sums))}, products={sorted(set(idx_prods))}")

# Try: the DIFFERENCE between QR pair indices
print("\n\nDifferences between QR pair indices:")
for label, bits_list in [("H=93467", class1), ("H=92411", class2)]:
    diffs_set = set()
    for bits in bits_list:
        S = bits_to_S(bits)
        qr_idx = sorted([j+1 for j in range(m) if S[j] in S_qr_set])
        diffs = []
        for i in range(len(qr_idx)):
            for j in range(i+1, len(qr_idx)):
                diffs.append(qr_idx[j] - qr_idx[i])
        diffs_set.add(frozenset(diffs) if diffs else frozenset())
    print(f"  {label}: diff sets = {sorted(tuple(sorted(d)) for d in diffs_set)}")

# The REAL invariant might be the QR character sum
# chi(S) = sum_{s in S} (s/p) where (s/p) is the Legendre symbol
print("\n\nLegendre character sum chi(S) = sum (s/p):")
for label, bits_list in [("H=93467", class1), ("H=92411", class2),
                          ("H=95095", [29]), ("H=93027", [0])]:
    vals = set()
    for bits in bits_list:
        S = bits_to_S(bits)
        chi = sum(1 if s in S_qr_set else -1 for s in S)
        vals.add(chi)
    print(f"  {label}: chi = {sorted(vals)}")

# Try: S ∩ (S+S) where S+S = {s1+s2 mod p : s1, s2 in S}
print("\n\nSumset |S+S| and |S cap (S+S)|:")
for label, bits_list in [("H=95095", [29]), ("H=93467", [4]),
                          ("H=93027", [0]), ("H=92411", [1])]:
    bits = bits_list[0]
    S = bits_to_S(bits)
    S_set = set(S)
    SS = set((s1+s2) % p for s1 in S for s2 in S)
    intersection = S_set & SS
    print(f"  {label}: S={sorted(S)}, |S+S|={len(SS)}, |S cap (S+S)|={len(intersection)}, "
          f"S+S={sorted(SS)}")

# Try: the PRODUCT set S*S = {s1*s2 mod p}
print("\n\nProduct set analysis:")
for label, bits_list in [("H=95095", [29]), ("H=93467", [4]),
                          ("H=93027", [0]), ("H=92411", [1])]:
    bits = bits_list[0]
    S = bits_to_S(bits)
    S_set = set(S)
    SS_prod = set((s1*s2) % p for s1 in S for s2 in S if (s1*s2) % p != 0)
    qr_in_prod = len(SS_prod & S_qr_set)
    nqr_in_prod = len(SS_prod & (set(range(1,p)) - S_qr_set))
    S_in_prod = len(SS_prod & S_set)
    print(f"  {label}: |S*S|={len(SS_prod)}, QR_in={qr_in_prod}, NQR_in={nqr_in_prod}, "
          f"|S cap (S*S)|={S_in_prod}")

# Try: Cayley table — how many products s1*s2 land BACK in S?
print("\n\nCayley self-return count (s1*s2 mod p in S):")
for label, bits_list in [("H=95095", [29]), ("H=93467", [4]),
                          ("H=93027", [0]), ("H=92411", [1])]:
    bits = bits_list[0]
    S = bits_to_S(bits)
    S_set = set(S)
    count = sum(1 for s1 in S for s2 in S if (s1*s2) % p in S_set)
    print(f"  {label}: {count}/{m*m}")

# The DEFINITIVE test: which elements of Z_11* map S to S?
# Aut(S) = {a in Z_11* : aS = S mod p}
print("\n\nStabilizer |Aut(S)| = |{a : aS = S}|:")
for label, bits_list in [("H=95095", [29]), ("H=93467", [4]),
                          ("H=93027", [0]), ("H=92411", [1])]:
    bits = bits_list[0]
    S = bits_to_S(bits)
    S_set = frozenset(S)
    aut = [a for a in range(1, p) if frozenset((a*s) % p for s in S) == S_set]
    print(f"  {label}: Aut = {aut}, |Aut| = {len(aut)}")

# Cross-class: which multipliers map class1 reps to class2 reps?
print("\n\nMultiplier mapping between classes:")
for a in range(1, p):
    S1 = frozenset(bits_to_S(4))  # class1 rep
    aS1 = frozenset((a*s) % p for s in S1)
    # Find which bits this corresponds to
    for bits in range(1 << m):
        if frozenset(bits_to_S(bits)) == aS1:
            if bits in class2:
                print(f"  a={a}: maps bits=4 (H=93467) -> bits={bits} (H=92411)")
            elif bits in class1:
                pass  # same class
            break

print("\nDONE.")
