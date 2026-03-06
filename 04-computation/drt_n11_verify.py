#!/usr/bin/env python3
"""Verify DRT properties at n=11 and analyze differences between the two classes."""
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

D_paley = {1, 3, 4, 5, 9}
D_other = {1, 2, 3, 5, 8}

T1 = build_drt(11, D_paley)
T2 = build_drt(11, D_other)

print(f"Paley DRT verified: {is_drt(T1, 11)}")
print(f"Other DRT verified: {is_drt(T2, 11)}")

print(f"\nPaley D: {sorted(D_paley)}, complement: {sorted(set(range(1,11)) - D_paley)}")
print(f"Other D: {sorted(D_other)}, complement: {sorted(set(range(1,11)) - D_other)}")

# Check if D_other = m*D_paley for any m
for m in range(1, 11):
    mD = frozenset((m*d) % 11 for d in D_paley)
    if mD == frozenset(D_other):
        print(f"\nD_other = {m} * D_paley mod 11")
        break
else:
    print("\nD_other is NOT a scalar multiple of D_paley -- genuinely different class")

# Check which multipliers preserve D_other
print("\nMultipliers preserving D_other:")
for m in range(1, 11):
    mD = frozenset((m*d) % 11 for d in D_other)
    if mD == frozenset(D_other):
        print(f"  {m}: {sorted(mD)}")

# This tells us the automorphism group
# For Paley: closed under all QR multipliers (5 of them), plus cyclic => |Aut| = 55
# For other: check

print(f"\nH(Paley) = 95095")
print(f"H(other) = 69311")
print(f"H(Paley)/H(other) = {95095/69311:.4f}")
print(f"H(Paley) - H(other) = {95095 - 69311}")
print(f"\nc3 ratio: 55/44 = {55/44:.4f}")
print(f"c5 ratio: 594/407 = {594/407:.4f}")

# Check: is {1,2,3,5,8} a known DRT construction?
# At n=11, the non-Paley DRT comes from the Szekeres-Whiteman construction
# or possibly a non-cyclic difference set.
# Let's check: is D_other related to D_paley by complementing a coset?
# The QR form an index-2 subgroup. The cosets of an index-2 subgroup of QR:
# QR = {1,3,4,5,9}. Index-2 subgroup of QR has order (p-1)/4 = 10/4 = 2.5
# Not integer! So QR doesn't have index-2 subgroups.
# For p=11: phi(p)=10, QR has order 5, prime => no proper subgroups.
# So Szekeres-Whiteman 2-block doesn't apply at p=11 (needs p=3 mod 8).
# 11 mod 8 = 3, so it should work... let me think.
# The 2-block SW construction for p=3 mod 8:
# Take C = QR. Partition C into C0, C1 where Ci = QR ∩ (i + NQR) somehow.
# Actually, the SW construction uses the factorization of (p-1)/2:
# (11-1)/2 = 5, prime. So no 2-block factorization.
# At p=11, the "other" DRT might come from a completely different construction.

# Let's see what the complement diff set looks like
D_other_comp = set(range(1, 11)) - D_other
print(f"\nComplement of D_other: {sorted(D_other_comp)}")
# {4, 6, 7, 9, 10}
# Check: is complement also a difference set?
diff_counts = {}
for d1 in D_other_comp:
    for d2 in D_other_comp:
        if d1 != d2:
            diff = (d1 - d2) % 11
            diff_counts[diff] = diff_counts.get(diff, 0) + 1
is_comp_diff_set = all(diff_counts.get(i, 0) == 2 for i in range(1, 11))
print(f"Complement is also a (11,5,2)-diff set: {is_comp_diff_set}")

# Check if complement gives isomorphic tournament
T_comp = build_drt(11, D_other_comp)
H_comp = hamiltonian_path_count(T_comp)
print(f"H(complement tournament) = {H_comp}")
print(f"This should equal H(D_other) since T_comp = T_other^op and H(T)=H(T^op)")

print("\nDone.")
