#!/usr/bin/env python3
"""
INV-029: Ballot sequence / Dyck path connection to C(L-2, 2k-1).

THM-007 states: For an L-cycle consecutively embedded in a Ham path P',
the number of internal signature patterns giving exactly k Type-II positions
within the embedding window is C(L-2, 2k-1).

This script:
1. Verifies C(L-2, 2k-1) computationally
2. Constructs an explicit bijection to ballot sequences
3. Tests whether the bijection is natural/canonical

A Type-II position at index j means s[j]=1, s[j+1]=0 (a "descent" in the
signature). The signature within an L-cycle window has:
  s[0] = 1 (first vertex is beaten by v, since v -> first vertex of cycle)
  s[L-2] = 0 (last vertex beats v, since last -> v in cycle)
  s[1],...,s[L-3] are free binary values

So we have L-2 signature values with s[0]=1 forced, s[L-2]=0 forced,
and L-4 free values in between (for L >= 5).

Type-II count = #{j : s[j]=1, s[j+1]=0} for j=0,...,L-3.

Claim: #{patterns with exactly k Type-II} = C(L-2, 2k-1).

Instance: opus-2026-03-05-S4b
"""
from math import comb
from itertools import product

def count_type_ii(sig):
    """Count Type-II positions (descents 1->0) in signature."""
    return sum(1 for j in range(len(sig)-1) if sig[j] == 1 and sig[j+1] == 0)

print("=" * 60)
print("VERIFICATION OF C(L-2, 2k-1) DISTRIBUTION")
print("=" * 60)

for L in range(3, 14, 2):  # odd cycle lengths
    n_free = L - 4  # number of free signature values
    if n_free < 0:
        # L=3: signature is (1, 0), exactly 1 Type-II, C(1,1)=1
        print(f"\nL={L}: sig=(1,0), TypeII=1, C({L-2},1)={comb(L-2,1)} -- {'MATCH' if comb(L-2,1)==1 else 'FAIL'}")
        continue

    # Enumerate all 2^{L-4} patterns
    counts = {}  # k -> count
    total = 0
    for bits in product([0, 1], repeat=n_free):
        sig = (1,) + bits + (0,)
        k = count_type_ii(sig)
        counts[k] = counts.get(k, 0) + 1
        total += 1

    print(f"\nL={L}: {total} patterns, {L-2} signature values")
    all_match = True
    for k in sorted(counts.keys()):
        expected = comb(L-2, 2*k-1) if 2*k-1 >= 0 else 0
        match = counts[k] == expected
        if not match:
            all_match = False
        print(f"  k={k}: count={counts[k]}, C({L-2},{2*k-1})={expected} {'MATCH' if match else 'FAIL'}")
    print(f"  Overall: {'ALL MATCH' if all_match else 'SOME FAIL'}")

print("\n" + "=" * 60)
print("BIJECTIVE PROOF VIA BALLOT SEQUENCES")
print("=" * 60)
print("""
Proof sketch for C(L-2, 2k-1):

Given: signature s = (s_0, s_1, ..., s_{L-3}) with s_0=1, s_{L-3}=0.
       (L-2 values total, L-4 free)

Define the "descent word" w_j = 1 if s_j != s_{j+1}, 0 if s_j == s_{j+1}
for j = 0, ..., L-4.

Then #TypeII = #descents of type (1,0).
But #descents(1->0) = #rises(0->1) + (s_0 - s_{L-3})/... no, that's not right.

Alternative: think of s as a lattice path.
  - s_j = 1 means "up" (v beats this vertex)
  - s_j = 0 means "down" (this vertex beats v)
  - Start at s_0=1 (up), end at s_{L-3}=0 (down)
  - A Type-II at position j is a step from 1 to 0 = "peak"

So #TypeII = #peaks in the binary string s.
Peaks in a binary string of length m starting with 1, ending with 0:
  s = 1 x_1 x_2 ... x_{m-2} 0
  Peak at position j: s_j=1, s_{j+1}=0, 0 <= j <= m-2

This is a well-known problem! The number of binary strings of length m
starting with 1, ending with 0, with exactly k peaks is C(m-1, 2k-1).
Wait, our m = L-2, so this gives C(L-3, 2k-1)?

Let me recheck...
""")

# Let's verify by a different encoding
print("Recheck: encoding transitions instead of values")
print()

for L in [5, 7, 9, 11]:
    m = L - 2  # signature length
    n_free = L - 4
    if n_free < 0:
        continue

    # Count by transitions: define t_j = |s_j - s_{j+1}| for j=0..m-2
    # t_j = 1 means transition, t_j = 0 means same
    # We have m-1 = L-3 transitions
    # Starting from s_0=1, the signature is determined by the transitions.
    # #TypeII = #(j where s_j=1, s_{j+1}=0) = #falling transitions from 1

    # Alternative: map to choosing positions of the 2k-1 transitions
    # From s_0=1, to get k peaks we need exactly 2k-1 transitions:
    #   1->0 (peak), 0->1 (valley), 1->0, 0->1, ..., ending with 1->0
    # since s_0=1 and s_{m-1}=0, we need an odd number of transitions
    # and specifically the pattern must start with a fall (1->0).
    #
    # Number of transitions = 2k-1 (must be odd since start=1, end=0).
    # Choose which of the m-1 = L-3 transition positions are "active":
    # C(L-3, 2k-1) transitions.
    #
    # But the claim is C(L-2, 2k-1), not C(L-3, 2k-1)!

    # Let me recount carefully at L=5
    # m = L-2 = 3, signature = (1, x, 0), x in {0,1}
    # x=0: sig=(1,0,0), TypeII at j=0 only (1->0). k=1.
    # x=1: sig=(1,1,0), TypeII at j=1 only (1->0). k=1.
    # So all 2 patterns have k=1. C(3,1)=3, C(2,1)=2.
    # Actual count for k=1 is 2. So it's C(L-3, 2k-1) = C(2,1) = 2. MATCH.
    # But C(L-2, 2k-1) = C(3,1) = 3. MISMATCH.

    counts = {}
    for bits in product([0, 1], repeat=n_free):
        sig = (1,) + bits + (0,)
        k = count_type_ii(sig)
        counts[k] = counts.get(k, 0) + 1

    print(f"L={L}, m={m}, #transitions={m-1}:")
    for k in sorted(counts.keys()):
        c_m1 = comb(m-1, 2*k-1) if 2*k-1 >= 0 else 0
        c_m = comb(m, 2*k-1) if 2*k-1 >= 0 else 0
        print(f"  k={k}: actual={counts[k]}, C(m-1,2k-1)=C({m-1},{2*k-1})={c_m1}, C(m,2k-1)=C({m},{2*k-1})={c_m}")

print("\n" + "=" * 60)
print("CORRECTED BIJECTION: TRANSITION ENCODING")
print("=" * 60)
print("""
The correct formula is:

  #{signatures with exactly k Type-II positions} = C(L-3, 2k-1)

where L is the cycle length, and the signature has L-2 values with
s_0=1 (forced), s_{L-3}=0 (forced), and L-4 free values.

PROOF (bijective):
The signature s = (1, x_1, ..., x_{L-4}, 0) has L-3 consecutive pairs.
Define transition indicators t_j = 1 if s_j != s_{j+1} for j=0,...,L-4.
This gives L-3 binary indicators.

Since s_0=1 and s_{L-3}=0, the total number of transitions must be ODD
(we must cross from 1 to 0 an odd net number of times).

The Type-II count k = #(j: s_j=1, s_{j+1}=0) equals the number of
FALLING transitions. Since transitions alternate fall-rise-fall-...-fall
(starting and ending with a fall because s_0=1, s_{L-3}=0), we have:
  #falls = (#transitions + 1) / 2
  #rises = (#transitions - 1) / 2

So k = (total_transitions + 1) / 2, meaning total_transitions = 2k - 1.

Choosing which 2k-1 of the L-3 transition positions are active gives
C(L-3, 2k-1) configurations.

QED.

NOTE: The original THM-007 states C(L-2, 2k-1). Let me check if the
indexing convention differs (L-2 vs L-3 depends on whether L counts
the full cycle length or the number of non-v vertices).
""")

# Final verification with BOTH conventions
print("Convention check:")
for L in [3, 5, 7, 9]:
    m_vertices = L - 1  # non-v vertices in the cycle
    m_sig = L - 2  # signature values (arcs between v and non-v vertices... wait)

    # Actually: in a cycle v -> a_1 -> a_2 -> ... -> a_{L-1} -> v
    # The signature s_j = T[v][a_j] for j = 1, ..., L-1
    # s_1 = 1 (v -> a_1), s_{L-1} = 0 (a_{L-1} -> v)
    # Free values: s_2, ..., s_{L-2}, which is L-3 values
    # Total signature length: L-1 values
    # Number of consecutive pairs: L-2

    # Wait, the TANGENT T001 says C(L-2, 2k-1), and the MASTER_FINDINGS says
    # "within the window s[j..j+L-2]" with "s[j]=1 and s[j+L-2]=0"
    # That's L-1 values in the window, L-2 transitions.
    # The window has L-1 vertices (the non-v vertices of the L-cycle)
    # and L-1 signature values.

    # Let me recount for L=5:
    # Cycle: v -> a -> b -> c -> d -> v (5 vertices including v)
    # Non-v vertices: a, b, c, d (4 = L-1 vertices)
    # Signature: s(a)=1, s(b)=?, s(c)=?, s(d)=0
    # Free: s(b), s(c) -- that's L-3 = 2 free values
    # Pairs: (s(a),s(b)), (s(b),s(c)), (s(c),s(d)) -- 3 = L-2 pairs
    # TypeII = #(pairs with 1->0)

    # For L=5: 4 patterns (2 free bits)
    # (1,0,0,0): peaks at j=0. k=1.
    # (1,0,1,0): peaks at j=0,2. k=2.
    # (1,1,0,0): peaks at j=1. k=1.
    # (1,1,1,0): peaks at j=2. k=1.
    # Counts: k=1: 3, k=2: 1
    # C(L-2, 2k-1): C(3,1)=3, C(3,3)=1. MATCH!

    if L == 3:
        # Cycle v -> a -> b -> v. Signature: s(a)=1, s(b)=0. 0 free values.
        # 1 pattern, k=1. C(1,1)=1. MATCH.
        print(f"  L={L}: 1 pattern, k=1, C({L-2},1)={comb(L-2,1)} MATCH")
        continue

    n_free = L - 3
    counts = {}
    for bits in product([0, 1], repeat=n_free):
        sig = (1,) + bits + (0,)
        k = count_type_ii(sig)
        counts[k] = counts.get(k, 0) + 1

    print(f"  L={L}: {2**n_free} patterns, sig length {L-1}")
    all_match = True
    for k in sorted(counts.keys()):
        expected = comb(L-2, 2*k-1) if 2*k-1 >= 0 else 0
        match = counts[k] == expected
        if not match: all_match = False
        print(f"    k={k}: count={counts[k]}, C({L-2},{2*k-1})={expected} {'OK' if match else 'FAIL'}")

print("""
RESOLUTION: THM-007 is CORRECT with C(L-2, 2k-1).

The L-cycle through v has L-1 non-v vertices, giving L-1 signature values.
Of these, s_1=1 and s_{L-1}=0 are forced, leaving L-3 free values.
There are L-2 consecutive pairs, and the number of pairs with
exactly 2k-1 transitions among L-2 pairs gives C(L-2, 2k-1).

BIJECTIVE PROOF:
The L-2 consecutive pairs form a binary string t of length L-2 where
t_j = 1 iff (s_j, s_{j+1}) is a transition (s_j != s_{j+1}).
Since s starts at 1 and ends at 0, the number of 1's in t must be odd.
A Type-II position (1->0 descent) occurs precisely at falling transitions.
Since the transitions alternate fall-rise-fall-..., having exactly k
Type-II positions requires exactly 2k-1 transitions total.

BUT the free bits are s_2, ..., s_{L-2} (L-3 of them), while the
transition string t has L-2 entries. These are NOT independent —
t_j is determined by s_j and s_{j+1}.

The correct argument: given s_0=1, choosing the transition string t
(which L-2 entries are transitions?) determines s uniquely. But t must
satisfy the constraint that the running parity of transitions equals
the forced endpoint. Since s_0=1 and s_{L-1}=0, we need sum(t) to be
odd. Among all C(L-2, m) choices of m transition positions, only those
with m odd are valid, and m = 2k-1 gives k Type-II positions.

The number of valid choices is C(L-2, 2k-1). QED.
""")
