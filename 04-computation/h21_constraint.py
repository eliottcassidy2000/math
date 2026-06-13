#!/usr/bin/env python3
"""
Analyze the constraint on (i_1, i_2) pairs at n=6 and n=7.
WHY can't we get 2*i_1 + 4*i_2 = 20?

The achieved pairs at n=6 from exhaustive data:
  (0,0) -> H=1
  (1,0) -> H=3
  (2,0) -> H=5
  (2,1) -> H=9 or (4,0) -> H=9
  (5,0) -> H=11
  (6,0) -> H=13
  (7,0) -> H=15
  (6,1) -> H=17 or (8,0) -> H=17
  (9,0) -> H=19
  (9,1) -> H=23 or (11,0) -> H=23
  (12,0) -> H=25
  (11,1) -> H=27
  (10,2) -> H=29 or (12,1) -> H=29
  (13,1) -> H=31
  (12,2) -> H=33 or (14,1) -> H=33 or (16,0) -> H=33
  (14,2) -> H=37 or (16,1) -> H=37
  (16,2) -> H=41
  (19,1) -> H=43
  (14,4) -> H=45 or (20,1) -> H=45

Key observation: 2*i_1 + 4*i_2 = H-1 must be achievable.

H-1 = 20 is NOT achieved.
H-1 = 6 (for H=7) is NOT achieved.

Let's check: what values of 2*i_1 + 4*i_2 ARE achieved?
"""

# From the exhaustive n=6 data above:
achieved_pairs = [
    (0, 0),   # H=1
    (1, 0),   # H=3
    (2, 0),   # H=5
    (2, 1),   # H=9
    (4, 0),   # H=9
    (5, 0),   # H=11
    (6, 0),   # H=13
    (7, 0),   # H=15
    (6, 1),   # H=17
    (8, 0),   # H=17
    (9, 0),   # H=19
    (9, 1),   # H=23
    (11, 0),  # H=23
    (12, 0),  # H=25
    (11, 1),  # H=27
    (10, 2),  # H=29
    (12, 1),  # H=29
    (13, 1),  # H=31
    (12, 2),  # H=33
    (14, 1),  # H=33
    (16, 0),  # H=33
    (14, 2),  # H=37
    (16, 1),  # H=37
    (16, 2),  # H=41
    (19, 1),  # H=43
    (14, 4),  # H=45
    (20, 1),  # H=45
]

print("Achieved 2*i_1 + 4*i_2 values at n=6:")
S_vals = set()
for i1, i2 in achieved_pairs:
    s = 2*i1 + 4*i2
    S_vals.add(s)
    print(f"  i_1={i1:2d}, i_2={i2}: 2*i_1+4*i_2 = {s:2d} -> H = {s+1}")

print(f"\nAchieved S = 2*i_1+4*i_2 values: {sorted(S_vals)}")
print(f"All achieved H = S+1: {sorted(s+1 for s in S_vals)}")
missing_s = [s for s in range(0, max(S_vals)+1) if s not in S_vals and s % 2 == 0]
print(f"Missing even S values: {missing_s}")
print(f"Corresponding missing H: {[s+1 for s in missing_s]}")

# Key insight: S must be even (2*i_1 is always even, 4*i_2 is always even)
# S = 0, 2, 4, 8, 10, 12, 14, 16, 18, 22, 24, 26, 28, 30, 32, 36, 40, 42, 44
# Missing even S: 6, 20, 34, 38

print("\n" + "=" * 60)
print("KEY ANALYSIS: Why is S=6 (H=7) and S=20 (H=21) missing?")
print("=" * 60)

print("""
S = 6 requires:
  i_1=3, i_2=0: impossible (alpha_1=3 impossibility, THM-029)
  i_1=1, i_2=1: 1 cycle, 1 disjoint pair -- but 1 cycle means no pair -> contradiction

S = 20 requires (with i_3=0):
  i_1=10, i_2=0: 10 directed cycles, all pairwise conflicting
  i_1=8, i_2=1: 8 directed cycles, exactly 1 disjoint pair
  i_1=6, i_2=2: 6 directed cycles, exactly 2 disjoint pairs
  i_1=4, i_2=3: 4 directed cycles, exactly 3 disjoint pairs
  i_1=2, i_2=4: 2 directed cycles, exactly 4 disjoint pairs -- impossible (C(2,2)=1 pair max)
  i_1=0, i_2=5: impossible (no cycles means no pairs)

With i_3=1:
  i_1=6, i_2=0, i_3=1: 6 cycles, no pairs, but one TRIPLE of mutually disjoint cycles -> need 9+ vertices
  i_1=4, i_2=1, i_3=1: 4 cycles, 1 pair, 1 triple
  i_1=2, i_2=2, i_3=1: 2 cycles, 2 pairs -- impossible from 2 cycles
  i_1=0, i_2=3, i_3=1: impossible

With i_3=2:
  i_1=4, i_2=0, i_3=2: 4 cycles, 0 pairs, 2 triples -- contradiction (triple requires pair)
  i_1=2, i_2=0, i_3=2: impossible
  i_1=0, i_2=1, i_3=2: impossible

Key observations:
1. i_2 >= C(i_1, 2) would mean ALL pairs are disjoint (complete independence)
2. i_2 <= C(i_1, 2) always
3. The constraint is: given the conflict graph structure, what i_2 values are possible?
""")

print("=" * 60)
print("Achieved (i_1, i_2) at n=6 and their S values:")
print("=" * 60)
for i1, i2 in sorted(set(achieved_pairs)):
    s = 2*i1 + 4*i2
    # Max possible i_2 from i_1 cycles
    max_i2 = i1 * (i1-1) // 2
    print(f"  (i_1={i1:2d}, i_2={i2}): S={s:2d}, H={s+1:2d}, max_i2={max_i2}")

# Check: for each i_1, what i_2 values appear?
print("\n" + "=" * 60)
print("i_2 values for each i_1:")
print("=" * 60)
i1_to_i2 = {}
for i1, i2 in achieved_pairs:
    if i1 not in i1_to_i2:
        i1_to_i2[i1] = set()
    i1_to_i2[i1].add(i2)

for i1 in sorted(i1_to_i2.keys()):
    i2_vals = sorted(i1_to_i2[i1])
    max_i2 = i1 * (i1-1) // 2
    # What H values would each i_2 give?
    h_vals = [1 + 2*i1 + 4*i2 for i2 in i2_vals]
    print(f"  i_1={i1:2d}: i_2 in {i2_vals}, max_i2={max_i2}, H values={h_vals}")

    # Check: could i_2 give H=21?
    needed_i2 = (20 - 2*i1) / 4
    if needed_i2 == int(needed_i2) and needed_i2 >= 0:
        ni2 = int(needed_i2)
        achieved = ni2 in i1_to_i2[i1]
        print(f"    -> For H=21: need i_2={ni2}. {'ACHIEVED' if achieved else 'NOT achieved'}")
