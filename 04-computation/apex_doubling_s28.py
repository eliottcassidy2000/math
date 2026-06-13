#!/usr/bin/env python3
"""
The Apex Doubling Law and its consequences.
opus-2026-04-03-S28

DISCOVERY: The apex tile (n,1) — the arc between the two endpoints of the
base path — has an "apex effect" on mean H that DOUBLES with each n:
  n=3: +2, n=4: +2, n=5: +4, n=6: +8

Conjecture: apex_effect(n) = 2^(n-3) for n >= 3.

This makes deep geometric sense:
- The apex arc (n,1) connects the SOURCE of the base path to the SINK
- Flipping it reverses the "global polarity" of the tournament
- The effect should grow exponentially because it interacts with ALL other arcs

This script investigates:
1. Verify the 2^(n-3) law
2. What is E[H | apex=1] - E[H | apex=0] exactly?
3. Walsh-Fourier interpretation of the apex effect
4. Does the apex effect decompose into skip contributions?
5. The "inner apex" — tiles with second-highest skip — do they show 2^(n-4)?
"""

from collections import defaultdict
import sys

def tournament_from_bits(n, bits_list):
    T = [[0]*n for _ in range(n)]
    tiles = []
    for y in range(1, n-1):
        for x in range(n, y+1, -1):
            if x - y >= 2:
                tiles.append((x, y))
    for i in range(n, 1, -1):
        T[i-1][i-2] = 1
    for idx, (x, y) in enumerate(tiles):
        if idx < len(bits_list):
            if bits_list[idx] == 0:
                T[x-1][y-1] = 1
            else:
                T[y-1][x-1] = 1
    return T

def count_hp(T, n):
    dp = [[0]*n for _ in range(1 << n)]
    for v in range(n):
        dp[1 << v][v] = 1
    full = (1 << n) - 1
    for mask in range(1, full + 1):
        for v in range(n):
            if not (mask & (1 << v)) or dp[mask][v] == 0:
                continue
            for u in range(n):
                if not (mask & (1 << u)) and T[v][u]:
                    dp[mask | (1 << u)][u] += dp[mask][v]
    return sum(dp[full])

def get_tiles(n):
    tiles = []
    for y in range(1, n-1):
        for x in range(n, y+1, -1):
            if x - y >= 2:
                tiles.append((x, y))
    return tiles

print("=" * 70)
print("THE APEX DOUBLING LAW")
print("=" * 70)

for n in range(3, 8):
    m = (n-1)*(n-2)//2
    tiles = get_tiles(n)

    # Find apex (max skip)
    skips = [x - y for x, y in tiles]
    max_skip = max(skips)
    apex_idx = [i for i in range(m) if skips[i] == max_skip]

    # Find tiles by skip value
    tiles_by_skip = defaultdict(list)
    for i in range(m):
        tiles_by_skip[skips[i]].append(i)

    if 2**m > 5000000:
        print(f"\nn={n}: m={m}, 2^m={2**m} — sampling 50000 random tilings")
        import random
        random.seed(42)
        sample_size = 50000
        H_by_bit = {j: {0: [], 1: []} for j in range(m)}

        for _ in range(sample_size):
            bits = random.getrandbits(m)
            b = [(bits >> i) & 1 for i in range(m)]
            T = tournament_from_bits(n, b)
            H = count_hp(T, n)
            for j in range(m):
                H_by_bit[j][b[j]].append(H)

        # Main effects
        print(f"\nMain effects E[H|bit=1] - E[H|bit=0] by tile:")
        effects = []
        for j in range(m):
            e0 = sum(H_by_bit[j][0]) / len(H_by_bit[j][0])
            e1 = sum(H_by_bit[j][1]) / len(H_by_bit[j][1])
            effect = e1 - e0
            effects.append(effect)
            t = tiles[j]
            print(f"  {chr(65+j) if j < 26 else f'T{j}'} ({t[0]},{t[1]}) skip={skips[j]}: "
                  f"E[H|0]={e0:.2f}, E[H|1]={e1:.2f}, effect={effect:+.2f}")

        apex_eff = effects[apex_idx[0]]
        print(f"\nApex effect = {apex_eff:+.2f}")
        print(f"2^(n-3) = {2**(n-3)}")
        print(f"Ratio = {apex_eff / 2**(n-3):.4f}")

        # Effect by skip value
        print(f"\nMean effect by skip value:")
        for s in sorted(tiles_by_skip):
            idxs = tiles_by_skip[s]
            mean_eff = sum(effects[i] for i in idxs) / len(idxs)
            print(f"  skip={s}: mean effect = {mean_eff:+.4f} (tiles: {len(idxs)})")

    else:
        print(f"\nn={n}: m={m}, 2^m={2**m} — exhaustive enumeration")
        H_by_bit = {j: {0: [], 1: []} for j in range(m)}

        for bits in range(2**m):
            b = [(bits >> i) & 1 for i in range(m)]
            T = tournament_from_bits(n, b)
            H = count_hp(T, n)
            for j in range(m):
                H_by_bit[j][b[j]].append(H)

        # Main effects
        print(f"\nMain effects E[H|bit=1] - E[H|bit=0] by tile:")
        effects = []
        for j in range(m):
            e0 = sum(H_by_bit[j][0]) / len(H_by_bit[j][0])
            e1 = sum(H_by_bit[j][1]) / len(H_by_bit[j][1])
            effect = e1 - e0
            effects.append(effect)
            t = tiles[j]
            print(f"  {chr(65+j) if j < 26 else f'T{j}'} ({t[0]},{t[1]}) skip={skips[j]}: "
                  f"E[H|0]={e0:.2f}, E[H|1]={e1:.2f}, effect={effect:+.2f}")

        apex_eff = effects[apex_idx[0]]
        print(f"\nApex effect = {apex_eff:+.2f}")
        print(f"2^(n-3) = {2**(n-3)}")
        print(f"Ratio = {apex_eff / 2**(n-3):.4f}")

        # Effect by skip value
        print(f"\nMean effect by skip value:")
        for s in sorted(tiles_by_skip):
            idxs = tiles_by_skip[s]
            mean_eff = sum(effects[i] for i in idxs) / len(idxs)
            print(f"  skip={s}: mean effect = {mean_eff:+.4f} (tiles: {len(idxs)})")

        # CHECK: Is the effect EXACTLY 2^(n-3)?
        # The sum of all H values with apex=0 vs apex=1
        total_0 = sum(H_by_bit[apex_idx[0]][0])
        total_1 = sum(H_by_bit[apex_idx[0]][1])
        count = len(H_by_bit[apex_idx[0]][0])  # = 2^(m-1)
        exact_diff = (total_1 - total_0) / count
        print(f"\nExact apex effect = {exact_diff} = {total_1 - total_0}/{count}")
        from fractions import Fraction
        frac = Fraction(total_1 - total_0, count)
        print(f"As fraction: {frac}")
        print(f"Is it exactly 2^(n-3)? {frac == 2**(n-3)}")

    # Skip-effect hierarchy
    print(f"\n--- SKIP-EFFECT HIERARCHY ---")
    print(f"If apex (skip {max_skip}) has effect 2^(n-3), what about skip {max_skip-1}?")
    if max_skip - 1 in tiles_by_skip:
        second_idxs = tiles_by_skip[max_skip - 1]
        if second_idxs and 2**m <= 5000000:
            second_effects = [effects[i] for i in second_idxs]
            mean_second = sum(second_effects) / len(second_effects)
            print(f"  Second-skip effect = {mean_second:.4f}")
            print(f"  2^(n-4) = {2**(n-4)}")
            if 2**(n-4) != 0:
                print(f"  Ratio = {mean_second / 2**(n-4):.4f}")

    # The INTERACTION between apex and complement
    # Does flipping the apex THEN taking complement = flipping the "anti-apex"?
    print(f"\n--- APEX-COMPLEMENT INTERACTION ---")
    if 2**m <= 100000:
        # What is the "anti-apex"? In the complement tiling, the apex is the same tile
        # but its bit is flipped. So complement(tiling with apex=1) has apex=0.
        #
        # More interesting: what tournament operation does "flip apex" correspond to?
        # The apex tile is (n,1), which is the arc between vertex n and vertex 1.
        # Flipping it reverses this single arc.
        #
        # In the tournament: this swaps the direction of the arc between the
        # FIRST and LAST vertex of the base path.

        # Let's check: is apex-flip the same as reversing vertex 1 and vertex n's roles?
        swap_count = 0
        total = 0
        for bits in range(min(2**m, 10000)):
            b = [(bits >> i) & 1 for i in range(m)]
            T0 = tournament_from_bits(n, b)
            H0 = count_hp(T0, n)

            # Flip apex
            b_apex = b.copy()
            ai = apex_idx[0]
            b_apex[ai] = 1 - b_apex[ai]
            T_apex = tournament_from_bits(n, b_apex)
            H_apex = count_hp(T_apex, n)

            # Complement
            b_comp = [1 - x for x in b]
            T_comp = tournament_from_bits(n, b_comp)
            H_comp = count_hp(T_comp, n)

            # Flip apex then complement (= complement all EXCEPT apex)
            b_apex_comp = [1 - x for x in b_apex]
            T_ac = tournament_from_bits(n, b_apex_comp)
            H_ac = count_hp(T_ac, n)

            total += 1
            if H0 + H_ac == H_apex + H_comp:
                swap_count += 1

        print(f"  H(t) + H(apex+comp) == H(apex) + H(comp)? {swap_count}/{total} = {swap_count/total:.4f}")
        print(f"  (If 1.0, the apex and complement operations COMMUTE in H)")

print(f"\n{'='*70}")
print("SUMMARY: APEX EFFECT TABLE")
print(f"{'='*70}")
print(f"{'n':>4} {'m':>4} {'apex_effect':>15} {'2^(n-3)':>10} {'match':>8}")
