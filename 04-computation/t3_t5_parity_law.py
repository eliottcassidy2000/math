#!/usr/bin/env python3
"""
Investigate the parity law: t_3 mod 4 determines t_5 parity at n=6.

Observation from exhaustive data:
  t_3=0: t_5 always even
  t_3=1: t_5 always even
  t_3=2: t_5 always even
  t_3=3: t_5 always ODD
  t_3=4: t_5 mixed
  t_3=5: t_5 always even
  t_3=6: t_5 mixed
  t_3=7: t_5 always ODD
  t_3=8: t_5 mixed

Pattern for t_3 where t_5 parity is FORCED:
  t_3=0: t_5 even, t_3 mod 4 = 0
  t_3=1: t_5 even, t_3 mod 4 = 1
  t_3=2: t_5 even, t_3 mod 4 = 2
  t_3=3: t_5 odd,  t_3 mod 4 = 3
  t_3=5: t_5 even, t_3 mod 4 = 1
  t_3=7: t_5 odd,  t_3 mod 4 = 3

So when t_3 is forced: t_3 ≡ 3 mod 4 → t_5 odd; otherwise → t_5 even.

Actually, the "forced" cases all have unique or very constrained score
sequences. Let me check:
  t_3=0: only transitive (one score seq)
  t_3=1: one score seq
  t_3=2: two score seqs (one with disjoint pair, one without)
  t_3=3: one score seq (must be regular 3,3,3,2,2,2 type)
  t_3=5: specific score seqs

The question: is there a GENERAL parity law connecting t_3 and t_5?

Actually, maybe I should look at this differently.
The key identity: alpha_1 = t_3 + t_5 (at n=6, since no 7-cycles).
And OCF: H = 1 + 2*alpha_1 + 4*alpha_2.

So H = 1 + 2*(t_3 + t_5) + 4*alpha_2.
H mod 4 = 1 + 2*(t_3 + t_5) mod 4 (since 4*alpha_2 ≡ 0 mod 4)
        = 1 + 2*(t_3 + t_5) mod 4

For H=21: 21 mod 4 = 1. So 2*(t_3+t_5) ≡ 0 mod 4, i.e., t_3+t_5 ≡ 0 mod 2.
This means t_3+t_5 must be even.

But the constraint is stronger: at n=6, when t_3=5 (odd), t_5 is always even,
so t_3+t_5 is always odd. This means H ≡ 1 + 2*(odd) = 3 mod 4.
So for t_3=5: H ≡ 3 mod 4, making H=21 (≡ 1 mod 4) impossible!

Wait, 21 mod 4 = 1. But t_3=5, t_5 even → t_3+t_5 odd → H ≡ 3 mod 4.
And 21 ≡ 1 mod 4. So H=21 with t_3=5 is impossible for PARITY REASONS!

But what about other t_3 values? For H=21 we need sw(1)+sw(2) = 20.
sw(1) = 2*t_3, so t_3 ≤ 10.

Let me check each t_3 case.

opus-2026-03-07-S39
"""
from collections import defaultdict

# From the exhaustive data at n=6:
# (t3, t5, p33) triples and their counts
data = {
    (0, 0, 0): 720,
    (1, 0, 0): 960,
    (2, 0, 0): 2160,
    (2, 0, 1): 80,
    (3, 1, 0): 2880,
    (4, 1, 0): 1440,
    (4, 2, 0): 1440,
    (4, 2, 1): 720,
    (4, 3, 0): 1920,
    (4, 4, 0): 720,
    (5, 2, 0): 288,
    (5, 4, 0): 1440,
    (5, 4, 1): 1440,
    (5, 6, 1): 480,
    (6, 4, 2): 720,
    (6, 5, 0): 1440,
    (6, 6, 0): 1440,
    (6, 6, 1): 2160,
    (6, 6, 2): 720,
    (6, 7, 1): 1440,
    (6, 8, 2): 720,
    (7, 7, 1): 1440,
    (7, 7, 2): 1440,
    (7, 9, 0): 480,
    (7, 9, 1): 1440,
    (8, 6, 4): 240,
    (8, 8, 2): 720,
    (8, 11, 1): 1440,
    (8, 12, 1): 240,
}

print("=== H=21 impossibility analysis at n=6 ===")
print("H = 1 + 2*t_3 + 2*t_5 + 4*p_33")
print("H = 21 requires: 2*t_3 + 2*t_5 + 4*p_33 = 20")
print("i.e., t_3 + t_5 + 2*p_33 = 10")
print()

for t3_target in range(11):
    # Find achievable t_5 values for this t_3
    t5_vals = set()
    p33_vals = set()
    for (t3, t5, p33), cnt in data.items():
        if t3 == t3_target:
            t5_vals.add(t5)
            p33_vals.add((t5, p33))

    if not t5_vals:
        print(f"  t_3={t3_target}: NOT ACHIEVABLE at n=6")
        continue

    # Check if any (t5, p33) satisfies t_3 + t_5 + 2*p_33 = 10
    solutions = []
    for (t5, p33) in p33_vals:
        if t3_target + t5 + 2*p33 == 10:
            solutions.append((t5, p33))

    near = [(t5, p33, t3_target + t5 + 2*p33) for (t5, p33) in p33_vals]
    near.sort(key=lambda x: abs(x[2] - 10))

    if solutions:
        print(f"  t_3={t3_target}: SOLUTION EXISTS! (t5, p33) = {solutions}")
    else:
        print(f"  t_3={t3_target}: no solution. achievable (t5,p33,target): {near[:3]}")
        # Explain why
        t5_parity = set(t5 % 2 for (t5, p33) in p33_vals)
        t3_plus_t5_mods = set((t3_target + t5) % 2 for t5 in t5_vals)
        need = 10 - t3_target
        print(f"    Need t_5 + 2*p_33 = {need}. t_5 parity: {'always even' if t5_parity == {0} else 'always odd' if t5_parity == {1} else 'mixed'}")
        if len(t5_parity) == 1:
            t5_p = list(t5_parity)[0]
            need_p = need % 2
            if t5_p != need_p:
                print(f"    PARITY OBSTRUCTION: t_5 is {'even' if t5_p==0 else 'odd'}, but need t_5 + 2*p_33 = {need} ({'even' if need_p==0 else 'odd'})")
                print(f"    Since 2*p_33 is always even, need t_5 ≡ {need_p} mod 2, but t_5 ≡ {t5_p} mod 2. IMPOSSIBLE!")

print()
print("=== SUMMARY ===")
print("H=21 requires t_3 + t_5 + 2*p_33 = 10.")
print()
print("For each achievable t_3:")
print("  t_3=0..2: target too large (need t_5+2p_33 >= 8, but t_5 is 0)")
print("  t_3=3: need t_5+2p_33=7. t_5=1 always → 1+2p_33=7 → p_33=3.")
print("         But p_33 max at t_3=3 is 0! (Only 3 3-cycles, need disjoint pair from 6 verts)")
print("  t_3=4: need t_5+2p_33=6. Achievable (t5,p33): (1,0),(2,0),(2,1),(3,0),(4,0)")
print("         t_5+2*p_33 values: 1,2,4,3,4. Max is 4, never reaches 6!")
print("  t_3=5: need t_5+2p_33=5. t_5 always EVEN → t_5+2p_33 always EVEN → can't be 5!")
print("         THIS IS THE PARITY OBSTRUCTION!")
print("  t_3=6: need t_5+2p_33=4. Achievable t_5+2p_33: 4,5,6,8,10,12. 4 IS achievable!")
print("         Wait... let me recheck.")

# Recheck t_3=6
print("\n  t_3=6 detail:")
for (t3, t5, p33), cnt in sorted(data.items()):
    if t3 == 6:
        target = t3 + t5 + 2*p33
        H = 1 + 2*t3 + 2*t5 + 4*p33
        print(f"    (t5={t5}, p33={p33}): target={target}, H={H}")

print()
print("  t_3=6, t_5=4, p_33=0: target=10? No, that's 6+4+0=10!")
print("  Wait, checking: (6, 4, 2): H=29, count=720")
print("  (6,4) exists but with p_33=2, giving target=6+4+4=14")
print("  (6,5,0): target=6+5+0=11")
print("  So t_3=6 never has target=10!")
print()
print("  The issue: when t_3=6, the MINIMUM t_5 achievable is 4 (not lower).")
print("  And with t_5=4, we'd need p_33=0. But (6,4,0) doesn't exist!")
print("  We only have (6,4,2): the t_5=4 comes WITH p_33=2.")
print()
print("This is NOT a parity obstruction for t_3=6.")
print("It's a VALUE obstruction: the constraints of tournament structure")
print("prevent t_5=4 with p_33=0 when t_3=6.")
