#!/usr/bin/env python3
"""
lrc_round_extension_s577.py   monad-researcher-2026-06-02-S577

Three tasks:
  1. Extend A000016 round-tournament verification to m=14..17, extending
     oracle S575's table (which stopped at m=13).
  2. Verify SC-round formula: SC-round(m) = 2^floor((m-1)/2) for m=3..17.
     Oracle S575 confirmed this through m=13. Extends to m=17.
  3. Test HYP-1998(C): at AP boundary (wall) times for n=14, are runner
     tournaments self-converse (and/or round)?
     - Open times: runner tournament is round (HYP-1998(A), confirmed)
     - Boundary times: does the runner tournament stay SC?

HYP-2093: SC-round(m) = 2^floor((m-1)/2) for ALL m (not just m<=13).
"""

import importlib.util
from math import gcd
from fractions import Fraction
from pathlib import Path

ROOT = Path(__file__).resolve().parents[1]
S574 = ROOT / "04-computation" / "lrc_round_count_m89_s574.py"


# ---------- load the S574 module ----------
def load_s574():
    spec = importlib.util.spec_from_file_location("s574_rounds", S574)
    mod = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(mod)
    return mod


# ---------- A000016 closed form ----------
def totient(d):
    return sum(1 for k in range(1, d + 1) if gcd(k, d) == 1)


def A000016(m):
    return sum(totient(d) * 2 ** (m // d)
               for d in range(1, m + 1) if m % d == 0 and d % 2 == 1) // (2 * m)


def opposite_adj(adj):
    m = len(adj)
    return [[0 if i == j else adj[j][i] for j in range(m)] for i in range(m)]


def do_row(rounds, m):
    """Compute round count, SC count, merged count for one m."""
    reps = {}
    for d in rounds.valid_dvectors(m):
        adj = rounds.build_adj(d, m)
        key = rounds.canon(adj, m)
        reps.setdefault(key, adj)

    sc_count = 0
    merged = set()
    for key, adj in reps.items():
        ckey = rounds.canon(opposite_adj(adj), m)
        if key == ckey:
            sc_count += 1
        merged.add(min(key, ckey))

    return len(reps), sc_count, len(merged)


# ============================================================
# TASK 1 & 2: Extend table to m=14..17
# ============================================================
print("=" * 70)
print("TASKS 1 & 2: Extend A000016 / SC-round / merged table to m=17")
print("=" * 70)

rounds = load_s574()

# Reproduce oracle S575 table first, then extend
# Oracle table: m=3..13
oracle_rows = [
    (3, 4, 2, 2, 2, 2, 2),
    (4, 8, 2, 2, 2, 2, 2),
    (5, 16, 4, 4, 4, 4, 4),
    (6, 32, 6, 6, 4, 4, 5),
    (7, 64, 10, 10, 8, 8, 9),
    (8, 128, 16, 16, 8, 8, 12),
    (9, 256, 30, 30, 16, 16, 23),
    (10, 512, 52, 52, 16, 16, 34),
    (11, 1024, 94, 94, 32, 32, 63),
    (12, 2048, 172, 172, 32, 32, 102),
    (13, 4096, 316, 316, 64, 64, 190),
]

import time

print(f"\n{'m':>3}  {'valid-d':>8}  {'round':>6}  {'A000016':>8}  "
      f"{'SC-round':>9}  {'2^fl':>5}  {'merged':>7}  {'match':>6}  {'t(s)':>6}")
print("-" * 75)

all_ok = True
rows_all = []

for m in range(3, 18):
    t0 = time.time()
    round_count, sc_count, merged_count = do_row(rounds, m)
    elapsed = time.time() - t0

    a16 = A000016(m)
    sc_formula = 2 ** ((m - 1) // 2)
    valid_d = 2 ** (m - 1)
    match = (round_count == a16) and (sc_count == sc_formula)
    if not match:
        all_ok = False

    print(f"{m:3d}  {valid_d:8d}  {round_count:6d}  {a16:8d}  "
          f"{sc_count:9d}  {sc_formula:5d}  {merged_count:7d}  "
          f"{'OK' if match else 'FAIL':>6}  {elapsed:6.2f}")
    rows_all.append((m, valid_d, round_count, a16, sc_count, sc_formula, merged_count, match))

print()
print(f"round == A000016 for ALL m=3..17: {all(r[7] for r in rows_all)}")
print(f"SC-round == 2^floor((m-1)/2) for ALL m=3..17: {all(r[7] for r in rows_all)}")

# Extract the new sequences
merged_seq = [r[6] for r in rows_all]
sc_seq = [r[4] for r in rows_all]
round_seq = [r[2] for r in rows_all]

print(f"\nround(m=3..17): {round_seq}")
print(f"SC-round(m=3..17): {sc_seq}")
print(f"merged(m=3..17): {merged_seq}")


# ============================================================
# TASK 3: HYP-1998(C) — boundary runner tournaments at n=14
# ============================================================
print("\n" + "=" * 70)
print("TASK 3: HYP-1998(C) — runner tournaments at AP boundary times, n=14")
print("=" * 70)
print("""
AP speed set: {1,2,...,13}, observer speed 0, threshold 1/14.
Half-turn runner tournament at time t: i→j if (v_j-v_i)*t mod 1 ∈ (0,1/2).
Wall times: t = k/(v_i + v_j) where k is an integer making ||(v_i+v_j)*k/(v_i+v_j)||=|k| small,
or more precisely: t such that some runner i is at ||v_i*t|| = 1/14.
For AP: wall times include t = k/(v_i * 14) for integer k, also pinch clocks t = j/(v_i+v_j).
""")

n = 14
m = n - 1  # 13 runners
AP = list(range(1, n))  # {1,...,13}

# Build the round class database for m=13
print(f"Building round class database for m={m}...")
reps13 = {}
for d in rounds.valid_dvectors(m):
    adj = rounds.build_adj(d, m)
    key = rounds.canon(adj, m)
    reps13.setdefault(key, adj)

sc_classes = set()
all_keys = set(reps13.keys())
for key, adj in reps13.items():
    ckey = rounds.canon(opposite_adj(adj), m)
    if key == ckey:
        sc_classes.add(key)

print(f"  {len(reps13)} round classes; {len(sc_classes)} are SC")


def half_turn_adj(speeds, t_num, t_den):
    """
    Build tournament adj matrix using half-turn comparator at time t = t_num/t_den.
    i→j if (v_j - v_i) * t mod 1 ∈ (0, 1/2).
    Use exact fraction arithmetic.
    """
    k = len(speeds)
    adj = [[0] * k for _ in range(k)]
    for i in range(k):
        for j in range(k):
            if i == j:
                continue
            # (v_j - v_i) * t_num / t_den mod 1
            diff = (speeds[j] - speeds[i]) * t_num
            # diff / t_den mod 1 ∈ (0, 1/2)?
            val = diff % t_den
            # val/t_den ∈ (0, 1/2) <=> val ∈ (0, t_den/2) <=> 0 < val < t_den/2
            if 0 < val < t_den / 2:
                adj[i][j] = 1
            elif val == 0 or val == t_den / 2:
                adj[i][j] = -1  # tie (exactly 0 or 1/2 turn)
    return adj


def is_tournament(adj):
    k = len(adj)
    for i in range(k):
        for j in range(i + 1, k):
            if adj[i][j] + adj[j][i] != 1:
                return False
    return True


def classify_adj(adj, reps_db, sc_set):
    """
    Classify adj (0/1 matrix, tournament) into:
    - is_round: is it in the round database?
    - is_sc: is the round class self-converse?
    - canon_key: the canonical key (None if not round)
    """
    k = len(adj)
    key = rounds.canon(adj, k)
    if key in reps_db:
        is_round = True
        is_sc = (key in sc_set)
    else:
        is_round = False
        is_sc = False
        # Check if SC even if not round
        opp = [[0 if i == j else adj[j][i] for j in range(k)] for i in range(k)]
        okey = rounds.canon(opp, k)
        if key == okey:
            is_sc = True
    return is_round, is_sc, key


# Test 1: Open time near t=0 (small t, all runners near 0)
# At t = 1/1000 (very small), check the runner tournament
print("\nTest 1: Open time t = 1/(14*1000) (very small, all runners near 0)")
t_num, t_den = 1, 14 * 1000
adj_open = half_turn_adj(AP, t_num, t_den)
if is_tournament(adj_open):
    ir, is_sc, key = classify_adj(adj_open, reps13, sc_classes)
    print(f"  is_round={ir}, is_sc={is_sc}")
else:
    # Count ties
    ties = sum(1 for i in range(m) for j in range(m) if adj_open[i][j] == -1)
    print(f"  NOT a tournament (has {ties} tied pairs)")

# Test 2: n-clock wall time t = 1/14
print("\nTest 2: n-clock wall time t = 1/14 (AP all-at-threshold)")
t_num, t_den = 1, 14
adj_wall = half_turn_adj(AP, t_num, t_den)
ties = sum(1 for i in range(m) for j in range(i + 1, m)
           if adj_wall[i][j] == -1 or adj_wall[j][i] == -1)
print(f"  Tied pairs at t=1/14: {ties}")
if is_tournament(adj_wall):
    ir, is_sc, key = classify_adj(adj_wall, reps13, sc_classes)
    print(f"  is_round={ir}, is_sc={is_sc}")

# The wall t=1/14 has ALL runners at threshold for AP (v_i/14 has ||v_i/14||=1/14 for v_i=1)
# Actually for v_1=1: v_1*1/14 = 1/14, ||1/14|| = 1/14 = threshold (tied)
# But for v_2=2: ||2/14|| = 2/14 > 1/14 (safe)
# So only v_1=1 and v_13=13 are at the threshold at t=1/14
# Let me recheck

print(f"\n  Checking ||v_i * 1/14|| for AP:")
for v in AP:
    pos = (v * 1) % 14  # numerator of v/14
    dist = min(pos, 14 - pos)  # distance to 0 (min(v/14, 1-v/14)*14)
    at_threshold = (dist == 1)
    print(f"    v={v:2d}: ||v/14||={dist}/14, at_threshold={at_threshold}")

# Test 3: Pinch clock t = 1/(v_i + v_j) for AP pair (1,13) - sum=14, so 1/(14)
# Already tested above. Try (1,6): sum=7, t=1/7
print("\nTest 3: Pinch clock t = 1/(v1+v6) = 1/7 (AP pair 1+6=7)")
t_num, t_den = 1, 7
adj_p1 = half_turn_adj(AP, t_num, t_den)
ties = sum(1 for i in range(m) for j in range(i + 1, m)
           if adj_p1[i][j] == -1 or adj_p1[j][i] == -1)
if ties > 0:
    print(f"  Tied pairs: {ties} — NOT a tournament (boundary state)")
    # Try t slightly before: t = (7-1)/(7*7) = slightly less than 1/7
    t_num2, t_den2 = 7 * 7 - 1, 7 * 7  # = 48/49 ≈ 0.9796... no that's > 1
    # t = 1/7 - 1/(7*100) = (100-1)/(7*100) = 99/700
    t_num2, t_den2 = 99, 700
    adj_before = half_turn_adj(AP, t_num2, t_den2)
    if is_tournament(adj_before):
        ir, is_sc, key = classify_adj(adj_before, reps13, sc_classes)
        print(f"  t just before 1/7: is_round={ir}, is_sc={is_sc}")
    t_num3, t_den3 = 101, 700
    adj_after = half_turn_adj(AP, t_num3, t_den3)
    if is_tournament(adj_after):
        ir, is_sc, key = classify_adj(adj_after, reps13, sc_classes)
        print(f"  t just after 1/7: is_round={ir}, is_sc={is_sc}")
else:
    ir, is_sc, key = classify_adj(adj_p1, reps13, sc_classes)
    print(f"  is_tournament=True, is_round={ir}, is_sc={is_sc}")

# Test 4: Wall time t = 1/(2*1) = 1/2 (runner 1 at 1/2 mark)
print("\nTest 4: Wall time t = 1/2 (runner 1 at position 1/2)")
t_num, t_den = 1, 2
adj_half = half_turn_adj(AP, t_num, t_den)
ties = sum(1 for i in range(m) for j in range(i + 1, m)
           if adj_half[i][j] == -1 or adj_half[j][i] == -1)
if ties > 0:
    print(f"  Tied pairs: {ties}")
else:
    ir, is_sc, key = classify_adj(adj_half, reps13, sc_classes)
    print(f"  is_tournament=True, is_round={ir}, is_sc={is_sc}")

# Systematic test: all wall times of the form t = a/b with b <= 30 (small denominators)
print("\nTest 5: Systematic — all open times t=a/b with b<=40, check round/SC class")
from fractions import Fraction

seen_classes = {}  # key -> (t_val, is_round, is_sc)
seen_nonround = []

for b in range(2, 41):
    for a in range(1, b):
        if gcd(a, b) > 1:
            continue
        t_val = Fraction(a, b)
        # Check if this is actually an open time (no ties) for AP at n=14
        adj = half_turn_adj(AP, a, b)
        ties = sum(1 for i in range(m) for j in range(i + 1, m)
                   if adj[i][j] == -1 or adj[j][i] == -1)
        if ties > 0:
            continue  # Wall time (has ties), skip
        if not is_tournament(adj):
            continue
        ir, is_sc, key = classify_adj(adj, reps13, sc_classes)
        if key not in seen_classes:
            seen_classes[key] = (t_val, ir, is_sc)
        if not ir:
            seen_nonround.append((t_val, ir, is_sc))

print(f"\n  Distinct open-time round classes seen (t=a/b, b<=40): {sum(1 for k,v in seen_classes.items() if v[1])}")
print(f"  Non-round classes seen: {len(seen_nonround)}")
if seen_nonround:
    print(f"  Non-round examples: {seen_nonround[:5]}")
else:
    print("  All open-time (tie-free) AP runner tournaments are ROUND. ✓")

print(f"\n  SC classes seen: {sum(1 for k,v in seen_classes.items() if v[2])}")
print(f"  Non-SC classes seen: {sum(1 for k,v in seen_classes.items() if not v[2])}")

# Also check: at wall times, are they SC?
print("\nTest 6: Wall times (tied pairs) for AP at n=14 — check runner tournament SC-ness")
wall_sc = []
wall_nonsc = []

for b in range(2, 41):
    for a in range(1, b):
        if gcd(a, b) > 1:
            continue
        adj_raw = half_turn_adj(AP, a, b)
        ties = [(i, j) for i in range(m) for j in range(i + 1, m)
                if adj_raw[i][j] == -1 or adj_raw[j][i] == -1]
        if not ties:
            continue

        # Build tournament by resolving ties (use both 0 and 1)
        # Actually: at boundary, the tournament may be undefined for tied pairs.
        # Check if removing tied pairs gives a partial tournament, then
        # check if the "definite" arcs could come from a self-converse extension.
        # More directly: check if a//b gives a self-converse partial tournament
        # by checking if swapping tied pairs still gives a consistent partial
        # tournament.

        # Simpler test: build the tournament with ties broken by circular order
        # (the "limit" tournament as t approaches a/b from left vs right)
        adj_before = half_turn_adj(AP, a * 1000 - 1, b * 1000)
        adj_after = half_turn_adj(AP, a * 1000 + 1, b * 1000)

        ok_before = is_tournament(adj_before)
        ok_after = is_tournament(adj_after)

        if ok_before and ok_after:
            ir_b, sc_b, key_b = classify_adj(adj_before, reps13, sc_classes)
            ir_a, sc_a, key_a = classify_adj(adj_after, reps13, sc_classes)
            wall_info = (Fraction(a, b), len(ties), ir_b, sc_b, key_b == key_a, ir_a, sc_a)
            if sc_b and sc_a:
                wall_sc.append(wall_info)
            else:
                wall_nonsc.append(wall_info)

print(f"\n  Wall times where both sides are SC: {len(wall_sc)}")
print(f"  Wall times where at least one side is non-SC: {len(wall_nonsc)}")
if wall_nonsc:
    print("  Non-SC wall examples:")
    for w in wall_nonsc[:5]:
        print(f"    t={w[0]}: {w[1]} ties, before:round={w[2]},sc={w[3]}, same_class={w[4]}, after:round={w[5]},sc={w[6]}")

# Summary
print("\n" + "=" * 70)
print("SUMMARY")
print("=" * 70)
print(f"\nHYP-2093: SC-round(m) = 2^floor((m-1)/2) for m=3..17")
print(f"  Verified: {all_ok}")
print(f"\nHYP-1998(B): Which round classes arise at open AP times (n=14)?")
print(f"  Distinct classes seen (t=a/b, b<=40): {len(seen_classes)}")
print(f"  All are round: {len(seen_nonround) == 0}")
print(f"\nHYP-1998(C): Boundary runner tournaments — are they SC?")
print(f"  Wall times with SC on both sides: {len(wall_sc)}")
print(f"  Wall times with non-SC on a side: {len(wall_nonsc)}")
