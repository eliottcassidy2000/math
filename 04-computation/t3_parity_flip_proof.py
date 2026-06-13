#!/usr/bin/env python3
"""
WHY does t3 parity flip under tiling flip?

Tiling flip: complement all non-backbone edges (j-i >= 2).
Backbone: edges (i, i+1) for i = 0,...,n-2, always A[i][i+1] = 1.

A 3-cycle on vertices {a,b,c} (a < b < c) exists iff one of:
  A[a][b]*A[b][c]*A[c][a] = 1  (clockwise: a->b->c->a)
  A[a][c]*A[c][b]*A[b][a] = 1  (counter: a->c->b->a)

Exactly one of these holds for each triple {a,b,c}.

Types of triples based on backbone edges:
1. Consecutive: {i, i+1, i+2} — TWO backbone edges (i->i+1 and i+1->i+2 fixed)
   Only non-backbone edge: (i, i+2). Flip changes it.
2. One backbone: {i, i+1, j} with j != i+2 — ONE backbone edge
   Two non-backbone edges. Flip changes both.
3. No backbone: {a, b, c} with no consecutive pair — ZERO backbone edges
   Three non-backbone edges. Flip changes all three.

For type 1 (consecutive): backbone gives i->i+1->i+2.
  3-cycle exists iff A[i+2][i] = 1 (completing the cycle i->i+1->i+2->i).
  After flip: A'[i+2][i] = 1 - A[i+2][i] = A[i][i+2].
  So the 3-cycle on {i,i+1,i+2} EXISTS before flip iff it DOESN'T after flip.
  Each consecutive triple contributes EXACTLY 1 to (t3_before + t3_after).
  There are n-2 consecutive triples.

For type 3 (no backbone): all three edges flip.
  Original: one of A[a][b]*A[b][c]*A[c][a] or A[a][c]*A[c][b]*A[b][a] = 1.
  After flip: A'[x][y] = 1 - A[x][y] for all pairs.
  Clockwise after flip: (1-A[a][b])*(1-A[b][c])*(1-A[c][a])
    This is the COUNTER-clockwise original! (Since if A[a][b]=0 then A[b][a]=1)
  Wait, A'[a][b] = 1 - A[a][b] means the edge REVERSES. So:
  A'[a][b]*A'[b][c]*A'[c][a] = A[b][a]*A[c][b]*A[a][c] = counter-clockwise original

  So for type 3: the clockwise 3-cycle after flip = counter-clockwise before flip.
  Total 3-cycles on this triple: same as before! (1 before, 1 after)
  Type 3 triples are INVARIANT under flip.

For type 2 (one backbone edge): say backbone edge is (i, i+1).
  The triple is {i, i+1, j} where j differs from i by >= 2 AND j differs from i+1 by >= 2.
  Wait: j != i+2 (else it's consecutive), so j - i >= 2 is possible for j > i+1.

  Actually, the edges are:
  - (i, i+1): backbone, A[i][i+1] = 1 fixed
  - (i, j): non-backbone, A[i][j] flips
  - (i+1, j): could be backbone if j = i+2 (but we excluded that). Otherwise non-backbone.

  Wait, I need to be more careful. If j > i+1, the edges of the triple {i, i+1, j} are:
  - (i, i+1): backbone (fixed, A[i][i+1] = 1)
  - (i, j) if j > i: non-backbone if j-i >= 2 (which is true since j != i+2 is not needed, just j > i+1 means j-i >= 2). So this IS non-backbone.
  - (i+1, j) if j > i+1: non-backbone if j-(i+1) >= 2, i.e., j >= i+3. But j could be i+2!

  Case 2a: j = i+2. This is the CONSECUTIVE triple, already handled as type 1.

  Case 2b: j > i+2. Then both (i,j) and (i+1,j) are non-backbone.
  Backbone contributes: A[i][i+1] = 1.
  Clockwise 3-cycle (i->i+1->j->i): needs A[i][i+1]*A[i+1][j]*A[j][i] = A[i+1][j]*A[j][i]
  Counter-clockwise (i->j->i+1->i): needs A[i][j]*A[j][i+1]*A[i+1][i]
    But A[i+1][i] = 0 (backbone says A[i][i+1]=1). So counter-CW is impossible!

  Hmm wait, I need to think about direction. The backbone says i->i+1 (A[i][i+1]=1, A[i+1][i]=0).

  For triple {i, i+1, j} with j > i+2:
  CW (i->i+1->j->i): A[i][i+1] * A[i+1][j] * A[j][i] = 1 * A[i+1][j] * (1-A[i][j])
  CCW (i->j->i+1->i): A[i][j] * A[j][i+1] * A[i+1][i] = A[i][j] * (1-A[i+1][j]) * 0 = 0

  Wait, A[i+1][i] = 0 because A[i][i+1] = 1 (tournament). So CCW through this triple
  is ALWAYS 0. Only CW is possible.

  CW exists iff A[i+1][j] = 1 AND A[i][j] = 0 (since A[j][i] = 1-A[i][j]).

  After flip: A'[i+1][j] = 1-A[i+1][j], A'[i][j] = 1-A[i][j].
  CW after flip exists iff (1-A[i+1][j]) = 1 AND (1-A[i][j]) = 0
  i.e., A[i+1][j] = 0 AND A[i][j] = 1.

  So CW before = (A[i+1][j]=1, A[i][j]=0), CW after = (A[i+1][j]=0, A[i][j]=1).
  These are DIFFERENT conditions! So the 3-cycle on type 2b may or may not be preserved.

  In fact: before + after = 1 iff A[i+1][j] != A[i][j].
  before + after = 0 iff A[i+1][j] = A[i][j] and both = 0 (no cycle before or after)
    or A[i+1][j] = A[i][j] = 1 (no cycle before: need A[i][j]=0; no cycle after: need A[i][j]=1).

  Wait, let me redo:
  before = 1 iff A[i+1][j]=1 AND A[i][j]=0
  after = 1 iff A[i+1][j]=0 AND A[i][j]=1

  If A[i+1][j]=A[i][j]=0: before=0, after=0 (sum=0)
  If A[i+1][j]=1, A[i][j]=0: before=1, after=0 (sum=1)
  If A[i+1][j]=0, A[i][j]=1: before=0, after=1 (sum=1)
  If A[i+1][j]=1, A[i][j]=1: before=0, after=0 (sum=0)

  So before + after = 1 iff A[i+1][j] XOR A[i][j] = 1, else = 0.

  Hmm, so type 2b triples can have before + after = 0 or 1.

For type 2 with j < i (if j < i, then the triple {j, i, i+1}):
  Similar analysis with backbone i->i+1.
  CW (j->i->i+1->j): A[j][i] * A[i][i+1] * A[i+1][j] = (1-A[i][j]) * 1 * A[i+1][j]
  CCW (j->i+1->i->j): A[j][i+1] * A[i+1][i] * A[i][j] = (1-A[i+1][j]) * 0 * A[i][j] = 0

  Same structure. before = 1 iff A[i][j]=0 AND A[i+1][j]=1.
  After flip: A'[i][j]=1-A[i][j], A'[i+1][j]=1-A[i+1][j].
  After: (1-A[i][j])=0 AND (1-A[i+1][j])=1, i.e., A[i][j]=1 AND A[i+1][j]=0.
  Same XOR condition.

OK so this is getting complicated for type 2. Let me verify computationally.

kind-pasteur-2026-03-06-S25h
"""
from itertools import combinations

def tournament_from_tiling(n, tiling_bits):
    A = [[0]*n for _ in range(n)]
    for i in range(n-1):
        A[i][i+1] = 1
    idx = 0
    for i in range(n):
        for j in range(i+2, n):
            if (tiling_bits >> idx) & 1:
                A[i][j] = 1
            else:
                A[j][i] = 1
            idx += 1
    return A

def num_tiling_bits(n):
    return n*(n-1)//2 - (n-1)

def count_3cycles(A, n):
    t3 = 0
    for i, j, k in combinations(range(n), 3):
        if A[i][j] and A[j][k] and A[k][i]:
            t3 += 1
        if A[i][k] and A[k][j] and A[j][i]:
            t3 += 1
    return t3

def flip_tiling(tiling_bits, m):
    return tiling_bits ^ ((1 << m) - 1)

def classify_triple(a, b, c, n):
    """How many backbone edges does this triple have?"""
    backbone = 0
    if b == a + 1:
        backbone += 1
    if c == b + 1:
        backbone += 1
    if c == a + 1:  # can't happen since a < b < c and b >= a+1
        pass
    return backbone

for n in [5, 7, 9]:
    m = num_tiling_bits(n)
    total = min(2**m, 10000)  # sample if too many

    print(f"\nn={n}, m={m}, sampling {total} tilings")

    # Count consecutive triples
    n_consec = n - 2
    print(f"  Consecutive triples: {n_consec}")

    t3_parity_always_flips = True
    t3_diff_stats = []

    import random
    for trial in range(total):
        if 2**m <= 10000:
            bits = trial
        else:
            bits = random.randint(0, 2**m - 1)

        A_before = tournament_from_tiling(n, bits)
        t3_before = count_3cycles(A_before, n)

        flipped = flip_tiling(bits, m)
        A_after = tournament_from_tiling(n, flipped)
        t3_after = count_3cycles(A_after, n)

        diff = t3_after - t3_before
        t3_diff_stats.append(diff)

        # Check parity
        if (t3_before + t3_after) % 2 != n_consec % 2:
            # For consecutive triples: each contributes 1 to (before + after)
            # Total contribution of consec triples to (before+after) = n-2
            # But type 2 triples can add 0 or 1 each
            pass

        if t3_before % 2 == t3_after % 2:
            t3_parity_always_flips = False

    print(f"  t3 parity ALWAYS flips under tiling flip: {t3_parity_always_flips}")

    # Distribution of t3_after - t3_before
    from collections import Counter
    diff_counts = Counter(t3_diff_stats)
    print(f"  t3_diff distribution: {dict(sorted(diff_counts.items()))}")

    # Check: t3_before + t3_after mod 2
    sum_mod2 = Counter((t3_diff_stats[i] + 2*count_3cycles(tournament_from_tiling(n, i if 2**m<=10000 else 0), n)) % 2
                       for i in range(min(100, total)))

    # Simpler: just check sum parities
    sum_parities = set()
    for trial in range(min(total, 1000)):
        if 2**m <= 10000:
            bits = trial
        else:
            bits = random.randint(0, 2**m - 1)
        A_b = tournament_from_tiling(n, bits)
        A_a = tournament_from_tiling(n, flip_tiling(bits, m))
        s = count_3cycles(A_b, n) + count_3cycles(A_a, n)
        sum_parities.add(s % 2)

    print(f"  (t3_before + t3_after) mod 2: {sum_parities}")
    print(f"  n-2 mod 2 = {(n-2) % 2}")

    # Detailed: count type 1 (consecutive) and type 2 contributions separately
    if 2**m <= 10000:
        for bits in range(min(5, 2**m)):
            A = tournament_from_tiling(n, bits)
            Af = tournament_from_tiling(n, flip_tiling(bits, m))

            # Type 1 (consecutive) 3-cycles
            t3_consec_before = 0
            t3_consec_after = 0
            for i in range(n-2):
                a, b, c = i, i+1, i+2
                if A[a][b] and A[b][c] and A[c][a]:
                    t3_consec_before += 1
                if A[a][c] and A[c][b] and A[b][a]:
                    t3_consec_before += 1
                if Af[a][b] and Af[b][c] and Af[c][a]:
                    t3_consec_after += 1
                if Af[a][c] and Af[c][b] and Af[b][a]:
                    t3_consec_after += 1

            t3_b = count_3cycles(A, n)
            t3_a = count_3cycles(Af, n)
            t3_other_b = t3_b - t3_consec_before
            t3_other_a = t3_a - t3_consec_after

            print(f"\n  bits={bits}: t3_before={t3_b}, t3_after={t3_a}")
            print(f"    consec: before={t3_consec_before}, after={t3_consec_after}, sum={t3_consec_before+t3_consec_after}")
            print(f"    other:  before={t3_other_b}, after={t3_other_a}, sum={t3_other_b+t3_other_a}")
