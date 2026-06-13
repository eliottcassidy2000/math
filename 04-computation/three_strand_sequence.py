"""
three_strand_sequence.py -- kind-pasteur-2026-03-14-S87
Investigate the sequence and its 3 subsequences.

THE SEQUENCE:
1, 1, 2, 3, 4, 6, 10, 15, 20, 35, 56, 70, 126, 210, 252, 462, 792, 924,
1716, 3003, 3432, 6435, 11440, 12870

FIRST: Identify what this is. These look like CENTRAL BINOMIAL COEFFICIENTS
or related to Pascal's triangle.

Let me check: 1, 1, 2, 3, 6, 10, 20, 35, 70, 126, 252, 462, 924, 1716, 3432, 6435, 12870
That's C(2n, n) / something... or maybe the sequence is the INTERLEAVING
of 3 sequences from Pascal's triangle.

Actually: let me look at the subsequences by position mod 3:
  Pos 0,3,6,9,...:  1, 3, 10, 35, 126, 462, 1716, 6435
  Pos 1,4,7,10,...: 1, 4, 15, 56, 210, 792, 3003, 11440
  Pos 2,5,8,11,...: 2, 6, 20, 70, 252, 924, 3432, 12870

These are C(2n,n-1), C(2n,n), C(2n,n+1) or similar...

Let me check:
  C(0,0)=1, C(2,1)=2, C(4,2)=6, C(6,3)=20, C(8,4)=70, C(10,5)=252, C(12,6)=924, C(14,7)=3432, C(16,8)=12870
  That's 1, 2, 6, 20, 70, 252, 924, 3432, 12870 = CENTRAL BINOMIAL COEFFICIENTS!

And: 1, 3, 10, 35, 126, 462, 1716, 6435 = C(2n+1, n)?
  C(1,0)=1, C(3,1)=3, C(5,2)=10, C(7,3)=35, C(9,4)=126, C(11,5)=462, C(13,6)=1716, C(15,7)=6435
  YES! These are C(2n+1, n).

And: 1, 4, 15, 56, 210, 792, 3003, 11440 = C(2n+2, n+1)/something?
  C(2,1)=2... no. Let me check: C(3,1)=3, C(4,1)=4, C(5,2)=10, C(6,2)=15, C(7,3)=35, C(8,3)=56
  Hmm: 1, 4, 15, 56, 210, 792, 3003, 11440
  4 = C(4,1), 15 = C(6,2), 56 = C(8,3), 210 = C(10,4), 792 = C(12,5), 3003 = C(14,6), 11440 = C(16,7)
  Pattern: C(2n+2, n)? C(4,1)=4, C(6,2)=15, C(8,3)=56, C(10,4)=210, C(12,5)=792, C(14,6)=3003, C(16,7)=11440
  YES! C(2n+2, n) starting from n=1.
  Actually for n=0: C(2,0)=1. Yes!

So the THREE STRANDS are:
  Strand A (pos mod 3 = 0): C(2n+1, n) = 1, 3, 10, 35, 126, 462, 1716, 6435, ...
  Strand B (pos mod 3 = 1): C(2n+2, n) = 1, 4, 15, 56, 210, 792, 3003, 11440, ...
  Strand C (pos mod 3 = 2): C(2n, n) = 1, 2, 6, 20, 70, 252, 924, 3432, 12870, ...
  (shifted: C(2n+2, n+1) = C(2n, n) for the central ones)

Actually let me just compute and verify.
"""

import sys, math
from collections import Counter, defaultdict

sys.stdout.reconfigure(encoding='utf-8')

def C(n, k):
    if k < 0 or k > n: return 0
    return math.comb(n, k)

def main():
    print("=" * 70)
    print("THREE-STRAND SEQUENCE INVESTIGATION")
    print("kind-pasteur-2026-03-14-S87")
    print("=" * 70)

    seq = [1, 1, 2, 3, 4, 6, 10, 15, 20, 35, 56, 70, 126, 210, 252,
           462, 792, 924, 1716, 3003, 3432, 6435, 11440, 12870]

    print(f"\n  The sequence ({len(seq)} terms):")
    for i, v in enumerate(seq):
        print(f"    a({i}) = {v}")

    # Split into 3 strands
    s0 = [seq[i] for i in range(0, len(seq), 3)]
    s1 = [seq[i] for i in range(1, len(seq), 3)]
    s2 = [seq[i] for i in range(2, len(seq), 3)]

    print(f"\n  Strand 0 (pos mod 3 = 0): {s0}")
    print(f"  Strand 1 (pos mod 3 = 1): {s1}")
    print(f"  Strand 2 (pos mod 3 = 2): {s2}")

    # ============================================================
    # IDENTIFY EACH STRAND
    # ============================================================
    print(f"\n{'='*70}")
    print("IDENTIFYING THE THREE STRANDS")
    print(f"{'='*70}")

    # Check strand 2: central binomial coefficients C(2n, n)?
    print(f"\n  Strand 2: {s2}")
    print(f"  Testing C(2n, n):")
    for i, v in enumerate(s2):
        candidate = C(2*i, i)
        print(f"    n={i}: C({2*i},{i}) = {candidate}, actual = {v}, match = {candidate == v}")

    # Check strand 0: C(2n+1, n)?
    print(f"\n  Strand 0: {s0}")
    print(f"  Testing C(2n+1, n):")
    for i, v in enumerate(s0):
        candidate = C(2*i+1, i)
        print(f"    n={i}: C({2*i+1},{i}) = {candidate}, actual = {v}, match = {candidate == v}")

    # Check strand 1: C(2n+2, n+1)?
    print(f"\n  Strand 1: {s1}")
    print(f"  Testing C(2n+2, n):")
    for i, v in enumerate(s1):
        candidate = C(2*i+2, i)
        print(f"    n={i}: C({2*i+2},{i}) = {candidate}, actual = {v}, match = {candidate == v}")

    # Actually let me try C(2n+2, n+1):
    print(f"\n  Testing C(2n+2, n+1):")
    for i, v in enumerate(s1):
        candidate = C(2*i+2, i+1)
        print(f"    n={i}: C({2*i+2},{i+1}) = {candidate}, actual = {v}, match = {candidate == v}")

    # ============================================================
    # THE UNIFIED FORMULA
    # ============================================================
    print(f"\n{'='*70}")
    print("THE UNIFIED FORMULA")
    print(f"{'='*70}")

    # The full sequence interleaves C(2n+1,n), C(2n+2,n+1), C(2n+2,n+1)... hmm
    # Let me check: position k in the sequence corresponds to what binomial?

    # Position 0: 1 = C(1,0)
    # Position 1: 1 = C(2,1)
    # Position 2: 2 = C(2,1)... or C(2,0)=1? No, C(2,1)=2. Hmm.
    # Wait: maybe it's ALL central-ish binomial coefficients from Pascal's triangle.

    # Let me look at the MIDDLE of Pascal's triangle row by row:
    # Row 0: C(0,0) = 1
    # Row 1: C(1,0)=1, C(1,1)=1
    # Row 2: C(2,1) = 2
    # Row 3: C(3,1)=3, C(3,2)=3
    # Row 4: C(4,2) = 6
    # Row 5: C(5,2)=10, C(5,3)=10
    # Row 6: C(6,3) = 20
    # ...

    # The central element(s) of each row of Pascal's triangle!
    # For even row 2k: one central element C(2k, k)
    # For odd row 2k+1: two central elements C(2k+1, k) = C(2k+1, k+1)

    print(f"\n  Central elements of Pascal's triangle, row by row:")
    for row in range(20):
        if row % 2 == 0:
            k = row // 2
            vals = [C(row, k)]
            print(f"    Row {row:2d}: C({row},{k}) = {vals[0]}")
        else:
            k = row // 2
            vals = [C(row, k), C(row, k+1)]
            print(f"    Row {row:2d}: C({row},{k}) = {vals[0]}, C({row},{k+1}) = {vals[1]}")

    # Read off: 1, 1, 1, 2, 3, 3, 6, 10, 10, 20, 35, 35, 70, 126, 126, 252, ...
    # This doesn't match the given sequence exactly. The given sequence has
    # 1, 1, 2, 3, 4, 6, 10, 15, 20, 35, 56, 70, ...

    # Hmm, let me reconsider. Maybe the 3 strands are:
    # Reading Pascal's triangle in a DIAGONAL way?

    # Or: the sequence is the SORTED central-ish binomials?
    # Let me try a different approach: look up in OEIS.

    # The sequence 1,1,2,3,4,6,10,15,20,35,56,70,126...
    # Let me check: this might be C(n, floor(n/3)) or similar.

    print(f"\n  Testing: a(k) = C(k, floor(k/3)):")
    for k in range(len(seq)):
        candidate = C(k, k//3)
        match = (candidate == seq[k])
        print(f"    k={k:2d}: C({k},{k//3}) = {candidate:6d}, actual = {seq[k]:6d}, "
              f"{'MATCH' if match else 'no'}")

    # That's not it either. Let me try C(k, floor(k/2)):
    # Actually, the strands were:
    # s0 = [1, 3, 10, 35, 126, 462, 1716, 6435]
    # s1 = [1, 4, 15, 56, 210, 792, 3003, 11440]
    # s2 = [2, 6, 20, 70, 252, 924, 3432, 12870]

    # s0[n] = C(2n+1, n): Catalan-related!
    # C(2n+1, n) = (2n+1)! / (n! * (n+1)!) = (2n+1) * C(2n, n) / (n+1)

    # s1[n] = C(2(n+1), n+1) for n >= 0?
    # C(2,1)=2, C(4,2)=6, C(6,3)=20, C(8,4)=70, C(10,5)=252...
    # That's strand 2, not strand 1!

    # Let me re-examine. Let me just check ALL binomial coefficients:
    print(f"\n  Checking all C(a,b) for small a,b:")
    target_set = set(seq)
    hits = defaultdict(list)
    for a in range(25):
        for b in range(a+1):
            v = C(a, b)
            if v in target_set:
                hits[v].append((a, b))

    for v in sorted(target_set):
        if v in hits:
            print(f"    {v}: {hits[v][:5]}")

    # ============================================================
    # REINTERPRET: Maybe the sequence is read along diagonals of Pascal
    # ============================================================
    print(f"\n{'='*70}")
    print("READING PASCAL'S TRIANGLE ALONG SPECIFIC DIAGONALS")
    print(f"{'='*70}")

    # The "shallow diagonal" of Pascal's triangle:
    # Reading C(n, floor(n/2)) for n = 0, 1, 2, 3, ...
    print(f"\n  C(n, floor(n/2)) for n=0..23:")
    diag = [C(n, n//2) for n in range(24)]
    print(f"  {diag}")
    print(f"  Given: {seq}")
    match_count = sum(1 for i in range(min(len(diag), len(seq))) if diag[i] == seq[i])
    print(f"  Matches: {match_count}/{min(len(diag), len(seq))}")

    # Try C(n, floor(n/3)):
    print(f"\n  C(n, floor(n/3)):")
    diag2 = [C(n, n//3) for n in range(24)]
    print(f"  {diag2}")

    # Try the middle third:
    print(f"\n  C(n, round(n/3)):")
    diag3 = [C(n, round(n/3)) for n in range(24)]
    print(f"  {diag3}")

    # ============================================================
    # CATALAN CONNECTION
    # ============================================================
    print(f"\n{'='*70}")
    print("CATALAN AND RELATED SEQUENCES")
    print(f"{'='*70}")

    # Catalan numbers: C_n = C(2n,n)/(n+1) = 1, 1, 2, 5, 14, 42, 132, 429, ...
    catalan = [C(2*n, n) // (n+1) for n in range(15)]
    print(f"  Catalan: {catalan}")

    # Central binomial: C(2n,n) = 1, 2, 6, 20, 70, 252, 924, 3432, 12870
    central = [C(2*n, n) for n in range(15)]
    print(f"  Central binomial C(2n,n): {central}")

    # The strand 2 IS central binomial starting from C(2,1)=2:
    # 2, 6, 20, 70, 252, 924, 3432, 12870 = C(2n,n) for n=1,2,3,...
    print(f"  Strand 2 = C(2n,n) for n=1,...: {[C(2*n,n) for n in range(1,9)]}")
    print(f"  Actual strand 2: {s2}")

    # Ratios between strands:
    print(f"\n  Ratios between strands:")
    for i in range(min(len(s0), len(s1), len(s2))):
        if s0[i] > 0 and s1[i] > 0 and s2[i] > 0:
            print(f"    n={i}: s0={s0[i]}, s1={s1[i]}, s2={s2[i]}, "
                  f"s1/s0={s1[i]/s0[i]:.4f}, s2/s1={s2[i]/s1[i]:.4f}, "
                  f"s2/s0={s2[i]/s0[i]:.4f}")

    # ============================================================
    # CONNECTION TO TOURNAMENTS
    # ============================================================
    print(f"\n{'='*70}")
    print("CONNECTION TO TOURNAMENTS")
    print(f"{'='*70}")

    # Do any of these numbers appear in tournament theory?
    print(f"\n  Tournament connections:")
    print(f"  C(n,2) arc counts: {[C(n,2) for n in range(3,12)]}")
    print(f"  C(n-1,2) tiling counts: {[C(n-1,2) for n in range(3,12)]}")
    print(f"  n! tournament labels: {[math.factorial(n) for n in range(1,10)]}")
    print(f"  2^C(n,2) total tournaments: {[2**C(n,2) for n in range(2,8)]}")

    # The GS degrees of freedom:
    print(f"  GS DOF: {[(C(n-1,2) + (n-1)//2)//2 for n in range(3,12)]}")

    # Do any strand values appear as tournament quantities?
    for v in sorted(target_set):
        for n in range(3, 10):
            if v == C(n, 2):
                print(f"    {v} = C({n},2) = number of arcs at n={n}")
            if v == math.factorial(n):
                print(f"    {v} = {n}! = number of labeled transitive tournaments")

    # The central binomials C(2n,n) have a natural role:
    # C(2n,n) = number of lattice paths from (0,0) to (n,n)
    # = number of ways to arrange n 0s and n 1s
    # In tournaments: a tournament on 2n vertices with equal scores (n-1 each)
    # relates to balanced binary sequences.

    print(f"\n  C(2n,n) counts lattice paths = balanced binary sequences.")
    print(f"  Tournament score sequences are like CONSTRAINED binary sequences.")
    print(f"  The strand interleaving mirrors how odd/even row parity")
    print(f"  determines whether Pascal's row has 1 or 2 central elements.")

    # ============================================================
    # THE THREE-FOLD SYMMETRY
    # ============================================================
    print(f"\n{'='*70}")
    print("THE THREE-FOLD SYMMETRY")
    print(f"{'='*70}")
    print(f"""
  The 3 strands correspond to residues mod 3 of the index.
  This THREE-FOLD splitting is like:

  1. TOURNAMENT: The 3-cycle (fundamental odd cycle)
     - 3-cycles are the basic building blocks of tournament complexity
     - The blue skeleton bipartition is by t3 PARITY (mod 2 of 3-cycles)

  2. FIBONACCI: The Tribonacci generalization
     - T(n) = T(n-1) + T(n-2) + T(n-3) has 3 strands
     - The characteristic polynomial x^3 - x^2 - x - 1 has 3 roots

  3. PASCAL: Three diagonals of Pascal's triangle
     - C(3n, n), C(3n+1, n), C(3n+2, n) form 3 strands
     - Related to the cube root of the generating function

  THE META-STRUCTURE:
  Just as the Fibonacci 2-strand interleaving corresponds to the
  ORDER-2 involution (T -> T^op), a 3-strand interleaving would
  correspond to an ORDER-3 SYMMETRY.

  For tournaments: is there a natural ORDER-3 operation?
  - The 3-cycle reversal (which preserves H at n=4 but not n>=5)
  - The Z/3 action on the 3 arcs of a triangle
  - The mod-3 structure of alpha_1 (which determines H mod 3)
""")

    # ============================================================
    # VERIFY: Is the sequence C(n, floor(n/2))?
    # ============================================================
    print(f"\n{'='*70}")
    print("FINAL IDENTIFICATION ATTEMPT")
    print(f"{'='*70}")

    # Let me try: the sequence might be the SORTED values of
    # {C(n+1, k) : 0 <= k <= n+1} read in a specific order.
    # Or it could be from a SPECIFIC diagonal.

    # Actually: 1,1,2,3,4,6,10,15,20,35,56,70,126,210,252,462,...
    # The Motzkin-like numbers? No.
    # Delannoy? No.
    # Narayana? No.

    # Let me check OEIS-style: 1,1,2,3,4,6,10,15,20,35
    # This doesn't match any single well-known sequence.
    # But the 3 strands DO match:
    #   C(2n+1, n): 1, 3, 10, 35, 126, 462, 1716, 6435
    #   ?: 1, 4, 15, 56, 210, 792, 3003, 11440
    #   C(2n, n): 1, 2, 6, 20, 70, 252, 924, 3432, 12870

    # For strand 1: 1, 4, 15, 56, 210, 792, 3003, 11440
    # Check: C(2n+2, n+1)?
    # C(2,1)=2, C(4,2)=6... no.
    # C(3,1)=3... no.
    # Let me try C(2n+1, n-1) for n=1,2,...:
    # C(3,0)=1, C(5,1)=5... no.
    # C(2n+2, n) for n=0,1,...:
    # C(2,0)=1, C(4,1)=4, C(6,2)=15, C(8,3)=56, C(10,4)=210, C(12,5)=792,
    # C(14,6)=3003, C(16,7)=11440
    # YES! Strand 1 = C(2n+2, n) for n=0,1,...

    print(f"\n  IDENTIFICATION:")
    print(f"  Strand 0 = C(2n+1, n): {[C(2*n+1,n) for n in range(8)]}")
    print(f"  Strand 1 = C(2n+2, n): {[C(2*n+2,n) for n in range(8)]}")
    print(f"  Strand 2 = C(2n+2, n+1) = C(2n, n) shifted:")
    print(f"    C(2n+2,n+1): {[C(2*n+2,n+1) for n in range(8)]}")
    print(f"    Actual s2: {s2}")

    # Hmm, C(2*0+2, 1) = C(2,1) = 2. s2[0] = 2. ✓
    # C(2*1+2, 2) = C(4,2) = 6. s2[1] = 6. ✓
    # C(2*2+2, 3) = C(6,3) = 20. s2[2] = 20. ✓
    # C(2*3+2, 4) = C(8,4) = 70. s2[3] = 70. ✓
    # YES!

    print(f"\n  *** IDENTIFICATION CONFIRMED ***")
    print(f"  The full sequence interleaves three binomial strands:")
    print(f"    a(3k)   = C(2k+1, k)     [odd rows, left of center]")
    print(f"    a(3k+1) = C(2k+2, k)     [even rows, left of center]")
    print(f"    a(3k+2) = C(2k+2, k+1)   [even rows, right of center]")
    print(f"")
    print(f"  This reads Pascal's triangle from BOTH sides of center,")
    print(f"  alternating: center of odd row, then left-of-center and")
    print(f"  right-of-center of even row.")

    # Verify the full sequence
    print(f"\n  Full verification:")
    reconstructed = []
    for k in range(8):
        reconstructed.append(C(2*k+1, k))
        reconstructed.append(C(2*k+2, k))
        reconstructed.append(C(2*k+2, k+1))

    print(f"  Reconstructed: {reconstructed}")
    print(f"  Original:      {seq}")
    print(f"  Match: {reconstructed[:len(seq)] == seq}")

    # ============================================================
    # GENERATING FUNCTION AND ASYMPTOTICS
    # ============================================================
    print(f"\n{'='*70}")
    print("GENERATING FUNCTIONS AND ASYMPTOTICS")
    print(f"{'='*70}")

    # Each strand has a known generating function:
    # sum C(2n+1, n) x^n = 1/sqrt(1-4x) * something
    # sum C(2n, n) x^n = 1/sqrt(1-4x)

    # Asymptotics: C(2n, n) ~ 4^n / sqrt(pi*n)
    print(f"\n  Asymptotic ratios:")
    for n in range(1, 9):
        s0_val = C(2*n+1, n)
        s1_val = C(2*n+2, n)
        s2_val = C(2*n+2, n+1)

        # Ratios between strands
        print(f"    n={n}: s0/s2_prev = {s0_val/C(2*n,n):.4f}, "
              f"s1/s0 = {s1_val/s0_val:.4f}, "
              f"s2/s1 = {s2_val/s1_val:.4f}")

    # Each consecutive ratio approaches 2:
    print(f"\n  Consecutive term ratios (approach 4^{1/3} ≈ 1.587?):")
    for i in range(1, len(seq)):
        ratio = seq[i] / seq[i-1]
        print(f"    a({i})/a({i-1}) = {seq[i]}/{seq[i-1]} = {ratio:.4f}")

    print(f"\n{'='*70}")
    print("DONE")
    print(f"{'='*70}")

if __name__ == '__main__':
    main()
