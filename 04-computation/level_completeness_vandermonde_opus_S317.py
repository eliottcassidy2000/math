# opus-2026-07-15-S317 -- HYP-6940: (A) THM-865 referee -- the leftmost-tie-
# split walk realizes EVERY level (verify n = 3..9: walk from the floor hits
# exactly the census levels; every split Landau-legal). (B) THE VANDERMONDE
# TRUNCATION LAW -- the polygonal/polyhedral triangles: row sums, skip sums,
# diagonal differences, the Pascal coefficient law.
from math import comb
from collections import defaultdict

# ---------- (A) THM-865 referee
print("(A) THM-865 -- leftmost-tie-split walk vs the level census:")
def landau(seq):
    s = sorted(seq)
    return (sum(s) == comb(len(s), 2) and
            all(sum(s[:k+1]) >= comb(k+1, 2) for k in range(len(s))))

for n in range(3, 10):
    # census of levels (exact, via score-sequence enumeration)
    levels = set()
    seqs = set()
    def rec(pre):
        k = len(pre)
        if k == n:
            if sum(pre) == comb(n, 2): seqs.add(tuple(pre))
            return
        lo = pre[-1] if pre else 0
        for s in range(lo, n):
            nxt = pre + [s]
            if sum(nxt) < comb(k+1, 2): continue
            if sum(nxt) + (n-k-1)*(n-1) < comb(n, 2): continue
            if sum(nxt) > comb(n, 2): continue
            rec(nxt)
    rec([])
    for q in seqs:
        levels.add(sum((2*s - (n-1))**2 for s in q))
    levels = sorted(levels)
    # the walk: start at the floor sequence (regular / near-regular), split
    # leftmost tie repeatedly; check every split legal and Delta x = +8
    if n % 2 == 1: cur = [ (n-1)//2 ]*n
    else: cur = [ (n-2)//2 ]*(n//2) + [ n//2 ]*(n//2)
    assert landau(cur)
    walk = [sum((2*s - (n-1))**2 for s in cur)]
    ok = True
    while True:
        s = sorted(cur)
        tie = next((i for i in range(n-1) if s[i] == s[i+1]), None)
        if tie is None: break                     # transitive
        s[tie] -= 1; s[tie+1] += 1                 # split leftmost tie
        if not landau(s): ok = False; break
        cur = s
        walk.append(sum((2*v - (n-1))**2 for v in cur))
    steps8 = all(b - a == 8 for a, b in zip(walk, walk[1:]))
    print(f"   n={n}: levels={len(levels)} walk-visits={len(walk)} "
          f"all-legal={ok} all-steps-8={steps8} "
          f"walk==census: {walk == levels}")

# ---------- (B) the two triangles
print("\n(B) THE VANDERMONDE TRUNCATION LAW:")
def pas(i, j):   # polyhedral: i-th entry (i>=1) of column j
    return comb(i + j - 1, j)
def pol(i, j):   # polygonal: i-th entry of column j: (j+1)-gonal; col0 = 1s
    if j == 0: return 1
    return i + (j - 1) * comb(i, 2)

# triangle rows (row n has columns j = 0..n): entry at (n, j) = value with
# index i = n - j + 1 (so column j starts at row j-? with 1). Row sums:
def rowsum(f, n):
    return sum(f(n - j + 1, j) for j in range(n + 1) if n - j + 1 >= 1)
print("   polyhedral row sums (2^n):", [rowsum(pas, n) for n in range(0, 9)])
print("   polygonal  row sums      :", [rowsum(pol, n) for n in range(0, 10)])
print("   A000127                  : 1 2 4 8 16 31 57 99 163 256")
def skipsum(f, n):
    return sum(f(n - 2*j + 1, j) for j in range(n // 2 + 1) if n - 2*j + 1 >= 1)
print("   polyhedral skip sums (Fib):", [skipsum(pas, n) for n in range(0, 12)])
print("   polygonal  skip sums      :", [skipsum(pol, n) for n in range(0, 14)])
print("   owner's sequence          : 1 1 2 3 5 8 13 21 33 51 76 111 157 218")

# diagonal differences: i-th entries, columns j
print("\n   diagonal differences Delta(i,j) = pas - pol, and the Pascal law:")
for i in range(3, 8):
    row = [pas(i, j) - pol(i, j) for j in range(2, 9)]
    law = [sum(comb(i, k+1) * comb(j-1, k) for k in range(2, i))
           for j in range(2, 9)]
    print(f"   i={i}: diffs={row}  law-match={row == law}")
# the closed second-diagonal form
sd = [pas(4, j) - pol(4, j) for j in range(3, 8)]
cf = [(j-1)*(j-2)*(j+9)//6 for j in range(3, 8)]
print(f"   i=4 closed form (j-1)(j-2)(j+9)/6: {sd == cf} ({sd})")

# the Moser analogy: A000127(n) = sum_{k<=4} C(n-1,k); the polygonal skip sum
# as a truncated-Fibonacci: derive its binomial form empirically
print("\n   polygonal skip sum in truncated-diagonal-binomial form:")
for n in range(0, 14):
    v = skipsum(pol, n)
    # Fibonacci diagonal: F(n+1) = sum_j C(n-j, j); try truncations with
    # Vandermonde-tail corrections: v =? sum_j [ (n-2j+1) + (j-1)C(n-2j+1,2) ]
    direct = sum((n-2*j+1) + (j-1)*comb(n-2*j+1, 2)
                 for j in range(1, n//2 + 1) if n-2*j+1 >= 1) + (1 if n >= 0 else 0)
    print(f"      n={n:2d}: skip={v:4d} direct-form={direct:4d} "
          f"F(n+1)={sum(comb(n-j, j) for j in range(n//2+1)):4d}")
