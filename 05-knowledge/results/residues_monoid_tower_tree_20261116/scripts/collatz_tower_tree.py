#!/usr/bin/env python3
"""collatz_tower_tree.py: (1) the mod-3 tower: rows, period 12, the halving and odd arrows of the accelerated map
T(n) = n/2 (even), (3n+1)/2 (odd) with their row lengths and families; (2) the pruned accelerated inverse tree: prune
the multiples of 3, contract unary nodes, check full binarity, edge lengths in {1,2,3,4}, ray/branch lengths, no two
consecutive sterile forks, and label reconstruction from (lengths, left/right)."""
def row(n): return 2*(n//3) if n % 3 == 0 else 2*((n-1)//3) + 1
def T(n): return n//2 if n % 2 == 0 else (3*n+1)//2
print("(1) rows: n -> row(n) for n = 1..14:", [(n, row(n)) for n in range(1, 15)])
print("    period 12 classes: multiples of 3 {0,3,6,9} (odd 3,9; even 0,6); even non-multiples {2,4,8,10} (4,8 = 0 mod 4; 2,10 = 2*odd); 6k+-1 {1,5,7,11}")
print("    halving arrows 2m -> m, m not divisible by 3, (m, row(2m) - row(m)):", [(m, row(2*m) - row(m)) for m in range(1, 30) if m % 3])
print("    law check: row(2m)-row(m) = 2*((2m-1)//3 - (m-1)//3) for m <= 10^5 not divisible by 3:", all(row(2*m)-row(m) == 2*((2*m-1)//3 - (m-1)//3) for m in range(1, 100001) if m % 3))
print("    growth: m -> m+3 adds 2 to the length (m <= 10^5):", all((row(2*(m+3))-row(m+3)) - (row(2*m)-row(m)) == 2 for m in range(1, 100001) if m % 3))
print("    short-row chain of 3: rows of 3, 6, 12, 24, 48, 96:", [row(3*2**k) for k in range(6)], "differences", [row(3*2**(k+1)) - row(3*2**k) for k in range(5)])
for res, name in [(1, "left (n = 1 mod 6)"), (3, "short (n = 3 mod 6)"), (5, "right (n = 5 mod 6)")]:
    fam = [(n, row(n), row(T(n))) for n in range(res, res + 6*6, 6)]
    print(f"    odd-arrow family {name}: (n, row(n), row(T(n))) =", fam, "; start rows step 4, lengths step 2:",
          all(fam[i+1][1] - fam[i][1] == 4 and (fam[i+1][2]-fam[i+1][1]) - (fam[i][2]-fam[i][1]) == 2 for i in range(len(fam)-1)))
# (2) the pruned tree
def children(n):   # accelerated inverse: 2n always; (2n-1)/3 when n = 2 mod 3
    c = [2*n]
    if n % 3 == 2: c.append((2*n - 1)//3)
    return c
def sterile(n): return n % 3 == 0
# forks: n = 2 mod 3; contracted tree on forks whose odd child is fertile
def contracted(n, depth):
    """return the contracted edges from fork n: list of (child_fork, length, kind) with kind 'ray' or 'branch'"""
    out = []
    # ray continuation: 2n (unary, 1 mod 3), 4n (fork). If 4n's odd child is sterile, continue to 16n (skip the sterile fork).
    m = 4*n; L = 2
    while (2*m - 1)//3 % 3 == 0:  # sterile odd child -> this fork is unary after pruning; contract
        m *= 4; L += 2
    out.append((m, L, 'ray'))
    k = (2*n - 1)//3            # odd child
    if k % 3 == 0: return out    # sterile branch (then n itself is unary: handled by the caller's contraction)
    # from k reach its first fork: k = 2 mod 3 -> k is a fork (length 1); k = 1 mod 3 -> 2k is a fork (length 2)
    if k % 3 == 2: m, L = k, 1
    else: m, L = 2*k, 2
    while (2*m - 1)//3 % 3 == 0: m *= 4; L += 2
    out.append((m, L, 'branch'))
    return out
root = 8  # 8 = 2 mod 3, odd child 5 (fertile): the first fork above the cycle {1, 2}
frontier = [root]; lengths = {}; binary = True; depth = 0; nodes = 1; labels = {root: ()}
max_depth = 18
while depth < max_depth:
    nxt = []
    for n in frontier:
        ch = contracted(n, depth)
        if len(ch) != 2: binary = False; print("   NOT BINARY at", n, ch)
        for (m, L, kind) in ch:
            lengths[(kind, L)] = lengths.get((kind, L), 0) + 1
            labels[m] = labels[n] + ((kind, L),); nxt.append(m)
    frontier = nxt; nodes += len(nxt); depth += 1
print("(2) pruned contracted tree from fork 8, depth", max_depth, ": nodes", nodes, "; full binary:", binary)
print("    edge (kind, length) counts:", dict(sorted(lengths.items())))
print("    all lengths in {1,2,3,4}:", all(L in (1,2,3,4) for (_, L) in lengths), "; ray lengths:", sorted(set(L for (k, L) in lengths if k == 'ray')), "; branch lengths:", sorted(set(L for (k, L) in lengths if k == 'branch')))
print("    labels distinct (reconstruction from lengths + left/right):", len(labels) == nodes)
print("    no two consecutive sterile forks on a ray (n = 5 mod 9 -> 4n = 2 mod 9), checked n <= 10^6:", all(not ((n % 9 == 5) and ((4*n) % 9 == 5)) for n in range(2, 10**6, 3)))
# which residues mod 9 of the fork give which lengths
print("    fork n mod 9 -> (ray length, branch length):", {r: (lambda n: tuple((L, k) for (m, L, k) in contracted(n, 0)))(r + 9*2) for r in (2, 5, 8)})
