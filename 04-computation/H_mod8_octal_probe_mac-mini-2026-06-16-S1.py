#!/usr/bin/env python3
"""
H_mod8_octal_probe — mac-mini-2026-06-16-S1

NEW probe motivated by the octal odd-square identity (2k+1)^2 = 8*T_k + 1:
odd squares are EXACTLY the residue-1 class mod 8, and stripping the octal '1'
gives a triangular number. The project's tile-count m = C(n-1,2) = T_{n-2} is
triangular; H(T) = #directed Hamiltonian paths (definitions.md) is ALWAYS ODD
(Redei), so H mod 8 in {1,3,5,7} -- the same odd-residue world odd squares live in.

Question: does H(T) mod 8 carry structure? In particular do the FORBIDDEN H-values
(7, 21, ...) and the realized H-spectrum cluster by residue mod 8?
  - 7  = 8*0 + 7  -> residue 7   (smallest forbidden)
  - 21 = 8*2 + 5  -> residue 5   (= 3*7, Tr(M_3^5))
H is computed EXACTLY by Held-Karp DP (O(2^n n^2)); fully brute force over all
2^C(n,2) labeled tournaments for n=3..6. Pure stdlib.
"""
from itertools import combinations

def ham_path_count(n, adj):
    """adj[v] = bitmask of out-neighbors of v. Count directed Hamiltonian paths."""
    full = (1 << n) - 1
    # dp[mask] = list over end-vertex v of #paths covering mask ending at v
    dp = [[0]*n for _ in range(1 << n)]
    for v in range(n):
        dp[1 << v][v] = 1
    for mask in range(1 << n):
        row = dp[mask]
        for v in range(n):
            c = row[v]
            if not c:
                continue
            outs = adj[v] & ~mask
            w = outs
            while w:
                u = (w & -w).bit_length() - 1
                dp[mask | (1 << u)][u] += c
                w &= w - 1
    return sum(dp[full])

def all_tournaments(n):
    """Yield adj (list of out-neighbor bitmasks) for every labeled tournament on n."""
    pairs = list(combinations(range(n), 2))
    P = len(pairs)
    for bits in range(1 << P):
        adj = [0]*n
        for k, (i, j) in enumerate(pairs):
            if bits & (1 << k):       # arc i->j
                adj[i] |= (1 << j)
            else:                      # arc j->i
                adj[j] |= (1 << i)
        yield adj

def line(t): print("\n" + "="*70 + "\n" + t + "\n" + "="*70)

line("H(T) = #directed Hamiltonian paths; spectrum + mod 8  (n=3..6 exhaustive)")
union = set()
spectra = {}
for n in range(3, 7):
    counts = {}
    for adj in all_tournaments(n):
        h = ham_path_count(n, adj)
        counts[h] = counts.get(h, 0) + 1
    realized = sorted(counts)
    spectra[n] = realized
    union |= set(realized)
    maxH = realized[-1]
    odd_gaps = [x for x in range(1, maxH+1, 2) if x not in counts]
    print(f"\nn={n}: realized H-values ({len(realized)}): {realized}")
    print(f"      all H odd: {all(h % 2 == 1 for h in realized)}  (Redei)")
    print(f"      ODD integers in [1,{maxH}] NOT realized (gaps): {odd_gaps}")
    # mod-8 split of realized values
    from collections import Counter
    rmod = Counter(h % 8 for h in realized)
    gmod = Counter(g % 8 for g in odd_gaps)
    print(f"      realized mod 8: {{r: rmod.get(r,0) for r in (1,3,5,7)}}  -> {dict((r,rmod.get(r,0)) for r in (1,3,5,7))}")
    print(f"      gaps     mod 8: {dict((r,gmod.get(r,0)) for r in (1,3,5,7))}")

line("Union H-spectrum n<=6 and the forbidden frontier through the mod-8 lens")
U = sorted(union)
maxU = U[-1]
forbidden = [x for x in range(1, maxU+1, 2) if x not in union]
print(f"realized (union n<=6): {U}")
print(f"odd integers in [1,{maxU}] NEVER realized at n<=6: {forbidden}")
print("Note the documented universal-forbidden values:")
for f in (7, 21):
    k = f // 8
    print(f"  {f:3d} = 8*{k}+{f-8*k}  -> residue {f%8} mod 8;  floor({f}/8)={k}"
          + ("  (triangular!)" if (lambda x:(lambda d:(int(d**0.5))**2==d)(8*x+1))(k) else ""))
print("\nInterpretation: H is odd (residues 1,3,5,7). The 'octal peel' floor(H/8)")
print("is the analogue of T_k for odd squares (which are exactly residue 1).")
print("If realized H clustered at residue 1 the analogy would be tight; the data")
print("above shows the ACTUAL residue distribution -- report it honestly.")

line("Cross-check vs odd squares: which realized H are themselves odd squares (res 1)?")
sq = [h for h in U if int(h**0.5)**2 == h]
print(f"realized H that are perfect squares: {sq}  (these are the residue-1 'odd-square-like' H)")
print("\nDONE.")
