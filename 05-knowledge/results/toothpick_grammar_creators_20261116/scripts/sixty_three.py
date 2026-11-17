#!/usr/bin/env python3
"""sixty_three.py: the two first creators of 63 (order 8): digit-0 vertices, out-degree sequences, converse pair."""
exec(open('nearjoin_check.py').read().split('# (1)')[0])
def rev(T): n = len(T); return [[T[j][i] for j in range(n)] for i in range(n)]
hs = [int(l.split()[1]) for l in open('counts.txt') if l.startswith('8 ')]
cl = classes(8); assert len(cl) == len(hs)
for i, T in enumerate(cl):
    if hs[i] != 63: continue
    n = 8; N = 1 << n; dp, H = ham_all(T); dpr, _ = ham_all(rev(T))
    B = [sum(dpr[m][t] * H[(N-1) & ~m] for m in range(1, N) if (m >> t) & 1) for t in range(n)]
    print(f"class {i}: H = {H[N-1]}, out-degrees {sorted(sum(r) for r in T)}, root responses {sorted(B)}, digit-0 vertices {sum(1 for b in B if b == 63)}, strong {strong(T)}, deletion counts {sorted(H[(N-1) & ~(1 << v)] for v in range(n))}")
