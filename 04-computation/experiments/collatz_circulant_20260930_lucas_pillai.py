"""Pillai gaps |2^K - 3^X| that are Lucas or Fibonacci numbers, or Fibonacci-group orders |O_j| = |L_j - 1 - (-1)^j|;
the free shapes they give (THM-4484: shape (K,X) free for 3y +- d iff (2^K - 3^X) | d)."""
from math import gcd
KMAX = 400; JMAX = 800
L = [2, 1]; F = [0, 1]
for j in range(2, JMAX+1): L.append(L[-1]+L[-2]); F.append(F[-1]+F[-2])
Lset = {L[j]: j for j in range(JMAX+1)}; Fset = {F[j]: j for j in range(2, JMAX+1)}
O = {abs(L[j] - 1 - (-1)**j): j for j in range(1, JMAX+1)}
hits = []
for K in range(1, KMAX+1):
    for X in range(0, K+1):
        g = 2**K - 3**X; a = abs(g)
        tags = []
        if a in Lset: tags.append(f"L_{Lset[a]}")
        if a in Fset: tags.append(f"F_{Fset[a]}")
        if a in O: tags.append(f"|O_{O[a]}|=|F(2,{O[a]})^ab|")
        if tags: hits.append((K, X, g, tags))
print("Pillai gaps 2^K - 3^X (K <= %d) equal to Lucas/Fibonacci numbers or Fibonacci-group orders:" % KMAX)
for K, X, g, tags in hits: print(f"  (K,X)=({K},{X}): 2^K-3^X = {g}  {tags}  -> free shape ({K},{X}) for 3y{'+' if g>0 else '-'}{abs(g)}  (density {X/K:.3f}, {'contracting' if g>0 else 'expanding'})")
print("Fibonacci-group orders |O_j| for j=1..12:", [abs(L[j]-1-(-1)**j) for j in range(1,13)])
