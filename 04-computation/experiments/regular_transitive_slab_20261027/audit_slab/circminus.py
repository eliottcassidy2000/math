"""Self-converse balanced test tournaments for Lemma B': for every circulant Z_N{S} (N odd, S ∩ -S = ∅, |S| = (N-1)/2)
delete vertex 0; the result is balanced on n = N-1 vertices and keeps the anti-automorphism x -> -x.  Output in
gentourng's upper-triangle format.  Usage: python3 circminus.py N
"""
import sys
N = int(sys.argv[1]); n = N - 1; m = (N - 1) // 2
for mask in range(1 << m):
    S = {d if mask >> (d - 1) & 1 else N - d for d in range(1, m + 1)}
    V = list(range(1, N))   # delete vertex 0
    A = [[1 if a != b and (b - a) % N in S else 0 for b in V] for a in V]
    print("".join("1" if A[i][j] else "0" for i in range(n) for j in range(i + 1, n)))
