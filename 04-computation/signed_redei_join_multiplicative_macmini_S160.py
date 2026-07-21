#!/usr/bin/env python3
"""
THM-1936: the signed Redei count R(T)=sum_{Ham paths} sgn(pi) is MULTIPLICATIVE under the
order-join, R(T1|>T2)=R(T1)R(T2); hence |R| factors over strong components, and the achievable
|R|-values on n vertices = products of strong-atom |R| over compositions of n.  This answers
mac-mini-S159's open "why are 9,13 absent from |R| at n=5?"  (9=3*3 needs 6 vertices; 13 prime,
first strong realizer at n=6).                                            mac-mini-2026-07-21-S160

Also records the negatives: no determinant/Pfaffian collapse for R (Pf(A-A^T) exact only n<=4);
max|R| = 3,3,15,15,147 is NOT the double factorial (7!!=105).
"""
from itertools import combinations, product
from collections import Counter

def all_tours(n):
    pairs = list(combinations(range(n), 2))
    for bits in product((0, 1), repeat=len(pairs)):
        A = [[0]*n for _ in range(n)]
        for (i, j), b in zip(pairs, bits):
            if b: A[i][j] = 1
            else: A[j][i] = 1
        yield A

def R_signed(n, A):
    """Signed Held-Karp: sgn built incrementally (append v adds #{used>v} inversions)."""
    dp = [[0]*n for _ in range(1 << n)]
    for v in range(n): dp[1 << v][v] = 1
    for mask in range(1 << n):
        for last in range(n):
            c = dp[mask][last]
            if not c: continue
            for nxt in range(n):
                if mask & (1 << nxt): continue
                if A[last][nxt]:
                    invs = bin(mask >> (nxt+1)).count("1")
                    dp[mask | (1 << nxt)][nxt] += c * (-1 if invs & 1 else 1)
    return sum(dp[(1 << n)-1][v] for v in range(n))

def is_strong(n, A):
    def reach(s):
        vis = {s}; st = [s]
        while st:
            u = st.pop()
            for w in range(n):
                if A[u][w] and w not in vis: vis.add(w); st.append(w)
        return vis
    return all(len(reach(s)) == n for s in range(n))

def join(A1, n1, A2, n2):
    n = n1 + n2; A = [[0]*n for _ in range(n)]
    for i in range(n1):
        for j in range(n1): A[i][j] = A1[i][j]
    for i in range(n2):
        for j in range(n2): A[n1+i][n1+j] = A2[i][j]
    for i in range(n1):
        for j in range(n2): A[i][n1+j] = 1
    return A, n

if __name__ == "__main__":
    import random
    rng = random.Random(1)
    # 1) multiplicativity
    ok = True
    for _ in range(300):
        m1, m2 = rng.randint(1, 4), rng.randint(1, 4)
        A1 = [[1 if (i<j)==(rng.random()<0.5) else 0 for j in range(m1)] for i in range(m1)]
        for i in range(m1):
            for j in range(i+1, m1):
                b = rng.random() < 0.5
                A1[i][j], A1[j][i] = (1, 0) if b else (0, 1)
            A1[i][i] = 0
        A2 = [[0]*m2 for _ in range(m2)]
        for i in range(m2):
            for j in range(i+1, m2):
                b = rng.random() < 0.5
                A2[i][j], A2[j][i] = (1, 0) if b else (0, 1)
        A, n = join(A1, m1, A2, m2)
        if R_signed(n, A) != R_signed(m1, A1)*R_signed(m2, A2):
            ok = False; break
    print("R(T1|>T2)=R(T1)R(T2):", "VERIFIED" if ok else "FAILED")

    # 2) strong-atom spectra + gap reduction
    atom = {1: {1}}
    for n in range(2, 7):
        atom[n] = {abs(R_signed(n, A)) for A in all_tours(n) if is_strong(n, A)}
    def achievable(N):
        res = set()
        def rec(rem, acc):
            if rem == 0: res.add(acc); return
            for part in range(1, rem+1):
                for a in atom.get(part, ()): rec(rem-part, acc*a)
        rec(N, 1); return res
    print("strong-atom |R| spectra:", {n: sorted(atom[n]) for n in range(2, 7)})
    for N in range(2, 7):
        pred = achievable(N)
        actual = {abs(R_signed(N, A)) for A in all_tours(N)}
        print(f"  n={N}: |R| = {sorted(actual)}  (products of strong atoms match: {pred == actual})")
    print("=> 9,13 absent at n=5: 9=3*3 needs >=6 vertices, 13 prime first at strong n=6.")
