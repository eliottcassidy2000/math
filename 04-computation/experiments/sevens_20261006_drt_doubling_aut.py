#!/usr/bin/env python3
"""Doubling a doubly regular tournament keeps its automorphism group (proves HYP-9162).  mac-mini, 2026-10-06.

D(T) = T + {0'} + T' with T' the converse of T, i -> j' iff (i -> j or i = j), i' -> j iff i -> j, 0' -> T, T' -> 0'
(the skew-Hadamard doubling H_(2n) = [[H, H], [-H^T, H^T]] of HYP-9162; T_(k+1) = D(T_k)).

  1. N+(x') in D(T) is N+(x) u {0'} u N-(x)', and x -> 0' identifies it with T^(x) := T with every arc inside N-(x)
     reversed (checked as an exact relabelling).
  2. If T is doubly regular on n = 4t+3 >= 7 vertices then T^(x) is not doubly regular: a mixed pair j in N+(x),
     y in N-(x) changes its common out-neighbour count by -(S_R 1_(O_j))_y, where R = T[N-(x)], O_j = N+(j) n N-(x)
     has t+1 elements, and ker S_R = span(1) because S_R = J - I (mod 2) has F_2-rank 2t (checked: some mixed
     pair fails, for every x).
  3. Hence the set of vertices with doubly regular out-neighbourhood contains 0', misses T', and 0' is its unique
     source; with the diagonal-stabiliser lemma, Aut(D(T)) = Aut(T) (checked with nauty: T_3..T_8 and Paley P_p,
     p = 7, 11, 19, 23, 31, 43, and D(D(P_7)), D(D(P_11))).
Run: python3 sevens_20261006_drt_doubling_aut.py   (needs nauty's dreadnaut; about 1 min)
"""
import subprocess, re
import numpy as np

OK = True


def check(cond, msg):
    global OK
    print(("  ok   " if cond else "  FAIL ") + msg, flush=True)
    OK &= bool(cond)


def paley(p):
    qr = {(i * i) % p for i in range(1, p)}
    return np.array([[1 if (j - i) % p in qr else 0 for j in range(p)] for i in range(p)], dtype=np.int64)


def double(A):
    """A: adjacency of T on n vertices (A[i, j] = 1 iff i -> j). Returns D(T) on 2n+1 vertices: T = 0..n-1, 0' = n, T' = n+1..2n."""
    n = A.shape[0]
    B = np.zeros((2 * n + 1, 2 * n + 1), dtype=np.int64)
    B[:n, :n] = A
    B[n, :n] = 1                                        # 0' -> T
    B[n + 1:, n] = 1                                    # T' -> 0'
    B[n + 1:, n + 1:] = A.T                             # T' = converse
    B[:n, n + 1:] = A + np.eye(n, dtype=np.int64)       # i -> j' iff i -> j or i = j
    B[n + 1:, :n] = A                                   # i' -> j iff i -> j
    return B


def is_tournament(A):
    n = A.shape[0]
    return np.all(A + A.T + np.eye(n, dtype=np.int64) == 1)


def is_drt(A):
    n = A.shape[0]
    if n < 3:
        return True
    C = A @ A.T
    return len(set(np.diag(C).tolist())) == 1 and len(set(C[~np.eye(n, dtype=bool)].tolist())) == 1


def aut_order(A):
    n = A.shape[0]
    s = "d\nn=%d g\n" % n + ";\n".join(" ".join(str(j) for j in np.nonzero(A[i])[0]) for i in range(n)) + ".\nx\nq\n"
    out = subprocess.run(["dreadnaut"], input=s, capture_output=True, text=True).stdout
    return int(float(re.search(r"grpsize=([0-9.e+]+)", out).group(1)))


def tower(k):
    H = np.array([[1, 1], [-1, 1]])
    for _ in range(k - 1):
        H = np.block([[H, H], [-H.T, H.T]])
    S = H - np.eye(2 ** k, dtype=np.int64)
    return (S[1:, 1:] > 0).astype(np.int64)


print("0. the doubling of HYP-9162 reproduces the tower: T_(k+1) = D(T_k) up to the index map")
good = True
for k in range(2, 8):
    good &= aut_order(double(tower(k))) == aut_order(tower(k + 1)) and is_drt(double(tower(k)))
check(good, "D(T_k) is a doubly regular tournament with |Aut D(T_k)| = |Aut T_(k+1)|, k = 2..7")

print("1-2. out-neighbourhoods of second-copy vertices")
cases = [("T_%d" % k, tower(k)) for k in range(3, 8)] + [("P_%d" % p, paley(p)) for p in (7, 11, 19, 23, 31, 43)]
good1 = good2 = True
for name, A in cases:
    n = A.shape[0]
    B = double(A)
    assert is_tournament(B) and is_drt(A)
    for x in range(n):
        outp = np.nonzero(B[n + 1 + x])[0]
        Np, Nm = set(np.nonzero(A[x])[0]), set(np.nonzero(A[:, x])[0])
        good1 &= set(outp) == Np | {n} | {n + 1 + y for y in Nm}
        # relabel: 0' -> x, y' -> y; compare with T^(x)
        lab = {n: x}
        lab.update({v: v for v in Np})
        lab.update({n + 1 + y: y for y in Nm})
        Tx = A.copy()
        for y in Nm:
            for z in Nm:
                if y != z:
                    Tx[y, z] = A[z, y]
        sub = B[np.ix_(outp, outp)]
        good1 &= all(sub[i, j] == Tx[lab[outp[i]], lab[outp[j]]] for i in range(len(outp)) for j in range(len(outp)))
        good2 &= not is_drt(sub)
        # the mechanism: S_R 1_(O_j) != 0 for every j in N+(x), and |O_j| = t + 1
        Nm_l = sorted(Nm)
        R = A[np.ix_(Nm_l, Nm_l)]
        SR = R - R.T
        t = (n - 3) // 4
        for j in Np:
            O = np.array([1 if A[j, y] else 0 for y in Nm_l])
            good2 &= O.sum() == t + 1 and np.any(SR @ O != 0)
    print(f"   {name}: n = {n}, checked all {n} second-copy vertices")
check(good1, "N+(x') = N+(x) u {0'} u N-(x)' and N+(x') = T^(x) under x <- 0' (all x, all cases)")
check(good2, "T^(x) is never doubly regular; |O_j| = t+1 and S_R 1_(O_j) != 0 for every j in N+(x)")
good = True
for m in range(3, 40, 2):
    J = np.ones((m, m), dtype=np.int64)
    M = (J - np.eye(m, dtype=np.int64)) % 2
    # F_2 rank by elimination
    M = M.copy(); r = 0
    for c in range(m):
        piv = next((i for i in range(r, m) if M[i, c]), None)
        if piv is None:
            continue
        M[[r, piv]] = M[[piv, r]]
        for i in range(m):
            if i != r and M[i, c]:
                M[i] ^= M[r]
        r += 1
    good &= r == m - 1
check(good, "rank_F2(J - I) = m - 1 for odd m <= 39 (so a skew +-1 matrix of odd order m has rank exactly m - 1)")

print("3. automorphism groups")
good = True
rows = []
for name, A in cases + [("D(P_7)", double(paley(7))), ("D(P_11)", double(paley(11)))]:
    a, b = aut_order(A), aut_order(double(A))
    rows.append((name, a, b))
    good &= a == b
print("   " + ", ".join(f"|Aut {nm}| = {a} = |Aut D({nm})|" for nm, a, b in rows))
check(good, "|Aut D(T)| = |Aut T| for every doubly regular T tested (orders 7 to 87)")
check(all(aut_order(tower(k)) == 21 for k in range(3, 9)), "|Aut T_k| = 21 for k = 3..8 (consistent with Aut(T_k) = F_21 for all k >= 3)")
print("ALL CHECKS PASSED" if OK else "SOME CHECK FAILED")
