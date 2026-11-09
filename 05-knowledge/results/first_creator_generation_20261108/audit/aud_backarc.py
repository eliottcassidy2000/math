"""aud_backarc.py - check of the one-back-arc formula quoted in AUDIT.md m9:
for T = (A => B) with the single arc a -> b reversed (a in A, b in B),
H(T) = H(A) H(B) - e_A(a) s_B(b) + U_A(a) U'_B(b), where e = paths ending at, s = paths starting at,
U_A(a) = ordered pairs (P1, P2) of disjoint paths covering A, P2 nonempty starting at a (P1 may be empty),
U'_B(b) = ordered pairs (Q1, Q2) covering B, Q1 nonempty ending at b (Q2 may be empty).  Brute force by permutations."""
import itertools, random
def is_path(A, seq): return all(A[seq[i]][seq[i + 1]] for i in range(len(seq) - 1))
def H(A, V): return sum(1 for p in itertools.permutations(V) if is_path(A, p)) if V else 1
def count_pairs(A, V, root, second_starts):
    c = 0
    for k in range(len(V) + 1):
        for P in itertools.permutations(V, k):
            rest = [v for v in V if v not in P]
            for Q in itertools.permutations(rest):
                if not (is_path(A, P) and is_path(A, Q)): continue
                if second_starts:
                    if Q and Q[0] == root: c += 1          # (P1, P2) with P2 starting at root
                else:
                    if P and P[-1] == root: c += 1         # (Q1, Q2) with Q1 ending at root
    return c
random.seed(1); bad = 0
for trial in range(300):
    na, nb = random.randint(1, 4), random.randint(1, 4); n = na + nb
    M = [[0] * n for _ in range(n)]
    for i in range(n):
        for j in range(i + 1, n):
            if i < na and j >= na: M[i][j] = 1
            elif random.random() < 0.5: M[i][j] = 1
            else: M[j][i] = 1
    a, b = random.randrange(na), na + random.randrange(nb)
    M[a][b], M[b][a] = 0, 1
    Av, Bv = list(range(na)), list(range(na, n))
    e = sum(1 for p in itertools.permutations(Av) if is_path(M, p) and p[-1] == a)
    s = sum(1 for p in itertools.permutations(Bv) if is_path(M, p) and p[0] == b)
    rhs = H(M, Av) * H(M, Bv) - e * s + count_pairs(M, Av, a, True) * count_pairs(M, Bv, b, False)
    if rhs != H(M, list(range(n))): bad += 1
print("one-back-arc formula: 300 random cases, failures", bad)
