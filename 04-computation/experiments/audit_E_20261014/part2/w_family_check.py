#!/usr/bin/env python3
"""Audit E part 2 (2026-10-15): THM-4655 statement 6, independent check by DFS enumeration of Hamiltonian paths.
W_n (n >= 5, t = n - 3): vertices 0..t-1 (transitive, i -> j for i < j), a, b, c with a -> b -> c -> a; every i -> {a,b,c}
except the single reversed arc b -> 0.  Checks for n = 5..14:
  strong; H = 3(2^(n-4)+1); the induced order on {a,b,c} uniform over abc, bca, cab (and never acb, bac, cba);
  p(c<a) = 2/3, p(c<b) = 1/3, p(c<0) = 2^(t-1)/H, p(c<i) = 2^(t-2)/H (1 <= i <= t-1); vertex balance of c exactly 1/3;
  W_n - b transitive in the order 0..t-1, c, a; W_n equals (as labelled tournaments, via u_i -> i-1, u_(t+1) -> c,
  u_(t+2) -> a, z -> b) the transitive-plus-one tournament with word 1 0^(m-3) 1 0 (m = n-1), whose reverse, relabelled
  u_i -> u_(m+1-i), is the one with the equality word 1 0 1^(m-3) 0; the vertex-wise constant of W_n is exactly 1/3."""
from fractions import Fraction as F

def build(n):
    t = n - 3; a, b, c = t, t + 1, t + 2
    A = [[0] * n for _ in range(n)]
    for i in range(t):
        for j in range(i + 1, t): A[i][j] = 1
        for x in (a, b, c): A[i][x] = 1
    A[a][b] = A[b][c] = A[c][a] = 1
    A[0][b] = 0; A[b][0] = 1
    return A, t, a, b, c

def paths(A):
    n = len(A); res = []
    def rec(p, used):
        if len(p) == n: res.append(tuple(p)); return
        for v in range(n):
            if not used >> v & 1 and A[p[-1]][v]:
                p.append(v); rec(p, used | 1 << v); p.pop()
    for s in range(n): rec([s], 1 << s)
    return res

def tpo(word):
    m = len(word); n = m + 1; Z = m
    A = [[0] * n for _ in range(n)]
    for i in range(m):
        for j in range(i + 1, m): A[i][j] = 1
        if word[i]: A[Z][i] = 1
        else: A[i][Z] = 1
    return A

allok = True
for n in range(5, 15):
    A, t, a, b, c = build(n)
    assert all(A[i][j] + A[j][i] == 1 for i in range(n) for j in range(n) if i != j)
    P = paths(A); H = len(P)
    # strong
    R = [[bool(A[i][j]) or i == j for j in range(n)] for i in range(n)]
    for k in range(n):
        for i in range(n):
            if R[i][k]:
                for j in range(n):
                    if R[k][j]: R[i][j] = True
    strong = all(all(r) for r in R)
    orders = {}
    for p in P:
        o = ''.join({a: 'a', b: 'b', c: 'c'}[x] for x in p if x in (a, b, c))
        orders[o] = orders.get(o, 0) + 1
    pos = [{v: k for k, v in enumerate(p)} for p in P]
    pr = lambda x, y: F(sum(1 for q in pos if q[x] < q[y]), H)
    ok = (strong and H == 3 * (2 ** (n - 4) + 1) and set(orders) == {'abc', 'bca', 'cab'} and len(set(orders.values())) == 1
          and pr(c, a) == F(2, 3) and pr(c, b) == F(1, 3) and pr(c, 0) == F(2 ** (t - 1), H)
          and all(pr(c, i) == F(2 ** (t - 2), H) for i in range(1, t)))
    bal = lambda x, y: min(pr(x, y), 1 - pr(x, y))
    cbal = max(bal(c, w) for w in range(n) if w != c)
    vw = min(max(bal(v, w) for w in range(n) if w != v) for v in range(n))
    # W_n - b transitive in the order 0..t-1, c, a
    order = list(range(t)) + [c, a]
    trans = all(A[order[i]][order[j]] for i in range(len(order)) for j in range(i + 1, len(order)))
    # identification with the transitive-plus-one word 1 0^(m-3) 1 0, m = n - 1 (labels: u_1..u_t = 0..t-1, u_(t+1) = c, u_(t+2) = a, z = b)
    m = n - 1; w1 = [1] + [0] * (m - 3) + [1, 0]
    T1 = tpo(w1); lab = list(range(t)) + [c, a] + [b]      # T1 vertex index -> W_n vertex
    same = all(T1[i][j] == A[lab[i]][lab[j]] for i in range(n) for j in range(n) if i != j)
    # reverse of the equality-word tournament, relabelled u_i -> u_(m+1-i), equals T1
    we = [1, 0] + [1] * (m - 3) + [0]
    Te = tpo(we); rev = [[Te[j][i] for j in range(n)] for i in range(n)]
    relab = [m - 1 - i for i in range(m)] + [m]           # u_i (index i-1) -> u_(m+1-i) (index m-i)
    rev_ok = all(rev[relab[i]][relab[j]] == T1[i][j] for i in range(n) for j in range(n) if i != j)
    allok &= ok and cbal == F(1, 3) and vw == F(1, 3) and trans and same and rev_ok
    print(f"n={n}: strong={strong} H={H} rotations={orders} p(c<a)={pr(c,a)} p(c<b)={pr(c,b)} p(c<0)={pr(c,0)} "
          f"p(c<1)={pr(c,1) if t > 1 else '-'} | c's balance={cbal} vertex-wise constant={vw} | W-b transitive={trans} "
          f"= tpo(1 0^(m-3) 1 0)={same} = reverse of tpo(1 0 1^(m-3) 0)={rev_ok} | all claims {'ok' if ok else 'FAIL'}")
print("ALL OK" if allok else "SOME CHECK FAILED")
