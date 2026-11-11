# audit-5129, step 18: items (iii), (iv), (v) of the brief WITHOUT any class generator: all labelled tournaments
# of orders 3..6 give every strong tournament with at most 23 paths (f(7) = 25 is read off all 2^21 labelled
# tournaments of order 7 here), hence the complete slope lists for offsets D <= 23, r(N) for odd N < 72 and the
# window counts for h <= 23.
import numpy as np
from numba import njit
import acore as A
import supply as S


@njit(cache=True)
def scan(n, Hout, strong, Pi2, Uarr):
    adj = np.zeros(n, np.int64)
    hp = np.zeros((1 << n, n), np.int64); sp = np.zeros((1 << n, n), np.int64); Hs = np.zeros(1 << n, np.int64)
    st = np.zeros(4, np.int64)
    for code in range(1 << (n * (n - 1) // 2)):
        A.adj_from_code(code, n, adj)
        strong[code] = 1 if A.is_strong(adj, n) else 0
        A.tables(adj, n, hp, sp, Hs)
        Hout[code] = Hs[(1 << n) - 1]
        Pi2[code] = A.pi2_from_tables(n, Hs)
        for y in range(n):
            A.root_stats(n, hp, sp, Hs, y, st)
            Uarr[code, y] = st[2] + st[3]


@njit(cache=True)
def min_strong_H(n):
    adj = np.zeros(n, np.int64)
    dp = np.zeros((1 << n, n), np.int64)
    best = 1 << 40
    cnt = 0
    for code in range(1 << (n * (n - 1) // 2)):
        A.adj_from_code(code, n, adj)
        if A.is_strong(adj, n):
            h = A.ham_count(adj, n, dp)
            if h < best:
                best = h
                cnt = 1
            elif h == best:
                cnt += 1
    return best, cnt


rows = []
for n in range(3, 7):
    m = 1 << (n * (n - 1) // 2)
    H = np.zeros(m, np.int64); st = np.zeros(m, np.int8); P = np.zeros(m, np.int64); U = np.zeros((m, n), np.int64)
    scan(n, H, st, P, U)
    seen = set()
    for c in np.nonzero(st)[0]:
        key = (int(H[c]), int(P[c]), tuple(sorted(int(x) for x in U[c])))
        if key not in seen:
            seen.add(key)
            rows.append((n, key[0], key[1], list(key[2])))
    print("order %d: %d labelled, %d strong labelled, least H of a strong one %d, distinct (H, Pi_2, multiset of U): %d"
          % (n, m, int(st.sum()), int(H[st == 1].min()), len(seen)))
f7, c7 = min_strong_H(7)
print("order 7 (all 2^21 labelled): least H of a strong tournament = %d (%d labelled copies)" % (f7, c7))
DMAX = f7 - 2
Root, PI, Aord, Suf, S0, minord = S.build(rows, DMAX)
DL = {D: S.delta(S0[D]) for D in S0}
print("complete for offsets D <= %d (strong components have order <= 6)" % DMAX)
print("attained odd numbers <= %d: %s" % (DMAX, sorted(Aord)))
for D in (3, 5, 9, 17):
    print("offset %d: %s" % (D, " ".join("%d:%d" % (d, DL[D][d]) for d in sorted(DL[D]))))
mine = {}
for line in open("s3_slope_lists.txt"):
    Dd, rest = line.split("|")
    mine[int(Dd)] = {int(t.split(":")[0]): int(t.split(":")[1]) for t in rest.split()}
print("all lists for D <= %d equal to the class-based lists of step 3: %s" %
      (DMAX, all(DL[D] == mine[D] for D in DL) and set(DL) == {D for D in mine if D <= DMAX}))
NM = 3 * DMAX + 2
r = {N: sum(1 for D in range(1, N // 3 + 1, 2) if D in DL and S.is_slope(DL[D], (N - D) // 2)) for N in range(3, NM + 1, 2)}
print("odd N in [3, %d] with r(N) = 0: %s" % (NM, [N for N in r if r[N] == 0]))
print("windows h <= %d: %s" % (DMAX, " ".join("%d:%d/%d" % (h, sum(1 for u in range(h, (3 * h) // 2 + 1) if S.is_slope(DL[h], u)),
                                                          (3 * h) // 2 - h + 1) for h in sorted(DL))))
