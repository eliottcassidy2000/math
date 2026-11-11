# audit-5129, step 6: literal recount of the level 7, 8, 9 sources printed at the end of THM-5129.
# Out-masks: bit j of row i set iff i -> j.
import numpy as np
import acore as A

SRC = [
    (7, (444, 273, 258, 6, 268, 286, 319, 382, 8), 8, (129, 144)),
    (8, (82, 36, 129, 135, 142, 157, 190, 3), 7, (257, 288)),
    (9, (402, 360, 3, 261, 270, 285, 317, 382, 4), 8, (513, 576)),
]
for s, masks, y, claim in SRC:
    n = len(masks)
    adj = np.array(masks, dtype=np.int64)
    # tournament check
    ok = True
    for i in range(n):
        if (adj[i] >> i) & 1 or adj[i] >> n:
            ok = False
        for j in range(i + 1, n):
            if ((adj[i] >> j) & 1) == ((adj[j] >> i) & 1):
                ok = False
    perms = A.perm_array(n)
    E = np.zeros(n, np.int64); Ue = np.zeros(n, np.int64); Uo = np.zeros(n, np.int64)
    Hlit = A.literal_stats(adj, n, perms, E, Ue, Uo)
    hp = np.zeros((1 << n, n), np.int64); sp = np.zeros((1 << n, n), np.int64); Hs = np.zeros(1 << n, np.int64)
    st = np.zeros(4, np.int64)
    A.tables(adj, n, hp, sp, Hs)
    A.root_stats(n, hp, sp, Hs, y, st)
    D = int(st[0]); U = int(st[2] + st[3])
    print("s = %d: order %d, valid tournament: %s, strong: %s, scores %s" %
          (s, n, ok, A.is_strong(adj, n), [A.popcount(int(a)) for a in adj]))
    print("   literal (all %d permutations): H(B) = %d, E(y) = %d, U_even = %d, U_odd = %d, U = %d" %
          (perms.shape[0], Hlit, E[y], Ue[y], Uo[y], Ue[y] + Uo[y]))
    print("   subset programme:            H(B) = %d, E(y) = %d, U = %d ; claimed (D, U) = %s ; match: %s" %
          (D, st[1], U, claim, (Hlit, int(Ue[y] + Uo[y])) == claim and (D, U) == claim))
    assert 2 ** s + 1 == D
    v = 0
    q = U
    while q % 2 == 0:
        q //= 2
        v += 1
    print("   U = %d * 2^%d ; grid: N(9, r, %d) = 1 + 2^%d + 9 * 2^(r+%d) needs t = r + %d - %d >= 1" %
          (q, v, s, s, s, s, v))
    cadj = np.zeros(16, np.int64)
    # head
    A.head_adj(adj, n, y, cadj)
    dp = np.zeros((1 << (n + 1), n + 1), np.int64)
    Hh = A.ham_count(cadj, n + 1, dp)
    print("   head H(B + z) = %d ; D - E + U = %d ; D + 2 U_even = %d" % (Hh, D - int(st[1]) + U, D + 2 * int(st[2])))
    for t in (1, 2, 3):
        m = n + t + 1
        A.tail_adj(adj, n, y, t, cadj)
        dp = np.zeros((1 << m, m), np.int64)
        Ht = A.ham_count(cadj, m, dp)
        extra = ""
        if m <= 10:
            extra = " ; literal over %d! permutations: %d" % (m, A.literal_H(cadj, m, A.perm_array(m)))
        print("   t = %d: H(closed tournament, order %d) = %d ; D + 2^t U = %d ; equal: %s%s" %
              (t, m, Ht, D + (1 << t) * U, Ht == D + (1 << t) * U, extra))
    # all pairs of this base
    allp = []
    for yy in range(n):
        A.root_stats(n, hp, sp, Hs, yy, st)
        allp.append(int(st[2] + st[3]))
    print("   slopes at all roots of this base:", allp)
