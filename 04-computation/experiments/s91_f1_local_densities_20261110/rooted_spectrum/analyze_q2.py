# analyze_q2.py - the rooted spectrum (question Q2): triples (x, s, Y) = (H(T), s(T, v), Yplus(T, v)).
#
# Input: rooted_le9.npz, agg10.npz, reps_le8.npz; klein's census data
# 05-knowledge/results/sole_exit_tower_tribonacci_20261104/data/triples_n.txt (n = 4..8) for a cross-check.
# Output: stdout (q2.out), q2_strong_pairs.txt (every strong rooted pair for odd x <= 123: complete at all orders,
# because a strong tournament of order >= 10 has H >= 125), q2_pairs_le8.txt (all pairs for x <= 61, orders <= 8).
import os
import sys
from collections import defaultdict
import numpy as np
import core

HERE = os.path.dirname(os.path.abspath(__file__))
ROOT = os.path.join(HERE, "..", "..", "..", "..")
FAIL = []


def check(name, ok, detail=""):
    print(("ok   " if ok else "FAIL ") + name + ((" : " + detail) if detail else ""), flush=True)
    if not ok:
        FAIL.append(name)


def main():
    z = np.load(os.path.join(HERE, "rooted_le9.npz"))
    a10 = np.load(os.path.join(HERE, "agg10.npz"))
    reps = np.load(os.path.join(HERE, "reps_le8.npz"))
    rows = {m: z["m%d" % m] for m in range(1, 10)}

    # ---- 1. cross-check against klein's census (h, e_t, B_t, strong) with multiplicities, orders 4..8
    kdir = os.path.join(ROOT, "05-knowledge", "results", "sole_exit_tower_tribonacci_20261104", "data")
    direct = {}
    for m in range(1, 9):
        buf = np.zeros((len(reps["n%d" % m]) * m, 10), np.int64)
        r = core.rooted_rows_direct(m, reps["n%d" % m], buf)
        direct[m] = buf[:r].copy()
    okk = True
    nk = 0
    for m in range(4, 9):
        fn = os.path.join(kdir, "triples_%d.txt" % m)
        if not os.path.exists(fn):
            okk = False
            break
        theirs = {}
        for line in open(fn):
            p = line.split()
            if len(p) == 5:
                theirs[(int(p[0]), int(p[1]), int(p[2]), int(p[3]))] = int(p[4])
        mine = defaultdict(int)
        for r in direct[m]:
            mine[(int(r[0]), int(r[2]), int(r[4]), int(r[6]))] += 1
        nk += len(theirs)
        if dict(mine) != theirs:
            okk = False
    check("rows (h, e, Yminus, strong) with multiplicities = klein's triples_n.txt, orders 4..8 (%d rows)" % nk, okk)
    tot = sum(len(np.unique(rows[m][:, [0, 2, 4, 6]], axis=0)) for m in range(4, 10))
    check("distinct (h, e, Yminus, strong) rows over orders 4..9 = 1,728,779 (klein's receipt coverage_9)",
          tot == 1728779, str(tot))

    # ---- 2. strong floors and the composition law
    SR = {}
    ALL = {}
    for m in range(1, 10):
        R = rows[m]
        ALL[m] = set(map(tuple, np.unique(R[:, [0, 1, 3]], axis=0).tolist()))
        SR[m] = set(map(tuple, np.unique(R[R[:, 6] == 1][:, [0, 1, 3]], axis=0).tolist()))
    floors = [min(t[0] for t in SR[m]) for m in range(3, 10)] + [int(a10["mx"][2])]
    check("least H of a strong tournament of order 3..10 = 3, 5, 9, 15, 25, 45, 75, 125",
          floors == [3, 5, 9, 15, 25, 45, 75, 125], str(floors))
    HP = {0: {(1, 1)}}
    HZ = {0: {1}}
    for m in range(1, 9):
        D = direct[m]
        HP[m] = set(zip(D[:, 0].tolist(), (2 * D[:, 0] + D[:, 8]).tolist()))
        HZ[m] = set(D[:, 0].tolist())
    okc = True
    for M in range(1, 10):
        G = set()
        for a in range(0, M):
            for m in range(1, M - a + 1):
                zz = M - a - m
                for (hA, pA) in HP[a]:
                    for hZ in HZ[zz]:
                        if a == 0:
                            G |= {(hA * x * hZ, s * hZ, pA * y * hZ) for (x, s, y) in SR[m]}
                        else:
                            G |= {(hA * x * hZ, 0, pA * y * hZ) for (x, s, y) in SR[m]}
        if G != ALL[M]:
            okc = False
            print("   composition law fails at order", M, len(G), len(ALL[M]))
    check("composition law: for T = A => C => Z (C the strong component of the root), "
          "(x, s, Y) = (H(A) x_C H(Z), [A empty] s_C H(Z), pi(A) Y_C H(Z)), pi = 2H + P2; "
          "generates exactly the rooted triples of every order 1..9", okc)
    pi_small = sorted(HP[3] | HP[4])
    check("pi(C3) = 12 and pi(strong 4-tournament) = 32", (3, 12) in HP[3] and (5, 32) in HP[4], str(pi_small))

    # ---- 3. shape of the strong spectrum
    strong = defaultdict(set)          # x -> {(s, Y)}
    first_order = {}
    for m in range(1, 10):
        for (x, s, y) in SR[m]:
            strong[x].add((s, y))
            first_order.setdefault((x, s, y), m)
    allp = defaultdict(set)
    for m in range(1, 9):
        for (x, s, y) in ALL[m]:
            if x <= 61:
                allp[x].add((s, y))
    with open(os.path.join(HERE, "q2_strong_pairs.txt"), "w") as fo:
        fo.write("# strong rooted pairs (s, Yplus) for every odd x = H(T) <= 123 (all orders: strong tournaments of "
                 "order >= 10 have H >= 125)\n")
        for x in sorted(strong):
            if x <= 123:
                fo.write("x=%d n=%d: %s\n" % (x, len(strong[x]), " ".join("(%d,%d)" % p for p in sorted(strong[x]))))
    with open(os.path.join(HERE, "q2_pairs_le8.txt"), "w") as fo:
        fo.write("# all rooted pairs (s, Yplus) for odd x = H(T) <= 61, rooted tournaments of order <= 8\n")
        for x in sorted(allp):
            fo.write("x=%d n=%d: %s\n" % (x, len(allp[x]), " ".join("(%d,%d)" % p for p in sorted(allp[x]))))
    print("strong rooted pairs (complete for x <= 123); D = Yplus - s (pendant increment); x: #pairs, s-range, "
          "Y-range, D-set")
    viol = 0
    for x in sorted(strong):
        if x > 61:
            break
        P = sorted(strong[x])
        D = sorted(set(y - s for s, y in P))
        print("  x=%2d: %3d pairs, s in [%d, %d], Y in [%d, %d], D = %s"
              % (x, len(P), min(p[0] for p in P), max(p[0] for p in P), min(p[1] for p in P),
                 max(p[1] for p in P), D if len(D) <= 14 else str(D[:12])[:-1] + ", ... %d values, max %d]" % (len(D), D[-1])))
    for x in sorted(strong):
        if x > 123 or x < 3:
            continue
        for (s, y) in strong[x]:
            if not (1 <= s <= (x - 1) // 2 and y >= x and (y - s) % 2 == 0 and 2 * (y - s) >= x + 1
                    and (x < 5 or 2 * (y - s) >= x + 3)):
                viol += 1
    check("strong rows, 3 <= x <= 123: 1 <= s <= (x-1)/2, Y >= x, Y = s mod 2, Y - s >= (x+1)/2 "
          "(>= (x+3)/2 for x >= 5)", viol == 0)
    eq = sorted(x for x in strong if 5 <= x <= 123 and any(2 * (y - s) == x + 3 for s, y in strong[x]))
    check("the floor Y - s = (x+3)/2 is attained exactly at x = 2^k + 1 (5, 9, 17, 33, 65), x <= 123",
          eq == [5, 9, 17, 33, 65], str(eq))
    smax = sorted(x for x in strong if 3 <= x <= 123 and any(2 * s == x - 1 for s, y in strong[x]))
    check("s = (x-1)/2 is attained exactly at x = 2^k + 1 (3, 5, 9, 17, 33, 65), x <= 123",
          smax == [3, 5, 9, 17, 33, 65], str(smax))
    nos1 = sorted(x for x in strong if x <= 123 and all(s != 1 for s, y in strong[x]))
    print("strong x <= 123 with no pair s = 1: %s" % nos1)
    ratio = max((y / x, x, s, y) for x in strong if x <= 123 for s, y in strong[x])
    print("largest Y/x among strong pairs with x <= 123: %.2f at (x, s, Y) = (%d, %d, %d)" % ratio)
    npairs = sum(len(strong[x]) for x in strong if x <= 123)
    print("number of strong pairs with x <= 123: %d; with x <= 61: %d" % (npairs, sum(len(strong[x]) for x in strong if x <= 61)))

    # ---- 4. how 7 and 21 look at the rooted level
    print()
    print("pendant extension: H(T + w) = x - s + Y.  Triples that would give 7: (1,0,6), (3,s,s+4), (5,s,s+2);")
    print("21: (1,0,20), (3,s,s+18), (5,s,s+16), (9,s,s+12), (11,s,s+10), (13,s,s+8), (15,s,s+6), (17,s,s+4), (19,s,s+2).")
    Dstr = {x: sorted(set(y - s for s, y in strong[x])) for x in (1, 3, 5, 9, 11, 13, 15, 17, 19)}
    for x in (3, 5, 9, 11, 13, 15, 17, 19):
        print("   strong x = %2d: D = %s; 7 needs %s, 21 needs %d" % (x, Dstr[x], (7 - x) if x < 7 else "-", 21 - x))
    need = {3: (4, 18), 5: (2, 16), 9: (None, 12), 11: (None, 10), 13: (None, 8), 15: (None, 6), 17: (None, 4),
            19: (None, 2)}
    ok7 = all((need[x][0] is None or need[x][0] not in Dstr[x]) and need[x][1] not in Dstr[x] for x in need)
    check("no strong rooted tournament has x - s + Y in {7, 21} (finite list, x <= 19)", ok7)
    Y3 = sorted(set(y for s, y in strong[3]))
    Y5 = sorted(set(y for s, y in strong[5]))
    check("Yplus on the 3-cycle is {3}; on the strong 4-tournament {5, 6, 7, 9}", Y3 == [3] and Y5 == [5, 6, 7, 9])
    # non-strong, root in the last component: N = x + 2^a pi-products ... ; exhaustive over orders <= 10
    pend_le9 = set()
    per_order = {}
    for m in range(1, 10):
        R = rows[m]
        v = set((R[:, 0] - R[:, 1] + R[:, 3]).tolist())
        pend_le9 |= v
        per_order[m] = set(pend_le9)
    p10 = set(np.nonzero(a10["pend"])[0].tolist())
    per_order[10] = pend_le9 | p10
    check("7 and 21 are not pendant values of any rooted tournament of order <= 10",
          7 not in per_order[10] and 21 not in per_order[10])
    print("pendant spectrum {x - s + Y} = {H(T') : T' has a vertex of in-degree one}, bases of order <= n:")
    tops = {}
    for m in (5, 6, 7, 8, 9, 10):
        v = 23
        while v in per_order[m]:
            v += 2
        tops[m] = v - 2
        small = [u for u in range(1, 23, 2) if u not in per_order[m]]
        print("   n <= %2d: every odd number of [23, %d] occurs; missing below 23: %s" % (m, v - 2, small))
    check("pendant spectrum: [1, 393] at order <= 7, [1, 1949] at <= 8, [1, 9151] at <= 9, [1, 48173] at <= 10, "
          "minus {7, 21}", [tops[m] for m in (7, 8, 9, 10)] == [393, 1949, 9151, 48173])
    spend = set()
    for m in range(1, 10):
        R = rows[m]
        R = R[R[:, 6] == 1]
        spend |= set((R[:, 0] - R[:, 1] + R[:, 3]).tolist())
    v = 23
    while v in spend:
        v += 2
    print("   strong rooted tournaments of order <= 9 only: every odd number of [23, %d]; missing below 23: %s"
          % (v - 2, [u for u in range(1, 23, 2) if u not in spend]))

    # ---- 5. the doubling gadget H + 2 Yminus (dually H + 2 Yplus)
    xy = defaultdict(set)
    for m in range(1, 10):
        R = rows[m]
        R = R[R[:, 0] <= 21]
        for x, y in set(zip(R[:, 0].tolist(), R[:, 3].tolist())):
            xy[x].add(y)
    print("doubling gadget H + 2Y: 7 needs (x, Y) in {(1,3)}; 21 needs (1,10), (3,9), (5,8) (Y >= x).")
    print("   Y-values (orders <= 9): x = 1: %s; x = 3: %s; x = 5: %s"
          % (sorted(xy[1]), sorted(xy[3]), sorted(xy[5])))
    check("none of (1,3), (1,10), (3,9), (5,8) is a pair (x, Y)",
          3 not in xy[1] and 10 not in xy[1] and 9 not in xy[3] and 8 not in xy[5])
    two = set()
    for m in range(1, 10):
        R = rows[m]
        two |= set((R[:, 0] + 2 * R[:, 3]).tolist())
    v = 69
    while v in two:
        v += 2
    miss = [u for u in range(1, 69, 2) if u not in two]
    check("two-step values x + 2Y, orders <= 9: every odd number of [69, 15449]; below 69 the misses are "
          "1, 7, 11, 13, 21, 47, 67 (klein's receipt)", v - 2 == 15449 and miss == [1, 7, 11, 13, 21, 47, 67],
          "%d %s" % (v - 2, miss))
    if FAIL:
        print("FAILED:", FAIL)
        sys.exit(1)
    print("ALL CHECKS PASSED")


if __name__ == "__main__":
    main()
