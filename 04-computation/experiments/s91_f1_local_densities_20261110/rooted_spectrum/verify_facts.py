# verify_facts.py - brute-force verification of the facts F1-F4 and of the kernels of core.py.
#
# Literal layer: every invariant of the base is obtained by enumerating all orderings of the vertices
# (no subset DP).  Children T + w of order <= 6 are also counted literally; children of order 7 by an
# independent Held-Karp count.
#   F1  s(T + w, w) = sum of s(T, u) over the out-neighbours u of w
#   F2  w loses only to b:  H(T + w) = H(T) - s(T, b) + Yplus(T, b)
#       (dual) w beats only a:  H(T + w) = H(T) - e(T, a) + Yminus(T, a)
#   F3  one-arc join A => B with the arc a -> b replaced by b -> a:
#       H = H(A) H(B) - e(A, a) s(B, b) + Yminus(A, a) Yplus(B, b)
#   F4  THM-4097 (4): H(T + x_S) = sum_{b in S} s_b + sum_{a notin S} e_a + sum_{a notin S, b in S} Q_ab
# Universe: every labelled tournament T of order 1..6 (1 + 2 + 8 + 64 + 1024 + 32768) and every cut;
# F3: every labelled pair (A, B) with |A| + |B| <= 7, every (a, b).
import sys
import numpy as np
from numba import njit
import core

FAIL = []


def check(name, ok, detail=""):
    print(("ok   " if ok else "FAIL ") + name + ((" : " + detail) if detail else ""))
    if not ok:
        FAIL.append(name)


@njit(cache=True)
def literal(n, out, perms, s, e, yp, ym, Q):
    # returns H; fills s, e, Yplus, Yminus, Q by enumeration of all orderings
    for v in range(n):
        s[v] = 0
        e[v] = 0
        yp[v] = 0
        ym[v] = 0
        for w in range(n):
            Q[v, w] = 0
    H = 0
    for p in range(perms.shape[0]):
        nbad = 0
        pos = -1
        for i in range(n - 1):
            if ((out[perms[p, i]] >> perms[p, i + 1]) & 1) == 0:
                nbad += 1
                pos = i
        if nbad == 0:
            H += 1
            s[perms[p, 0]] += 1
            e[perms[p, n - 1]] += 1
            # every split is admissible
            for k in range(1, n + 1):          # P1 = first k vertices, P2 possibly empty
                yp[perms[p, k - 1]] += 1
            for k in range(0, n):              # P2 = last n - k vertices, P1 possibly empty
                ym[perms[p, k]] += 1
            for k in range(1, n):              # exposed slot between positions k-1, k
                Q[perms[p, k - 1], perms[p, k]] += 1
        elif nbad == 1:
            k = pos + 1                        # the only admissible split
            yp[perms[p, k - 1]] += 1
            ym[perms[p, k]] += 1
            Q[perms[p, k - 1], perms[p, k]] += 1
    return H


@njit(cache=True)
def literal_h_start(n, out, perms, w):
    H = 0
    sw = 0
    for p in range(perms.shape[0]):
        ok = True
        for i in range(n - 1):
            if ((out[perms[p, i]] >> perms[p, i + 1]) & 1) == 0:
                ok = False
                break
        if ok:
            H += 1
            if perms[p, 0] == w:
                sw += 1
    return H, sw


@njit(cache=True)
def hk(n, out):
    # independent Held-Karp count: number of Hamiltonian paths, and of those starting at the last vertex
    size = 1 << n
    dp = np.zeros((size, n), np.int64)
    for v in range(n):
        dp[1 << v, v] = 1
    for S in range(1, size):
        for v in range(n):
            c = dp[S, v]
            if c == 0:
                continue
            m = out[v] & (~S) & (size - 1)
            for w in range(n):
                if (m >> w) & 1:
                    dp[S | (1 << w), w] += c
    H = 0
    for v in range(n):
        H += dp[size - 1, v]
    # paths starting at vertex n-1: count by reversing the roles
    dq = np.zeros((size, n), np.int64)
    dq[1 << (n - 1), n - 1] = 1
    for S in range(1, size):
        for v in range(n):
            c = dq[S, v]
            if c == 0:
                continue
            m = out[v] & (~S) & (size - 1)
            for w in range(n):
                if (m >> w) & 1:
                    dq[S | (1 << w), w] += c
    sw = 0
    for v in range(n):
        sw += dq[size - 1, v]
    return H, sw


@njit(cache=True)
def run_f124(n, perms_n, perms_c, use_literal_child, cnt):
    # all labelled bases of order n, all cuts; cnt = (bases, cuts, F1 fails, F2 fails, F2dual fails, F4 fails,
    #                                                 kernel fails)
    npairs = n * (n - 1) // 2
    out = np.zeros(n, np.int64)
    outc = np.zeros(n + 1, np.int64)
    s = np.zeros(n, np.int64)
    e = np.zeros(n, np.int64)
    yp = np.zeros(n, np.int64)
    ym = np.zeros(n, np.int64)
    Q = np.zeros((n, n), np.int64)
    size = 1 << n
    f = np.zeros((size, n), np.int64)
    g = np.zeros((size, n), np.int64)
    hS = np.zeros(size, np.int64)
    res = np.zeros((n, 8), np.int64)
    s2 = np.zeros(n, np.int64)
    e2 = np.zeros(n, np.int64)
    yp2 = np.zeros(n, np.int64)
    ym2 = np.zeros(n, np.int64)
    Q2 = np.zeros((n, n), np.int64)
    full = size - 1
    for code in range(1 << npairs):
        core.out_masks(code, n, out)
        H = literal(n, out, perms_n, s, e, yp, ym, Q)
        cnt[0] += 1
        # kernels of core.py against the literal layer
        core.dp_tables(n, out, f, g, hS)
        core.rooted_direct(n, f, g, hS, res)
        core.base_data(n, f, g, hS, s2, e2, yp2, ym2, Q2)
        bad = 0
        if hS[full] != H:
            bad += 1
        for v in range(n):
            if res[v, 1] != s[v] or res[v, 2] != e[v] or res[v, 3] != yp[v] or res[v, 4] != ym[v]:
                bad += 1
            if s2[v] != s[v] or e2[v] != e[v] or yp2[v] != yp[v] or ym2[v] != ym[v]:
                bad += 1
            for w in range(n):
                if Q2[v, w] != Q[v, w]:
                    bad += 1
        cnt[6] += bad
        for S in range(size):
            cnt[1] += 1
            for v in range(n):
                outc[v] = out[v]
                if ((S >> v) & 1) == 0:
                    outc[v] |= (1 << n)
            outc[n] = S
            if use_literal_child:
                Hc, sw = literal_h_start(n + 1, outc, perms_c, n)
            else:
                Hc, sw = hk(n + 1, outc)
            # F1
            t = 0
            for v in range(n):
                if (S >> v) & 1:
                    t += s[v]
            if t != sw:
                cnt[2] += 1
            # F4
            t4 = 0
            for a in range(n):
                if (S >> a) & 1:
                    t4 += s[a]
                else:
                    t4 += e[a]
                    for b in range(n):
                        if (S >> b) & 1:
                            t4 += Q[a, b]
            if t4 != Hc:
                cnt[5] += 1
            # F2: w loses only to b  <=>  S = full minus b
            for b in range(n):
                if S == (full ^ (1 << b)):
                    if Hc != H - s[b] + yp[b]:
                        cnt[3] += 1
                if S == (1 << b):
                    if Hc != H - e[b] + ym[b]:
                        cnt[4] += 1


@njit(cache=True)
def run_f3(na, nb, perms_a, perms_b, perms_c, use_literal, cnt):
    pa = na * (na - 1) // 2
    pb = nb * (nb - 1) // 2
    n = na + nb
    outa = np.zeros(na, np.int64)
    outb = np.zeros(nb, np.int64)
    outc = np.zeros(n, np.int64)
    sa = np.zeros(na, np.int64)
    ea = np.zeros(na, np.int64)
    ypa = np.zeros(na, np.int64)
    yma = np.zeros(na, np.int64)
    Qa = np.zeros((na, na), np.int64)
    sb = np.zeros(nb, np.int64)
    eb = np.zeros(nb, np.int64)
    ypb = np.zeros(nb, np.int64)
    ymb = np.zeros(nb, np.int64)
    Qb = np.zeros((nb, nb), np.int64)
    fullb = (1 << nb) - 1
    for ca in range(1 << pa):
        core.out_masks(ca, na, outa)
        HA = literal(na, outa, perms_a, sa, ea, ypa, yma, Qa)
        for cb in range(1 << pb):
            core.out_masks(cb, nb, outb)
            HB = literal(nb, outb, perms_b, sb, eb, ypb, ymb, Qb)
            for a in range(na):
                for b in range(nb):
                    # A => B, all arcs from A to B, except b -> a
                    for v in range(na):
                        outc[v] = outa[v] | (fullb << na)
                    for v in range(nb):
                        outc[na + v] = outb[v] << na
                    outc[a] &= ~(1 << (na + b))
                    outc[na + b] |= (1 << a)
                    if use_literal:
                        Hc, sw = literal_h_start(n, outc, perms_c, 0)
                    else:
                        Hc, sw = hk(n, outc)
                    cnt[0] += 1
                    if Hc != HA * HB - ea[a] * sb[b] + yma[a] * ypb[b]:
                        cnt[1] += 1


def main():
    perms = {n: core.perm_table(n) for n in range(1, 8)}
    tot = np.zeros(7, np.int64)
    for n in range(1, 7):
        cnt = np.zeros(7, np.int64)
        lit = (n + 1 <= 6)
        run_f124(n, perms[n], perms[n + 1], lit, cnt)
        print("order %d bases: %d labelled tournaments, %d cuts, child counted %s; failures F1 %d  F2 %d  "
              "F2-dual %d  F4 %d  kernels %d" % (n, cnt[0], cnt[1], "literally" if lit else "by Held-Karp",
                                                cnt[2], cnt[3], cnt[4], cnt[5], cnt[6]))
        tot += cnt
    check("F1 on %d (T, cut) pairs, orders 1-6" % tot[1], tot[2] == 0)
    check("F2 (w loses only to b) and its dual on every (T, b), orders 1-6", tot[3] == 0 and tot[4] == 0)
    check("F4 = THM-4097 (4) on %d (T, cut) pairs, orders 1-6" % tot[1], tot[5] == 0)
    check("core.py kernels (dp_tables, rooted_direct, base_data) agree with literal enumeration on %d tournaments"
          % tot[0], tot[6] == 0)
    expected = sum((1 << (n * (n - 1) // 2)) for n in range(1, 7))
    check("universe size %d labelled tournaments" % expected, tot[0] == expected)

    f3 = np.zeros(2, np.int64)
    for na in range(1, 7):
        for nb in range(1, 7):
            if na + nb > 7:
                continue
            cnt = np.zeros(2, np.int64)
            lit = (na + nb <= 6)
            run_f3(na, nb, perms[na], perms[nb], perms[na + nb], lit, cnt)
            print("F3 |A| = %d, |B| = %d: %d (A, B, a, b), joined count %s, failures %d"
                  % (na, nb, cnt[0], "literal" if lit else "Held-Karp", cnt[1]))
            f3 += cnt
    check("F3 on %d quadruples (A, B, a, b), |A| + |B| <= 7" % f3[0], f3[1] == 0)

    # negative control: a deliberately wrong version of F2 must fail somewhere
    n = 4
    out = np.zeros(n, np.int64)
    s = np.zeros(n, np.int64); e = np.zeros(n, np.int64); yp = np.zeros(n, np.int64); ym = np.zeros(n, np.int64)
    Q = np.zeros((n, n), np.int64)
    wrong = 0
    for code in range(64):
        core.out_masks(code, n, out)
        H = literal(n, out, perms[n], s, e, yp, ym, Q)
        for b in range(n):
            outc = np.zeros(n + 1, np.int64)
            S = 15 ^ (1 << b)
            for v in range(n):
                outc[v] = out[v] | ((1 << n) if not (S >> v) & 1 else 0)
            outc[n] = S
            Hc, _ = literal_h_start(n + 1, outc, perms[n + 1], n)
            if Hc != H - e[b] + yp[b]:
                wrong += 1
    check("negative control: the wrong law H - e + Yplus fails (%d failures at order 4)" % wrong, wrong > 0)

    if FAIL:
        print("FAILED:", FAIL)
        sys.exit(1)
    print("ALL CHECKS PASSED")


if __name__ == "__main__":
    main()
