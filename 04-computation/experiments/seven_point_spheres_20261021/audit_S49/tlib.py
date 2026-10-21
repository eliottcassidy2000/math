"""Shared numba tournament tools for the THM-4743 audit (own code; nothing imported from the repo)."""
import numpy as np
from numba import njit


@njit(cache=True)
def popc(x):
    c = 0
    while x:
        x &= x - 1
        c += 1
    return c


@njit(cache=True)
def wl_col(out, n, init):
    """invariant colour refinement starting from init colours (canonical hashing)."""
    col = init.copy()
    ncol = -1
    for it in range(n + 2):
        h = np.zeros(n, np.int64)
        for u in range(n):
            cnt = np.zeros(n + 1, np.int64)
            for v in range(n):
                if (out[u] >> v) & 1:
                    cnt[col[v]] += 1
            x = (col[u] * 1000003 + 7) % 2305843009213693951
            for c in range(n + 1):
                x = (x * 1000033 + cnt[c] + 1) % 2305843009213693951
            h[u] = x
        order = np.argsort(h)
        newcol = np.zeros(n, np.int64)
        k = 0
        for i in range(n):
            if i > 0 and h[order[i]] != h[order[i - 1]]:
                k += 1
            newcol[order[i]] = k
        col = newcol
        if k + 1 == ncol:
            break
        ncol = k + 1
        if ncol == n:
            break
    return col


@njit(cache=True)
def all_isos(outA, outB, colA, colB, n, cap):
    """all colour-respecting isomorphisms A -> B (u -> v in A iff img u -> img v in B), up to cap; returns (arr, count)."""
    res = np.zeros((cap, n), np.int64)
    cnt = 0
    img = -np.ones(n, np.int64)
    used = np.zeros(n, np.bool_)
    cand = np.zeros(n, np.int64)
    depth = 0
    while depth >= 0:
        if depth == n:
            if cnt < cap:
                res[cnt, :] = img
            cnt += 1
            if cnt >= cap:
                break
            depth -= 1
            used[img[depth]] = False
            img[depth] = -1
            cand[depth] += 1
            continue
        u = depth
        found = False
        while cand[u] < n:
            x = cand[u]
            if (not used[x]) and colA[u] == colB[x]:
                ok = True
                for w in range(u):
                    if ((outA[u] >> w) & 1) != ((outB[x] >> img[w]) & 1):
                        ok = False
                        break
                if ok:
                    found = True
                    break
            cand[u] += 1
        if found:
            img[u] = cand[u]; used[cand[u]] = True
            depth += 1
            if depth < n:
                cand[depth] = 0
        else:
            cand[u] = 0
            depth -= 1
            if depth >= 0:
                used[img[depth]] = False
                img[depth] = -1
                cand[depth] += 1
    return res[:min(cnt, cap)], cnt


@njit(cache=True)
def converse(out, n):
    full = (1 << n) - 1
    c = np.zeros(n, np.int64)
    for u in range(n):
        c[u] = full ^ out[u] ^ (1 << u)
    return c


@njit(cache=True)
def joint_colours(outA, outB, n, initA, initB):
    """WL on the disjoint union so that colours are comparable across A and B."""
    N = 2 * n
    U = np.zeros(N, np.int64)
    for u in range(n):
        U[u] = outA[u]
        U[n + u] = outB[u] << n
    init = np.zeros(N, np.int64)
    for u in range(n):
        init[u] = initA[u]; init[n + u] = initB[u]
    # canonical relabel of init
    order = np.argsort(init)
    lab = np.zeros(N, np.int64); k = 0
    for i in range(N):
        if i > 0 and init[order[i]] != init[order[i - 1]]:
            k += 1
        lab[order[i]] = k
    col = wl_col(U, N, lab)
    return col[:n].copy(), col[n:].copy()


@njit(cache=True)
def is_regular_pm(out, n, S):
    """True iff no nonidentity automorphism and no anti-automorphism of T maps S to itself."""
    conv = converse(out, n)
    initA = np.zeros(n, np.int64); initB = np.zeros(n, np.int64)
    for u in range(n):
        initA[u] = ((S >> u) & 1) * 1000 + popc(out[u])
        initB[u] = ((S >> u) & 1) * 1000 + popc(conv[u])
    cA, cB = joint_colours(out, conv, n, initA, initB)
    # anti: any iso (T,S) -> (T^conv,S)
    isos, c = all_isos(out, conv, cA, cB, n, 1)
    if c > 0:
        return False
    cA2, cA3 = joint_colours(out, out, n, initA, initA)
    auts, c = all_isos(out, out, cA2, cA3, n, 2)
    return c == 1


@njit(cache=True)
def image(S, perm, n):
    r = 0
    for i in range(n):
        if (S >> i) & 1:
            r |= 1 << perm[i]
    return r


@njit(cache=True)
def group_pm(out, n, cap):
    """Aut and anti (full lists, up to cap) by backtracking with WL colours."""
    conv = converse(out, n)
    initA = np.zeros(n, np.int64); initB = np.zeros(n, np.int64)
    for u in range(n):
        initA[u] = popc(out[u]); initB[u] = popc(conv[u])
    cA, cB = joint_colours(out, conv, n, initA, initB)
    cA2, cA3 = joint_colours(out, out, n, initA, initA)
    aut, ca = all_isos(out, out, cA2, cA3, n, cap)
    anti, cn = all_isos(out, conv, cA, cB, n, cap)
    return aut, ca, anti, cn


@njit(cache=True)
def count_regular(n, els, k):
    """number of subsets fixed by none of the k permutations els[0..k-1]."""
    tot = 0
    for S in range(1 << n):
        reg = True
        for e in range(k):
            if image(S, els[e], n) == S:
                reg = False
                break
        if reg:
            tot += 1
    return tot


def nonid(arr):
    n = arr.shape[1]
    idt = np.arange(n)
    return np.array([g for g in arr if not np.array_equal(g, idt)], dtype=np.int64).reshape(-1, n)


def from_rows(rows):
    return np.array(rows, dtype=np.int64)


def circ(n, S):
    return from_rows([sum(1 << ((x + d) % n) for d in S) for x in range(n)])


def lexprod(A, B):
    A = [int(x) for x in A]; B = [int(x) for x in B]
    na, nb = len(A), len(B); T = [0] * (na * nb)
    for a in range(na):
        for b in range(nb):
            for a2 in range(na):
                for b2 in range(nb):
                    if (a != a2 and (A[a] >> a2) & 1) or (a == a2 and (B[b] >> b2) & 1):
                        T[a * nb + b] |= 1 << (a2 * nb + b2)
    return from_rows(T)


def subst(A, Bs):
    """A[B_1, ..., B_m] with unequal blocks."""
    A = [int(x) for x in A]; Bs = [[int(x) for x in B] for B in Bs]
    off = [0]
    for B in Bs:
        off.append(off[-1] + len(B))
    n = off[-1]; T = [0] * n
    for a in range(len(A)):
        for b in range(len(Bs[a])):
            u = off[a] + b
            for a2 in range(len(A)):
                for b2 in range(len(Bs[a2])):
                    v = off[a2] + b2
                    if (a != a2 and (A[a] >> a2) & 1) or (a == a2 and (Bs[a][b] >> b2) & 1):
                        T[u] |= 1 << v
    return from_rows(T)


def is_tournament(T):
    n = len(T)
    for u in range(n):
        if (int(T[u]) >> u) & 1:
            return False
        for v in range(u + 1, n):
            if ((int(T[u]) >> v) & 1) == ((int(T[v]) >> u) & 1):
                return False
    return True


def iso_exists(A, B):
    n = len(A)
    if n != len(B):
        return False
    initA = np.array([popc(int(x)) for x in A], dtype=np.int64); initB = np.array([popc(int(x)) for x in B], dtype=np.int64)
    cA, cB = joint_colours(A, B, n, initA, initB)
    _, c = all_isos(A, B, cA, cB, n, 1)
    return c > 0


def random_regular_search(T, tries, rng):
    n = len(T)
    for t in range(tries):
        S = int(rng.integers(0, 1 << n))
        if is_regular_pm(T, n, S):
            return S, t + 1
    return None, tries
