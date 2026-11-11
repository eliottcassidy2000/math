#!/usr/bin/env python3
"""S91 helper (census_k2): exact by-value census of the root-change words of THM-4916 / HYP-9506, with the
confinement classification of every counted word.  Own code (numba); the enumeration scheme (prefix states, one tail
loop per root-change seed) is that of rc_lf.c / decomp.c (kind-pasteur S15), the classification is new here.

Model: states (D, E, U, d), root t = (1, 0, 2, 1), letters
    O: (D, 0, 2U, d)        Z: (D - E + U, d, 2U - d, d)        J: (D - E + U, D - E, D + U, D).
Counted words end in Z or J, value = first coordinate of the final state; M_K(N) = number of counted words with at
most K letters J and value N.  Every counted word is (prefix w) J (tail m): seed (h, e, R, d) = J(state of w), value
    H(m) = h + R m - d Q(m) - e (m mod 2),        Q(1) = 0, Q(2m) = Q(m), Q(2m+1) = Q(m) + m + (m mod 2),
where the J-free words are the tails of the root seed (1, 0, 2, 1) = J(virtual state (D, E, U) = (1, 1, 1)).

Dictionary (dict_k1.py): with the pre-J state (D, E, U) and the tail m <-> (c, h3) (m even: c = 0, h3 = m/2; m odd:
c = v2(m+1), h3 = m >> (c+1)), x3 = x(h3), y3 = y(h3):
        value = D x3 + 2^c U y3 - [c = 0] E.
Classification of a counted word, for an odd prime p:
    S   seed-confined at the last root change:  p | gcd(D, U) of the pre-J state        (THM-4916 (1a), cascade)
    T   tail-confined:                          p | gcd(x3, y3)                         (THM-4916 (1b) corollary)
    C1  (two-J words) the FIRST root change is confined: p | gcd(D1, U1)
    shift type st: 0 = odd tail (c >= 1), 1 = even tail and p | E, 2 = even tail and p does not divide E.
A word with S or T has value = -[c = 0] E (mod p) exactly (checked for every such word: err must stay 0); its
distinguished class is 0 for st in {0, 1} and -E != 0 for st = 2.
Flags of a two-J word A J B J C:  k1 = 1 when A ends in O (or is empty) or B starts with Z ("sheet s >= 1 at the
first root change"), l2 = 1 when B ends in O or C starts with Z ("sheet s >= 1 at the second root change").
"""
import numpy as np
from numba import njit, prange

NJ = 3      # number of J's: 0, 1, 2


@njit(cache=True)
def build_Q(size):
    Q = np.zeros(size, np.int32)
    for m in range(1, size):
        half = m >> 1
        if m & 1:
            Q[m] = Q[half] + half + (half & 1)
        else:
            Q[m] = Q[half]
    return Q


@njit(cache=True)
def build_tailmask(Q, primes):
    size = Q.shape[0]
    tm = np.zeros(size, np.uint8)
    for j in range(primes.shape[0]):
        p = primes[j]
        for m in range(p - 1, size, p):          # p | m + 1  <=>  p | y3
            if m & 1:
                if (Q[m] + 1) % p == 0:           # x3 = m - Q(m) = -(Q(m) + 1) mod p
                    tm[m] |= (1 << j)
            else:
                if Q[m] % p == 0:                 # x3 = m + 1 - Q(m) = -Q(m) mod p
                    tm[m] |= (1 << j)
    return tm


@njit(cache=True)
def gcd_mask(D, U, primes):
    mk = 0
    for j in range(primes.shape[0]):
        p = primes[j]
        if D % p == 0 and U % p == 0:
            mk |= (1 << j)
    return mk


@njit(cache=True)
def run_seed(D, E, U, d, nJ, c1mask, k1, LO, X, Q, tm, primes, tot, confhist, cnt, allcnt, allp, cnt2, err):
    """all counted words (pre-J state (D, E, U, d)) J (tail m) with value in [LO, X]."""
    h = D - E + U
    if h > X:
        return
    R = D + U
    e = D - E
    dd = D
    npr = primes.shape[0]
    smask = 0
    ezmask = 0
    for j in range(npr):
        p = primes[j]
        if D % p == 0 and U % p == 0:
            smask |= (1 << j)
        if E % p == 0:
            ezmask |= (1 << j)
    z = 1 if E == 0 else 0
    base = smask | c1mask
    if h >= LO:
        mlo = 0
    else:
        mlo = (LO - h + R - 1) // R
    mhi = (X - h) // (R - dd)
    n_odd = 0
    n_all = 0
    for m in range(mlo, mhi + 1):
        par = m & 1
        v = h + R * m - dd * Q[m] - e * par
        if v < LO or v > X:
            continue
        tot[(v - LO) >> 1] += 1
        n_all += 1
        n_odd += par
        tmm = tm[m]
        mk = base | tmm
        if mk != 0:
            j = 0
            while mk != 0:
                if mk & 1:
                    p = primes[j]
                    S = (smask >> j) & 1
                    T = (tmm >> j) & 1
                    C1 = (c1mask >> j) & 1
                    if par == 1:
                        st = 0
                    elif (ezmask >> j) & 1:
                        st = 1
                    else:
                        st = 2
                    cnt[j, nJ, S, T, C1, st] += 1
                    if S | T:
                        confhist[j, v % p] += 1
                        sh = 0 if par == 1 else E
                        if (v + sh) % p != 0:
                            err[0] += 1
                mk >>= 1
                j += 1
    n_even = n_all - n_odd
    allcnt[nJ, 1, z] += n_odd
    allcnt[nJ, 0, z] += n_even
    for j in range(npr):
        allp[j, nJ, 0] += n_odd
        if (ezmask >> j) & 1:
            allp[j, nJ, 1] += n_even
        else:
            allp[j, nJ, 2] += n_even
    if nJ == 2:
        cnt2[k1, 1] += n_odd
        cnt2[k1, z] += n_even
    err[1] += n_all


@njit(cache=True)
def run_root_seed(LO, X, Q, primes, tot, confhist, cnt, allcnt, allp, err):
    """the J-free words: tails of the root seed (h, e, R, d) = (1, 0, 2, 1) = J(virtual pre-state (1, 1, 1)).
    Q and the tail mask are computed on the fly beyond the table (m can reach X - 1)."""
    npr = primes.shape[0]
    size = Q.shape[0]
    n_odd = 0
    n_all = 0
    mlo = 0
    if LO > 1:
        mlo = (LO - 1 + 1) // 2
    for m in range(mlo, X):
        if m < size:
            q = Q[m]
        else:
            half = m >> 1
            q = Q[half]
            if m & 1:
                q += half + (half & 1)
        par = m & 1
        v = 1 + 2 * m - q
        if v < LO or v > X:
            continue
        tot[(v - LO) >> 1] += 1
        n_all += 1
        n_odd += par
        for j in range(npr):
            p = primes[j]
            if (m + 1) % p == 0:
                if (par == 1 and (q + 1) % p == 0) or (par == 0 and q % p == 0):
                    st = 0 if par == 1 else 2          # E = 1: never 0 mod p
                    cnt[j, 0, 0, 1, 0, st] += 1
                    confhist[j, v % p] += 1
                    sh = 0 if par == 1 else 1
                    if (v + sh) % p != 0:
                        err[0] += 1
    n_even = n_all - n_odd
    allcnt[0, 1, 0] += n_odd
    allcnt[0, 0, 0] += n_even
    for j in range(npr):
        allp[j, 0, 0] += n_odd
        allp[j, 0, 2] += n_even
    err[1] += n_all


@njit(cache=True)
def run_subtree(D0, E0, U0, d0, j0, c1m0, e10, fb0, K, LO, X, Q, tm, primes, tot, confhist, cnt, allcnt, allp, cnt2, err):
    """depth-first over the prefixes below (and including) the prefix state given; a seed for every J applied."""
    stack = np.empty((1024, 8), np.int64)
    stack[0, 0] = D0
    stack[0, 1] = E0
    stack[0, 2] = U0
    stack[0, 3] = d0
    stack[0, 4] = j0
    stack[0, 5] = c1m0
    stack[0, 6] = e10
    stack[0, 7] = fb0
    sp = 1
    while sp > 0:
        sp -= 1
        D = stack[sp, 0]
        E = stack[sp, 1]
        U = stack[sp, 2]
        d = stack[sp, 3]
        j = stack[sp, 4]
        c1m = stack[sp, 5]
        e1 = stack[sp, 6]
        fb = stack[sp, 7]
        if j >= K:
            continue
        hJ = D - E + U
        if hJ <= X:
            if j == 0:
                run_seed(D, E, U, d, 1, 0, 0, LO, X, Q, tm, primes, tot, confhist, cnt, allcnt, allp, cnt2, err)
            else:
                k1 = 1 if (e1 == 0 or fb == 2) else 0
                run_seed(D, E, U, d, 2, c1m, k1, LO, X, Q, tm, primes, tot, confhist, cnt, allcnt, allp, cnt2, err)
            if j + 1 < K:
                stack[sp, 0] = hJ
                stack[sp, 1] = D - E
                stack[sp, 2] = D + U
                stack[sp, 3] = D
                stack[sp, 4] = j + 1
                stack[sp, 5] = gcd_mask(D, U, primes)
                stack[sp, 6] = E
                stack[sp, 7] = 0
                sp += 1
        if 2 * U <= X:                              # letter O
            stack[sp, 0] = D
            stack[sp, 1] = 0
            stack[sp, 2] = 2 * U
            stack[sp, 3] = d
            stack[sp, 4] = j
            stack[sp, 5] = c1m
            stack[sp, 6] = e1
            stack[sp, 7] = fb if (j == 0 or fb != 0) else 1
            sp += 1
        if 2 * U - d <= X:                          # letter Z
            stack[sp, 0] = D - E + U
            stack[sp, 1] = d
            stack[sp, 2] = 2 * U - d
            stack[sp, 3] = d
            stack[sp, 4] = j
            stack[sp, 5] = c1m
            stack[sp, 6] = e1
            stack[sp, 7] = fb if (j == 0 or fb != 0) else 2
            sp += 1


@njit(cache=True)
def gen_tasks(K, X, usplit, primes, maxtasks):
    """split the prefix tree: a prefix with U >= usplit becomes a subtree task (type 1); a shallower prefix gives a
    seed task (type 0: apply J there and run its tails) and is expanded further."""
    tasks = np.empty((maxtasks, 9), np.int64)
    nt = 0
    stack = np.empty((1024, 8), np.int64)
    stack[0, 0] = 1
    stack[0, 1] = 0
    stack[0, 2] = 2
    stack[0, 3] = 1
    stack[0, 4] = 0
    stack[0, 5] = 0
    stack[0, 6] = 0
    stack[0, 7] = 0
    sp = 1
    while sp > 0:
        sp -= 1
        D = stack[sp, 0]
        E = stack[sp, 1]
        U = stack[sp, 2]
        d = stack[sp, 3]
        j = stack[sp, 4]
        c1m = stack[sp, 5]
        e1 = stack[sp, 6]
        fb = stack[sp, 7]
        if j >= K:
            continue
        if U >= usplit:
            tasks[nt, 0] = 1
            for i in range(8):
                tasks[nt, 1 + i] = stack[sp, i]
            nt += 1
            continue
        hJ = D - E + U
        if hJ <= X:
            tasks[nt, 0] = 0
            for i in range(8):
                tasks[nt, 1 + i] = stack[sp, i]
            nt += 1
            if j + 1 < K:
                stack[sp, 0] = hJ
                stack[sp, 1] = D - E
                stack[sp, 2] = D + U
                stack[sp, 3] = D
                stack[sp, 4] = j + 1
                stack[sp, 5] = gcd_mask(D, U, primes)
                stack[sp, 6] = E
                stack[sp, 7] = 0
                sp += 1
        if 2 * U <= X:
            stack[sp, 0] = D
            stack[sp, 1] = 0
            stack[sp, 2] = 2 * U
            stack[sp, 3] = d
            stack[sp, 4] = j
            stack[sp, 5] = c1m
            stack[sp, 6] = e1
            stack[sp, 7] = fb if (j == 0 or fb != 0) else 1
            sp += 1
        if 2 * U - d <= X:
            stack[sp, 0] = D - E + U
            stack[sp, 1] = d
            stack[sp, 2] = 2 * U - d
            stack[sp, 3] = d
            stack[sp, 4] = j
            stack[sp, 5] = c1m
            stack[sp, 6] = e1
            stack[sp, 7] = fb if (j == 0 or fb != 0) else 2
            sp += 1
    return tasks[:nt]


@njit(parallel=True, cache=True)
def run_tasks(tasks, K, LO, X, Q, tm, primes, nthreads, tot, confhist, cnt, allcnt, allp, cnt2, err):
    ntasks = tasks.shape[0]
    for t in prange(nthreads):
        for i in range(t, ntasks, nthreads):
            typ = tasks[i, 0]
            D = tasks[i, 1]
            E = tasks[i, 2]
            U = tasks[i, 3]
            d = tasks[i, 4]
            j = tasks[i, 5]
            c1m = tasks[i, 6]
            e1 = tasks[i, 7]
            fb = tasks[i, 8]
            if typ == 0:
                if j == 0:
                    run_seed(D, E, U, d, 1, 0, 0, LO, X, Q, tm, primes, tot[t], confhist[t], cnt[t], allcnt[t], allp[t],
                             cnt2[t], err[t])
                else:
                    k1 = 1 if (e1 == 0 or fb == 2) else 0
                    run_seed(D, E, U, d, 2, c1m, k1, LO, X, Q, tm, primes, tot[t], confhist[t], cnt[t], allcnt[t],
                             allp[t], cnt2[t], err[t])
            else:
                run_subtree(D, E, U, d, j, c1m, e1, fb, K, LO, X, Q, tm, primes, tot[t], confhist[t], cnt[t],
                            allcnt[t], allp[t], cnt2[t], err[t])


def census(K, LO, HI, primes, nthreads=3, usplit=512, seed=12345):
    """exact census of the counted words with at most K letters J and value in [LO, HI) (LO, HI even or powers of 2).
    Returns a dict: tot (int64 per odd N of the window), confhist, cnt, allcnt, allp, cnt2, nwords."""
    import numba
    numba.set_num_threads(nthreads)
    primes = np.array(primes, dtype=np.int64)
    assert primes.shape[0] <= 8
    X = HI - 1
    Q = build_Q(HI // 2 + 4)
    tm = build_tailmask(Q, primes)
    nodd = (HI - LO) // 2
    maxp = int(primes.max())
    tot = np.zeros((nthreads, nodd), np.uint16)
    confhist = np.zeros((nthreads, primes.shape[0], maxp), np.int64)
    cnt = np.zeros((nthreads, primes.shape[0], NJ, 2, 2, 2, 3), np.int64)
    allcnt = np.zeros((nthreads, NJ, 2, 2), np.int64)
    allp = np.zeros((nthreads, primes.shape[0], NJ, 3), np.int64)
    cnt2 = np.zeros((nthreads, 2, 2), np.int64)
    err = np.zeros((nthreads, 2), np.int64)
    run_root_seed(LO, X, Q, primes, tot[0], confhist[0], cnt[0], allcnt[0], allp[0], err[0])
    if K >= 1:
        tasks = gen_tasks(K, X, usplit, primes, 4000000)
        rng = np.random.default_rng(seed)
        tasks = tasks[rng.permutation(tasks.shape[0])]
        run_tasks(tasks, K, LO, X, Q, tm, primes, nthreads, tot, confhist, cnt, allcnt, allp, cnt2, err)
        ntasks = int(tasks.shape[0])
    else:
        ntasks = 0
    tot64 = tot.astype(np.int64).sum(axis=0)
    out = dict(K=K, LO=LO, HI=HI, primes=[int(p) for p in primes], tot=tot64, confhist=confhist.sum(axis=0),
               cnt=cnt.sum(axis=0), allcnt=allcnt.sum(axis=0), allp=allp.sum(axis=0), cnt2=cnt2.sum(axis=0),
               nwords=int(err[:, 1].sum()), lemma_violations=int(err[:, 0].sum()), ntasks=ntasks,
               maxM=int(tot64.max()), overflow_ok=bool(int(tot64.sum()) == int(err[:, 1].sum())
                                                         and int(tot.max()) < 65535))
    return out


def class_means(tot, LO, HI, p):
    """mean of M over the odd N of [LO, HI) in each class mod p, and the number of such N."""
    N = np.arange(LO + 1, HI, 2, dtype=np.int64)
    res = N % p
    sums = np.bincount(res, weights=tot.astype(np.float64), minlength=p)
    n = np.bincount(res, minlength=p)
    return sums / n, n
