"""Independent census for THM-4743 (6) / HYP-9353: exceptional tournaments (no subset with trivial stabilizer in
Aut+-(T)) among ALL self-converse tournaments on n vertices.

Completeness: every self-converse tournament has an anti-automorphism that is an involution with <= 1 fixed point
(take an odd power of any anti-automorphism to get one of 2-power order; a cycle of length 2^j >= 4 would give an
automorphism sigma^(2^(j-1)) swapping two vertices; two fixed points of an anti map would reverse their arc).
So every self-converse tournament is isomorphic to one for which sigma0 = (0 1)(2 3)... (+ fixed n-1 if n odd) is an
anti-automorphism; we enumerate ALL such labelled tournaments (2^bits of them).
Non-self-converse tournaments have odd |Aut+-| = |Aut| and are never exceptional (Gluck Cor. 1).
For each labelled T: WL colour refinement (invariant); discrete => Aut = 1 => Aut+- = {1, sigma0} => not exceptional.
Otherwise all automorphisms by backtracking, Aut+- = Aut u Aut.sigma0, and an exhaustive subset scan."""
import sys, time
import numpy as np
from numba import njit, prange


def orbit_tables(n):
    sig = list(range(n))
    for i in range(0, n - 1, 2):
        sig[i], sig[i + 1] = i + 1, i
    seen = set(); orbs = []
    for u in range(n):
        for v in range(u + 1, n):
            if (u, v) in seen:
                continue
            # bit = 1 means u -> v ; anti-invariance then forces sig(v) -> sig(u)
            a = (u, v); b = (sig[v], sig[u])
            arcs = [a]
            key_b = (min(b), max(b))
            if key_b != (u, v):
                arcs.append(b)
            seen.add((u, v)); seen.add(key_b)
            orbs.append(arcs)
    m = len(orbs)
    tab = np.full((m, 2, 2), -1, dtype=np.int64)
    for j, arcs in enumerate(orbs):
        for t, (x, y) in enumerate(arcs):
            tab[j, t, 0] = x; tab[j, t, 1] = y
    return np.array(sig, dtype=np.int64), tab


@njit(cache=True)
def build(mask, n, tab):
    out = np.zeros(n, np.int64)
    m = tab.shape[0]
    for j in range(m):
        b = (mask >> j) & 1
        for t in range(2):
            x = tab[j, t, 0]
            if x < 0:
                continue
            y = tab[j, t, 1]
            if b:
                out[x] |= 1 << y
            else:
                out[y] |= 1 << x
    return out


@njit(cache=True)
def popc(x):
    c = 0
    while x:
        x &= x - 1
        c += 1
    return c


@njit(cache=True)
def wl(out, n):
    col = np.zeros(n, np.int64)
    for u in range(n):
        col[u] = popc(out[u])
    ncol = -1
    for it in range(n + 1):
        h = np.zeros(n, np.int64)
        for u in range(n):
            cnt = np.zeros(n + 1, np.int64)
            for v in range(n):
                if (out[u] >> v) & 1:
                    cnt[col[v]] += 1
            x = col[u] * 1000003 + 7
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
    return col, ncol


@njit(cache=True)
def automorphisms(out, n, col, cap):
    """all automorphisms preserving the (invariant) colouring col; returns array (cnt, n)."""
    res = np.zeros((cap, n), np.int64)
    cnt = 0
    img = -np.ones(n, np.int64)
    used = np.zeros(n, np.bool_)
    cand = np.zeros(n, np.int64)          # current candidate index per depth
    depth = 0
    # vertex order: 0..n-1
    while depth >= 0:
        if depth == n:
            if cnt < cap:
                res[cnt, :] = img
            cnt += 1
            depth -= 1
            if depth >= 0:
                used[img[depth]] = False
                img[depth] = -1
                cand[depth] += 1
            continue
        u = depth
        found = False
        while cand[u] < n:
            x = cand[u]
            if (not used[x]) and col[x] == col[u]:
                ok = True
                for w in range(u):
                    if ((out[u] >> w) & 1) != ((out[x] >> img[w]) & 1):
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
def image(S, perm, n):
    r = 0
    for i in range(n):
        if (S >> i) & 1:
            r |= 1 << perm[i]
    return r


@njit(cache=True)
def has_regular(n, aut, sig):
    na = aut.shape[0]
    els = np.zeros((2 * na, n), np.int64)
    k = 0
    for g in range(na):
        isid = True
        for i in range(n):
            if aut[g, i] != i:
                isid = False
        if not isid:
            els[k, :] = aut[g, :]; k += 1
        for i in range(n):
            els[k, i] = aut[g, sig[i]]
        k += 1
    for S in range(1 << n):
        reg = True
        for e in range(k):
            if image(S, els[e], n) == S:
                reg = False
                break
        if reg:
            return True
    return False


@njit(parallel=True, cache=True)
def census(n, tab, sig, nbits, chunk, maxkeep):
    total = 1 << nbits
    nchunks = (total + chunk - 1) // chunk
    keep = np.zeros((nchunks, maxkeep), np.int64)
    nkeep = np.zeros(nchunks, np.int64)
    nontriv = np.zeros(nchunks, np.int64)
    nondisc = np.zeros(nchunks, np.int64)
    for c in prange(nchunks):
        lo = c * chunk; hi = min(total, lo + chunk)
        for mask in range(lo, hi):
            out = build(mask, n, tab)
            col, ncol = wl(out, n)
            if ncol == n:
                continue
            nondisc[c] += 1
            aut, cnt = automorphisms(out, n, col, 4096)
            if cnt > 1:
                nontriv[c] += 1
            if cnt == 1:
                continue          # Aut+- = {1, sigma0}; sigma0 moves a vertex -> singleton regular
            if not has_regular(n, aut, sig):
                if nkeep[c] < maxkeep:
                    keep[c, nkeep[c]] = mask
                nkeep[c] += 1
    return keep, nkeep, nontriv.sum(), nondisc.sum()


def py_iso(A, B, n):
    """all isomorphisms A -> B (bitmask rows), python backtracking"""
    sa = [bin(x).count('1') for x in A]; sb = [bin(x).count('1') for x in B]
    res = []
    img = [-1] * n; used = [False] * n
    def bt(u):
        if u == n:
            res.append(tuple(img)); return
        for x in range(n):
            if used[x] or sa[u] != sb[x]:
                continue
            if all(((A[u] >> w) & 1) == ((B[x] >> img[w]) & 1) for w in range(u)):
                img[u] = x; used[x] = True
                bt(u + 1)
                img[u] = -1; used[x] = False
    bt(0)
    return res


def circ(n, S):
    return [sum(1 << ((x + d) % n) for d in S) for x in range(n)]


def lexprod(A, B):
    na, nb = len(A), len(B); T = [0] * (na * nb)
    for a in range(na):
        for b in range(nb):
            u = a * nb + b
            for a2 in range(na):
                for b2 in range(nb):
                    v = a2 * nb + b2
                    if (a != a2 and (A[a] >> a2) & 1) or (a == a2 and (B[b] >> b2) & 1):
                        T[u] |= 1 << v
    return T


if __name__ == '__main__':
    ns = [int(x) for x in sys.argv[1:]]
    fam = {3: circ(3, {1}), 5: circ(5, {1, 2}), 7: circ(7, {1, 2, 4}), 9: lexprod(circ(3, {1}), circ(3, {1}))}
    for n in ns:
        t0 = time.time()
        sig, tab = orbit_tables(n)
        nbits = tab.shape[0]
        # sanity: sigma0 is an anti-automorphism of a few random members
        rng = np.random.default_rng(n)
        for _ in range(50):
            mk = int(rng.integers(0, 1 << nbits))
            out = build(mk, n, tab)
            for u in range(n):
                assert (out[u] >> u) & 1 == 0
                for v in range(n):
                    if u != v:
                        assert ((out[u] >> v) & 1) ^ ((out[v] >> u) & 1) == 1
                        assert ((out[u] >> v) & 1) == ((out[sig[v]] >> sig[u]) & 1)
        chunk = 1 << max(0, min(16, nbits - 4))
        keep, nkeep, nontriv, nondisc = census(n, tab, sig, nbits, chunk, 64)
        exc = [int(keep[c, i]) for c in range(keep.shape[0]) for i in range(min(nkeep[c], 64))]
        total_exc = int(nkeep.sum())
        # identify
        ident = {}
        for mk in exc:
            T = [int(x) for x in build(mk, n, tab)]
            nm = 'NEW'
            if n in fam and py_iso(T, fam[n], n):
                nm = {3: 'C3', 5: 'T5', 7: 'P7', 9: 'C3[C3]'}[n]
            ident[nm] = ident.get(nm, 0) + 1
        print(f'n={n}: labelled sigma0-anti tournaments 2^{nbits} = {1 << nbits}; WL non-discrete {nondisc}; |Aut|>1 {nontriv}; '
              f'exceptional labelled {total_exc}; identified (of kept) {ident}  [{time.time() - t0:.0f}s]', flush=True)
