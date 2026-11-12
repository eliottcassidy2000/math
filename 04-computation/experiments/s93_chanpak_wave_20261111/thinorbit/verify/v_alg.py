#!/usr/bin/env python3
"""Verifier (lane thinorbit, S93): independent algebra checks on MY OWN block matrices (verify/v_blocks2_k*.txt).
Sections: (a) rank/determinant census, (b) parity shape, (c) corollary P3, (d) Proposition Z facts, (e) algebra
dimension, (f) SL_3(F_p), (g) 2-adic images, (h) orbit mod q, (i) (h,R)-closed blocks, (j) invariant planes, (k) C3.
"""
import sys
import math
from fractions import Fraction
from itertools import product
import numpy as np

V = "C:/tmp/wt-s93/04-computation/experiments/s93_chanpak_wave_20261111/thinorbit/verify/"


def load(k):
    out = []
    with open(V + f"v_blocks2_k{k}.txt") as f:
        for line in f:
            p = line.split()
            if len(p) == 9:
                out.append(tuple(int(x) for x in p))
    return out


def det(m):
    return (m[0] * (m[4] * m[8] - m[5] * m[7]) - m[1] * (m[3] * m[8] - m[5] * m[6]) + m[2] * (m[3] * m[7] - m[4] * m[6]))


def c2(m):
    return (m[0] * m[4] - m[1] * m[3]) + (m[0] * m[8] - m[2] * m[6]) + (m[4] * m[8] - m[5] * m[7])


def tr(m):
    return m[0] + m[4] + m[8]


def rank(m):
    A = np.array(m, dtype=object).reshape(3, 3)
    rows = [[Fraction(int(x)) for x in r] for r in A]
    r = 0
    for c in range(3):
        piv = None
        for i in range(r, 3):
            if rows[i][c] != 0:
                piv = i
                break
        if piv is None:
            continue
        rows[r], rows[piv] = rows[piv], rows[r]
        for i in range(3):
            if i != r and rows[i][c] != 0:
                fct = rows[i][c] / rows[r][c]
                rows[i] = [a - fct * b for a, b in zip(rows[i], rows[r])]
        r += 1
    return r


def mul(a, b):
    return tuple(sum(a[3 * i + k] * b[3 * k + j] for k in range(3)) for i in range(3) for j in range(3))


def section_a(lib):
    for k in (3, 4, 5, 6):
        ms = lib[k]
        dets = [det(m) for m in ms]
        rk = [rank(m) for m in ms]
        uni = [m for m, d in zip(ms, dets) if abs(d) == 1]
        print(f"(a) order <= {k}: {len(ms)} matrices; rank1 {rk.count(1)}, rank2 {rk.count(2)}, invertible {rk.count(3)}, "
              f"unimodular {len(uni)} (det +1: {sum(1 for m in uni if det(m) == 1)}, det -1: {sum(1 for m in uni if det(m) == -1)}); "
              f"singular {sum(1 for d in dets if d == 0)}")
    d4 = sorted(set(det(m) for m in lib[4] if det(m) != 0))
    print("    determinants of the invertible generators of order <= 4:", d4)


def section_b(lib):
    for k in (4, 5, 6):
        bad = 0
        for m in lib[k]:
            mm = [x % 2 for x in m]
            ok = (mm[0], mm[1], mm[2]) == (1, 0, 0) and mm[4] == mm[8] and mm[5] == mm[7]
            if not ok:
                bad += 1
        print(f"(b) parity shape [[1,0,0],[p,x,y],[q,y,x]] mod 2, order <= {k}: violations {bad} / {len(lib[k])}")


def section_c(lib):
    for k in (4, 5, 6):
        bmin = None
        bad = 0
        n2 = 0
        for m in lib[k]:
            if det(m) != 0:
                continue
            a, b = tr(m), c2(m)
            if rank(m) == 2:
                n2 += 1
                if a % 2 != 1 or b % 2 != 0:
                    bad += 1
                if bmin is None or abs(b) < bmin:
                    bmin = abs(b)
            else:
                if a % 2 != 1 or b % 2 != 0:
                    bad += 1
        print(f"(c) singular blocks of order <= {k}: charpoly t(t^2 - a t + b) with a odd, b even: violations {bad}; "
              f"rank-2 blocks {n2}, min |b| over rank-2 blocks = {bmin}")


def section_d(lib):
    T0 = (1, 2, 0, 0, 0, 1, 1, 1, 0)
    g = (1, 4, 4, 2, 3, 2, 2, 4, 3)
    print(f"(d) T0 in order<=3 list: {T0 in lib[3]}; g in order<=3 list: {g in lib[3]}; det T0 = {det(T0)}, det g = {det(g)}")
    print(f"    charpoly T0: t^3 - {tr(T0)} t^2 + ({c2(T0)}) t - ({det(T0)}); values at +-1: {1 - tr(T0) + c2(T0) - det(T0)}, {-1 - tr(T0) - c2(T0) - det(T0)}")
    # discriminant of t^3 + a t^2 + b t + c
    a, b, c = -tr(T0), c2(T0), -det(T0)
    disc = 18 * a * b * c - 4 * a ** 3 * c + a * a * b * b - 4 * b ** 3 - 27 * c * c
    print(f"    discriminant = {disc}")
    # g^{-1} (integer since det 1): adjugate
    def adj(m):
        return (m[4] * m[8] - m[5] * m[7], m[2] * m[7] - m[1] * m[8], m[1] * m[5] - m[2] * m[4],
                m[5] * m[6] - m[3] * m[8], m[0] * m[8] - m[2] * m[6], m[2] * m[3] - m[0] * m[5],
                m[3] * m[7] - m[4] * m[6], m[1] * m[6] - m[0] * m[7], m[0] * m[4] - m[1] * m[3])
    gi = adj(g)
    assert mul(g, gi) == (1, 0, 0, 0, 1, 0, 0, 0, 1)
    c1 = mul(mul(g, T0), gi)
    print(f"    g T0 g^-1 = {c1}; commutes with T0: {mul(c1, T0) == mul(T0, c1)}")
    uni4 = [m for m in lib[4] if abs(det(m)) == 1]
    print(f"    the 8 unimodular generators of order <= 4 (all det +1: {all(det(m) == 1 for m in uni4)}):")
    for m in uni4:
        print("     ", m, "det", det(m), "charpoly (tr, c2):", tr(m), c2(m))


def span_dim(gens, maxlen=6):
    # dimension of the linear span of all words of length <= maxlen (exact rational elimination)
    basis = []

    def reduce(v):
        v = [Fraction(x) for x in v]
        for b, piv in basis:
            if v[piv] != 0:
                f = v[piv] / b[piv]
                v = [x - f * y for x, y in zip(v, b)]
        return v

    def add(v):
        v = reduce(v)
        for i, x in enumerate(v):
            if x != 0:
                basis.append((v, i))
                return True
        return False

    frontier = []
    for g in gens:
        if add(g):
            frontier.append(g)
    for _ in range(maxlen):
        nxt = []
        for w in frontier:
            for g in gens:
                p = mul(w, g)
                if add(p):
                    nxt.append(p)
        if not nxt:
            break
        frontier = nxt
    return len(basis)


def closure_size(gens, m, vec=None, chunk=200000):
    """size of the multiplicative closure mod m of the matrices (vec None) or of the orbit of vec mod m."""
    G = np.array(gens, dtype=np.int64).reshape(-1, 3, 3) % m
    if vec is None:
        nd = 9
    else:
        nd = 3
    pw = (m ** np.arange(nd - 1, -1, -1)).astype(np.int64)

    def enc(A):
        return (A.reshape(-1, nd) * pw).sum(axis=1)

    def dec(codes):
        out = np.zeros((len(codes), nd), np.int64)
        c = codes.copy()
        for j in range(nd - 1, -1, -1):
            out[:, j] = c % m
            c //= m
        return out

    use_bitmap = m ** nd <= 200_000_000
    if use_bitmap:
        seen = np.zeros(m ** nd, dtype=np.bool_)
    else:
        seen_arr = np.zeros(0, np.int64)
    if vec is None:
        start = np.unique(enc(G))
    else:
        start = np.array([int((np.array(vec) % m * pw).sum())], dtype=np.int64)
    if use_bitmap:
        seen[start] = True
    else:
        seen_arr = start.copy()
    frontier = start
    total = len(start)
    while len(frontier):
        newcodes = []
        for i in range(0, len(frontier), chunk):
            F = dec(frontier[i:i + chunk])
            if vec is None:
                F = F.reshape(-1, 3, 3)
                Pr = np.einsum('nij,gjk->ngik', F, G) % m
                codes = np.unique(enc(Pr.reshape(-1, 9)))
            else:
                Pr = np.einsum('gij,nj->ngi', G, F) % m
                codes = np.unique(enc(Pr.reshape(-1, 3)))
            if use_bitmap:
                codes = codes[~seen[codes]]
                seen[codes] = True
            newcodes.append(codes)
        new = np.unique(np.concatenate(newcodes)) if newcodes else np.zeros(0, np.int64)
        if not use_bitmap:
            new = new[~np.isin(new, seen_arr, assume_unique=True)]
            seen_arr = np.union1d(seen_arr, new)
        total += len(new)
        frontier = new
    if vec is None:
        return total
    # also number of residues of h
    if use_bitmap:
        codes = np.nonzero(seen)[0]
    else:
        codes = seen_arr
    hres = np.unique(codes // (m * m))
    return total, len(hres)


def section_efg(lib):
    uni4 = [m for m in lib[4] if abs(det(m)) == 1]
    uni6 = [m for m in lib[6] if abs(det(m)) == 1]
    print(f"(e) dimension of the algebra spanned by words in the 8 unimodular generators: {span_dim(uni4)}; "
          f"by T0 and g alone: {span_dim([(1, 2, 0, 0, 0, 1, 1, 1, 0), (1, 4, 4, 2, 3, 2, 2, 4, 3)])}")
    for p in (2, 3, 5, 7):
        n = closure_size(uni4, p)
        full = p ** 3 * (p ** 2 - 1) * (p ** 3 - 1)
        print(f"(f) group generated mod {p} by the 8 unimodular generators: {n}; |SL_3(F_{p})| = {full}; equal {n == full}", flush=True)
    prev = None
    for e in (1, 2, 3, 4):
        n = closure_size(uni4, 2 ** e)
        print(f"(g) 8 unimodular generators (order <= 4) mod 2^{e}: {n} = 2^{math.log2(n):.3f}" + (f"; ratio {n // prev}" if prev else ""), flush=True)
        prev = n
    for e in (1, 2, 3):
        n = closure_size(uni6, 2 ** e)
        print(f"(g) 46 unimodular generators (order <= 6) mod 2^{e}: {n} = 2^{math.log2(n):.3f}", flush=True)


def section_h(lib):
    all4 = lib[4]
    inv4 = [m for m in all4 if det(m) != 0]
    okall = True
    for q in (2, 4, 8, 16, 32, 64, 3, 9, 27, 5, 25, 7, 49, 11, 13, 15, 21, 35, 45, 63, 105, 17, 19, 23, 33, 55, 77, 81, 121, 125, 128):
        r1 = closure_size(all4, q, vec=(1, 0, 1))
        r2 = closure_size(inv4, q, vec=(1, 0, 1))
        exp_h = q // 2 if q % 2 == 0 else q
        full_states = q ** 3 // 2 if q % 2 == 0 else None
        ok = r1[1] == exp_h and r2[1] == exp_h and (full_states is None or (r1[0] == full_states and r2[0] == full_states))
        okall = okall and ok
        print(f"(h) q = {q}: all 56 -> {r1}, invertible 21 -> {r2}; expected h-residues {exp_h}"
              + (f", expected states {full_states}" if full_states else "") + f"; ok {ok}", flush=True)
    print("(h) all moduli complete:", okall)


def section_i(lib):
    for k in (4, 5, 6):
        hb = [m for m in lib[k] if m[1] == m[2] and m[4] == m[5] and m[7] == m[8]]
        evens = sum(1 for m in hb if (m[0] * (m[4] + m[7]) - m[1] * (m[3] + m[6])) % 2 == 0)
        # induced 2x2 on (h, R = E + O): [[m0, m1], [m3 + m6, m4 + m7]]
        print(f"(i) blocks of order <= {k} killing (0,1,-1): {len(hb)}; induced 2x2 on (h,R) with even determinant: {evens}; "
              f"singular among them: {sum(1 for m in hb if det(m) == 0)}")


def int_roots(m):
    t, c, d = tr(m), c2(m), det(m)
    roots = set()
    if d == 0:
        roots.add(0)
        disc = t * t - 4 * c
        if disc >= 0:
            s = math.isqrt(disc)
            if s * s == disc and (t + s) % 2 == 0:
                roots.add((t + s) // 2)
                roots.add((t - s) // 2)
    else:
        for x in range(1, abs(d) + 1):
            if d % x == 0:
                for lam in (x, -x):
                    if lam ** 3 - t * lam ** 2 + c * lam - d == 0:
                        roots.add(lam)
    return roots


def cross(a, b):
    return (a[1] * b[2] - a[2] * b[1], a[2] * b[0] - a[0] * b[2], a[0] * b[1] - a[1] * b[0])


def prim(v):
    g = math.gcd(math.gcd(abs(v[0]), abs(v[1])), abs(v[2]))
    if g == 0:
        return None
    v = tuple(x // g for x in v)
    for x in v:
        if x != 0:
            if x < 0:
                v = tuple(-y for y in v)
            break
    return v


def left_eigs(m):
    """list of (lam, n) with n primitive integer, n^T N = lam n^T (rank-2 case), and pencils (lam, None, c)."""
    res = []
    for lam in int_roots(m):
        A = [m[0] - lam, m[1], m[2], m[3], m[4] - lam, m[5], m[6], m[7], m[8] - lam]
        cols = [(A[0], A[3], A[6]), (A[1], A[4], A[7]), (A[2], A[5], A[8])]   # columns of N - lam: n must be orthogonal to all
        found = None
        for i in range(3):
            for j in range(i + 1, 3):
                c = cross(cols[i], cols[j])
                if c != (0, 0, 0):
                    found = prim(c)
        if found is not None:
            assert all(sum(found[i] * (m[3 * i + j]) for i in range(3)) == lam * found[j] for j in range(3))
            res.append((lam, found, None))
        else:
            nz = [c for c in cols if c != (0, 0, 0)]
            res.append((lam, None, prim(nz[0]) if nz else None))
    return res


def small_orbit(gens, lim):
    seen = {(1, 0, 1)}
    stack = [(1, 0, 1)]
    while stack:
        h, e, o = stack.pop()
        for m in gens:
            h2 = m[0] * h + m[1] * e + m[2] * o
            if h2 > lim:
                continue
            s = (h2, m[3] * h + m[4] * e + m[5] * o, m[6] * h + m[7] * e + m[8] * o)
            if s not in seen:
                seen.add(s)
                stack.append(s)
    return seen


def section_j(lib):
    ms = lib[6]
    inv = [m for m in ms if det(m) != 0]
    sing = [m for m in ms if det(m) == 0]
    s0 = (1, 0, 1)
    planes = {}
    pencils = []
    for m in inv:
        for lam, n, c in left_eigs(m):
            if n is None:
                pencils.append((m, lam, c))
                continue
            planes.setdefault(n, []).append((m, lam))
    print(f"(j) invertible blocks of order <= 6: {len(inv)}; with a rational invariant plane: "
          f"{len(set(m for v in planes.values() for m, _ in v))}; distinct rational invariant planes: {len(planes)}; "
          f"pencil cases (N - lam of rank 1): {len(pencils)}")
    thru = {n: v for n, v in planes.items() if sum(a * b for a, b in zip(n, s0)) == 0}
    print(f"    planes through (1,0,1): {len(thru)}; shared by >= 2 invertible blocks: "
          f"{[(n, len(v)) for n, v in thru.items() if len(v) >= 2]}")
    shared = {n: v for n, v in planes.items() if len(v) >= 2}
    print(f"    ALL rational invariant planes shared by >= 2 invertible blocks (through the start state or not): {len(shared)}")
    orb = small_orbit(ms, 300)
    print(f"    (orbit states of the full order-<=6 semigroup with h <= 300: {len(orb)})")
    for n, v in sorted(shared.items(), key=lambda kv: -len(kv[1])):
        on = [s for s in orb if sum(a * b for a, b in zip(n, s)) == 0]
        # induced determinant on the plane = det / lam
        info = sorted(set((lam, det(m), Fraction(det(m), lam)) for m, lam in v))
        print(f"      normal {n}: {len(v)} invertible blocks; (lam, det N, det on plane) = {[(a, b, str(c)) for a, b, c in info]}; "
              f"orbit states (h <= 300) on the plane: {len(on)}" + (f", e.g. {sorted(on)[:3]}" if on else ""))
    # the plane h = O
    n0 = (1, 0, -1)
    pres = [(m, lam) for m, lam in planes.get(n0, [])]
    into = [m for m in sing if (m[0] - m[6], m[1] - m[7], m[2] - m[8]) == (0, 0, 0)]
    pres_sing = [m for m in sing if any(n == n0 for lam, n, c in left_eigs(m) if n is not None)]
    print(f"    plane h = O: preserved by {len(pres)} invertible blocks; singular blocks with image inside it (row h = row O): {len(into)}; "
          f"singular blocks having (1,0,-1) as a left eigenvector (any lam): {len(pres_sing)}")
    for m, lam in pres:
        # induced matrix on (h, E) with O = h:  h' = (m0 + m2) h + m1 E ; E' = (m3 + m5) h + m4 E
        A = (m[0] + m[2], m[1], m[3] + m[5], m[4])
        print(f"      block {m}: lam {lam}, det N {det(m)}, induced on (h,E): {A}, det {A[0] * A[3] - A[1] * A[2]}")
    ind = set()
    for m in into:
        ind.add((m[0] + m[2], m[1], m[3] + m[5], m[4]))
    dets = sorted(set(A[0] * A[3] - A[1] * A[2] for A in ind))
    print(f"      singular blocks into the plane: {len(ind)} distinct induced matrices; determinants all even or 0: "
          f"{all(d % 2 == 0 for d in dets)}; odd determinants: {[d for d in dets if d % 2]}")
    # unimodular blocks with eigenvalue +-1 planes (cyclic unimodular 2x2 systems anywhere)
    cnt = 0
    for n, v in planes.items():
        for m, lam in v:
            if abs(det(m)) == 1 and abs(lam) == 1:
                cnt += 1
                on = [s for s in orb if sum(a * b for a, b in zip(n, s)) == 0]
                print(f"      unimodular block {m} preserves plane normal {n} (lam {lam}); orbit states with h<=300 on it: {len(on)}")
    print(f"    unimodular blocks with a rational invariant plane (lam = +-1): {cnt}")
    # products of two blocks of order <= 4 (invertible ones): shared planes through any orbit point?
    inv4 = [m for m in lib[4] if det(m) != 0]
    words = set(inv4)
    for a in inv4:
        for b in inv4:
            words.add(mul(a, b))
    for a in inv4:
        for b in inv4:
            for c in inv4:
                words.add(mul(mul(a, b), c))
    pl2 = {}
    for m in words:
        for lam, n, c in left_eigs(m):
            if n is not None:
                pl2.setdefault(n, []).append((m, lam))
    orb4 = small_orbit(lib[4], 2000)
    sh = {n: v for n, v in pl2.items() if len(v) >= 2}
    print(f"    words of length <= 3 in the 21 invertible order-<=4 blocks: {len(words)} distinct matrices; rational invariant "
          f"planes shared by >= 2 of them: {len(sh)}")
    for n, v in sorted(sh.items(), key=lambda kv: -len(kv[1]))[:12]:
        on = [s for s in orb4 if sum(a * b for a, b in zip(n, s)) == 0]
        uni = [(m, lam) for m, lam in v if abs(det(m)) == abs(lam)]
        # do the plane-unimodular elements commute pairwise?
        noncomm = any(mul(a[0], b[0]) != mul(b[0], a[0]) for a in uni for b in uni)
        print(f"      normal {n}: {len(v)} words, {len(uni)} act unimodularly on the plane (|det| = |lam|), non-commuting pair among them: {noncomm}; "
              f"order-<=4 orbit states (h <= 2000) on the plane: {len(on)}")


def section_k(lib):
    inv = [m for m in lib[6] if det(m) != 0]
    cnt = 0
    for m in inv:
        # N (0,1,-1)^T = lam (0,1,-1)^T
        v = (m[1] - m[2], m[4] - m[5], m[7] - m[8])
        if v[0] == 0 and v[1] == -v[2] and v[1] in (1, -1):
            cnt += 1
    eig = sum(1 for m in inv if (m[1] - m[2]) == 0 and (m[4] - m[5]) == -(m[7] - m[8]))
    print(f"(k) invertible blocks of order <= 6 with (0,1,-1) as an eigenvector: {eig}; with eigenvalue +-1: {cnt}")


def main():
    lib = {k: load(k) for k in (3, 4, 5, 6)}
    which = sys.argv[1] if len(sys.argv) > 1 else "abcdefghijk"
    if "a" in which:
        section_a(lib)
    if "b" in which:
        section_b(lib)
    if "c" in which:
        section_c(lib)
    if "d" in which:
        section_d(lib)
    if "e" in which:
        section_efg(lib)
    if "h" in which:
        section_h(lib)
    if "i" in which:
        section_i(lib)
    if "j" in which:
        section_j(lib)
    if "k" in which:
        section_k(lib)


if __name__ == "__main__":
    main()
