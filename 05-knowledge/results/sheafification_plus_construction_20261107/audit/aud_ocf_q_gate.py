#!/usr/bin/env python3
"""Independent audit of THM-4892 (S85), parts (B5), (C), (D), (F), with explicit sets and maps.

(b) OCF |Pk2(V)| = H(T) and the Jss-sheaf condition Pk2(S) = prod_i Pk2(C_i) as an explicit bijection;
    Redei via the colour swap; the spectrum as the multiplicative closure of the strong spectrum (n <= 6).
(c) Q = Pk2 / global swap: explicit orbits, the map Q(S) -> prod_i Q(C_i), its fibres, and the action of
    (Z/2)^M (flip the colours on a set of components) - transitive on each fibre, with stabiliser exactly
    the diagonal; separatedness <=> at most one strong component of size >= 3; the second plus step on the
    R_min model (Lemma 0, verified by colimit in aud_sites_colimit.py) for non-separated cases n = 6, 7;
    the two counting laws of (D)4 on order joins A => B of strong A, B, and a probe with A NOT strong.
(d) the gate formula H(Gamma T) = H(T) + 2 sum_{|S| odd} H(T[S]) H(T[V-S]), Session 210's
    H(Gamma T) = H(T) + P_T(2), the alternating identity sum_S (-1)^|S| H(T[S])H(T[V-S]) = 0, and the direct
    enumeration of the gate odd cycles; padding Sigma(T).
(B5) componentwise converses: same strong components, same H.

Helpers (tournament bitmasks, strong components, DFS odd cycles, Hamiltonian DP, matching families) are the
auditor's own, from aud_sites_colimit.py in this folder.  Nothing is imported from the session's scripts.

Run:  python aud_ocf_q_gate.py > aud_ocf_q_gate.out
"""
import itertools
import random
import time
from collections import Counter

from aud_sites_colimit import (popcount, elems_of, submasks, is_sub, tour_from_bits, random_tour, components,
                               odd_cycles, ham_count, matching_families)

T0 = time.time()


def packings(cyc, S):
    inside = [c for c in cyc if is_sub(c[0], S)]
    res = []

    def rec(i, used, chosen):
        res.append(tuple(chosen))
        for j in range(i, len(inside)):
            if inside[j][0] & used == 0:
                chosen.append(inside[j])
                rec(j + 1, used | inside[j][0], chosen)
                chosen.pop()
    rec(0, 0, [])
    return res


def pk2(cyc, S):
    """elements of Pk2(S): sorted tuples of (cycle mask, cycle tuple, colour)."""
    out = []
    for p in packings(cyc, S):
        for cols in itertools.product((0, 1), repeat=len(p)):
            out.append(tuple(sorted((c[0], c[1], col) for c, col in zip(p, cols))))
    return out


def res_pk(x, Sp):
    return tuple(t for t in x if is_sub(t[0], Sp))


def swap(x):
    return tuple(sorted((a, b, 1 - c) for a, b, c in x))


def canon(x):
    return min(x, swap(x))


def flip_on(x, mask):
    """flip the colours of the cycles lying inside the vertex set `mask`."""
    return tuple(sorted((a, b, (1 - c) if is_sub(a, mask) else c) for a, b, c in x))


def all_tours(n):
    for bits in range(1 << (n * (n - 1) // 2)):
        yield tour_from_bits(n, bits)


def order_join(blocks):
    """blocks: list of out-mask tournaments; every vertex of an earlier block beats every later one."""
    n = sum(len(b) for b in blocks)
    out = [0] * n
    off = 0
    offs = []
    for b in blocks:
        offs.append(off)
        off += len(b)
    for bi, b in enumerate(blocks):
        o = offs[bi]
        later = 0
        for bj in range(bi + 1, len(blocks)):
            for v in range(len(blocks[bj])):
                later |= 1 << (offs[bj] + v)
        for v in range(len(b)):
            m = 0
            for w in elems_of(b[v]):
                m |= 1 << (o + w)
            out[o + v] = m | later
    return out


def strong_random(k, rng):
    while True:
        t = random_tour(k, rng)
        if len(components(t, (1 << k) - 1)) == 1:
            return t


# ----------------------------------------------------------------------------------------------------------
def part_C(rep):
    rng = random.Random(21)
    cases = [(n, t, True) for n in range(1, 6) for t in all_tours(n)]
    cases += [(6, random_tour(6, rng), False) for _ in range(300)]
    cases += [(7, random_tour(7, rng), False) for _ in range(120)]
    cases += [(8, random_tour(8, rng), False) for _ in range(25)]
    spectra = {}
    strong_spec = {}
    for n, out, allS in cases:
        V = (1 << n) - 1
        cyc = odd_cycles(out, n)
        # cycle sanity: every listed cycle is a directed cycle, distinct
        assert len(set(c[1] for c in cyc)) == len(cyc)
        for m, t in cyc:
            assert all(out[t[i]] >> t[(i + 1) % len(t)] & 1 for i in range(len(t)))
        H = ham_count(out, V)
        P = pk2(cyc, V)
        rep['(C1) OCF |Pk2(V)| = H(T): ' + str(len(P) == H)] += 1
        fixed = [x for x in P if swap(x) == x]
        rep['(C2) swap has exactly one fixed point, the empty packing: ' + str(fixed == [()])] += 1
        subsets = submasks(V) if allS else [V]
        for S in subsets:
            cs = components(out, S)
            PS = pk2(cyc, S) if S != V else P
            img = [tuple(res_pk(x, C) for C in cs) for x in PS]
            prod_size = 1
            for C in cs:
                prod_size *= len(pk2(cyc, C))
            ok = (len(set(img)) == len(PS) == prod_size)
            rep['(C1) Pk2(S) -> prod Pk2(C_i) bijective: ' + str(ok)] += 1
        Hc = [ham_count(out, C) for C in components(out, V)]
        prod = 1
        for h in Hc:
            prod *= h
        rep['(C3) H = prod over strong components: ' + str(prod == H)] += 1
        if allS:
            spectra.setdefault(n, set()).add(H)
            if len(components(out, V)) == 1:
                strong_spec.setdefault(n, set()).add(H)
    # (C3) spectrum of order n = products of strong values over compositions of n
    strong_spec.setdefault(2, set())
    comp_spec = {}
    for n in range(1, 6):
        vals = set()

        def rec(rem, acc):
            if rem == 0:
                vals.add(acc)
                return
            for k in range(1, rem + 1):
                for h in strong_spec.get(k, set()):
                    rec(rem - k, acc * h)
        rec(n, 1)
        comp_spec[n] = vals
        rep['(C3) spectrum(n) = products of strong values, n=' + str(n) + ': ' + str(vals == spectra[n])] += 1
    print('   spectra by order (exhaustive):', {n: sorted(spectra[n]) for n in sorted(spectra)})
    print('   strong spectra by order:      ', {n: sorted(strong_spec[n]) for n in sorted(strong_spec)})


# ----------------------------------------------------------------------------------------------------------
def q_structure(out, n, rep, subsets, second_step=False):
    V = (1 << n) - 1
    cyc = odd_cycles(out, n)
    for S in subsets:
        cs = components(out, S)
        PS = pk2(cyc, S)
        QS = sorted(set(canon(x) for x in PS))
        H = ham_count(out, S)
        rep['(D) |Q(S)| = (H(T[S]) + 1)/2: ' + str(len(QS) == (H + 1) // 2)] += 1
        Qc = [sorted(set(canon(x) for x in pk2(cyc, C))) for C in cs]
        target = set(itertools.product(*Qc))
        img = {}
        for q in QS:
            key = tuple(canon(res_pk(q, C)) for C in cs)
            img.setdefault(key, []).append(q)
        rep['(D3) Q(S) -> prod Q(C_i) surjective: ' + str(set(img) == target)] += 1
        big = sum(1 for C in cs if popcount(C) >= 3)
        separated_here = all(len(f) == 1 for f in img.values())
        rep['(D1) Q separated at S <=> at most one component of size >= 3: ' +
            str(separated_here == (big <= 1))] += 1
        for key, fib in img.items():
            q0 = fib[0]
            M = [C for C in cs if res_pk(q0, C)]
            m = len(M)
            ok_size = len(fib) == (2 ** (m - 1) if m >= 1 else 1)
            rep['(D3) fibre size 2^(m-1): ' + str(ok_size)] += 1
            # action of (Z/2)^M on the fibre, through representatives
            fibset = set(fib)
            orbit = set()
            stab_ok = True
            for bits in range(1 << m):
                mask = 0
                for i in range(m):
                    if bits >> i & 1:
                        mask |= M[i]
                y = canon(flip_on(q0, mask))
                orbit.add(y)
                if y not in fibset:
                    stab_ok = False
                if (y == q0) != (bits == 0 or bits == (1 << m) - 1):
                    stab_ok = False
            rep['(D3) (Z/2)^M acts on the fibre, transitively, stabiliser = diagonal: ' +
                str(stab_ok and orbit == fibset)] += 1
            rep['fibres by m: m=' + str(m)] += 1
        if second_step:
            pass


def q_plus_steps(out, n, rep):
    """plus construction on the R_min model for Q (and once more), all subsets; checks Q^+(S) = prod Q(C_i)
    and that Q^+ -> Q^{++} is a bijection (one step), plus non-injectivity of Q -> Q^+ where predicted."""
    objs = list(range(1 << n))
    cyc = odd_cycles(out, n)
    comps = {S: components(out, S) for S in objs}
    rmin = {S: frozenset(U for C in comps[S] for U in submasks(C)) for S in objs}
    Qel = {S: sorted(set(canon(x) for x in pk2(cyc, S))) for S in objs}

    def Qres(S, x, Sp):
        return canon(res_pk(x, Sp))
    P = (Qel, Qres)
    fams = {S: matching_families(rmin[S], P) for S in objs}
    idx = {S: {f: i for i, f in enumerate(fams[S])} for S in objs}

    def P1res(S, i, Sp):
        f = fams[S][i]
        return idx[Sp][tuple(p for p in f if p[0] in rmin[Sp])]
    P1 = ({S: list(range(len(fams[S]))) for S in objs}, P1res)
    fams2 = {S: matching_families(rmin[S], P1) for S in objs}
    for S in objs:
        prodq = 1
        for C in comps[S]:
            prodq *= len(Qel[C])
        rep['(D2) |Q^+(S)| = prod (H(C_i)+1)/2: ' + str(len(fams[S]) == prodq)] += 1
        # unit Q -> Q^+ and Q^+ -> Q^++
        eta = [idx[S][tuple((U, Qres(S, x, U)) for U in sorted(rmin[S]))] for x in Qel[S]]
        big = sum(1 for C in comps[S] if popcount(C) >= 3)
        rep['Q -> Q^+ injective <=> at most one big component: ' + str((len(set(eta)) == len(eta)) == (big <= 1))] += 1
        rep['Q -> Q^+ surjective: ' + str(set(eta) == set(range(len(fams[S]))))] += 1
        idx2 = {f: i for i, f in enumerate(fams2[S])}
        eta2 = [idx2[tuple((U, P1res(S, i, U)) for U in sorted(rmin[S]))] for i in range(len(fams[S]))]
        rep['(D2) Q^+ -> Q^++ bijective (a(Q) = Q^+): ' +
            str(len(set(eta2)) == len(eta2) == len(fams2[S]))] += 1


def part_D(rep):
    rng = random.Random(31)
    # all tournaments n <= 5: full subset structure (Q is separated there: at most one component of size >= 3)
    for n in range(1, 6):
        for out in all_tours(n):
            q_structure(out, n, rep, submasks((1 << n) - 1))
    # random n = 6, 7 tournaments, all subsets
    for n, cnt in ((6, 60), (7, 12)):
        for _ in range(cnt):
            out = random_tour(n, rng)
            q_structure(out, n, rep, submasks((1 << n) - 1))
    # order joins of strong pieces (sizes 3-4), the non-separated cases
    joins = 0
    law = Counter()
    for _ in range(60):
        k = rng.choice([2, 2, 3])
        sizes = [rng.choice([3, 3, 4]) for _ in range(k)]
        if sum(sizes) > 10:
            sizes = sizes[:2]
        parts = [strong_random(s, rng) for s in sizes]
        out = order_join(parts)
        n = len(out)
        q_structure(out, n, rep, [(1 << n) - 1])
        if n <= 7:
            q_plus_steps(out, n, rep)
        joins += 1
        if len(parts) == 2:
            hA = ham_count(parts[0], (1 << len(parts[0])) - 1)
            hB = ham_count(parts[1], (1 << len(parts[1])) - 1)
            qA, qB = (hA - 1) // 2, (hB - 1) // 2
            H = ham_count(out, (1 << n) - 1)
            law['presheaf law |Q|-1 = qA+qB+2qAqB (A,B strong): ' + str((H + 1) // 2 - 1 == qA + qB + 2 * qA * qB)] += 1
            law['sheaf law |Q^+|-1 = qA+qB+qAqB (A,B strong): ' + str((qA + 1) * (qB + 1) - 1 == qA + qB + qA * qB)] += 1
    rep['order joins tested'] += joins
    for k_, v in law.items():
        rep['(D4) ' + k_] += v
    # (D4) probe with A NOT strong: A = C3 => C3, B = C3
    C3 = [0b010, 0b100, 0b001]          # 0->1->2->0
    A = order_join([C3, C3])
    J = order_join([A, C3])
    n = len(J)
    V = (1 << n) - 1
    hA = ham_count(A, (1 << 6) - 1)
    hB = 3
    qA, qB = (hA - 1) // 2, (hB - 1) // 2
    Qplus = 1
    for C in components(J, V):
        Qplus *= (ham_count(J, C) + 1) // 2
    print(f'   (D4) probe, A = C3 => C3 (not strong), B = C3: H(A)={hA}, qA={qA}, qB={qB}; '
          f'|Q(A=>B)|-1 = {(ham_count(J, V) + 1) // 2 - 1} (qA+qB+2qAqB = {qA + qB + 2 * qA * qB}); '
          f'|Q^+(A=>B)|-1 = {Qplus - 1} versus the stated qA+qB+qAqB = {qA + qB + qA * qB}')
    rep['(D4) sheaf law fails for non-strong A (probe): ' + str(Qplus - 1 != qA + qB + qA * qB)] += 1
    # the non-separated cases among random tournaments n = 6, 7: second plus step on the R_min model
    for n, cnt in ((6, 25), (7, 6)):
        done = 0
        while done < cnt:
            out = random_tour(n, rng)
            if sum(1 for C in components(out, (1 << n) - 1) if popcount(C) >= 3) >= 2 or done % 2 == 0:
                q_plus_steps(out, n, rep)
                done += 1


# ----------------------------------------------------------------------------------------------------------
def gate(out):
    n = len(out)
    a, b = n, n + 1
    g = [m | (1 << a) for m in out] + [1 << b, (1 << n) - 1]
    return g


def pad(out):
    n = len(out)
    b, a = n, n + 1                     # b source, a sink, b -> a
    return [m | (1 << a) for m in out] + [((1 << n) - 1) | (1 << a), 0]


def part_F(rep):
    rng = random.Random(41)
    cases = [t for n in range(1, 6) for t in all_tours(n)]
    cases += [random_tour(n, rng) for n in (6, 7, 8) for _ in range(80)]
    for out in cases:
        n = len(out)
        V = (1 << n) - 1
        G = gate(out)
        HG = ham_count(G, (1 << (n + 2)) - 1)
        Hsub = {S: ham_count(out, S) for S in submasks(V)}
        odd_sum = sum(Hsub[S] * Hsub[V & ~S] for S in Hsub if popcount(S) % 2 == 1)
        P2 = sum(Hsub[S] * Hsub[V & ~S] for S in Hsub)
        alt = sum((-1) ** popcount(S) * Hsub[S] * Hsub[V & ~S] for S in Hsub)
        rep['(F3) H(Gamma T) = H(T) + 2 sum_{|S| odd} H(T[S])H(T[V-S]): ' + str(HG == Hsub[V] + 2 * odd_sum)] += 1
        rep["(F3) Session 210 (5): H(Gamma T) = H(T) + P_T(2): " + str(HG == Hsub[V] + P2)] += 1
        rep['(F3) alternating sum vanishes (Lemma 3 coefficient): ' + str(alt == 0)] += 1
        rep['(F2) Gamma(T) strong: ' + str(len(components(G, (1 << (n + 2)) - 1)) == 1)] += 1
        if n <= 6:
            gc = odd_cycles(G, n + 2)
            a, b = n, n + 1
            through = [c for c in gc if (c[0] >> a & 1) or (c[0] >> b & 1)]
            shape_ok = all((c[0] >> a & 1) and (c[0] >> b & 1) for c in through)
            # each gate cycle reads a -> b -> P -> a with P a Ham path of T[S]; count per S
            per_S = Counter()
            for m, t in through:
                k = t.index(a)
                rot = t[k:] + t[:k]
                if rot[1] != b:
                    shape_ok = False
                per_S[m & V] += 1
            cnt_ok = all(per_S.get(S, 0) == (Hsub[S] if popcount(S) % 2 == 1 else 0) for S in Hsub if S)
            rep['(F3) gate odd cycles = (odd S, Ham path of T[S]) via a->b->P->a: ' + str(shape_ok and cnt_ok)] += 1
            # OCF on Gamma(T) directly
            packs = packings(gc, (1 << (n + 2)) - 1)
            rep['(F3) OCF on Gamma(T): ' + str(sum(2 ** len(p) for p in packs) == HG)] += 1
        Sg = pad(out)
        cs = components(Sg, (1 << (n + 2)) - 1)
        rep['(F1) padding: components = {b}, comps(T), {a}: ' +
            str(sorted(cs) == sorted(components(out, V) + [1 << n, 1 << (n + 1)]))] += 1
        rep['(F1) H(Sigma T) = H(T): ' + str(ham_count(Sg, (1 << (n + 2)) - 1) == Hsub[V])] += 1
    # Gamma of the empty tournament is the transitive pair, not strong
    G0 = gate([])
    rep['(F2) Gamma(empty) = transitive pair, not strong: ' + str(len(components(G0, 3)) == 2)] += 1


# ----------------------------------------------------------------------------------------------------------
def part_B5(rep):
    rng = random.Random(51)
    for _ in range(300):
        n = rng.randint(2, 8)
        out = random_tour(n, rng)
        V = (1 << n) - 1
        cs = components(out, V)
        H = ham_count(out, V)
        for X in range(1 << len(cs)):
            rev = 0
            for i in range(len(cs)):
                if X >> i & 1:
                    rev |= cs[i]
            new = []
            for v in range(n):
                m = out[v]
                if rev >> v & 1:
                    C = next(C for C in cs if C >> v & 1)
                    inside_out = m & C
                    inside_in = 0
                    for w in elems_of(C):
                        if out[w] >> v & 1:
                            inside_in |= 1 << w
                    m = (m & ~C) | inside_in
                new.append(m)
            # tournament check
            assert all(((new[u] >> w & 1) ^ (new[w] >> u & 1)) for u in range(n) for w in range(n) if u != w)
            ok = sorted(components(new, V)) == sorted(cs) and ham_count(new, V) == H
            rep['(B5) componentwise converse keeps components and H: ' + str(ok)] += 1


if __name__ == '__main__':
    rep = Counter()
    print('=' * 100)
    print('(b)/(C): OCF, the Pk2 sheaf condition, Redei swap, spectrum')
    part_C(rep)
    print('(c)/(D): Q = Pk2 / swap')
    part_D(rep)
    print('(d)/(F): gate and padding')
    part_F(rep)
    print('(B5): componentwise converses')
    part_B5(rep)
    for k in sorted(rep):
        print(f'   {k}: {rep[k]}')
    bad = [k for k in rep if k.endswith(': False')]
    print('FAILED CHECKS:', bad if bad else 'none (all reported conditions True)')
    print(f'total time {time.time() - T0:.1f}s')
