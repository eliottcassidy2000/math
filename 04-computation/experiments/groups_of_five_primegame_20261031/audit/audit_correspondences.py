#!/usr/bin/env python3
"""Independent audit of the correspondence tests (P3-P6) and of section 4 / 6 claims of
groups_of_five_primegame_merged_metagraph_20261031.md.  Own code (imports only the sibling audit module).

Run from the repository root:
    python3 04-computation/experiments/groups_of_five_primegame_20261031/audit/audit_correspondences.py
"""
import itertools
import os
import random
import sys
from fractions import Fraction

HERE = os.path.dirname(os.path.abspath(__file__))
sys.dont_write_bytecode = True
sys.path.insert(0, HERE)
import audit_metagraph as AM  # noqa: E402  (own audit code)

FAILS = []


def check(name, ok, info=''):
    print(f"  [{'OK' if ok else 'FAIL'}] {name}" + (f"  -- {info}" if info else ''))
    if not ok:
        FAILS.append(name)


PR = [2, 3, 5, 7, 11, 13, 17, 19, 23, 29]
REG, STA = PR[:4], PR[4:]
HEAD = [(17, 91), (78, 85), (19, 51), (23, 38), (29, 33), (77, 29), (95, 23), (77, 19), (1, 17), (11, 13), (13, 11)]
LG = HEAD + [(15, 14), (15, 2), (55, 1)]
LO = HEAD + [(15, 2), (1, 7), (55, 1)]


def primes_of(x):
    return [p for p in PR if x % p == 0]


# ------------------------------------------------------------------ C(a) circulants and palindromic paths
def circ(q, S):
    return [[1 if (j - i) % q in S else 0 for j in range(q)] for i in range(q)]


def is_palindromic(A, path):
    n = len(path)
    rev = {path[i]: path[n - 1 - i] for i in range(n)}
    return all(A[u][w] == A[rev[w]][rev[u]] for u in range(n) for w in range(n) if u != w)


def all_ham_paths(A):
    n = len(A)
    out = []

    def ext(path, used):
        if len(path) == n:
            out.append(tuple(path))
            return
        for w in range(n):
            if not used >> w & 1 and A[path[-1]][w]:
                path.append(w)
                ext(path, used | 1 << w)
                path.pop()
    for v in range(n):
        ext([v], 1 << v)
    return out


def find_nonpal(A, tries=200, seed=1):
    """random DFS for a Hamiltonian path that is not palindromic; returns one or None."""
    n = len(A)
    rnd = random.Random(seed)
    for _ in range(tries):
        path = [rnd.randrange(n)]
        used = 1 << path[0]
        stack = [list(range(n))]
        rnd.shuffle(stack[0])
        while path:
            if len(path) == n:
                if not is_palindromic(A, path):
                    return path
                break
            cand = stack[-1]
            moved = False
            while cand:
                w = cand.pop()
                if not used >> w & 1 and A[path[-1]][w]:
                    path.append(w)
                    used |= 1 << w
                    nxt = list(range(n))
                    rnd.shuffle(nxt)
                    stack.append(nxt)
                    moved = True
                    break
            if not moved:
                stack.pop()
                used &= ~(1 << path.pop())
    return None


def conn_sets(q):
    half = [s for s in range(1, q) if s <= q - s] if q > 1 else []
    out = []
    for choice in itertools.product(*[(s, q - s) for s in range(1, (q - 1) // 2 + 1)]):
        out.append(frozenset(choice))
    return out if q > 1 else [frozenset()]


def Ca():
    print('=' * 100)
    print('C(a)  circulant tournaments: orders 1, 3, 5 all pure-blue, 7 none; order 17 (Fermat) none')
    print('=' * 100)
    for q in (1, 3, 5, 7):
        res = []
        for S in conn_sets(q):
            A = circ(q, S)
            paths = all_ham_paths(A)
            res.append(all(is_palindromic(A, p) for p in paths))
        if q <= 5:
            check(f'order {q}: all {len(res)} circulant tournaments pure-blue (every Hamiltonian path palindromic)',
                  all(res))
        else:
            check(f'order 7: none of the {len(res)} circulant tournaments is pure-blue', not any(res))
    found = 0
    sets17 = conn_sets(17)
    for S in sets17:
        if find_nonpal(circ(17, S)) is not None:
            found += 1
    check(f'order 17: each of the {len(sets17)} circulant tournaments has a non-palindromic Hamiltonian path '
          f'(none is pure-blue, although all are rigid: THM-871 / THM-4644)', found == len(sets17),
          f'{found}/{len(sets17)}')
    # multiplier orbits on connection sets (= isomorphism classes at prime order, Elspas-Turner)
    orbits = {}
    for q in (3, 5, 17):
        seen, cnt = set(), 0
        for S in conn_sets(q):
            if S in seen:
                continue
            cnt += 1
            for a in range(1, q):
                seen.add(frozenset(a * s % q for s in S))
        orbits[q] = cnt
    check('rotational classes (multiplier orbits): 1 at order 3, 1 at order 5, 16 at order 17 -- uniqueness at 3, 5 '
          'is not a Fermat-order property', orbits == {3: 1, 5: 1, 17: 16}, str(orbits))


# ------------------------------------------------------------------ C(b) Terras cycles
def Cb():
    print('=' * 100)
    print('C(b)  rational Collatz (Terras) cycles of the circulant parity words, q = 1, 3, 5, 7')
    print('=' * 100)
    found = {}
    for q in (1, 3, 5, 7):
        k = (q - 1) // 2
        Dq = 2 ** q - 3 ** k
        necks = set()
        for S in conn_sets(q):
            w = tuple(1 if i in S else 0 for i in range(q))
            necks.add(min(w[i:] + w[:i] for i in range(q)))
        cyc = []
        for w in sorted(necks):
            # solve x = T_w(x) with T(x) = x/2 (bit 0), (3x+1)/2 (bit 1): x = c / (2^q - 3^k)
            c = Fraction(0)
            a = Fraction(1)
            for b in w:   # T_w(x) = a x + c, composed left to right
                if b:
                    a, c = a * Fraction(3, 2), c * Fraction(3, 2) + Fraction(1, 2)
                else:
                    a, c = a / 2, c / 2
            x = c / (1 - a)
            orbit = []
            y = x
            for b in w:
                orbit.append(y)
                y = y / 2 if b == 0 else (3 * y + 1) / 2
            assert y == x
            ints = [int(v * Dq) for v in orbit]
            assert all(v * Dq == int(v * Dq) for v in orbit)
            # check directly: integer cycle of z -> z/2 (even), (3z + Dq)/2 (odd), parities = w
            z = ints[0]
            okc = True
            for i, b in enumerate(w):
                if z % 2 != b:
                    okc = False
                z = z // 2 if z % 2 == 0 else (3 * z + Dq) // 2
            okc &= (z == ints[0])
            cyc.append((sorted(ints), okc))
        found[q] = (Dq, cyc)
        print(f'     q = {q}: D = {Dq}; cycles {[c for c, _ in cyc]}; all integer cycles with the right parities: '
              f'{all(o for _, o in cyc)}')
    check('q=1: {0} with D=1; q=3: {1,2,4} with D=5; q=5: {5,10,19,20,40} and {7,11,14,22,28} with D=23; '
          'q=7: five cycles with D=101',
          found[1] == (1, [([0], True)]) and found[3] == (5, [([1, 2, 4], True)])
          and found[5][0] == 23 and sorted(c for c, _ in found[5][1]) == [[5, 10, 19, 20, 40], [7, 11, 14, 22, 28]]
          and found[7][0] == 101 and len(found[7][1]) == 5 and all(o for _, o in found[7][1]))


# ------------------------------------------------------------------ C(c) relabellings
def relabel(prog, mr, ms):
    mp = dict(zip(REG, mr))
    mp.update(zip(STA, ms))

    def img(x):
        y = 1
        for p in PR:
            if x % p == 0:
                y *= mp[p]
        return y
    return [(img(a), img(b)) for a, b in prog]


def Cc():
    print('=' * 100)
    print('C(c)  relabellings of registers {2,3,5,7} and states {11..29} of the owner list (4!*6! = 17280)')
    print('=' * 100)
    tot = ev_a = ev_b = ev_full = 0
    src = {}
    for mr in itertools.permutations(REG):
        for ms in itertools.permutations(STA):
            pr = relabel(LO, mr, ms)
            tot += 1
            hits = [j for j, (a, b) in enumerate(pr) if b == 23 and a in (95, 77)]
            if hits:
                ev_a += 1
                for j in hits:
                    src[LO[j]] = src.get(LO[j], 0) + 1
                j = hits[0]
                other = 77 if pr[j][0] == 95 else 95
                if any(a == other for k, (a, _) in enumerate(pr) if k != j):
                    ev_full += 1
            nums = [a for a, _ in pr]
            if 95 in nums and 77 in nums:
                ev_b += 1
    print(f'     total {tot}; event (a) {ev_a} = {Fraction(ev_a, tot)}; sources {src}; '
          f'event (b) {ev_b} = {Fraction(ev_b, tot)}; full coincidence {ev_full} = {Fraction(ev_full, tot)}')
    check('864 (= 1/20, exactly 5%) for "some fraction is 95/23 or 77/23"; three disjoint sources of 288 each',
          tot == 17280 and ev_a == 864 and sorted(src.values()) == [288, 288, 288])
    check('96 (0.556%) for "95 and 77 both numerators"; 48 (= 1/360) for the full coincidence',
          ev_b == 96 and ev_full == 48)
    check('P4(a) as pre-registered ("fewer than 5%") fails: 864/17280 is exactly 5%', not (ev_a / tot < 0.05))


# ------------------------------------------------------------------ C(d) look-elsewhere calibration
def scc_count(A):
    n = len(A)
    R = [[A[i][j] or i == j for j in range(n)] for i in range(n)]
    for k in range(n):
        for i in range(n):
            if R[i][k]:
                for j in range(n):
                    if R[k][j]:
                        R[i][j] = True
    comps = {frozenset(j for j in range(n) if R[i][j] and R[j][i]) for i in range(n)}
    return len(comps)


def cycles_k(A, k):
    n = len(A)
    c = 0
    for S in itertools.combinations(range(n), k):
        for rest in itertools.permutations(S[1:]):
            cy = (S[0],) + rest
            if all(A[cy[i]][cy[(i + 1) % k]] for i in range(k)):
                c += 1
    return c


def battery(A):
    n = len(A)
    sc = [sum(r) for r in A]
    H = AM.ham_paths_dp(A)
    aut = AM.aut_count(A)
    c3, c4, c5 = cycles_k(A, 3), cycles_k(A, 4), cycles_k(A, 5)
    return {'n': n, 'arcs': n * (n - 1) // 2, 'H': H, 'Aut': aut, 'tilings': H // aut, 'c3': c3, 'c4': c4, 'c5': c5,
            'all cycles': c3 + c4 + c5, 'transitive triples': 10 - c3, 'strong components': scc_count(A),
            'score variance*n': sum((s - 2) ** 2 for s in sc), 'max score': max(sc), 'min score': min(sc),
            'sum sq scores': sum(s * s for s in sc), 'Ham cycles': c5}


def Cd():
    print('=' * 100)
    print('C(d)  look-elsewhere calibration on (TT5, 1->C3->1, T5)')
    print('=' * 100)
    TT5 = [[1 if i < j else 0 for j in range(5)] for i in range(5)]
    M = [[1 if i < j else 0 for j in range(5)] for i in range(5)]
    for u, w in ((1, 2), (2, 3), (3, 1)):
        M[u][w], M[w][u] = 1, 0
    T5 = circ(5, {1, 2})
    pts = [TT5, M, T5]
    bats = [battery(A) for A in pts]
    names = list(bats[0])
    vals = {nm: tuple(b[nm] for b in bats) for nm in names}
    for nm in names:
        print(f'     {nm:20} {vals[nm]}')
    EXPECT = {'n': (5, 5, 5), 'arcs': (10, 10, 10), 'H': (1, 3, 15), 'Aut': (1, 3, 5), 'tilings': (1, 1, 3),
              'c3': (0, 1, 5), 'c4': (0, 0, 5), 'c5': (0, 0, 2), 'all cycles': (0, 1, 12),
              'transitive triples': (10, 9, 5), 'strong components': (5, 3, 1), 'score variance*n': (10, 8, 0),
              'max score': (4, 4, 2), 'min score': (0, 0, 2), 'sum sq scores': (30, 28, 20), 'Ham cycles': (0, 0, 2)}
    check('the 16 battery values listed in correspondence_tests.out are right', vals == EXPECT)
    mults = (Fraction(1, 2), Fraction(1), Fraction(2), Fraction(3))
    offs = (-1, 0, 1)

    def hits(t):
        return [(nm, m, o) for nm in names for m in mults for o in offs
                if all(m * v + o == x for v, x in zip(vals[nm], t))]
    for tn, t in (('F', (0, 2, 10)), ('S', (0, 1, 7)), ('U', (0, 1, 3)), ('tau', (2, 4, 12))):
        print(f'     target {tn} = {t}: {hits(t)}')
    hF, hS, hU, ht = hits((0, 2, 10)), hits((0, 1, 7)), hits((0, 1, 3)), hits((2, 4, 12))
    check('F = (0,2,10) matched only by 2*c3; S, U, tau matched by nothing',
          [(nm, m, o) for nm, m, o in hF] == [('c3', 2, 0)] and not hS and not hU and not ht)
    matched = [(a, b) for a in range(16) for b in range(a, 16) if hits((0, a, b))]
    print(f'     matched (a, b): {matched}')
    check('13 of the 136 targets (0,a,b), 0<=a<=b<=15, are matched (9.6%): P5 (predicted >= 20%) fails',
          len(matched) == 13)


# ------------------------------------------------------------------ C(e) graphs
def mono_exists(G, H, n=10):
    """is there a bijection f with every edge of G mapped to an edge of H?"""
    Ga = [set() for _ in range(n)]
    Ha = [set() for _ in range(n)]
    for a, b in G:
        Ga[a].add(b)
        Ga[b].add(a)
    for a, b in H:
        Ha[a].add(b)
        Ha[b].add(a)
    order = sorted(range(n), key=lambda v: -len(Ga[v]))
    f = {}
    used = set()

    def rec(i):
        if i == n:
            return True
        v = order[i]
        for w in range(n):
            if w in used or len(Ha[w]) < len(Ga[v]):
                continue
            if all((u not in f) or (f[u] in Ha[w]) for u in Ga[v]):
                f[v] = w
                used.add(w)
                if rec(i + 1):
                    return True
                del f[v]
                used.discard(w)
        return False
    return rec(0)


def Ce():
    print('=' * 100)
    print('C(e)  prime graphs vs the merged-metagraph support')
    print('=' * 100)
    idx = {p: i for i, p in enumerate(PR)}
    arcs = [(p, q) for a, b in LO for p in primes_of(b) for q in primes_of(a)]
    pairs = {frozenset(e) for e in arcs if e[0] != e[1]}
    loops = [e for e in arcs if e[0] == e[1]]
    check('owner-list register-flow graph: 24 arcs, 19 unordered pairs, no loops',
          len(arcs) == 24 and len(pairs) == 19 and not loops)
    co = set()
    for a, b in LO:
        sup = sorted(set(primes_of(a)) | set(primes_of(b)))
        for p, q in itertools.combinations(sup, 2):
            co.add((p, q))
    deg = {p: 0 for p in PR}
    for p, q in co:
        deg[p] += 1
        deg[q] += 1
    check('co-occurrence graph: 28 edges, degree sequence (3,3,5,6,6,6,6,7,7,7)',
          len(co) == 28 and sorted(deg.values()) == [3, 3, 5, 6, 6, 6, 6, 7, 7, 7])
    R = AM.analyse(5, labels=True)
    sup = set()
    lp = 0
    for (a, b, ca, cb, blue) in R['lines']:
        x, y = sorted((R['node'][ca], R['node'][cb]))
        if x == y:
            continue
        sup.add((x, y))
    loopsM = {R['node'][ca] for (a, b, ca, cb, blue) in R['lines'] if R['node'][ca] == R['node'][cb]}
    check('metagraph support: 18 edges plus loops at 2 nodes', len(sup) == 18 and len(loopsM) == 2)
    flowE = [tuple(idx[p] for p in e) for e in pairs]
    coE = [(idx[p], idx[q]) for p, q in co]
    mE = list(sup)
    print(f'     edge counts: flow support 19, co-occurrence 28, metagraph 18 -> no isomorphism (trivially)')
    print(f'     one-way adjacency-preserving bijections: metagraph->flow {mono_exists(mE, flowE)}, '
          f'flow->metagraph {mono_exists(flowE, mE)}, metagraph->co-occurrence {mono_exists(mE, coE)}')


# ------------------------------------------------------------------ section 4 densities and F = U + S
def densities(prog):
    res = [Fraction(0)] * len(prog)
    for bits in range(1 << 10):
        S = [PR[i] for i in range(10) if bits >> i & 1]
        w = Fraction(1)
        for i, p in enumerate(PR):
            w *= Fraction(1, p) if bits >> i & 1 else Fraction(p - 1, p)
        for k, (a, b) in enumerate(prog):
            if all(p in S for p in primes_of(b)):
                res[k] += w
                break
    return res


def FUS_brute(N):
    divs = [d for d in range(2, N) if N % d == 0]
    sqf = lambda d: all(d % (p * p) for p in range(2, int(d ** 0.5) + 1))
    isp = lambda d: d > 1 and all(d % p for p in range(2, int(d ** 0.5) + 1))
    return len(divs), sum(1 for d in divs if isp(d)), sum(1 for d in divs if sqf(d))


def FUS(N):
    """F = #{1<d<N, d|N}, U = #prime d in that set, S = #squarefree d in that set, from the factorization."""
    e = []
    x = N
    p = 2
    while p * p <= x:
        k = 0
        while x % p == 0:
            x //= p
            k += 1
        if k:
            e.append(k)
        p += 1
    if x > 1:
        e.append(1)
    tau = 1
    for k in e:
        tau *= k + 1
    w = len(e)
    prime = (e == [1])
    sqfree = all(k == 1 for k in e)
    return tau - 2, w - prime, 2 ** w - 1 - sqfree


def Cf():
    print('=' * 100)
    print('C(f)  section 4 densities (exact, by divisibility patterns) and F = U + S facts')
    print('=' * 100)
    dg, do = densities(LG), densities(LO)
    q = Fraction(1)
    for p in STA:
        q *= Fraction(p - 1, p)
    check('first eleven densities identical; 17/91 = 1/91, 13/11 = 129024/2800733', dg[:11] == do[:11]
          and dg[0] == Fraction(1, 91) and dg[10] == Fraction(129024, 2800733))
    check('15/14, 15/2, 55/1 = (138240, 829440, 967680)/2800733 and 15/2, 1/7, 55/1 = (967680, 138240, 829440)/2800733',
          dg[11:] == [Fraction(x, 2800733) for x in (138240, 829440, 967680)]
          and do[11:] == [Fraction(x, 2800733) for x in (967680, 138240, 829440)])
    check('q = prod_{11<=p<=29}(1-1/p) = 21288960/30808063; L_G = (q/14, 3q/7, q/2), L_O = (q/2, q/14, 3q/7)',
          q == Fraction(21288960, 30808063) and dg[11:] == [q / 14, 3 * q / 7, q / 2]
          and do[11:] == [q / 2, q / 14, 3 * q / 7] and sum(dg) == 1 == sum(do))
    check('fast F, U, S (from the factorization) = the definition by divisor enumeration, N <= 3000',
          all(FUS(N) == FUS_brute(N) for N in range(2, 3001)))
    bal = [N for N in range(2, 20001) if FUS(N)[0] == FUS(N)[1] + FUS(N)[2]]

    def shape(N):
        f = {}
        x = N
        p = 2
        while p * p <= x:
            while x % p == 0:
                f[p] = f.get(p, 0) + 1
                x //= p
            p += 1
        if x > 1:
            f[x] = f.get(x, 0) + 1
        return tuple(sorted(f.values(), reverse=True))
    check('F = U + S exactly at shapes p, p^3, p^2qr (N <= 20000)', {shape(N) for N in bal} == {(1,), (3,), (2, 1, 1)}
          and all(shape(N) in {(1,), (3,), (2, 1, 1)} for N in bal)
          and all((FUS(N)[0] == FUS(N)[1] + FUS(N)[2]) for N in range(2, 20001) if shape(N) in {(1,), (3,), (2, 1, 1)}))
    check('7, 17 balanced; 55 has F=2, U+S=4; among {15,16,17,23,25} only 17, 23 balanced; p^2qr: (F,U,S)=(10,3,7)',
          FUS(7) == (0, 0, 0) and FUS(17) == (0, 0, 0) and FUS(55) == (2, 2, 2)
          and [N for N in (15, 16, 17, 23, 25) if FUS(N)[0] == FUS(N)[1] + FUS(N)[2]] == [17, 23]
          and FUS(60) == (10, 3, 7))


def main():
    Ca()
    Cb()
    Cc()
    Cd()
    Ce()
    Cf()
    print()
    if FAILS:
        print('FAILED:', FAILS)
        return 1
    print('ALL AUDIT CHECKS PASSED')
    return 0


if __name__ == '__main__':
    sys.exit(main())
