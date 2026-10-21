"""Gluck's trivial set-stabilizers with the converse symmetry added (session opus-2026-10-21-S49).
Note: 05-knowledge/results/seven_point_pencil_fano_and_converse_gluck_20261021.md; THM-4743, HYP-9353.

Gluck (Can. J. Math. 35 (1983) 59-67), Corollary 1: a permutation group of odd order has a subset with trivial setwise
stabilizer; so every tournament T has one for Aut(T) (D(T) <= 2; Albertson-Collins, Meslem-Sopena, Lozano).  Here the
group is Aut^+-(T) = automorphisms and anti-automorphisms (isomorphisms T -> T^converse), of order 2|Aut| when T is
self-converse.  Gluck's corollary no longer applies.

  (L)  PROVED (Lemma): Aut is normal of index 2 in Aut^+-; an anti-automorphism sigma permutes the regular Aut-orbits on
       vertex subsets as an involution, and T has an Aut^+--regular subset iff sigma moves at least one regular Aut-orbit;
       the number of regular Aut^+--orbits is (number of regular Aut-orbits moved by sigma)/2.  Checked on every class.
  (G1) Paley tournaments P_q, q prime = 3 mod 4: closed form for the subsets with trivial stabilizer in Aut(P_q) (centre
       stabiliser C = squares) and in Aut^+-(P_q) = AGL(1, q) (C = F_q^*):
           #regular = 2^q - 2 - q sum_{d | |C|, d > 1} (-mu(d)) (2^(1 + (q-1)/d) - 2)
       (two different centres generate translations, so the families are disjoint); cross-checked by brute force for
       q <= 19.  Aut^+-(P_q) has a regular subset iff q not in {3, 7} (THM-4743; Gluck Lemma 6 / Theorem 1).
  (G2) P_27: Aut^+-(P_27) = AGammaL(1, 27) (order 2106); Gluck's set GF(3) u {x}, Tr(x) != 0, has trivial stabilizer.
  (G3) FINITE-EXACT census: every isomorphism class of tournaments with n <= 9 (classes from the S36 generator,
       04-computation/experiments/s36_audit_20261014); for each self-converse class the full group Aut^+- by
       backtracking, every subset tested.  Exceptions (no Aut^+--regular subset), the minimum size rho+- of a regular
       subset, and the distribution of regular-orbit counts.  Result: exactly C_3, Z_5{1,2}, P_7 and C_3[C_3].
  (G4) FINITE-EXACT, the substitution rule at n = 15: Z_5{1,2}[C_3] has no Aut^+--regular subset, C_3[Z_5{1,2}] has some;
       the numbers r(B) of regular Aut(B)-orbits on subsets are 2, 6, 4, 2 for B = C_3, Z_5{1,2}, P_7, C_3[C_3], and
       |Aut(A[B])| = |Aut(B)|^|A| |Aut(A)| for the products tested (the wreath structure assumed by the substitution lemma).
  (G5) FINITE-EXACT: every circulant tournament on Z_n, n odd, 3 <= n <= 15 (all 2^((n-1)/2) connection sets): the
       exceptional ones are, up to isomorphism, exactly C_3, Z_5{1,2}, P_7, C_3[C_3] = Z_9{1,3,4,7} and Z_5{1,2}[C_3];
       none at n = 11, 13; C_3[Z_5{1,2}] is not exceptional.
Prints ALL CHECKS PASSED.  Runtime about 1 minute."""
import sys, os, time, itertools
from math import factorial
import numpy as np

HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, os.path.join(HERE, '..', 's36_audit_20261014'))
from audit_tournaments_classes import pair_index, code_to_adj           # noqa: E402
from audit_tournaments_classes_n9 import canon_refined                  # noqa: E402

FAIL = []


def check(c, m):
    print(('  ok   ' if c else '  FAIL ') + m, flush=True)
    if not c:
        FAIL.append(m)


def mobius(n):
    res, d, m = 1, 2, n
    while d * d <= m:
        if m % d == 0:
            m //= d
            if m % d == 0:
                return 0
            res = -res
        d += 1
    if m > 1:
        res = -res
    return res


def paley_regular_count(q, C):
    """#subsets of Z_q with trivial stabilizer in Z_q x| C (C cyclic, order m), q prime."""
    m = len(C)
    N = sum(-mobius(d) * (2 ** (1 + (q - 1) // d) - 2) for d in range(2, m + 1) if m % d == 0)
    return 2 ** q - 2 - q * N


def brute_regular_count(q, mults):
    """Union of fixed subsets of prime-order elements x -> a x + b (a in mults) and translations."""
    bad = set()
    elems = [(a, b) for a in mults for b in range(q) if not (a == 1 and b == 0)]
    for a, b in elems:
        g = [(a * x + b) % q for x in range(q)]
        # cycles of g
        seen = [False] * q; cyc = []
        for x in range(q):
            if not seen[x]:
                c = 0; y = x
                while not seen[y]:
                    seen[y] = True; c |= 1 << y; y = g[y]
                cyc.append(c)
        for mask in range(1 << len(cyc)):
            S = 0
            for i, c in enumerate(cyc):
                if (mask >> i) & 1:
                    S |= c
            bad.add(S)
    return 2 ** q - len(bad)


def section_G1():
    print('(G1) Paley tournaments P_q, q prime = 3 mod 4')
    rows = []; ok = True
    for q in (3, 7, 11, 19, 23, 31, 43, 47, 59, 67, 71, 79, 83):
        sq = sorted({x * x % q for x in range(1, q)})
        allu = list(range(1, q))
        r_aut = paley_regular_count(q, sq)
        r_pm = paley_regular_count(q, allu)
        if q <= 19:
            ok &= r_aut == brute_regular_count(q, sq) and r_pm == brute_regular_count(q, allu)
        ok &= r_aut % (q * len(sq)) == 0 and r_pm % (q * (q - 1)) == 0
        rows.append((q, r_aut // (q * len(sq)), r_pm // (q * (q - 1))))
    print(f'     (q, regular Aut-orbits, regular Aut+--orbits on vertex subsets): {rows[:6]} ...')
    ex = [q for q, a, b in rows if b == 0]
    check(ok and ex == [3, 7] and all(a > 0 for q, a, b in rows),
          'Aut(P_q) always has regular subsets (Gluck); Aut+-(P_q) = AGL(1,q) has none exactly for q = 3, 7 among primes q <= 83; '
          'closed form = brute force for q <= 19 (P_7: 4 and 0 regular orbits) (FINITE-EXACT; THM-4743 by Gluck Lemma 6)')


def section_G2():
    print('(G2) P_27 and AGammaL(1, 27)')
    # F_27 = F_3[x]/(x^3 - x - 1); elements as 3-tuples (c0, c1, c2)
    def fmul(a, b):
        r = [0] * 5
        for i in range(3):
            for j in range(3):
                r[i + j] += a[i] * b[j]
        for k in (4, 3):                    # x^3 = x + 1
            c = r[k]; r[k] = 0
            r[k - 3] += c; r[k - 2] += c
        return tuple(x % 3 for x in r[:3])
    def fadd(a, b): return tuple((x + y) % 3 for x, y in zip(a, b))
    els = [(a, b, c) for a in range(3) for b in range(3) for c in range(3)]
    zero = (0, 0, 0); one = (1, 0, 0)
    nz = [e for e in els if e != zero]
    def fpow(a, k):
        r = one
        for _ in range(k):
            r = fmul(r, a)
        return r
    frob = lambda a: fpow(a, 3)
    squares = {fmul(a, a) for a in nz}
    idx = {e: i for i, e in enumerate(els)}
    def tr(a): return fadd(fadd(a, frob(a)), frob(frob(a)))
    x = next(e for e in nz if tr(e) != zero)
    S = {idx[e] for e in [zero, one, (2, 0, 0), x]}
    def stab_size(mults):
        cnt = 0
        for a in mults:
            for k in range(3):
                for b in els:
                    img = set()
                    for e in els:
                        if idx[e] in S:
                            y = e
                            for _ in range(k):
                                y = frob(y)
                            img.add(idx[fadd(fmul(a, y), b)])
                    cnt += img == S
        return cnt
    st_pm = stab_size(nz); st_aut = stab_size(sorted(squares))
    check(st_pm == 1 and st_aut == 1 and len(squares) == 13,
          'P_27: the set GF(3) u {x} with Tr(x) != 0 has trivial stabilizer in AGammaL(1,27) = Aut+-(P_27) (order 2106) and in '
          'Aut(P_27) (FINITE-EXACT; Gluck Lemma 5)')


def classes_upto(nmax):
    fact = np.array([factorial(k) for k in range(10)], dtype=np.int64)
    permoff = np.zeros(10, dtype=np.int64); chunks = []; off = 0
    for k in range(10):
        permoff[k] = off
        P = np.array(list(itertools.permutations(range(k))), dtype=np.int64).reshape(-1) if k > 0 else np.zeros(0, dtype=np.int64)
        chunks.append(P); off += len(P)
    permtab = np.concatenate(chunks)
    out = {1: ([[0]], np.array([0]), np.array([1]))}
    classes = [[0]]
    for n in range(2, nmax + 1):
        idxp = pair_index(n)
        cand = []
        for adj in classes:
            for S in range(1 << (n - 1)):
                new = list(adj) + [S]
                for v in range(n - 1):
                    if not (S >> v) & 1:
                        new[v] |= 1 << (n - 1)
                cand.append(new)
        cand = np.array(cand, dtype=np.int64)
        codes, au = canon_refined(n, cand, idxp, permtab, permoff, fact)
        uniq, first = np.unique(codes, return_index=True)
        classes = [code_to_adj(n, int(c)) for c in uniq]
        A = np.array(classes, dtype=np.int64)
        full = (1 << n) - 1
        conv = np.array([[full ^ int(a) ^ (1 << v) for v, a in enumerate(row)] for row in classes], dtype=np.int64)
        ccodes, _ = canon_refined(n, conv, idxp, permtab, permoff, fact)
        out[n] = (classes, uniq, au[first], ccodes)
    return out


def isomorphisms(n, src, dst):
    """All permutations phi with u -> v in src iff phi(u) -> phi(v) in dst (bitmask rows)."""
    sc_s = [bin(r).count('1') for r in src]; sc_d = [bin(r).count('1') for r in dst]
    order = sorted(range(n), key=lambda v: sc_s[v])
    res = []; phi = [-1] * n; used = [False] * n
    def bt(i):
        if i == n:
            res.append(tuple(phi)); return
        u = order[i]
        for x in range(n):
            if used[x] or sc_d[x] != sc_s[u]:
                continue
            ok = True
            for j in range(i):
                w = order[j]
                if ((src[u] >> w) & 1) != ((dst[x] >> phi[w]) & 1):
                    ok = False; break
            if ok:
                phi[u] = x; used[x] = True
                bt(i + 1)
                used[x] = False; phi[u] = -1
    bt(0)
    return res


def regular_mask(n, group):
    """Boolean array over all 2^n subsets: True where no non-identity element of group fixes the subset."""
    masks = np.arange(1 << n, dtype=np.int64)
    fixed = np.zeros(1 << n, dtype=bool)
    ident = tuple(range(n))
    for g in group:
        if tuple(g) == ident:
            continue
        img = np.zeros(1 << n, dtype=np.int64)
        for i in range(n):
            img |= ((masks >> i) & 1) << g[i]
        fixed |= img == masks
    return ~fixed


def lex(A, B):
    """Lexicographic product A[B] of tournaments given as bitmask rows: vertex (a, b) -> a |B| + b."""
    na, nb = len(A), len(B); n = na * nb
    T = [0] * n
    for a in range(na):
        for b in range(nb):
            u = a * nb + b
            for a2 in range(na):
                for b2 in range(nb):
                    v = a2 * nb + b2
                    if (a != a2 and (A[a] >> a2) & 1) or (a == a2 and (B[b] >> b2) & 1):
                        T[u] |= 1 << v
    return T


def circulant(n, S):
    return [sum(1 << ((x + d) % n) for d in S) for x in range(n)]


def act(phi, S):
    r = 0
    while S:
        b = S & -S; i = b.bit_length() - 1
        r |= 1 << phi[i]; S ^= b
    return r


def section_G3(nmax=9):
    print(f'(G3) all tournament classes n <= {nmax}: the converse-extended Gluck problem')
    t0 = time.time()
    data = classes_upto(nmax)
    exceptions = []; lemma_ok = True; summary = []
    for n in range(1, nmax + 1):
        if n == 1:
            continue
        classes, codes, auts, ccodes = data[n]
        sc_classes = [i for i in range(len(classes)) if codes[i] == ccodes[i]]
        full = (1 << n) - 1
        rho_hist = {}; nregs = []
        for i in sc_classes:
            T = classes[i]
            Tc = [full ^ a ^ (1 << v) for v, a in enumerate(T)]
            aut = isomorphisms(n, T, T); anti = isomorphisms(n, T, Tc)
            lemma_ok &= len(aut) == auts[i] and len(anti) == len(aut)
            regm = regular_mask(n, aut + anti)
            reg = [int(S) for S in np.nonzero(regm)[0]]
            # lemma: regular Aut-orbits and the action of one anti-automorphism on them
            regA = [int(S) for S in np.nonzero(regular_mask(n, aut))[0]]
            orbitA = {}
            for S in regA:
                key = min(act(g, S) for g in aut)
                orbitA.setdefault(key, S)
            sigma = anti[0]
            moved = sum(1 for key in orbitA if min(act(g, act(sigma, key)) for g in aut) != key)
            lemma_ok &= len(reg) == moved * len(aut)          # moved orbits pair up: each Aut+- orbit = two Aut-orbits
            nreg = len(reg) // (2 * len(aut))
            nregs.append(nreg)
            if not reg:
                exceptions.append((n, T, len(aut)))
                continue
            rho = min(bin(S).count('1') for S in reg)
            rho_hist[rho] = rho_hist.get(rho, 0) + 1
        summary.append((n, len(classes), len(sc_classes), dict(sorted(rho_hist.items())), min(nregs) if nregs else None))
        print(f'     n = {n}: {len(classes)} classes, {len(sc_classes)} self-converse; min size of an Aut+--regular subset '
              f'(histogram over self-converse classes): {dict(sorted(rho_hist.items()))}; fewest regular Aut+--orbits: '
              f'{min(nregs) if nregs else None}  [{time.time() - t0:.0f}s]', flush=True)
    counts_ok = [s[1] for s in summary] == [1, 2, 4, 12, 56, 456, 6880, 191536][:nmax - 1]
    sc_ok = [s[2] for s in summary] == [1, 2, 2, 8, 12, 88, 176, 2752][:nmax - 1]
    exc = [(n, T, a) for n, T, a in exceptions]
    print(f'     exceptions (n, adjacency rows, |Aut|): {exc}')
    def is_circ(n, T, S):
        return all(((T[x] >> ((x + d) % n)) & 1) == (d in S) for x in range(n) for d in range(1, n))
    names = []
    for n, T, a in exc:
        found = None
        if n > 7:
            names.append(None)
            continue
        for perm in itertools.permutations(range(n)):
            Tp = [0] * n
            for u in range(n):
                for w in range(n):
                    if (T[u] >> w) & 1:
                        Tp[perm[u]] |= 1 << perm[w]
            for S in ({1}, {1, 2}, {1, 2, 4}):
                if is_circ(n, Tp, S):
                    found = f'circulant Z_{n}{sorted(S)}'; break
            if found:
                break
        names.append(found)
    print(f'     identified: {names}')
    check(counts_ok and sc_ok and lemma_ok,
          'class counts 1, 2, 4, 12, 56, 456, 6880, 191536 (A000568) and self-converse counts 1, 2, 2, 8, 12, 88, 176, 2752 '
          '(A002785) for n = 2..9; |anti| = |Aut| for every self-converse class; the orbit-pairing lemma holds on every class')
    c33 = None
    if len(exc) == 4 and exc[3][0] == 9:
        C3 = circulant(3, {1}); L = lex(C3, C3)
        c33 = len(isomorphisms(9, exc[3][1], L)) > 0
    check(names[:3] == ['circulant Z_3[1]', 'circulant Z_5[1, 2]', 'circulant Z_7[1, 2, 4]'] and len(exc) == 4 and c33,
          'for n <= 9 the only tournaments whose automorphism/anti-automorphism group has no regular orbit on vertex subsets are '
          'C_3, the circulant Z_5{1,2}, the Paley tournament P_7 and the lexicographic product C_3[C_3] (FINITE-EXACT; '
          'evidence for HYP-9353)')


def section_G4():
    print('(G4) the substitution rule at n = 15')
    C3 = circulant(3, {1}); T5 = circulant(5, {1, 2}); P7 = circulant(7, {1, 2, 4}); C33 = lex(C3, C3)
    def info(T):
        n = len(T); full = (1 << n) - 1
        Tc = [full ^ a ^ (1 << v) for v, a in enumerate(T)]
        aut = isomorphisms(n, T, T); anti = isomorphisms(n, T, Tc)
        rA = int(regular_mask(n, aut).sum()) // len(aut)
        rpm = int(regular_mask(n, aut + anti).sum()) if anti else None
        return len(aut), len(anti), rA, rpm
    rows = {}
    for name, T in (('C3', C3), ('Z5{1,2}', T5), ('P7', P7), ('C3[C3]', C33), ('Z5{1,2}[C3]', lex(T5, C3)),
                    ('C3[Z5{1,2}]', lex(C3, T5))):
        rows[name] = info(T)
        print(f'     {name}: |Aut| = {rows[name][0]}, |anti| = {rows[name][1]}, regular Aut-orbits r = {rows[name][2]}, '
              f'Aut+--regular subsets = {rows[name][3]}', flush=True)
    ok = [rows[x][2] for x in ('C3', 'Z5{1,2}', 'P7', 'C3[C3]')] == [2, 6, 4, 2]
    ok &= rows['Z5{1,2}[C3]'][3] == 0 and rows['C3[Z5{1,2}]'][3] > 0 and rows['C3[C3]'][3] == 0
    ok &= rows['Z5{1,2}[C3]'][0] == 3 ** 5 * 5 and rows['C3[Z5{1,2}]'][0] == 5 ** 3 * 3 and rows['C3[C3]'][0] == 3 ** 4
    check(ok, 'r = 2, 6, 4, 2 for C_3, Z_5{1,2}, P_7, C_3[C_3]; Z_5{1,2}[C_3] (n = 15) has NO Aut+--regular subset while '
              'C_3[Z_5{1,2}] has some; the automorphism groups are the wreath products (FINITE-EXACT; the substitution rule '
              'predicts both: 2 colours cannot distinguish D_10, 6 colours can distinguish S_3)')


def section_G5():
    print('(G5) all circulant tournaments on Z_n, n odd, n <= 15')
    C3 = circulant(3, {1}); T5 = circulant(5, {1, 2})
    family = {3: [circulant(3, {1})], 5: [T5], 7: [circulant(7, {1, 2, 4})], 9: [lex(C3, C3)], 15: [lex(T5, C3)]}
    ok = True; table = []
    for n in range(3, 16, 2):
        full = (1 << n) - 1
        exc_sets = []
        for bits in range(1 << ((n - 1) // 2)):
            S = {d if (bits >> (d - 1)) & 1 else n - d for d in range(1, (n - 1) // 2 + 1)}
            T = circulant(n, S)
            Tc = [full ^ a ^ (1 << v) for v, a in enumerate(T)]
            aut = isomorphisms(n, T, T); anti = isomorphisms(n, T, Tc)
            if not regular_mask(n, aut + anti).any():
                exc_sets.append(sorted(S))
        fam = family.get(n, [])
        all_in_family = all(any(isomorphisms(n, circulant(n, set(S)), F) for F in fam) for S in exc_sets)
        ok &= all_in_family and (bool(exc_sets) == bool(fam))
        table.append((n, len(exc_sets), exc_sets[:2]))
        print(f'     n = {n}: {1 << ((n - 1) // 2)} connection sets, {len(exc_sets)} exceptional, e.g. {exc_sets[:2]}; '
              f'all isomorphic to the predicted family member: {all_in_family}', flush=True)
    check(ok, 'among circulant tournaments on Z_n (n odd, n <= 15) the exceptional ones are exactly C_3, Z_5{1,2}, P_7, '
              'C_3[C_3] and Z_5{1,2}[C_3] up to isomorphism; none at n = 11, 13 (FINITE-EXACT; evidence for HYP-9353)')


if __name__ == '__main__':
    t0 = time.time()
    section_G1(); section_G2(); section_G3(9); section_G4(); section_G5()
    print(f'\n{"ALL CHECKS PASSED" if not FAIL else "FAILURES: " + str(FAIL)}  ({time.time() - t0:.0f}s)')
