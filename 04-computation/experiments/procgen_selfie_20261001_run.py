#!/usr/bin/env python3
"""procgen_selfie_20261001_run.py -- selfie tournaments, arc-HP counts, shaving, A049313 partner.

Sections
  A  gauge theorem for selfie tournaments; loop-Walsh expansion of H
  B  loops as odd 1-cycles: the fixed-point OCF polynomial H(T;x) and its loop-set Redei analogue
  C  the owner's conjecture: arcs on no HP / arc-HP-count parities (classes N<=9, N=10 sample)
  D  shaved tournaments: HP-core, HP-blocking number, Redei shaving
  E  OPEN-Q-060: A049313 = Euler graphs with orientation-even automorphism group

Run:  nice python3 -u 04-computation/experiments/procgen_selfie_20261001_run.py
Prints to stdout only; ends with ALL CHECKS PASSED.
"""
import itertools
import math
import os
import sys
import time
from collections import Counter, defaultdict

import numpy as np

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
import procgen_selfie_20261001_lib as L  # noqa: E402

NCHECK = 0
T0 = time.time()


def check(cond, msg):
    global NCHECK
    if not cond:
        print('CHECK FAILED:', msg)
        sys.exit(1)
    NCHECK += 1


def hdr(s):
    print()
    print('=' * 78)
    print(s)
    print('=' * 78, flush=True)


def H_all_labeled(n):
    """H(T) for every labeled tournament on n vertices (index convention of the lib)."""
    if n == 1:
        return np.array([1], dtype=np.int64)
    out = L.run_c(['hlab', str(n)])
    arr = np.array(out.split(), dtype=np.int64)
    check(arr.size == 1 << (n * (n - 1) // 2), 'hlab size n=%d' % n)
    return arr


def gauge_data(n):
    P = L.pairs(n)
    bit = {p: b for b, p in enumerate(P)}
    tiles = [(i, j) for (i, j) in P if j >= i + 2]
    m = len(tiles)
    idx = np.arange(1 << m, dtype=np.int64)
    T_of_t = np.zeros(1 << m, dtype=np.int64)
    for k, (i, j) in enumerate(tiles):
        T_of_t |= ((idx >> k) & 1) << bit[(i, j)]
    cut = np.zeros(1 << n, dtype=np.int64)
    for Ls in range(1 << n):
        c = 0
        for b, (i, j) in enumerate(P):
            if ((Ls >> i) & 1) != ((Ls >> j) & 1):
                c |= 1 << b
        cut[Ls] = c
    return P, bit, tiles, T_of_t, cut


def order_at_minus1(p):
    """Multiplicity of y = -1 as a root of the integer polynomial sum p[k] y^k (p != 0)."""
    order = 0
    while len(p) > 1 and sum(c * (-1) ** k for k, c in enumerate(p)) == 0:
        d = len(p) - 1
        q = [0] * d
        q[d - 1] = p[d]
        for k in range(d - 1, 0, -1):
            q[k - 1] = p[k] - q[k]
        p = q
        order += 1
    return order


# --------------------------------------------------------------------------------------
def section_A():
    hdr('A. Selfie gauge theorem: (tiling t, loop set L) -> switch_L(T_t)')
    fact = math.factorial
    for n in range(2, 8):
        P, bit, tiles, T_of_t, cut = gauge_data(n)
        ne, m, full = len(P), len(tiles), (1 << n) - 1
        check(m + n == ne + 1, 'dimension identity C(n-1,2)+n = C(n,2)+1 at n=%d' % n)
        img = T_of_t[:, None] ^ cut[None, :]
        counts = np.bincount(img.ravel(), minlength=1 << ne)
        check(counts.size == 1 << ne and np.all(counts == 2), '2-to-1 onto labeled tournaments n=%d' % n)
        comp = full ^ np.arange(1 << n)
        check(np.array_equal(img, img[:, comp]), 'fibre {L, L^c} n=%d' % n)
        Ls = np.arange(1 << n, dtype=np.int64)
        for i in range(n - 1):
            b = bit[(i, i + 1)]
            rev = (img >> b) & 1
            der = ((Ls >> i) ^ (Ls >> (i + 1))) & 1
            check(np.array_equal(rev, np.broadcast_to(der, rev.shape)), 'base-path arc = derivative n=%d i=%d' % (n, i))
        H = H_all_labeled(n)
        hv = H[img]
        check(np.all(hv.sum(axis=1) == 2 * fact(n)), 'class sum of H = 2 n! over loop sets, n=%d' % n)
        W = L.wht(hv)  # W[t, A] = sum_L H(switch_L T_t) (-1)^{|L cap A|}
        pc = np.array([bin(a).count('1') for a in range(1 << n)])
        check(np.all(W[:, pc % 2 == 1] == 0), 'odd loop-Walsh coefficients vanish n=%d' % n)
        check(np.all(W[:, 0] == 2 * fact(n)), 'constant loop-Walsh coefficient n=%d' % n)
        dmax = 2 * (n // 3)
        check(np.all(W[:, pc > 2 * n / 3] == 0), 'loop-Walsh degree <= 2 floor(n/3), n=%d' % n)
        attained = np.mean(np.any(W[:, pc == dmax] != 0, axis=1)) if dmax > 0 else 1.0
        alt = (hv * ((-1) ** pc)[None, :]).sum(axis=1)
        check(np.all(alt == 0), 'sum_L (-1)^|L| H = 0, n=%d' % n)
        # loop enumerator Lambda_t(y) = sum_L H y^|L| is divisible by (1+y)^(n - 2 floor(n/3))
        r = n - dmax
        lam = np.zeros((hv.shape[0], n + 1), dtype=np.int64)
        for k in range(n + 1):
            lam[:, k] = hv[:, pc == k].sum(axis=1)
        exact_orders = Counter(order_at_minus1([int(v) for v in row]) for row in lam)
        check(min(exact_orders) >= r, 'loop enumerator divisible by (1+y)^%d at n=%d' % (r, n))
        print('n=%d: tilings=%d, selfies=%d = 2 x %d labeled; class sum 2n!=%d; loop-Walsh degree <= %d '
              '(attained by %.3f of tilings); (1+y)-order of loop enumerator: %s'
              % (n, 1 << m, (1 << m) * (1 << n), 1 << ne, 2 * fact(n), dmax, attained, dict(sorted(exact_orders.items()))))
    # constant-H switching classes: need N!/2^(N-1) odd, i.e. N a power of 2; none at N = 8
    rows = [ln.split() for ln in L.run_c(['list', '8'], feed=L.gentourng(8)).splitlines()]
    h315 = [r[0] for r in rows if r[1] == '315']
    P8 = L.pairs(8)
    feed = []
    for srep in h315:
        for Ls in range(1 << 7):
            feed.append(''.join(str(int(ch) ^ (((Ls >> i) & 1) ^ ((Ls >> j) & 1))) for ch, (i, j) in zip(srep, P8)))
    Hs = [int(ln.split()[1]) for ln in L.run_c(['list', '8'], feed=['printf', '%s\\n'] + feed).splitlines()]
    const = [srep for k, srep in enumerate(h315) if all(h == 315 for h in Hs[128 * k:128 * (k + 1)])]
    check(len(Hs) == 128 * len(h315) and not const, 'no constant-H switching class at n=8')
    print('n=8: %d iso classes have H = 315 = 8!/2^7; none has H constant on its switching class' % len(h315))
    # A3: every ordering is gauge-fixable by exactly two loop sets; forward-edge distribution
    for n in range(2, 6):
        P, bit, tiles, T_of_t, cut = gauge_data(n)
        ne = len(P)
        X = np.arange(1 << ne, dtype=np.int64)
        ak = np.zeros((1 << ne, n), dtype=np.int64)
        for pi in itertools.permutations(range(n)):
            mask = val = 0
            fwd = np.zeros(1 << ne, dtype=np.int64)
            for k in range(n - 1):
                a, b = pi[k], pi[k + 1]
                bb = bit[(min(a, b), max(a, b))]
                need = 1 if a < b else 0
                mask |= 1 << bb
                val |= need << bb
                fwd += (((X >> bb) & 1) == need).astype(np.int64)
            ok = ((X[:, None] ^ cut[None, :]) & mask) == val
            check(np.all(ok.sum(axis=1) == 2), 'ordering fixed by exactly two loop sets n=%d' % n)
            ak[X, fwd] += 1
        # sum over loop sets of a_k(switch_L T) = 2 n! C(n-1,k)
        target = np.array([2 * math.factorial(n) * math.comb(n - 1, k) for k in range(n)])
        for x in range(1 << ne):
            s = ak[x ^ cut].sum(axis=0)
            check(np.array_equal(s, target), 'switching sum of a_k binomial n=%d' % n)
        euler = [sum((-1) ** j * math.comb(n + 1, j) * (k + 1 - j) ** n for j in range(k + 2)) for k in range(n)]
        check(np.array_equal(ak[0][::-1], np.array(euler)) or np.array_equal(ak[0], np.array(euler)),
              'transitive a_k = Eulerian n=%d' % n)
        print('n=%d: every ordering is a HP of exactly 2 selfies per class; sum over the 2^n loop sets of a_k = 2 n! C(n-1,k), '
              'i.e. per switching class n! C(n-1,k) = %s; per-class sum of the deformation a_k - A(n,k): %s'
              % (n, [int(v) for v in target // 2], [math.factorial(n) * math.comb(n - 1, k) - 2 ** (n - 1) * euler[k] for k in range(n)]))
    # A6: segment formula and the Ising form for n <= 5
    for n in (3, 4, 5):
        P, bit, tiles, T_of_t, cut = gauge_data(n)
        H = H_all_labeled(n)
        perms = list(itertools.permutations(range(n)))
        sample = T_of_t if n == 5 else np.arange(1 << len(P))
        for x in sample:
            x = int(x)
            Wx = L.wht(H[x ^ cut])
            A = L.adj_from_index(x, n)

            def s(a, b):
                return 1 if A[a][b] else -1

            for Aset in range(1 << n):
                if bin(Aset).count('1') % 2:
                    continue
                tot = 0
                for pi in perms:
                    pos = [k for k in range(n) if (Aset >> pi[k]) & 1]
                    prod = 1
                    for q in range(0, len(pos), 2):
                        for k in range(pos[q], pos[q + 1]):
                            prod *= s(pi[k], pi[k + 1])
                    tot += prod
                check(2 * tot == Wx[Aset], 'segment formula n=%d x=%d A=%d' % (n, x, Aset))
            for a in range(n):
                for b in range(a + 1, n):
                    others = [w for w in range(n) if w not in (a, b)]
                    J = 0
                    for k in range(1, n - 1, 2):
                        for seq in itertools.permutations(others, k):
                            walk = (a,) + seq + (b,)
                            pr = 1
                            for q in range(len(walk) - 1):
                                pr *= s(walk[q], walk[q + 1])
                            J += math.factorial(n - k - 1) * pr
                    # J_ab = 2^{2-n} * J ;  W[{a,b}] = 2^n J_ab = 4 J
                    check(4 * J == Wx[(1 << a) | (1 << b)], 'Ising coupling formula n=%d' % n)
        print('n=%d: segment formula for every loop-Walsh coefficient and the Ising coupling formula verified (%d tournaments)'
              % (n, len(sample)))


# --------------------------------------------------------------------------------------
def section_B():
    hdr('B. Loops as odd 1-cycles: H(T;x) = sum_S 2^|S| x^(n-|V(S)|) = sum_U (x-1)^(n-|U|) H(T[U])')
    for n in range(1, 8):
        reps = L.classes(n) if n >= 2 else ['']
        vals3, valsm1, valsm1_sign, val0 = Counter(), Counter(), Counter(), Counter()
        for srep in reps:
            A = L.adj_from_string(srep, n) if n >= 2 else [[0]]
            cyc = L.odd_cycles(A)
            cov = L.indep_sets_by_cover(cyc)
            Pdir = [0] * (n + 1)
            for vs, sizes in cov.items():
                for sz in sizes:
                    Pdir[n - len(vs)] += 2 ** sz
            hs = L.all_sub_H(A)
            Hsub = {frozenset(v for v in range(n) if (S >> v) & 1): hs[S] for S in range(1 << n)}
            Pind = [0]
            for U, h in Hsub.items():
                Pind = L.poly_add(Pind, [h * c for c in L.poly_pow([-1, 1], n - len(U))])
            check(L.trim(Pind) == L.trim(Pdir), 'fixed-point OCF identity n=%d %s' % (n, srep))
            H = Hsub[frozenset(range(n))]
            check(sum(Pdir) == H, 'H(T;1) = H(T)')
            # Grinberg-Stanley normalisation: U_T(1,1) = sum over ordered 2-block covers
            u2 = sum(Hsub[U] * Hsub[frozenset(range(n)) - U] for U in Hsub)
            gs = sum(4 ** sz * 2 ** (n - len(vs)) for vs, sizes in cov.items() for sz in sizes)
            check(u2 == gs, 'U_T(1^2) = sum_S (2k)^|S| k^(n-|V(S)|) at k=2')
            npaths = sum(Hsub.values())
            check(npaths == sum(c * 2 ** k for k, c in enumerate(Pdir)), 'H(T;2) = number of paths incl. empty')
            # loop-set Redei analogue
            for Ls in range(1 << n):
                Lset = [v for v in range(n) if (Ls >> v) & 1]
                Z = 0
                for r in range(len(Lset) + 1):
                    for M in itertools.combinations(Lset, r):
                        Z += 2 ** r * Hsub[frozenset(range(n)) - frozenset(M)]
                direct = sum(2 ** sz * 3 ** len(set(Lset) - vs) for vs, sizes in cov.items() for sz in sizes)
                check(Z == direct, 'loop OCF = I(Omega + loops, 2)')
                check((Z - H - 2 * len(Lset)) % 4 == 0, 'Z_L = H + 2|L| mod 4')
            v3 = sum(c * 3 ** k for k, c in enumerate(Pdir))
            vm1 = sum(c * (-1) ** k for k, c in enumerate(Pdir))
            Im2 = sum((-2) ** sz for vs, sizes in cov.items() for sz in sizes)
            check(vm1 == (-1) ** n * Im2, 'H(T;-1) = (-1)^n I(Omega,-2)')
            vals3[v3] += 1
            valsm1[Im2] += 1
            val0[Pdir[0]] += 1
        print('n=%d (%d classes): H(T;3)=I(Omega+all loops,2) values %s' % (n, len(reps), sorted(vals3)[:12]))
        print('       I(Omega,-2) = (-1)^n H(T;-1) value range [%d, %d], sign counts %s; H(T;0) (odd-cycle partitions) values %s'
              % (min(valsm1), max(valsm1), dict(Counter(int(np.sign(v)) for v in valsm1.elements())), sorted(val0)))
        if n <= 4:
            for srep in reps:
                A = L.adj_from_string(srep, n) if n >= 2 else [[0]]
                cov = L.indep_sets_by_cover(L.odd_cycles(A))
                Pdir = [0] * (n + 1)
                for vs, sizes in cov.items():
                    for sz in sizes:
                        Pdir[n - len(vs)] += 2 ** sz
                print('       T=%s  H(T;x) coefficients (x^0..x^n) = %s' % (srep or '-', Pdir))



# --------------------------------------------------------------------------------------
def parse_summary(text):
    """Parse the SUMMARY block of the C engine into a dict."""
    d = {}
    for line in text.splitlines():
        line = line.strip()
        if line.startswith('SUMMARY'):
            for tok in line.split()[1:]:
                k, v = tok.split('=')
                d[k] = float(v)
        elif 'classes=' in line and 'labeled=' in line and not line.startswith(('ALL', 'STRONG')):
            name = line.split()[0]
            vals = dict(t.split('=') for t in line.split()[1:])
            d[name] = (float(vals['classes']), float(vals['labeled']))
        elif line.startswith(('zero-arc', 'odd-arc histogram', 'odd-arc-graph components', 'odd-arc-graph isolated')):
            key = line.split('(')[0].strip()
            hist = {}
            for tok in line.split('): ', 1)[1].split():
                k, rest = tok.split(':')
                c, lab = rest.split('/')
                hist[int(k)] = (float(c), float(lab))
            d[key] = hist
        elif line.startswith('max_H='):
            head, tail = line.split(':', 1)
            kv = dict(t.split('=') for t in head.split()[:2])
            d['max_H'] = int(kv['max_H'])
            d['maximizer_odd'] = sorted(int(t) for t in tail.split())
        elif line.startswith('identity_failures'):
            for tok in line.split():
                k, v = tok.split('=')
                d[k] = int(v)
    return d


ALLPOS_LAB = {3: 2, 4: 40, 5: 664, 6: 26048}


def section_C():
    hdr('C. The owner\'s conjecture: arcs on no HP, and arc-HP-count parities')
    summaries = {}
    for n in range(3, 10):
        t = time.time()
        out = L.run_c(['stats', str(n)], feed=L.gentourng(n))
        d = parse_summary(out)
        summaries[n] = d
        ne = n * (n - 1) // 2
        check(d['labeled'] == 2 ** ne, 'Aut-weighted labeled total n=%d' % n)
        check(d['identity_failures'] == 0 and d['lemma_failures'] == 0, 'per-vertex identities n=%d' % n)
        if n in ALLPOS_LAB:
            check(d['all_arcs_on_some_HP'][1] == ALLPOS_LAB[n], 'orchestrator all-arcs-on-HP count n=%d' % n)
        odd = d['all_arcs_odd']
        even = d['all_arcs_even']
        if n % 4 in (0, 3):
            check(odd[0] == 0, 'all-odd impossible n=0,3 mod 4 (n=%d)' % n)
        if n % 2 == 0:
            check(even[0] == 0, 'all-even impossible for even n (n=%d)' % n)
        hist = d['odd-arc histogram']
        check(all(k % 2 == (n - 1) % 2 for k in hist), '#odd arcs = n-1 mod 2 (n=%d)' % n)
        print('n=%d: classes=%d labeled=2^%d | all arcs on some HP: %d classes / %d labeled (%.4f) | '
              'strong & some dead arc: %d classes | all-odd: %d/%d | all-even: %d/%d | all-equal: %d | '
              'min/max #odd arcs: %d/%d of %d | odd-arc graph components: %s  [%.1fs]'
              % (n, d['classes'], ne, d['all_arcs_on_some_HP'][0], d['all_arcs_on_some_HP'][1],
                 d['all_arcs_on_some_HP'][1] / 2 ** ne,
                 d['strong'][0] - d['strong_and_all_on_HP'][0], odd[0], odd[1], even[0], even[1],
                 d['all_arcs_equal'][0], min(hist), max(hist), ne,
                 {k: int(v[0]) for k, v in d['odd-arc-graph components'].items()}, time.time() - t), flush=True)
        if n % 2 == 0:
            check(min(hist) == n - 1, 'even n: minimum #odd arcs = n-1 (n=%d)' % n)
        print('       H-maximizers: max H = %d, %d class(es), #odd arcs of each: %s'
              % (d['max_H'], len(d['maximizer_odd']), d['maximizer_odd']))
    for n in range(3, 10):
        line = [ln for ln in L.run_c(['oddfilter', str(n), '-1'], feed=L.gentourng(n)).splitlines()
                if ln.startswith('ODDFILTER')][0]
        hist2 = {int(k): int(v) for k, v in (tok.split(':') for tok in line.split('odd_hist:')[1].split())}
        hist1 = {k: int(v[0]) for k, v in summaries[n]['odd-arc histogram'].items()}
        check(hist1 == hist2 and 'H_even_failures=0' in line, 'independent mod-2 DP reproduces the #odd histogram n=%d' % n)
    print('independent mod-2 bitmask DP reproduces every class #odd-arc histogram, n = 3..9')
    check(summaries[6]['all_arcs_odd'] == (1, 240), 'n=6 all-odd = 1 class / 240 labeled')
    check(summaries[5]['all_arcs_odd'][0] == 0 and summaries[9]['all_arcs_odd'][0] == 0, 'no all-odd at n=5,9')
    print('all-even classes n=3,5,7,9:', [int(summaries[n]['all_arcs_even'][0]) for n in (3, 5, 7, 9)])
    # labeled cross-check n<=6
    for n in range(3, 7):
        d = parse_summary(L.run_c(['lab', str(n)]))
        for key in ('all_arcs_on_some_HP', 'all_arcs_odd', 'all_arcs_even', 'strong', 'some_arc_on_every_HP'):
            check(d[key][1] == summaries[n][key][1], 'labeled enumeration agrees with class weights n=%d %s' % (n, key))
    print('labeled direct enumeration n=3..6 agrees with Aut-weighted class counts')
    # structural lemma on dead arcs, n<=7
    for n in range(3, 8):
        for srep in L.classes(n):
            A = L.adj_from_string(srep, n)
            H, c, st, en = L.arc_counts(A)
            z = sum(1 for v in c.values() if v == 0)
            comps = L.strong_components_ordered(A)
            pred = 0
            for i in range(len(comps)):
                for j in range(i + 2, len(comps)):
                    pred += len(comps[i]) * len(comps[j])
            for C in comps:
                if len(C) >= 3:
                    Hc, cc, _, _ = L.arc_counts(L.sub_adj(A, sorted(C)))
                    pred += sum(1 for v in cc.values() if v == 0)
            check(z == pred, 'dead-arc formula n=%d %s' % (n, srep))
    print('dead-arc formula z(T) = sum_{j>=i+2}|C_i||C_j| + sum_i z(T[C_i]) verified for all classes n<=7')
    # Cayley theorem: every Cayley tournament on an odd abelian group is all-even
    fams = []
    for n in (3, 5, 7, 9, 11, 13):
        fams.append(('circulants Z_%d' % n, n, L.circulant_strings(n)))
    fams.append(('Cayley Z_3xZ_3', 9, L.cayley_z3z3_strings()))
    for name, n, strs in fams:
        out = L.run_c(['stats', str(n), 'noaut'], feed=['printf', '%s\\n'] + strs)
        d = parse_summary(out)
        check(d['all_arcs_even'][0] == len(strs), 'Cayley all-even: %s' % name)
        print('%s: %d Cayley tournaments, all have every arc on an even number of HPs' % (name, len(strs)))
    # Paley minus a vertex: all-odd for p = 7, 11, 19
    for p in (7, 11, 19):
        n = p - 1
        d = parse_summary(L.run_c(['stats', str(n), 'noaut'], feed=['printf', '%s\\n', L.paley_string(p, 0)]))
        check(d['all_arcs_odd'][0] == 1, 'QR%d minus a vertex is all-odd' % p)
        print('QR_%d minus a vertex (N=%d): every one of the %d arcs lies on an odd number of HPs' % (p, n, n * (n - 1) // 2))
    for q in (7, 11, 19, 23, 27):
        t = time.time()
        line = L.run_c([str(q)], binary=L.BIN_PALEY).strip()
        kv = dict(tok.split('=') for tok in line.split()[1:])
        check(kv['all_odd'] == '1' and kv['H_parity'] == '1', 'mod-2 DP: Paley q=%d minus a vertex all-odd' % q)
        print('   mod-2 DP: %s  [%.1fs]' % (line, time.time() - t), flush=True)
    out = L.run_c(['stats', '10', 'noaut'], feed=['printf', '%s\\n', L.paley_string(11, 0)])
    hline = [ln for ln in out.splitlines() if ln.startswith('ALLODD')][0]
    print('   QR11 - v:', hline[:60], '... H =', hline.split('H=')[1].split()[0])
    # controls
    T15 = L.drt15_doubled()
    strs = [L.string_from_adj(L.sub_adj(T15, [v for v in range(15) if v != z])) for z in range(15)]
    d = parse_summary(L.run_c(['stats', '14', 'noaut'], feed=['printf', '%s\\n'] + strs))
    check(d['all_arcs_odd'][0] == 0, 'DRT15 (doubling) minus a vertex is never all-odd')
    d15 = parse_summary(L.run_c(['stats', '15', 'noaut'], feed=['printf', '%s\\n', L.string_from_adj(T15)]))
    print('control: doubled DRT(15): all-even=%d; its 15 vertex-deletions: all-odd=%d, #odd arcs %s of 91'
          % (d15['all_arcs_even'][0], d['all_arcs_odd'][0], sorted(d['odd-arc histogram'])))
    # N = 14 (= 2 mod 4, no Paley tournament of order 15): circulants on Z_15 minus a vertex
    strs = []
    for choice in itertools.product([0, 1], repeat=7):
        S = {r if c == 0 else 15 - r for r, c in zip(range(1, 8), choice)}
        V = list(range(1, 15))
        strs.append(''.join('1' if ((V[b] - V[a]) % 15) in S else '0' for a in range(14) for b in range(a + 1, 14)))
    d = parse_summary(L.run_c(['stats', '14', 'noaut'], feed=['printf', '%s\\n'] + strs))
    print('N=14: none of the 128 circulants on Z_15 minus a vertex is all-odd (all-odd count %d; max #odd arcs %d of 91)'
          % (d['all_arcs_odd'][0], max(d['odd-arc histogram'])))
    # n=10 extra all-odd example from the random-restart search
    ex = '110101111011111101011100000001101000100000001'
    d = parse_summary(L.run_c(['stats', '10', 'noaut'], feed=['printf', '%s\\n', ex]))
    check(d['all_arcs_odd'][0] == 1, 'search example at n=10 is all-odd')
    print('n=10 search example %s is all-odd (H=3929): all-odd classes at n=10 are not unique' % ex)
    # n=10 uniform random sample
    d = parse_summary(L.run_c(['rnd', '10', '20000', '20261001']))
    hist = d['odd-arc histogram']
    print('n=10 random sample (20000 labeled): all arcs on some HP %.4f, all-odd %d, min/max #odd arcs %d/%d, '
          'strong %.4f' % (d['all_arcs_on_some_HP'][1] / 20000, d['all_arcs_odd'][1], min(hist), max(hist), d['strong'][1] / 20000))
    # probes for conjecture C4 (even N: at least N-1 odd arcs) beyond exhaustive range
    for n, restarts in ((12, 40),):
        t = time.time()
        line = [ln for ln in L.run_c(['minodd', str(n), str(restarts), '20261001']).splitlines() if ln.startswith('MINODD')][0]
        best = int(line.split('min_odd_found=')[1].split()[0])
        print('C4 probe n=%d (local search, %d near-transitive restarts): fewest odd arcs found = %d (bound n-1 = %d)  [%.1fs]'
              % (n, restarts, best, n - 1, time.time() - t), flush=True)
    # n=10 exhaustive census (separate run, same engine)
    cens = os.path.join(L.REPO, '05-knowledge', 'results', 'procgen_selfie_20261001_n10_census.out')
    if os.path.exists(cens):
        d = parse_summary(open(cens).read())
        if 'labeled' in d:
            check(d['labeled'] == 2 ** 45, 'n=10 census labeled total')
            check(d['identity_failures'] == 0, 'n=10 census identities')
            check(d['all_arcs_even'][0] == 0, 'n=10 no all-even')
            hist = d['odd-arc histogram']
            print('Conjecture C4 at n=10: minimum #odd arcs = %d (%s), attained by %d classes'
                  % (min(hist), 'holds' if min(hist) >= 9 else 'FAILS: REFUTED', int(hist[min(hist)][0])))
            odds = [ln for ln in open(cens) if ln.startswith('ALLODD')]
            print('   all-odd classes at n=10: %s' % ['H=' + ln.split('H=')[1].split()[0] + ' |Aut|=%d'
                  % round(math.factorial(10) / float(ln.split('weight=')[1].split()[0])) for ln in odds])
            check(len(odds) == int(d['all_arcs_odd'][0]), 'all-odd lines = all-odd count at n=10')
            filt = os.path.join(L.REPO, '05-knowledge', 'results', 'procgen_selfie_20261001_n10_oddfilter.out')
            if os.path.exists(filt):
                txt = open(filt).read()
                line = [ln for ln in txt.splitlines() if ln.startswith('ODDFILTER')][0]
                hist2 = {int(k): int(v) for k, v in (tok.split(':') for tok in line.split('odd_hist:')[1].split())}
                check(hist2 == {k: int(v[0]) for k, v in hist.items()}, 'n=10 census histogram = independent mod-2 histogram')
                few = [ln.split()[2] for ln in txt.splitlines() if ln.startswith('FEWODD') and ln.endswith('odd=5')]
                print('   independent mod-2 census (%s) agrees; the %d classes with 5 odd arcs: %s'
                      % (os.path.basename(filt), len(few), [t[2:] for t in few]))
                for t in few:
                    A = L.adj_from_string(t[2:], 10)
                    H, c, st, en = L.arc_counts(A)
                    check(sum(1 for v in c.values() if v % 2) == 5, 'witness has 5 odd arcs')
                print('   (each witness re-verified with exact counts in Python)')
            print('n=10 EXHAUSTIVE (all %d classes; %s): all arcs on some HP %d classes / %.6f of labeled; '
                  'all-odd %d classes / %d labeled; min #odd arcs %d; odd-arc graph components %s'
                  % (d['classes'], os.path.basename(cens), d['all_arcs_on_some_HP'][0],
                     d['all_arcs_on_some_HP'][1] / 2 ** 45, d['all_arcs_odd'][0], d['all_arcs_odd'][1], min(hist),
                     {k: int(v[0]) for k, v in d['odd-arc-graph components'].items()}))
    return summaries


# --------------------------------------------------------------------------------------
def section_D():
    hdr('D. Shaved tournaments: HP-core, HP-blocking number beta, parity-break rho, Redei shaving')
    for n in range(3, 8):
        t = time.time()
        out = L.run_c(['shave', str(n)], feed=L.gentourng(n))
        rows = []
        for line in out.splitlines():
            if not line.startswith('SHAVE'):
                continue
            kv = dict(tok.split('=') for tok in line.split()[1:])
            rows.append(kv)
        check(len(rows) == len(L.classes(n)), 'shave rows n=%d' % n)
        ne = n * (n - 1) // 2
        check(all(r['beta'] == r['hall'] for r in rows), 'beta(T) = Hall-deficiency bound for all classes n=%d' % n)
        check(all(int(r['hall']) <= int(r['twosrc']) for r in rows), 'Hall bound <= two-source bound n=%d' % n)
        below = [r['T'] for r in rows if int(r['beta']) < int(r['twosrc'])]
        stuck = [r for r in rows if r['redei_shavable'] == '0']
        allodd = [r for r in rows if int(r['oddarcs']) == ne]
        check(sorted(r['T'] for r in stuck) == sorted(r['T'] for r in allodd), 'Redei-stuck iff all-odd n=%d' % n)
        rho2 = [r for r in rows if r['rho'] != '1']
        check(all(int(r['oddarcs']) == 0 for r in rho2) and all(r['rho'] == '2' for r in rho2),
              'rho = 1 unless all-even, then rho = 2 (n=%d)' % n)
        core = Counter(int(r['core']) for r in rows)
        betas = Counter(int(r['beta']) for r in rows)
        print('n=%d: beta distribution %s; beta = Hall bound in every class; beta < two-source bound in %d class(es) %s; '
              'Redei-shavable %d/%d (stuck = the all-odd classes: %d); rho distribution %s; HP-core size range [%d, %d]  [%.1fs]'
              % (n, dict(sorted(betas.items())), len(below), below, len(rows) - len(stuck), len(rows), len(allodd),
                 dict(Counter(r['rho'] for r in rows)), min(core), max(core), time.time() - t), flush=True)
        if n % 2 == 1:
            check(max(betas) == n - 1, 'max beta = n-1 for odd n=%d' % n)
        else:
            check(max(betas) == n - 2, 'max beta = n-2 for even n=%d' % n)
    for n in (8, 9):
        t = time.time()
        out = L.run_c(['beta', str(n)], feed=L.gentourng(n))
        line = [ln for ln in out.splitlines() if ln.startswith('BETA N=')][0]
        kv = dict(tok.split('=') for tok in line.split(' bound_hist')[0].split()[1:])
        check(int(kv['beta_below_hall_bound']) == 0, 'branch and bound: beta = Hall bound for all classes n=%d' % n)
        print('n=%d (branch and bound): %s  [%.1fs]' % (n, line, time.time() - t))


# --------------------------------------------------------------------------------------
def a049313_burnside(n):
    """THM-479 / Babai-Cameron Thm 7.2 closed form."""
    from fractions import Fraction
    total = Fraction(0)

    def partitions(m, maxp):
        if m == 0:
            yield []
            return
        for k in range(min(m, maxp), 0, -1):
            for rest in partitions(m - k, k):
                yield [k] + rest

    def v2(x):
        c = 0
        while x % 2 == 0:
            x //= 2
            c += 1
        return c

    for mu in partitions(n, n):
        cnt = Counter(mu)
        z = 1
        for part, mult in cnt.items():
            z *= part ** mult * math.factorial(mult)
        o2 = sum(l // 2 for l in mu) + sum(math.gcd(mu[i], mu[j]) for i in range(len(mu)) for j in range(i + 1, len(mu)))
        k = len(mu)
        if all(l % 2 for l in mu):
            total += Fraction(2 ** (o2 - k + 1), z)
        elif all(l % 2 == 0 for l in mu) and len({v2(l) for l in mu}) == 1:
            total += Fraction(2 ** (o2 - k), z)
    check(total.denominator == 1, 'A049313 Burnside integral n=%d' % n)
    return int(total)


def section_E():
    hdr('E. OPEN-Q-060: A049313 counts the orientation-parity-preserving Euler graphs')
    a = {n: a049313_burnside(n) for n in range(1, 13)}
    print('A049313 from the Babai-Cameron/THM-479 Burnside formula, n=1..12:', [a[n] for n in range(1, 13)])
    check([a[n] for n in range(1, 10)] == [1, 1, 1, 2, 2, 6, 12, 79, 792], 'A049313 known terms')
    # (1) direct count of untwisted Euler graphs
    for n in range(1, 11):
        out = L.run_c([str(n)], feed=['geng', '-q', str(n)], binary=L.BIN_EULER)
        line = [ln for ln in out.splitlines() if ln.startswith('EULER')][0]
        kv = dict(tok.split('=') for tok in line.split()[1:])
        check(int(kv['cycle_form_failures']) == 0, 'cycle form of the twist character n=%d' % n)
        check(int(kv['untwisted']) == a[n], 'untwisted Euler graphs = A049313(%d)' % n)
        print('n=%d: Euler graphs %s, orientation-parity-preserving %s = A049313(%d) = %d' % (n, kv['euler_graphs'], kv['untwisted'], n, a[n]))
    # (2) per-permutation twisted Brauer identity: all permutations n <= 6, one per cycle type n = 7, 8.
    #     Euler graphs are generated from the triangle basis {0,i,j} of the cycle space (2^C(n-1,2) of them).
    for n in range(2, 9):
        P, bit, tiles, T_of_t, cut = gauge_data(n)
        ne = len(P)
        euler = np.zeros(1, dtype=np.int64)
        for i in range(1, n):
            for j in range(i + 1, n):
                b = (1 << bit[(0, i)]) | (1 << bit[(0, j)]) | (1 << bit[(i, j)])
                euler = np.concatenate([euler, euler ^ b])
        check(euler.size == 2 ** ((n - 1) * (n - 2) // 2) and np.unique(euler).size == euler.size,
              'cycle space size n=%d' % n)
        seen_types = set()
        nperm = 0
        for sigma in itertools.permutations(range(n)):
            ctype = tuple(sorted(cycle_lengths(sigma)))
            if n >= 7 and ctype in seen_types:
                continue
            seen_types.add(ctype)
            nperm += 1
            img_bit = np.zeros(ne, dtype=np.int64)
            flip = np.zeros(ne, dtype=np.int64)
            for b, (i, j) in enumerate(P):
                a1, b1 = sigma[i], sigma[j]
                img_bit[b] = bit[(min(a1, b1), max(a1, b1))]
                flip[b] = 1 if a1 > b1 else 0

            def act(x, twist):
                y = np.zeros_like(x)
                for b in range(ne):
                    vb = (x >> b) & 1
                    if twist and flip[b]:
                        vb ^= 1
                    y |= vb << int(img_bit[b])
                return y

            def rep(x):
                Lv = np.zeros_like(x)
                lcur = np.zeros_like(x)
                for k in range(n - 1, 0, -1):
                    lnext = lcur ^ ((x >> bit[(k - 1, k)]) & 1)
                    Lv |= lnext << (k - 1)
                    lcur = lnext
                return x ^ cut[Lv]

            fixed_classes = int(np.sum(rep(act(T_of_t, True)) == T_of_t))
            fixedF = euler[act(euler, False) == euler]
            inv_mask = 0
            for b in range(ne):
                if flip[b]:
                    inv_mask |= 1 << b
            par = fixedF & inv_mask
            for sh in (32, 16, 8, 4, 2, 1):
                par = par ^ (par >> sh)
            rhs = int(np.sum(1 - 2 * (par & 1)))
            check(fixed_classes == rhs, 'twisted Brauer identity n=%d sigma=%s' % (n, sigma))
        print('n=%d: #switching classes fixed by s = sum over s-invariant Euler graphs F of eps_F(s), checked for %d %s'
              % (n, nperm, 'cycle types' if n >= 7 else 'permutations (all)'), flush=True)
    print('THEOREM (twisted Mallows-Sloane): A049313(n) = number of unlabeled Euler graphs on n vertices whose automorphisms '
          'all reverse an even number of edges of a fixed orientation.')


def cycle_lengths(sigma):
    n = len(sigma)
    seen = [False] * n
    out = []
    for s in range(n):
        if not seen[s]:
            l = 0
            x = s
            while not seen[x]:
                seen[x] = True
                x = sigma[x]
                l += 1
            out.append(l)
    return out


if __name__ == '__main__':
    L.build()
    section_A()
    section_B()
    section_C()
    section_D()
    section_E()
    print()
    print('total checks: %d, elapsed %.1f s' % (NCHECK, time.time() - T0))
    print('ALL CHECKS PASSED')
