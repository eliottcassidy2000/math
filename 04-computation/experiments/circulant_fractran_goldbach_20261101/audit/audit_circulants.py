#!/usr/bin/env python3
"""Independent audit of THM-4855 (circulant tournaments as S-sequencings).

Written from scratch for the audit; imports nothing from the session's scripts.
Every palindromicity verdict below uses the LITERAL arc test: the positional reversal
rho(v_i) = v_{q+1-i} of a Hamiltonian path is an anti-automorphism iff for all u != v,
u -> v  <=>  rho(v) -> rho(u).  The session's segment-sum criterion is tested AGAINST this.

usage (from the repo root):
  python3 04-computation/experiments/circulant_fractran_goldbach_20261101/audit/audit_circulants.py [QMAX_WITNESS] [QMAX_B1]
"""
import itertools
import random
import sys
import time

FAILS = []
MISMATCH = []


def report(name, ok, detail=''):
    print(f'  [{"PASS" if ok else "FAIL"}] {name}' + (f'  -- {detail}' if detail else ''))
    if not ok:
        FAILS.append(name)


def claim(name, ok, detail=''):
    """a number or statement printed by the session, compared with the audit's own value."""
    print(f'  [{"CLAIM OK" if ok else "CLAIM WRONG"}] {name}' + (f'  -- {detail}' if detail else ''))
    if not ok:
        MISMATCH.append(name)


# ----------------------------------------------------------------------------- basics
def connection_sets(q):
    """all S with 1 in S, exactly one of {x, -x} for each x != 0 (q odd)."""
    k = (q - 1) // 2
    for bits in itertools.product((0, 1), repeat=k - 1):
        S = {1}
        for x, b in zip(range(2, k + 1), bits):
            S.add(x if b == 0 else q - x)
        yield frozenset(S)


def member_table(q, S):
    return [(d % q) in S for d in range(q)]


def anti_aut_violation(q, inS, path):
    """literal arc test; returns None if rho is an anti-automorphism, else a violated ordered pair (u, v)."""
    pos = [0] * q
    for i, v in enumerate(path):
        pos[v] = i
    rho = [path[q - 1 - pos[v]] for v in range(q)]
    for u in range(q):
        ru = rho[u]
        for v in range(u + 1, q):
            # u -> v  iff  v - u in S ;  required: rho(v) -> rho(u) iff  rho(u) - rho(v) in S
            if inS[(v - u) % q] != inS[(ru - rho[v]) % q]:
                return (u, v) if inS[(v - u) % q] else (v, u)
    return None


def path_from_word(q, w):
    p = [0]
    for d in w:
        p.append((p[-1] + d) % q)
    return p


def is_ham_path(q, inS, path):
    return len(path) == q and len(set(path)) == q and all(inS[(path[i + 1] - path[i]) % q] for i in range(q - 1))


def segment_criterion(q, inS, w):
    """the session's criterion (1): all segment sums of w and rev(w) lie in S together."""
    r = w[::-1]
    n = len(w)
    for i in range(n):
        a = b = 0
        for j in range(i, n):
            a += w[j]
            b += r[j]
            if inS[a % q] != inS[b % q]:
                return False
    return True


def sequencings(q, S):
    """all words over S of length q-1 with distinct partial sums mod q (DFS)."""
    letters = sorted(S)
    out = []
    seen = [False] * q
    seen[0] = True
    word = []

    def rec(cur):
        if len(word) == q - 1:
            out.append(tuple(word))
            return
        for d in letters:
            nx = (cur + d) % q
            if not seen[nx]:
                seen[nx] = True
                word.append(d)
                rec(nx)
                word.pop()
                seen[nx] = False
    rec(0)
    return out


def multiplier_class_reps(q):
    """one S per orbit of the unit group acting by S -> aS, renormalised by x -> -x so that 1 in S."""
    from math import gcd
    units = [a for a in range(1, q) if gcd(a, q) == 1]
    seen = set()
    reps = []
    for S in connection_sets(q):
        if S in seen:
            continue
        orbit = set()
        for a in units:
            T = frozenset((a * x) % q for x in S)
            if 1 not in T:
                T = frozenset((-x) % q for x in T)
            orbit.add(T)
        seen |= orbit
        reps.append((min(orbit, key=lambda T: sorted(T)), len(orbit)))
    return reps


# ----------------------------------------------------------------------------- part 1: criterion (1) and claim (2)
def part1():
    print('PART 1: criterion (1) vs literal arc test; prime word test (2); counts')
    t0 = time.time()
    total = 0
    disagree = 0
    for q in (3, 5, 7, 9, 11):
        prime = all(q % p for p in range(2, q))
        for S in connection_sets(q):
            inS = member_table(q, S)
            words = sequencings(q, S)
            npal_word = npal_path = 0
            for w in words:
                path = path_from_word(q, w)
                lit = anti_aut_violation(q, inS, path) is None
                crit = segment_criterion(q, inS, w)
                total += 1
                if lit != crit:
                    disagree += 1
                if w == w[::-1]:
                    npal_word += 1
                    if not lit:
                        disagree += 1      # a palindromic word must give a palindromic path
                if lit:
                    npal_path += 1
            if prime:
                report(f'(2) q={q} S={sorted(S)}: palindromic paths = palindromic words',
                       npal_word == npal_path, f'{len(words)} words, {npal_word} palindromic words, '
                                                f'{npal_path} palindromic paths, H = {q * len(words)}')
            else:
                print(f'     q={q} S={sorted(S)}: {len(words)} words, {npal_word} palindromic words, '
                      f'{npal_path} palindromic paths (from 0), H = {q * len(words)}')
    # q = 13: one connection set per multiplier class (the criterion is multiplier-equivariant)
    q = 13
    words13 = 0
    for S, size in multiplier_class_reps(q):
        inS = member_table(q, S)
        words = sequencings(q, S)
        words13 += len(words)
        npal_word = npal_path = 0
        for w in words:
            path = path_from_word(q, w)
            lit = anti_aut_violation(q, inS, path) is None
            crit = segment_criterion(q, inS, w)
            total += 1
            if lit != crit:
                disagree += 1
            if w == w[::-1]:
                npal_word += 1
            if lit:
                npal_path += 1
        report(f'(2) q=13 S={sorted(S)} (class size {size}): palindromic paths = palindromic words',
               npal_word == npal_path, f'{len(words)} words, {npal_word} palindromic')
    report('criterion (1) = literal arc test on every enumerated word (q<=11 all S; q=13 class reps)',
           disagree == 0, f'{total} words checked, {disagree} disagreements, {time.time() - t0:.0f}s')
    # the session's count 1,741,634 = words from 0 for one S per multiplier class, q = 3..13
    tot_reps = words13
    for q in (3, 5, 7, 9, 11):
        for S, size in multiplier_class_reps(q):
            tot_reps += len(sequencings(q, S))
    claim('"1,741,634 Hamiltonian paths" = paths from vertex 0 of one S per multiplier class, q = 3..13',
          tot_reps == 1741634, f'audit count {tot_reps}; the number of ALL Hamiltonian paths of all these circulants '
                               f'is far larger (q = 13 alone: 32 sets x ~3.7e6)')
    # Z_9 hostile
    q, S = 9, frozenset({1, 3, 4, 7})
    inS = member_table(q, S)
    words = sequencings(q, S)
    bad = [w for w in words if w != w[::-1] and anti_aut_violation(q, inS, path_from_word(q, w)) is None]
    w0 = (1, 1, 1, 4, 7, 1, 7, 4)
    report('Z_9 hostile {1,3,4,7}: 66 palindromic paths (from 0) with non-palindromic words; '
           '(1,1,1,4,7,1,7,4) among them', len(bad) == 66 and w0 in bad, f'{len(bad)} found')


# ----------------------------------------------------------------------------- part 2: q = 3, 5 and THM-4849 data
def terras_cycle(word, c):
    """rational cycle of x -> x/2 (even), (3x+c)/2 (odd) with the given q-periodic parity word, as integers."""
    from fractions import Fraction
    # solve x = T^q(x) with parity word fixed: x_{j+1} = (3^{w_j} x_j + w_j c) / 2
    a, b = Fraction(1), Fraction(0)       # x_j = a x_0 + b
    for wj in word:
        a, b = a * (3 if wj else 1) / 2, (b * (3 if wj else 1) + (c if wj else 0)) / 2
    x0 = b / (1 - a)
    pts = [x0]
    for wj in word[:-1]:
        x = pts[-1]
        pts.append((3 * x + c) / 2 if wj else x / 2)
    # check parities agree with the word and values are integers
    ok = all(p.denominator == 1 and (p.numerator % 2 == wj) for p, wj in zip(pts, word))
    return [int(p) for p in pts], ok


def part2():
    print('PART 2: q = 3, 5 enumeration; THM-4849 cycle data')
    w3 = sequencings(3, frozenset({1}))
    report('q=3: only S-sequencing is (1,1)', w3 == [(1, 1)], str(w3))
    w5 = sorted(sequencings(5, frozenset({1, 2})))
    report('q=5, S={1,2}: S-sequencings exactly 1111, 2112, 2222, all palindromes, H = 15',
           w5 == [(1, 1, 1, 1), (2, 1, 1, 2), (2, 2, 2, 2)] and all(w == w[::-1] for w in w5), str(w5))
    w5b = sorted(sequencings(5, frozenset({1, 3})))
    print(f'     q=5, S={{1,3}}: {w5b}')
    inS = member_table(5, frozenset({1, 2}))
    allpaths = [[(p + t) % 5 for p in path_from_word(5, w)] for w in w5 for t in range(5)]
    report('T5: all 15 Hamiltonian paths palindromic by the literal test',
           len(allpaths) == 15 and all(anti_aut_violation(5, inS, p) is None for p in allpaths))
    # |Aut(T5)| by brute force over S_5
    aut = [p for p in itertools.permutations(range(5))
           if all(inS[(p[v] - p[u]) % 5] == inS[(v - u) % 5] for u in range(5) for v in range(5) if u != v)]
    report('|Aut(T5)| = 5, so tilings = H/|Aut| = 3 = number of words from 0', len(aut) == 5, f'|Aut| = {len(aut)}')
    # THM-4849: parity word u_j = [j in S], c = 2^q - 3^|S|
    for q, S, claimed in ((3, {1}, {1, 2, 4}), (5, {1, 2}, {5, 19, 40, 20, 10}), (5, {1, 3}, {7, 11, 22, 28, 14})):
        word = [1 if j in S else 0 for j in range(q)]
        c = 2 ** q - 3 ** len(S)
        pts, ok = terras_cycle(word, c)
        report(f'THM-4849 cycle q={q} S={sorted(S)}: 3x+{c} cycle = {sorted(claimed)}',
               ok and set(pts) == claimed, f'cycle from base point: {pts}')
    # T5 word 2112 visits 10, 19, 40, 20, 5 when position j = T^j(10)
    word = [1 if j in {1, 2} else 0 for j in range(5)]
    pts, _ = terras_cycle(word, 23)
    for w, claimed in (((1, 1, 1, 1), [10, 5, 19, 40, 20]), ((2, 1, 1, 2), [10, 19, 40, 20, 5]),
                       ((2, 2, 2, 2), [10, 19, 20, 5, 40])):
        visit = [pts[s] for s in path_from_word(5, w)]
        report(f'word {w} visits {claimed}', visit == claimed, str(visit))


# ----------------------------------------------------------------------------- part 3: witness families (3)
def case_of(q, S):
    """THM-4855 (3) case rule; returns list of (case, word, named segment (a, b) 1-based inclusive, claimed w-in-S,
    claimed rev-in-S). For A1 every admissible s is returned, for B every admissible rise u."""
    k = (q - 1) // 2
    out = []
    if 2 in S:
        ss = [s for s in S if 4 <= s <= k + 1 and (s - 1) not in S]
        for s in ss:
            w = (1,) * (s - 2) + (2,) + (1,) * (q - 1 - s) + (s,)
            out.append(('A1', w, (s, q - 1), False, True))
        if not ss:
            if S != frozenset(range(1, k + 1)):
                return [('A2-DERIVATION-FAILS', None, None, None, None)]
            w = (1, 2) + (1,) * (q - 4) + (3,)
            out.append(('A2', w, (1, k - 1), True, False))
    else:
        rises = [u for u in range(1, q - 1) if u not in S and (u + 1) in S]
        if not rises:
            return [('NO-RISE', None, None, None, None)]
        us = [u for u in rises if u != k]
        for u in us:
            w = (1,) * (u - 1) + (q - u, u + 1) + (1,) * (q - 2 - u)
            out.append(('B', w, (u, u + 1), True, False))
        if not us:
            if S != frozenset([1] + list(range(k + 1, q - 1))):
                return [('D-DERIVATION-FAILS', None, None, None, None)]
            w = (1,) + (q - 3, 1) * (k - 1) + (q - 2,)
            out.append(('D', w, (1, 2), True, False))
    return out


def part3(qmax):
    print(f'PART 3: explicit witnesses A1/A2/B/D for every S with 1 in S, odd 7 <= q <= {qmax} (literal arc test)')
    print('     every admissible s (A1) and every admissible rise u (B) is tested, not one choice per S;')
    print('     each witness: Hamiltonian path + named segment = literally violated arc (v_a, v_{b+1});')
    print('     one witness per S additionally gets the full O(q^2) anti-automorphism scan.')
    grand = 0
    t0 = time.time()
    for q in range(7, qmax + 1, 2):
        nsets = 0
        cnt = {'A1': 0, 'A2': 0, 'B': 0, 'D': 0}
        nwords = 0
        bad = []
        for S in connection_sets(q):
            nsets += 1
            inS = member_table(q, S)
            cases = case_of(q, S)
            kinds = set()
            full_done = False
            for kind, w, seg, cw, cr in cases:
                if w is None:
                    bad.append((sorted(S), kind))
                    continue
                kinds.add(kind)
                nwords += 1
                path = path_from_word(q, w)
                if not is_ham_path(q, inS, path):
                    bad.append((sorted(S), kind, w, 'not a Hamiltonian path'))
                    continue
                if not full_done:
                    full_done = True
                    if anti_aut_violation(q, inS, path) is None:
                        bad.append((sorted(S), kind, w, 'palindromic by the full literal scan'))
                        continue
                a, b = seg
                sw = sum(w[a - 1:b]) % q
                sr = sum(w[::-1][a - 1:b]) % q
                if (sw in S) != cw or (sr in S) != cr:
                    bad.append((sorted(S), kind, w, f'named segment sums {sw},{sr} not as claimed'))
                    continue
                pos = [0] * q
                for i, v in enumerate(path):
                    pos[v] = i
                va, vb = path[a - 1], path[b]                     # v_a, v_{b+1}
                rva, rvb = path[q - 1 - pos[va]], path[q - 1 - pos[vb]]
                if inS[(vb - va) % q] == inS[(rva - rvb) % q]:
                    bad.append((sorted(S), kind, w, 'named segment is not a violated arc'))
                if kind == 'A1' and not (sw == q - 1 and sr == (1 - w[-1]) % q):
                    bad.append((sorted(S), kind, w, 'A1 values'))
                if kind == 'A2' and not (sw == (q - 1) // 2 and sr == (q + 1) // 2):
                    bad.append((sorted(S), kind, w, 'A2 values'))
                if kind == 'B' and not (sw == 1 and sr == 2):
                    bad.append((sorted(S), kind, w, 'B values'))
                if kind == 'D' and not (sw == q - 2 and sr == q - 1):
                    bad.append((sorted(S), kind, w, 'D values'))
            if len(kinds) != 1:
                bad.append((sorted(S), 'case split not a partition', kinds))
            for kd in kinds:
                cnt[kd] += 1
        grand += nsets
        report(f'q={q}: {nsets} sets, every listed witness valid and failing at its named segment',
               not bad, f'cases {cnt}, {nwords} witness words checked' + (f'; first problems {bad[:3]}' if bad else ''))
    print(f'     total connection sets with 1 in S over odd q = 7..{qmax}: {grand}  ({time.time() - t0:.0f}s)')
    if qmax == 35:
        claim('"131,071 connection sets ... at every odd q from 7 to 35"', grand == 131071,
              f'actual total {grand} = 2^17 - 4 (131,071 = 2^17 - 1 would also count q = 3, 5)')


# ----------------------------------------------------------------------------- part 4: THM-4730 (B1) in general
def three_interval_path(q, a, b):
    x, y = a - 1, b - 1
    return list(range(0, x)) + list(range(x + y, q)) + list(range(x, x + y))


def part4(qmax_b1):
    print(f'PART 4: THM-4730 (B1) for ALL admissible (a, b) (not only b = 2), composite odd q <= {qmax_b1}')
    for q in range(9, qmax_b1 + 1, 2):
        if all(q % p for p in range(2, q)):
            continue   # at prime q (2) + non-palindromic word settles (B1)
        tested = pal = 0
        examples = []
        for S in connection_sets(q):
            inS = member_table(q, S)
            rest = sorted(S - {1})
            for a in rest:
                for b in rest:
                    if a == b or a + b > q + 1:
                        continue
                    path = three_interval_path(q, a, b)
                    assert is_ham_path(q, inS, path), (q, S, a, b)
                    tested += 1
                    if anti_aut_violation(q, inS, path) is None:
                        pal += 1
                        if len(examples) < 3:
                            examples.append((sorted(S), a, b))
        report(f'(B1) q={q}: every three-interval path P(a-1,b-1) non-palindromic', pal == 0,
               f'{tested} paths tested, {pal} palindromic' + (f', e.g. {examples}' if examples else ''))


# ----------------------------------------------------------------------------- part 5: (2) sampled at larger primes
def random_sequencing(q, S, rng, tries=2000):
    letters = sorted(S)
    for _ in range(tries):
        seen = {0}
        cur = 0
        w = []
        while len(w) < q - 1:
            opts = [d for d in letters if (cur + d) % q not in seen]
            if not opts:
                break
            d = rng.choice(opts)
            cur = (cur + d) % q
            seen.add(cur)
            w.append(d)
        if len(w) == q - 1:
            return tuple(w)
    return None


def descending_pairs_path(q):
    p = [0, 1]
    for j in range(1, (q - 3) // 2 + 1):
        p += [q - 2 * j, q - 2 * j + 1]
    return p + [2]


def part4b(qmax):
    print(f'PART 4b: are the THM-4855 witnesses members of THM-4730 (B)\'s own families? (odd 7 <= q <= {qmax})')
    ok = {'A1': True, 'A2': True, 'B': True, 'D': True}
    n = {'A1': 0, 'A2': 0, 'B': 0, 'D': 0}
    for q in range(7, qmax + 1, 2):
        for S in connection_sets(q):
            for kind, w, seg, cw, cr in case_of(q, S):
                path = path_from_word(q, w)
                n[kind] += 1
                if kind == 'A1':
                    s = w[-1]
                    ok['A1'] &= path == three_interval_path(q, s, 2)          # P(s-1, 1): (a, b) = (s, 2)
                elif kind == 'A2':
                    ok['A2'] &= path == three_interval_path(q, 3, 2)          # P(2, 1): (a, b) = (3, 2)
                elif kind == 'B':
                    u = seg[0]
                    ok['B'] &= path == three_interval_path(q, u + 1, q - u)   # P(u, q-1-u): a + b = q + 1
                else:
                    ok['D'] &= path == descending_pairs_path(q)
    report('A1 = three-interval P(s-1, 1), (a, b) = (s, 2)', ok['A1'], f'{n["A1"]} witnesses')
    report('A2 = three-interval P(2, 1), (a, b) = (3, 2)', ok['A2'], f'{n["A2"]} witnesses')
    report('B = three-interval P(u, q-1-u), (a, b) = (u+1, q-u), a + b = q + 1 (THM-4730 (B)\'s instance for one '
           'small element a = u+1, d = a-1): family B is NOT new', ok['B'], f'{n["B"]} witnesses')
    report('D = THM-4730 descending-pairs path', ok['D'], f'{n["D"]} witnesses')


def part5():
    print('PART 5: (2) sampled: random S-sequencings at q = 17, 19, 23 (non-palindromic word => non-palindromic path)')
    rng = random.Random(4855)
    for q in (17, 19, 23):
        sets = list(connection_sets(q))
        n = bad = 0
        for _ in range(3000):
            S = rng.choice(sets)
            w = random_sequencing(q, S, rng)
            if w is None or w == w[::-1]:
                continue
            n += 1
            if anti_aut_violation(q, member_table(q, S), path_from_word(q, w)) is None:
                bad += 1
        report(f'q={q}: sampled non-palindromic words give non-palindromic paths', bad == 0, f'{n} samples, {bad} bad')


def main():
    qmax_w = int(sys.argv[1]) if len(sys.argv) > 1 else 35
    qmax_b1 = int(sys.argv[2]) if len(sys.argv) > 2 else 27
    part1()
    part2()
    part3(qmax_w)
    part4(qmax_b1)
    part4b(qmax_w)
    part5()
    print()
    print('SESSION CLAIMS FOUND WRONG: ' + ('; '.join(MISMATCH) if MISMATCH else 'none'))
    print('AUDIT CHECK FAILURES: ' + ('; '.join(FAILS) if FAILS else 'none'))
    return 1 if FAILS else 0


if __name__ == '__main__':
    sys.exit(main())
