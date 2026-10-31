#!/usr/bin/env python3
"""mac-mini-2026-10-31-fivegroups: pre-registered tests of the owner's PRIMEGAME / metagraph correspondences.

PRE-REGISTRATION (written and committed before the first run of this file; any failure is reported as such).

  P1 (OEIS attribution). Counting the start 2 as step 0, the step at which 2^(p_k) first appears is
     - OEIS A007546 ("fast version") for Guy's list (..., 15/14, 15/2, 55/1), for all 34 listed terms;
     - OEIS A007547 ("slow version") for the owner's list (..., 15/2, 1/7, 55/1), for all 34 listed terms.
     (Exploratory before this file: both agreed for the first 10 terms, primegame_decode.out.)
  P2 (step-difference law). For every k checked,
         A007547(k) - A007546(k) = sum_{n=2}^{p_k - 1} (D(n) - 1),   D(n) = n / spf(n) (largest proper divisor),
     for k <= 34 against the OEIS data, and for all 168 primes p < 1000 with an accelerated simulator that is first
     validated configuration-by-configuration against the literal FRACTRAN simulator (first 10 primes).
  P3 (points of symmetry = short tournament cycles). Every circulant tournament of order q in {1, 3, 5} is pure-blue
     (all Hamiltonian paths palindromic, THM-4730 sense), no circulant tournament of order 7 is; under the THM-4849
     dictionary the q <= 5 circulants give exactly the Terras cycles {0}, (1,2,4)/5 and the two cycles with denominator
     23 whose integer forms (cycles of x -> x/2, (3x+23)/2) are {5,19,40,20,10} and {7,11,22,28,14}.
  P4 (label dependence of the 95/23 coincidence). PRIMEGAME contains 95/23 = (5*19)/23 and the numerator 77 = 7*11
     (exploratory observation before this file). Over all 4! * 6! = 17280 relabellings (registers among {2,3,5,7},
     states among {11,...,29}) of the owner's list, the event "some fraction is (odd-point product of a T5 cycle)/23"
     holds for fewer than 5% of relabellings, and "both 95 and 77 occur as numerators" for fewer than 10%.
     Outcome as predicted => NUMEROLOGY (not invariant under the program's relabelling symmetry).
  P5 (look-elsewhere calibration for F = U + S). The observation 2*c3 = F = (0, 2, 10) on (TT5, 1->C3->1, T5) was
     made before this file (exploratory). Prediction: over all target triples (0, a, b) with 0 <= a <= b <= 15, the
     fraction matched by at least one battery invariant under the allowed transformations (multiplier in
     {1/2, 1, 2, 3}, offset in {-1, 0, +1}) is at least 20%, so one match carries no evidential weight.
  P6 (cross-check of opus-2026-10-31-S70 THM-4852 (4)). The register-flow graph (den prime -> num prime) has 24 arcs,
     19 unordered pairs and no loops; the merged metagraph's simple non-loop support has 18 edges; hence no
     adjacency-preserving bijection. Also the co-occurrence graph (primes adjacent when they share a fraction) is not
     isomorphic to the support graph.
  P7 (exact step formula; derived by hand from the machine table, checked by hand only for N = 2..5 before this
     file). The phase that tests N = n + 1 (from the state-free configuration 2^n 7^(D(n)-1)) takes
         cost_G(N) = 7N + 1 + 2 spf(N) + (N - 1 - D(N)) (6N + 2) + 2 * sum_{D(N) < d < N} floor(N/d)
     steps in Guy's list and cost_O(N) = cost_G(N) + D(N-1) - 1 in the owner's list, so
         A007546(k) = sum_{N=2}^{p_k} cost_G(N),   A007547(k) = sum_{N=2}^{p_k} cost_O(N).
     Predicted to hold for all 34 OEIS terms of each and for all 168 primes < 1000 (jump simulator).

usage: python3 correspondence_tests.py   (needs merged_metagraph_n5.json from merged_metagraph_n5.py; networkx)
"""
import itertools
import json
import math
import os
import sys
from fractions import Fraction

import networkx as nx

HERE = os.path.dirname(os.path.abspath(__file__))
PR = [2, 3, 5, 7, 11, 13, 17, 19, 23, 29]
REG, STA = PR[:4], PR[4:]
HEAD = [(17, 91), (78, 85), (19, 51), (23, 38), (29, 33), (77, 29), (95, 23), (77, 19), (1, 17), (11, 13), (13, 11)]
GUY = HEAD + [(15, 14), (15, 2), (55, 1)]
OWNER = HEAD + [(15, 2), (1, 7), (55, 1)]

# OEIS data retrieved 2026-10-31 through the JSON API (https://oeis.org/search?fmt=json), all listed terms.
A007546 = [19, 69, 280, 707, 2363, 3876, 8068, 11319, 19201, 36866, 45551, 75224, 101112, 117831, 152025, 215384,
           293375, 327020, 428553, 507519, 555694, 700063, 808331, 989526, 1273490, 1434366, 1530213, 1710923,
           1818254, 2019962, 2833089, 3104685, 3546320, 3720785]
A007547 = [19, 69, 281, 710, 2375, 3893, 8102, 11361, 19268, 36981, 45680, 75417, 101354, 118093, 152344, 215797,
           293897, 327571, 429229, 508284, 556494, 701008, 809381, 990746, 1274952, 1435957, 1531854, 1712701,
           1820085, 2021938, 2835628, 3107393, 3549288, 3723821]

FAIL = []


def check(name, ok, detail=''):
    print(f'  [{"PASS" if ok else "FAIL"}] {name}' + (f'  {detail}' if detail else ''))
    if not ok:
        FAIL.append(name)


def section(t):
    print()
    print('=' * 100)
    print(t)
    print('=' * 100)


def fac(x):
    v = [0] * 10
    for i, p in enumerate(PR):
        while x % p == 0:
            x //= p
            v[i] += 1
    assert x == 1
    return v


def primes_upto(n):
    s = bytearray([1]) * (n + 1)
    s[0:2] = b'\x00\x00'
    for i in range(2, int(n ** 0.5) + 1):
        if s[i]:
            s[i * i::i] = bytearray(len(s[i * i::i]))
    return [i for i in range(n + 1) if s[i]]


def spf(n):
    for p in range(2, int(n ** 0.5) + 1):
        if n % p == 0:
            return p
    return n


def D(n):
    return 1 if n == 1 else n // spf(n)


# ---------------------------------------------------------------- literal FRACTRAN simulator
def literal_run(prog, nemit, record=False):
    rules = [(fac(a), fac(b)) for a, b in prog]
    v = fac(2)
    steps, out, trace = 0, [], []
    while len(out) < nemit:
        for num, den in rules:
            if all(v[i] >= den[i] for i in range(10)):
                v = [v[i] - den[i] + num[i] for i in range(10)]
                break
        steps += 1
        if record:
            trace.append(tuple(v))
        if v[0] > 0 and not any(v[1:]):
            out.append((v[0], steps))
    return out, trace


# ---------------------------------------------------------------- state-dispatch simulator with 2-cycle jumps
def machine(prog):
    """state -> ordered list of (dec vector over registers, inc vector, next state); valid by the one-hot lemma
    (primegame_decode.py section 2)."""
    tab = {}
    for a, b in prog:
        na, nb = fac(a), fac(b)
        sin = next((STA[i - 4] for i in range(4, 10) if nb[i]), 1)
        sout = next((STA[i - 4] for i in range(4, 10) if na[i]), 1)
        tab.setdefault(sin, []).append((tuple(nb[:4]), tuple(na[:4]), sout))
    return tab


def fast_run(prog, nemit, validate_trace=None):
    """Exact step counts. Executes single transitions, but jumps over runs of a self-loop or of a 2-cycle of
    transitions s -> t -> s whose guards stay true: the number of whole iterations is computed exactly."""
    tab = machine(prog)
    r = [1, 0, 0, 0]
    s = 1
    steps = 0
    out = []
    tpos = 0

    def pick(state, regs):
        for dec, inc, nxt in tab[state]:
            if all(regs[i] >= dec[i] for i in range(4)):
                return dec, inc, nxt
        raise RuntimeError('no rule')

    def full(regs, state):
        v = list(regs) + [0] * 6
        if state != 1:
            v[4 + STA.index(state)] = 1
        return tuple(v)

    while len(out) < nemit:
        dec, inc, nxt = pick(s, r)
        if validate_trace is None:
            # try to jump: self-loop (s == nxt) or 2-cycle s -> nxt -> s
            if nxt == s:
                delta = [inc[i] - dec[i] for i in range(4)]
                # repeat while guard holds: guard needs r[i] >= dec[i]; r changes by delta each time
                k = None
                for i in range(4):
                    if dec[i] and delta[i] < 0:
                        kk = (r[i] - dec[i]) // (-delta[i]) + 1
                        k = kk if k is None else min(k, kk)
                if k is not None and k > 1:
                    # but an earlier rule in the list must stay false: earlier rules of state s have stricter guards;
                    # check by recomputing pick at the last iterate is unnecessary only if no earlier rule can wake up.
                    # Safe fallback: verify earlier rules stay false along the run (registers move monotonically).
                    ok = True
                    for dec2, inc2, nxt2 in tab[s]:
                        if (dec2, inc2, nxt2) == (dec, inc, nxt):
                            break
                        # an earlier rule fails now; it can become true only if a register it needs increases
                        if any(dec2[i] and delta[i] > 0 for i in range(4)):
                            ok = False
                    if ok:
                        r = [r[i] + k * delta[i] for i in range(4)]
                        steps += k
                        if s == 1 and r[0] > 0 and not any(r[1:]):
                            pass  # emissions are checked below at single steps; a jump never ends at a pure power
                        continue
            else:
                dec2, inc2, nxt2 = pick(nxt, [r[i] - dec[i] + inc[i] for i in range(4)])
                if nxt2 == s and nxt != 1 and s != 1:
                    # 2-cycle s -> nxt -> s; iterate while the same two rules are chosen
                    d = [inc[i] - dec[i] + inc2[i] - dec2[i] for i in range(4)]
                    # guards: at s need r >= dec (and earlier rules false); at nxt need r' >= dec2 (earlier false)
                    # compute the max k such that all k iterations keep both guards; fall back to stepping if any
                    # register needed by an earlier rule of either state increases.
                    def earlier_ok(state, rule):
                        for rr in tab[state]:
                            if rr == rule:
                                return True
                            if any(rr[0][i] and d[i] > 0 for i in range(4)):
                                return False
                        return True
                    if earlier_ok(s, (dec, inc, nxt)) and earlier_ok(nxt, (dec2, inc2, nxt2)):
                        k = None
                        mid = [inc[i] - dec[i] for i in range(4)]
                        for i in range(4):
                            if d[i] < 0:
                                # guard at s: r_i + j d_i >= dec_i ; guard at nxt: r_i + j d_i + mid_i >= dec2_i
                                need = max(dec[i], dec2[i] - mid[i])
                                kk = (r[i] - need) // (-d[i]) + 1
                                k = kk if k is None else min(k, kk)
                        if k is not None and k > 1:
                            r = [r[i] + k * d[i] for i in range(4)]
                            steps += 2 * k
                            continue
        # single transition
        r = [r[i] - dec[i] + inc[i] for i in range(4)]
        s = nxt
        steps += 1
        if validate_trace is not None:
            if full(r, s) != validate_trace[tpos]:
                raise RuntimeError(f'validation mismatch at step {steps}')
            tpos += 1
        if s == 1 and r[0] > 0 and not any(r[1:]):
            out.append((r[0], steps))
    return out


def P1_P2():
    section('P1/P2: OEIS attribution of the two lists and the step-difference law')
    # validate the dispatch simulator (no jumps) against the literal simulator, configuration by configuration
    for name, prog in (('GUY', GUY), ('OWNER', OWNER)):
        lit, trace = literal_run(prog, 6, record=True)
        disp = fast_run(prog, 6, validate_trace=trace)
        check(f'{name}: state-dispatch simulator = literal FRACTRAN, every configuration, first 6 emissions',
              disp == lit, f'{len(trace)} configurations')
        jump = fast_run(prog, 10)
        lit10, _ = literal_run(prog, 10)
        check(f'{name}: jump simulator = literal FRACTRAN emissions and steps, first 10 primes', jump == lit10)
    g34 = fast_run(GUY, 34)
    o34 = fast_run(OWNER, 34)
    P = primes_upto(200)
    check('GUY list emits the first 34 primes in order', [p for p, _ in g34] == P[:34])
    check('OWNER list emits the first 34 primes in order', [p for p, _ in o34] == P[:34])
    check('P1: GUY steps = OEIS A007546 (34 terms)', [t for _, t in g34] == A007546)
    check('P1: OWNER steps = OEIS A007547 (34 terms)', [t for _, t in o34] == A007547)
    law = [sum(D(n) - 1 for n in range(2, p)) for p in P[:34]]
    diff = [b - a for a, b in zip(A007546, A007547)]
    check('P2: A007547 - A007546 = sum_{n=2}^{p_k-1} (D(n)-1), k <= 34 (OEIS data)', diff == law,
          f'first differences {diff[:10]}')
    P1000 = primes_upto(1000)
    g = fast_run(GUY, len(P1000))
    o = fast_run(OWNER, len(P1000))
    ok_emit = [p for p, _ in g] == P1000 and [p for p, _ in o] == P1000
    check('both lists emit the 168 primes < 1000 in order (jump simulator)', ok_emit)
    law1000 = [sum(D(n) - 1 for n in range(2, p)) for p in P1000]
    check('P2: step-difference law for all 168 primes < 1000', [b - a for (_, a), (_, b) in zip(g, o)] == law1000,
          f'k=168: steps GUY {g[-1][1]}, OWNER {o[-1][1]}, difference {o[-1][1] - g[-1][1]} = law {law1000[-1]}')
    # asymptotic reading: sum_{n<=x} D(n) ~ c x^2 / 2 with c = sum_p p^-2 prod_{q<p} (1-1/q)
    c = 0.0
    prod = 1.0
    for p in primes_upto(10 ** 6):
        c += prod / (p * p)
        prod *= 1 - 1 / p
    x = P1000[-1]
    print(f'  growth: law(997) = {law1000[-1]}; c x^2/2 with c = sum_p p^-2 prod_(q<p)(1-1/q) = {c:.6f}: '
          f'{c * x * x / 2:.0f} (ratio {law1000[-1] / (c * x * x / 2):.4f})')
    print(f'  share of the slow list spent in the extra 1/7 steps at k=168: {law1000[-1] / o[-1][1]:.2e}')
    # P7: exact closed forms
    def cost_G(N):
        d0 = D(N)
        return (7 * N + 1 + 2 * spf(N) + (N - 1 - d0) * (6 * N + 2)
                + 2 * sum(N // d for d in range(d0 + 1, N)))

    def cost_O(N):
        return cost_G(N) + D(N - 1) - 1

    def closed(costf, plist):
        out, acc, N = [], 0, 2
        for p in plist:
            while N <= p:
                acc += costf(N)
                N += 1
            out.append(acc)
        return out
    check('P7: A007546(k) = sum_{N=2}^{p_k} cost_G(N), 34 OEIS terms', closed(cost_G, P[:34]) == A007546)
    check('P7: A007547(k) = sum_{N=2}^{p_k} cost_O(N), 34 OEIS terms', closed(cost_O, P[:34]) == A007547)
    check('P7: closed forms = simulated steps for all 168 primes < 1000 (both lists)',
          closed(cost_G, P1000) == [t for _, t in g] and closed(cost_O, P1000) == [t for _, t in o])
    lead = [closed(cost_G, [p])[0] / (2 * (1 - c) * p ** 3) for p in (101, 499, 997)]
    print(f'  growth: A007546-type count / (2 (1 - c) p^3) at p = 101, 499, 997: {[round(x, 4) for x in lead]}')
    big = primes_upto(3000)
    vals = closed(cost_G, [big[-1]])
    print(f'  closed form at p = {big[-1]}: {vals[0]}; ratio to 2(1-c)p^3 = {vals[0] / (2 * (1 - c) * big[-1] ** 3):.5f}')


# ---------------------------------------------------------------- P3
def ham_paths(A):
    n = len(A)
    for p in itertools.permutations(range(n)):
        if all(A[p[i]][p[i + 1]] for i in range(n - 1)):
            yield p


def palindromic(A, path):
    n = len(path)
    # positional reversal v_i -> v_{n+1-i} must be an anti-automorphism: A[u][w] = A[rev u][rev w] reversed
    rev = {path[i]: path[n - 1 - i] for i in range(n)}
    return all(A[u][w] == A[rev[w]][rev[u]] for u in range(n) for w in range(n) if u != w)


def circulant(q, S):
    return [[1 if (j - i) % q in S else 0 for j in range(q)] for i in range(q)]


def terras_cycle(word):
    """rational periodic point of T(x) = x/2 (even), (3x+1)/2 (odd) with the given parity word; returns (x, D)"""
    L, k = len(word), sum(word)
    c = 0
    for j, b in enumerate(word):
        if b:
            c += 3 ** sum(word[j + 1:]) * 2 ** j
    Dq = 2 ** L - 3 ** k
    x = Fraction(c, Dq)
    orbit = []
    y = x
    for b in word:
        orbit.append(y)
        y = y / 2 if b == 0 else (3 * y + 1) / 2
    assert y == x
    return orbit, Dq


def P3():
    section('P3: the strong points of symmetry are the circulants of order <= 5, i.e. the short tournament cycles')
    for q in (1, 3, 5, 7):
        sets = [S for r in range(q) for S in itertools.combinations(range(1, q), r)
                if all(((-s) % q) not in S for s in S) and len(S) == (q - 1) // 2]
        allpb = []
        for S in sets:
            A = circulant(q, set(S))
            pb = all(palindromic(A, p) for p in ham_paths(A))
            allpb.append(pb)
        print(f'  q={q}: {len(sets)} connection sets; pure-blue: {allpb}')
        if q <= 5:
            check(f'P3: every circulant tournament of order {q} is pure-blue', all(allpb))
        else:
            check('P3: no circulant tournament of order 7 is pure-blue', not any(allpb))
        if q >= 3:
            necks = {}
            for S in sets:
                w = tuple(1 if i in S else 0 for i in range(q))
                key = min(w[i:] + w[:i] for i in range(q))
                necks.setdefault(key, []).append(S)
            for key, Ss in sorted(necks.items()):
                orbit, Dq = terras_cycle(list(key))
                ints = [int(y * Dq) for y in orbit]
                odd = sorted(v for v, b in zip(ints, key) if b)
                print(f'     necklace {key}: connection sets {Ss}; rational cycle with denominator D = 2^{q} - 3^'
                      f'{(q - 1) // 2} = {Dq}; integer cycle of x -> x/2, (3x+{Dq})/2: {ints}, odd points {odd}')
    o5 = []
    for key in ((0, 0, 0, 1, 1), (0, 0, 1, 0, 1)):
        orbit, Dq = terras_cycle(list(key))
        o5.append((sorted(int(y * Dq) for y in orbit), Dq))
    check('P3: the two q = 5 cycles are {5,10,19,20,40} and {7,11,14,22,28} with D = 23',
          sorted(o5) == sorted([([5, 10, 19, 20, 40], 23), ([7, 11, 14, 22, 28], 23)]), str(o5))


# ---------------------------------------------------------------- P4
def relabel(prog, mr, ms):
    m = dict(zip(REG, mr))
    m.update(dict(zip(STA, ms)))

    def img(x):
        y = 1
        for p in PR:
            while x % p == 0:
                x //= p
                y *= m[p]
        return y
    return [(img(a), img(b)) for a, b in prog]


def P4():
    section('P4: label dependence of the 95/23 and 77 coincidence')
    odd_products = {5 * 19, 7 * 11}
    n_tot = n_a = n_b = 0
    for mr in itertools.permutations(REG):
        for ms in itertools.permutations(STA):
            pr = relabel(OWNER, mr, ms)
            n_tot += 1
            if any(a in odd_products and b == 23 for a, b in pr):
                n_a += 1
            nums = {a for a, _ in pr}
            if 95 in nums and 77 in nums:
                n_b += 1
    fa, fb = n_a / n_tot, n_b / n_tot
    print(f'  relabellings: {n_tot}; "(odd-point product)/23 occurs": {n_a} ({fa:.4%}); "95 and 77 both numerators": '
          f'{n_b} ({fb:.4%})')
    check('P4: (odd-point product)/23 occurs for < 5% of relabellings', fa < 0.05)
    check('P4: 95 and 77 both numerators for < 10% of relabellings', fb < 0.10)


# ---------------------------------------------------------------- P5
def tournament_invariants(A):
    n = len(A)
    out = {}
    sc = sorted(sum(r) for r in A)
    paths = list(ham_paths(A))
    H = len(paths)
    aut = sum(1 for p in itertools.permutations(range(n))
              if all(A[p[i]][p[j]] == A[i][j] for i in range(n) for j in range(n)))
    def cyc(k):
        c = 0
        for S in itertools.combinations(range(n), k):
            for rest in itertools.permutations(S[1:]):
                cy = (S[0],) + rest
                if all(A[cy[i]][cy[(i + 1) % k]] for i in range(k)):
                    c += 1
        return c
    out['n'] = n
    out['arcs'] = n * (n - 1) // 2
    out['H'] = H
    out['|Aut|'] = aut
    out['tilings H/|Aut|'] = H // aut
    out['c3'] = cyc(3) if n >= 3 else 0
    out['c4'] = cyc(4) if n >= 4 else 0
    out['c5'] = cyc(5) if n >= 5 else 0
    out['all cycles'] = out['c3'] + out['c4'] + out['c5']
    out['transitive triples'] = math.comb(n, 3) - out['c3']
    out['strong components'] = None  # filled by caller
    out['score variance*n'] = sum((s - (n - 1) / 2) ** 2 for s in sc)
    out['max score'] = max(sc)
    out['min score'] = min(sc)
    out['sum of squares of scores'] = sum(s * s for s in sc)
    out['Ham cycles'] = out['c5'] if n == 5 else (out['c3'] if n == 3 else 0)
    return out


def P5():
    section('P5: look-elsewhere calibration for the three points of symmetry against F, S, U')
    TT5 = [[1 if i < j else 0 for j in range(5)] for i in range(5)]
    C31 = [[0] * 5 for _ in range(5)]
    for i in range(5):
        for j in range(5):
            if i != j:
                C31[i][j] = 1 if i < j else 0
    # 1 -> C3 -> 1 : vertex 0 beats all, vertex 4 loses to all, middle {1,2,3} a 3-cycle
    for (u, w) in ((1, 2), (2, 3), (3, 1)):
        C31[u][w], C31[w][u] = 1, 0
    T5 = circulant(5, {1, 2})
    pts = {'TT5': TT5, '1->C3->1': C31, 'T5': T5}
    inv = {k: tournament_invariants(A) for k, A in pts.items()}
    comps = {'TT5': 5, '1->C3->1': 3, 'T5': 1}
    for k in inv:
        inv[k]['strong components'] = comps[k]
    names = [k for k in inv['TT5'] if inv['TT5'][k] is not None]
    print('  battery on (TT5, 1->C3->1, T5):')
    for nm in names:
        print(f'    {nm:28} {[inv[k][nm] for k in pts]}')
    targets = {'F = d-2': (0, 2, 10), 'S (squarefree in D)': (0, 1, 7), 'U (prime in D)': (0, 1, 3),
               'tau': (2, 4, 12)}
    mults = (Fraction(1, 2), Fraction(1), Fraction(2), Fraction(3))
    offs = (-1, 0, 1)

    def matches(tgt):
        hits = []
        for nm in names:
            vals = [Fraction(inv[k][nm]) for k in pts]
            for m in mults:
                for o in offs:
                    if all(m * v + o == t for v, t in zip(vals, tgt)):
                        hits.append(f'{m}*{nm}{"+" if o >= 0 else ""}{o if o else ""}'.replace('+', '+') if o else
                                    f'{m}*{nm}')
        return hits
    for tn, tgt in targets.items():
        print(f'  target {tn} = {tgt}: matched by {matches(tgt)}')
    tot = hit = 0
    for a in range(0, 16):
        for b in range(a, 16):
            tot += 1
            if matches((0, a, b)):
                hit += 1
    rate = hit / tot
    print(f'  calibration: {hit} of {tot} target triples (0, a, b), 0 <= a <= b <= 15, are matched ({rate:.1%})')
    check('P5: calibration hit rate >= 20%', rate >= 0.20)
    print('  pairing test on the atoms K1, C3, T5 (orders 1, 3, 5) against (p, p^3, p^2qr) (Omega 1, 3, 4; tau 2, 4, 12):')
    print('    order (1,3,5) vs Omega (1,3,4): differ at T5; |Aut| (1,3,5) vs tau (2,4,12): no affine relation')


# ---------------------------------------------------------------- P6
def P6():
    section('P6: cross-check of S70 THM-4852 (4) and a second natural prime graph')
    with open(os.path.join(HERE, 'merged_metagraph_n5.json')) as f:
        M = json.load(f)
    G = nx.Graph()
    G.add_nodes_from(range(len(M['merged_nodes'])))
    loops = 0
    for e in M['merged_edges']:
        if e['u'] == e['v']:
            loops += 1
        else:
            G.add_edge(e['u'], e['v'])
    flow = nx.DiGraph()
    flow.add_nodes_from(PR)
    arcs = 0
    for a, b in OWNER:
        for p in PR:
            for q in PR:
                if b % p == 0 and a % q == 0:
                    arcs += 1
                    flow.add_edge(p, q)
    U = nx.Graph()
    U.add_nodes_from(PR)
    for (p, q) in flow.edges():
        if p != q:
            U.add_edge(p, q)
    selfloops = sum(1 for p, q in flow.edges() if p == q)
    check('P6: register-flow graph has 24 arcs (with multiplicity), 19 unordered pairs, no loops',
          arcs == 24 and U.number_of_edges() == 19 and selfloops == 0,
          f'arcs {arcs}, pairs {U.number_of_edges()}, loops {selfloops}')
    check('P6: merged metagraph simple support has 18 edges and 2 loops', G.number_of_edges() == 18 and loops == 2)
    check('P6: no isomorphism flow-support vs metagraph-support', not nx.is_isomorphic(U, G))
    CO = nx.Graph()
    CO.add_nodes_from(PR)
    for a, b in OWNER:
        sup = [p for p in PR if (a * b) % p == 0]
        for p, q in itertools.combinations(sup, 2):
            CO.add_edge(p, q)
    check('P6: no isomorphism co-occurrence graph vs metagraph-support', not nx.is_isomorphic(CO, G),
          f'co-occurrence edges {CO.number_of_edges()}, degree sequence {sorted(d for _, d in CO.degree())}; '
          f'metagraph degree sequence {sorted(d for _, d in G.degree())}')
    B = nx.Graph()
    B.add_nodes_from(range(len(M['merged_nodes'])))
    for e in M['merged_edges']:
        if e['colour'] == 'blue':
            B.add_edge(e['u'], e['v'])
    print(f'  blue graph on merged nodes: {sorted(B.edges())}; cycle basis {nx.cycle_basis(B)}')


def main():
    P1_P2()
    P3()
    P4()
    P5()
    P6()
    print()
    if FAIL:
        print('FAILED CHECKS:', FAIL)
        return 1
    print('ALL PRE-REGISTERED CHECKS PASSED')
    return 0


if __name__ == '__main__':
    sys.exit(main())
