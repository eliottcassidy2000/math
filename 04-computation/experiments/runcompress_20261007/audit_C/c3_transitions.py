#!/usr/bin/env python3
"""audit_C item 3: THM-4603 (3) long-run transitions, independent recomputation.

Own code for: cycle enumeration (primitive U-words, Sum <= 8, length <= 4, plus the halving point 0), the exact limit pair
iteration, and a STRICTER classification than the session's transitions.py:
  ABSORB  : u = v and k = 0 at some time (before the pair becomes periodic)
  ANCHOR  : jointly periodic, u and v on the same cycle and in phase (u = v), k constant
  SHIFT   : same cycle, out of phase (k bounded, periodic)
  DRIFT   : different cycles; rate = dens(u-cycle) - dens(v-cycle) (checked against the measured rate);
            different cycles with EQUAL density are reported as 'DRIFT0' (the session's code would call them SHIFT)
Realisability at the universal state: at the end of the two-run, X - 1 = 27 (Y - 1) = 3^J 2^r u0 with r in {1, 2}, so a run
starting right there must have v2(c - 1) in {1, 2} at its cycle point (child: c'; source: c), i.e. its first letter is 1 or
>= 3. Words starting with the letter 2 (and the halving point 0) are NOT realisable post-run.
Actual-integer checks (residual sources constructed with a prescribed post-run class, K, t, J random):
  child runs of letters 3, 4 (re-anchoring, debts 4, 9), child run of letter 5 (drift 7/15), child (1,1,4) absorption at 23,
  source (1,2,5) absorption at 20, child ones-run (-5 at 12), source ones-run (+1/2 after 8).
"""
import itertools, random
from fractions import Fraction as Fr
from math import gcd

rnd = random.Random(31337)


def par(q):
    return q.numerator & 1


def T(q):
    return (3 * q + 1) / 2 if par(q) else q / 2


def v2q(q):
    n = q.numerator
    return (n & -n).bit_length() - 1 if n else 10 ** 9


def cpoint(w):
    m, S, B, Q = len(w), sum(w), 0, 1
    for a in w:
        B = 3 * B + Q
        Q *= 2 ** a
    return Fr(B, 2 ** S - 3 ** m)


def tcycle(c):
    cyc = [c]
    x = T(c)
    while x != c:
        cyc.append(x)
        x = T(x)
        assert len(cyc) < 10 ** 4
    return cyc


cycles = []
seen = set()
for L in range(1, 5):
    for w in itertools.product(range(1, 9), repeat=L):
        if sum(w) > 8:
            continue
        if any(L % p == 0 and w == w[:p] * (L // p) for p in range(1, L)):
            continue
        c = cpoint(w)
        cy = tcycle(c)
        key = frozenset(cy)
        if key in seen:
            continue
        seen.add(key)
        cycles.append((w, c, cy))
cycles.append(((0,), Fr(0), [Fr(0)]))
cycle_of = {}
for idx, (w, c, cy) in enumerate(cycles):
    for x in cy:
        cycle_of[x] = idx


def dens(cy):
    return Fr(sum(par(x) for x in cy), len(cy))


def limit(k, e, which, c, maxit=100000):
    if which == 'child':
        v = c
        u = Fr(3) ** k * v + e
    else:
        u = c
        v = (u - e) / Fr(3) ** k
    seen = {}
    for s in range(maxit):
        if u == v and k == 0:
            return dict(kind='ABSORB', time=s)
        key = (u, v)
        if key in seen:
            s0, k0 = seen[key]
            per = s - s0
            rate = Fr(k - k0, per)
            # cycles of u and v (u, v are periodic points now)
            cu, cv = tcycle(u), tcycle(v)
            same = frozenset(cu) == frozenset(cv)
            if same:
                kind = 'ANCHOR' if u == v else 'SHIFT'
                assert rate == 0
            else:
                rho = dens(cu) - dens(cv)
                assert rho == rate, (rho, rate)
                kind = 'DRIFT' if rho != 0 else 'DRIFT0'
            return dict(kind=kind, entry=s0, k_entry=k0, period=per, rate=rate, u=u, v=v)
        seen[key] = (s, k)
        k += par(u) - par(v)
        u, v = T(u), T(v)
    return dict(kind='UNRESOLVED')


states = {
    'universal (3,-26)': (3, Fr(-26)),
    'J=2 (4,10)': (4, Fr(10)),
    'J=1 D=3 (3,13)': (3, Fr(13)),
    're-anchored (-5, 3^-5 - 1)': (-5, Fr(1, 243) - 1),
    'reset (1,1)': (1, Fr(1)),
}
print(f"{len(cycles)} cycles (primitive U-words with Sum <= 8, length <= 4, plus 0)")
allres = {}
for name, (k, e) in states.items():
    for w, c, cy in cycles:
        for which in ('child', 'source'):
            allres[(name, w, which)] = limit(k, e, which, c)
kinds = {}
for key, r in allres.items():
    kinds[r['kind']] = kinds.get(r['kind'], 0) + 1
print("outcome counts over 5 states x 55 cycles x 2:", kinds)

# ---- compare with the session's transitions.out ----
import re
sess = {}
state = None
names = {'universal (3,1-27)': 'universal (3,-26)', 'J=2 (4,10)': 'J=2 (4,10)', 'J=1 D=3 (3,13)': 'J=1 D=3 (3,13)',
         're-anchored (-5, 3^-5 - 1)': 're-anchored (-5, 3^-5 - 1)', 'reset (1,1)': 'reset (1,1)'}
for line in open('../transitions.out'):
    m = re.match(r'== state (.*): u = ', line)
    if m:
        state = names[m.group(1)]
        continue
    m = re.match(r'\s+(child|source)\s+run on (\(.*?\)): (\w+)(.*)', line)
    if m and state:
        which, w, kind, rest = m.groups()
        w = tuple(int(x) for x in w.strip('()').split(',') if x.strip())
        d = dict(kind=kind)
        if kind == 'ABSORB':
            d['time'] = int(re.search(r'time (\d+)', rest).group(1))
        else:
            d['entry'] = int(re.search(r't=\s*(\d+)', rest).group(1))
            d['k_entry'] = int(re.search(r'k_entry=\s*(-?\d+)', rest).group(1))
            d['rate'] = Fr(re.search(r'drift/step=(\S+)', rest).group(1))
        sess[(state, w, which)] = d
agree = disagree = 0
relabel = []
for key, d in sess.items():
    r = allres[key]
    if d['kind'] == 'ABSORB':
        ok = r['kind'] == 'ABSORB' and r['time'] == d['time']
    else:
        ok = r['kind'] != 'ABSORB' and r['entry'] == d['entry'] and r['k_entry'] == d['k_entry'] and r['rate'] == d['rate']
        if ok and r['kind'] != d['kind']:
            relabel.append((key, d['kind'], r['kind']))
    agree += ok
    disagree += not ok
    if not ok:
        print("DISAGREE", key, d, r)
print(f"session table lines compared: {len(sess)}; agree (time / entry / k_entry / rate): {agree}; disagree: {disagree}")
print(f"label differences (session label -> audit label): {len(relabel)}")
for key, a, b in relabel:
    print("   ", key, a, '->', b)


def realisable(which, c, k, e):
    """post-run realisability at the universal state: the run's starting point needs v2(. - 1) in {1, 2}"""
    pt = c if which == 'source' else c
    return v2q(pt - 1) in (1, 2)


# ---- universal state: complete lists ----
name = 'universal (3,-26)'
k, e = states[name]
print(f"\n== complete universal-state lists (R = realisable right at the end of the two-run) ==")
cats = {'ABSORB': [], 'ANCHOR': [], 'DRIFT+': [], 'DRIFT-': [], 'DRIFT0': [], 'SHIFT': []}
for w, c, cy in cycles:
    for which in ('child', 'source'):
        r = allres[(name, w, which)]
        R = 'R' if realisable(which, c, k, e) else '-'
        if r['kind'] == 'ABSORB':
            cats['ABSORB'].append(f"{which} {w} t={r['time']} {R}")
        elif r['kind'] == 'DRIFT':
            cats['DRIFT+' if r['rate'] > 0 else 'DRIFT-'].append(f"{which} {w} {float(r['rate']):+.3f} {R}")
        elif r['kind'] == 'ANCHOR':
            cats['ANCHOR'].append(f"{which} {w} k={r['k_entry']} t={r['entry']} {R}")
        else:
            cats[r['kind']].append(f"{which} {w} {R}")
for cat, lst in cats.items():
    print(f"{cat} ({len(lst)}):")
    for s in lst:
        print("    ", s)

# ---- spot re-derivations ----
r = allres[(name, (1,), 'child')]
print("\nchild ones-run:", r['kind'], 'entry', r['entry'], 'k', r['k_entry'], 'limit pair start', Fr(3) ** 3 * (-1) - 26, -1)
r = allres[(name, (1,), 'source')]
print("source ones-run:", r['kind'], 'entry', r['entry'], 'k', r['k_entry'], 'rate', r['rate'], 'child limit', (Fr(-1) + 26) / 27)
for w in ((3,), (4,)):
    r = allres[(name, w, 'child')]
    print(f"child run {w}:", r['kind'], 'entry', r['entry'], 'k', r['k_entry'])


# ---- actual integers ----
def make_source(K, J, target_Y, M):
    """residual source n = 2^K t - 1 with two-run length J (>= 3) whose post-run child Y = T^(2J)(y_3) satisfies
    Y = target_Y mod 2^M (target_Y a 2-adic rational); r is forced by v2(target_Y - 1)."""
    r = v2q(target_Y - 1)
    assert r in (1, 2)
    j = 2 * J + r
    Mod = 1 << (M + j + 8)
    w = (target_Y - 1) / 2 ** r                         # odd 2-adic unit
    wmod = (w.numerator * pow(w.denominator, -1, Mod)) % Mod
    u0 = (wmod * pow(3, -(J - 3), Mod)) % Mod           # Y - 1 = 3^(J-3) 2^r u0
    x = (1 + (1 << j) * u0) % Mod
    t = ((x + 1) // 2 * pow(3, -(K - 1), Mod)) % Mod
    if t % 2 == 0:
        raise ValueError
    t += Mod * rnd.getrandbits(64)
    return t


def chain_after_tworun(K, t):
    n = (t << K) - 1
    x = 2 * 3 ** (K - 1) * t - 1
    j = ((x - 1) & -(x - 1)).bit_length() - 1
    J = (j - 1) // 2
    y = (x + 1) // 27 - 1
    u, v, k = x, y, 3
    for s in range(2 * J):
        k += (u & 1) - (v & 1)
        u = (3 * u + 1) >> 1 if u & 1 else u >> 1
        v = (3 * v + 1) >> 1 if v & 1 else v >> 1
    assert k == 3 and u - 1 == 27 * (v - 1)
    return u, v, k


def step(u, v, k):
    k += (u & 1) - (v & 1)
    u = (3 * u + 1) >> 1 if u & 1 else u >> 1
    v = (3 * v + 1) >> 1 if v & 1 else v >> 1
    return u, v, k


def check(desc, which, w, ntrials, test):
    c = cpoint(w)
    target_Y = c if which == 'child' else 1 + (c - 1) / 27
    ok = 0
    for _ in range(ntrials):
        K, J = rnd.randint(6, 40), rnd.randint(3, 9)
        t = make_source(K, J, target_Y, 260)
        u, v, k = chain_after_tworun(K, t)
        assert test(u, v, k), (desc, K, J)
        ok += 1
    print(f"actual integers: {desc}: {ok}/{ntrials} sources pass")


def t_letter3(u, v, k):
    # re-anchoring: at time 15 the debt is 4 and u, v are 2-adically near the same cycle point p of the 1/5-cycle
    for s in range(15):
        u, v, k = step(u, v, k)
    if k != 4:
        return False
    # the limit point at time 15
    r = allres[(name, (3,), 'child')]
    p = r['u']
    # u - p = 3^4 (v - p) exactly (anchored at the periodic point p), as rationals
    return Fr(u) - p == 81 * (Fr(v) - p) and all(step(u, v, k)[2] == 4 for _ in range(1))


def t_letter4(u, v, k):
    for s in range(25):
        u, v, k = step(u, v, k)
    p = allres[(name, (4,), 'child')]['u']
    return k == 9 and Fr(u) - p == 3 ** 9 * (Fr(v) - p)


def t_letter5(u, v, k):
    ks = []
    for s in range(9 + 60 * 3):
        ks.append(k)
        u, v, k = step(u, v, k)
    ks.append(k)
    return ks[9] == 6 and ks[69] == 6 + 28 and ks[129] == 6 + 56 and ks[189] == 6 + 84


def t_absorb(time):
    def f(u, v, k):
        for s in range(time):
            if u == v and k == 0:
                return False          # earlier than claimed
            u, v, k = step(u, v, k)
        return u == v and k == 0
    return f


def t_child_ones(u, v, k):
    for s in range(12):
        u, v, k = step(u, v, k)
    return k == -5 and 243 * (u + 1) == v + 1


def t_source_ones(u, v, k):
    ks = []
    for s in range(200):
        u, v, k = step(u, v, k)
        ks.append(k)
    return ks[7] == 6 and all(ks[7 + 2 * i] == 6 + i for i in range(90))


check("child run letter 3: debt 4 at step 15, anchored at the limit cycle point", 'child', (3,), 25, t_letter3)
check("child run letter 4: debt 9 at step 25, anchored", 'child', (4,), 25, t_letter4)
check("child run letter 5: debt 6 at step 9, +28 per 60 steps (rate 7/15)", 'child', (5,), 25, t_letter5)
check("child run (1,1,4): absorbed exactly at step 23", 'child', (1, 1, 4), 25, t_absorb(23))
check("source run (1,2,5): absorbed exactly at step 20", 'source', (1, 2, 5), 25, t_absorb(20))
check("child ones-run: debt -5 at step 12, anchored at -1", 'child', (1,), 25, t_child_ones)
check("source ones-run: debt 6 at step 8, then +1 per 2 steps (200 steps)", 'source', (1,), 25, t_source_ones)
# unrealisable: source run of (2,6) needs v2(X - 1) = 4
c = cpoint((2, 6))
print("source run (2,6): v2(c - 1) =", v2q(c - 1), "but v2(X - 1) = r in {1,2} at the end of every two-run: not realisable")
print("ALL CHECKS PASSED")
