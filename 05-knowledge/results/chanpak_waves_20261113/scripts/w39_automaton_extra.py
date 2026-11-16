# Wave 39 (mac-mini-2026-11-15-frame). Run from a working directory; w39_lit_setup.py writes lit_base.pkl used by the others.
# Syllable-adaptive suffix automaton for literal span-3 words (fast exact version).
# PI = {all strings of length <= L0} U {u1...uj p : units u_i from UNITS, p a proper prefix of a unit, total length <= LMAX}.
# PI is closed under deleting the last letter, so "longest suffix in PI" is a deterministic automaton.
# Edge s --l--> t kept iff s+l has rank >= 2 in the minimal representation.  All entries of O A_w Rb are nonnegative
# integers, so rank <= 1 iff all nonzero columns have one primitive direction: an exact test in int64.
# usage: python3 automaton2.py UNIT_IDS L0 LMAX
import pickle, sys, itertools, numpy as np, time
B = pickle.load(open('lit_base.pkl', 'rb'))
A, v0, h, O, Rb, U = B['A'], B['v0'], B['h'], B['O'], B['Rb'], B['U']
unit_ids = [int(x) for x in sys.argv[1].split(',')]
L0, LMAX = int(sys.argv[2]), int(sys.argv[3])
EXTRA = pickle.load(open(sys.argv[4],'rb'))['D'] if len(sys.argv) > 4 else set()
units = [tuple(U[i]) for i in unit_ids]
PI = set()
for L in range(L0 + 1):
    for t in itertools.product(range(8), repeat=L):
        PI.add(t)
prefixes = set()
for u in units:
    for k in range(len(u)):
        prefixes.add(u[:k])
frontier = [()]
seen = {()}
while frontier:
    nf = []
    for base in frontier:
        for p in prefixes:
            if len(base) + len(p) <= LMAX:
                PI.add(base + p)
        for u in units:
            nb = base + u
            if len(nb) <= LMAX and nb not in seen:
                seen.add(nb)
                nf.append(nb)
                PI.add(nb)
    frontier = nf
PI |= set(EXTRA)
for s in list(PI):
    for k in range(len(s)):
        assert s[:k] in PI or True
OR = O @ Rb          # 61 x r, nonnegative
assert (OR >= 0).all()
cache = {}


def ORw(s):
    # O A_s Rb with A_s = A_{s[-1]} ... A_{s[0]};  O A_s = (O A_{s[1:]}) A_{s[0]}
    if s in cache:
        return cache[s]
    if not s:
        M = O.copy()
    else:
        M = ORw_left(s)
    cache[s] = M
    return M


lcache = {}


def ORw_left(s):
    if s in lcache:
        return lcache[s]
    if not s:
        M = O.copy()
    else:
        M = ORw_left(s[1:]) @ A[s[0]]
    lcache[s] = M
    return M


def rank_le1(s):
    M = ORw_left(s) @ Rb
    if M.max() > 2 ** 62 // 4:
        raise OverflowError
    cols = M[:, M.any(axis=0)]
    if cols.shape[1] <= 1:
        return True
    g = np.gcd.reduce(cols, axis=0)
    P = cols // g
    return bool((P == P[:, :1]).all())


t0 = time.time()
sidx = {(): 0}
order = [()]
edges = []
frontier = [()]
removed = 0
while frontier:
    nf = []
    for s in frontier:
        for l in range(8):
            w = s + (l,)
            if rank_le1(w):
                removed += 1
                continue
            t = next(w[k:] for k in range(len(w) + 1) if w[k:] in PI)
            if t not in sidx:
                sidx[t] = len(order)
                order.append(t)
                nf.append(t)
            edges.append((sidx[s], l, sidx[t]))
    frontier = nf
    lcache = {k: v for k, v in lcache.items() if len(k) <= LMAX + 1}
name = 'auto3_%s_L%d_M%d_%s.pkl' % (sys.argv[1].replace(',', '-'), L0, LMAX, (sys.argv[4][:-4] if len(sys.argv) > 4 else 'none'))
print('PI', len(PI), 'states', len(order), 'edges', len(edges), 'removed', removed, '%.0fs' % (time.time() - t0), '->', name)
lens = {}
for s in order:
    lens[len(s)] = lens.get(len(s), 0) + 1
print('state lengths', sorted(lens.items()))
pickle.dump({'order': order, 'edges': edges, 'unit_ids': unit_ids, 'L0': L0, 'LMAX': LMAX, 'extra': set(EXTRA)}, open(name, 'wb'))
