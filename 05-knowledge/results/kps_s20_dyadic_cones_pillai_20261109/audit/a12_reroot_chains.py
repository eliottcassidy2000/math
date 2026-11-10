# a12: cone chains when the next dyadic port may be ANY vertex of the output (re-rooting allowed), or the source of a
# new source above the output (Q(N) a power of two).  Breadth-first from the dyadic rooted tournaments of order <= NB,
# cone steps = one-arc joins with TT_m (m <= MMAX) at entry p >= 2 (value h + 2^(p-1) R); a state is a rooted
# tournament (N, w) with R_w(N) a power of two.  Deduplication by an isomorphism-invariant key (approximate: the
# multiset of vertex triples plus the root triple), which can only merge states, never invent chains: every chain
# printed is re-verified from scratch by explicit construction.
import sys
from a_lib import *
NB = int(sys.argv[1]); MMAX = int(sys.argv[2]); DEPTH = int(sys.argv[3]); NMAX = int(sys.argv[4])
S19 = sys.argv[5]
def key(r, w):
    deck = tuple(sorted(zip(r['e'], r['R'])))
    return (r['H'], r['e'][w], r['R'][w], deck)
level = []; seen = set()
for n in range(1, NB + 1):
    codes = [int(l.split()[0]) for l in open('%s/cls_%d.txt' % (S19, n)) if l.strip()]
    tours = [from_code(c, n) for c in codes]
    for A, r in zip(tours, run(tours)):
        for w in range(n):
            if ispow2(r['R'][w]):
                k = key(r, w)
                if k not in seen: seen.add(k); level.append((A, w, [(r['H'], r['e'][w], r['R'][w])]))
    # source above A (Q dyadic) is already a rooted tournament of order n + 1 <= NB + 1: covered by the census of order n+1
print('depth 0: %d dyadic rooted states (order <= %d)' % (len(level), NB))
best = []
for d in range(1, DEPTH + 1):
    nxt = []
    for (T, v, hist) in level:
        n = len(T)
        jobs = []
        for m in range(2, MMAX + 1):
            if n + m + 1 > NMAX: continue
            for p in range(2, m + 1): jobs.append(one_arc_join(T, v, TT(m), p - 1))
        if not jobs: continue
        for N, r in zip(jobs, run(jobs)):
            for w in range(len(N)):
                if ispow2(r['R'][w]):
                    k = key(r, w)
                    if k not in seen: seen.add(k); nxt.append((N, w, hist + [(r['H'], r['e'][w], r['R'][w])]))
            if ispow2(r['Q']) and len(N) + 1 <= NMAX:
                S = ordinal_sum(TT(1), N); rs = run([S])[0]
                k = key(rs, 0)
                if k not in seen: seen.add(k); nxt.append((S, 0, hist + [('src', rs['H'], rs['e'][0], rs['R'][0])]))
    level = nxt
    print('depth %d: %d new dyadic rooted states (orders up to %d); values: %s' % (d, len(level), NMAX, sorted(set(s[2][-1][0] if s[2][-1][0] != 'src' else s[2][-1][1] for s in level))[:30]))
    for s in level[:6]: print('    chain:', s[2])
    if not level: break
