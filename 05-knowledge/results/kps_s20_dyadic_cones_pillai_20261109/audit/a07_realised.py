# a07: realised dyadic cone steps (P4 re-done from the audit's own census a03_triples_le8.txt, not from S19 ports).
#  (a) formula search: every distinct dyadic triple of order <= 8, TT_m with m <= 20, p >= 2, every exit q;
#      second step from every output with m' <= 20;
#  (b) brute force: for one realising tournament of each dyadic triple, explicit joins with TT_m, m <= 5, all p >= 2,
#      all exits q; report every dyadic exit found (must agree with (a) restricted to m <= 5);
#  (c) explicit witness of 33 -> 289 and of the two 1 -> 9 steps, with the second-step test on the outputs.
import sys, collections
from a_lib import *
rows = [tuple(map(int, l.split())) for l in open('a03_triples_le8.txt')]
witness = {}
for n, c, v, h, e, R in rows:
    if ispow2(R) and (h, e, R) not in witness: witness[(h, e, R)] = (n, c, v)
print('distinct dyadic triples of order <= 8:', len(witness))

def steps(h, e, R, MM=20):
    out = []
    for m in range(2, MM + 1):
        for p in range(2, m + 1):
            hp = h + 2 ** (p - 1) * R
            for q in range(1, m + 1):
                if q >= p:
                    Rp = 2 ** (m - q) * hp
                    ep = hp if q == m else 0      # q = m > p: e' = h + 2^(p-1) R ; q = m = p: e' = h
                    if q == m == p: ep = h
                elif q >= 2:
                    Rp = 2 ** (m - q) * h + 3 ** (p - q - 1) * 2 ** (m - p + q) * R
                    ep = 2 ** (q - 1) * R if p == m else 0
                else:
                    Rp = 2 ** (m - 1) * h - 2 ** (m - p) * e + 3 ** (p - 2) * 2 ** (m - p + 1) * R
                    ep = R if p == m else 0
                if ispow2(Rp): out.append(((m, p, q), (hp, ep, Rp)))
    return out
one = {}; two = []
for b in sorted(witness):
    s = steps(*b)
    if s: one[b] = s
    for g, o in s:
        for g2, o2 in steps(*o):
            two.append((b, g, o, g2, o2))
print('(a) bases with a dyadic one-step (m <= 20):', len(one))
for b, s in one.items():
    print('    ', b, '->', sorted(set((g[1], g[2], o[0]) for g, o in s)), ' (m values: %d..%d)' % (min(g[0] for g, o in s), max(g[0] for g, o in s)))
print('    two-step chains:', len(two))

# (b) brute force with actual tournaments
found_bf = collections.defaultdict(set); joins = 0
for b, (n, c, v) in sorted(witness.items()):
    T = from_code(c, n)
    jobs = []; keys = []
    for m in range(2, 6):
        for p in range(2, m + 1):
            jobs.append(one_arc_join(T, v, TT(m), p - 1)); keys.append((m, p))
    res = run(jobs); joins += len(jobs)
    for (m, p), r in zip(keys, res):
        for q in range(1, m + 1):
            z = n + q - 1
            if ispow2(r['R'][z]): found_bf[b].add((m, p, q, r['H'], r['e'][z], r['R'][z]))
formula_m5 = {b: set((g[0], g[1], g[2], o[0], o[1], o[2]) for g, o in s if g[0] <= 5) for b, s in one.items()}
formula_m5 = {b: s for b, s in formula_m5.items() if s}
print('(b) brute-force joins with TT_m, m <= 5:', joins, ' bases with a dyadic exit:', len(found_bf),
      ' agreement with formula search (m <= 5):', dict(found_bf) == formula_m5)
for b, s in sorted(found_bf.items()): print('    ', b, sorted(s))

# (c) explicit witnesses
def show(b, m, p, q):
    n, c, v = witness[b]; T = from_code(c, n)
    N = one_arc_join(T, v, TT(m), p - 1); r = run([N])[0]; z = n + q - 1
    out = (r['H'], r['e'][z], r['R'][z])
    print('(c) base %s realised by class code %d (order %d, root %d); TT_%d entry p=%d exit q=%d -> (h,e,R) = %s; '
          'second dyadic step from it (m <= 20): %s' % (b, c, n, v, m, p, q, out, bool(steps(*out))))
show((33, 4, 64), 3, 3, 1)
show((1, 1, 1), 4, 4, 2)
show((1, 0, 2), 3, 3, 1)
