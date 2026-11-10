# a11: scope of "no two dyadic cone steps chain" (THM-5109 (E3)) vs "pure cone inductions are impossible" (README,
# THM-5109 consequence). (E3) only covers new dyadic ports INSIDE the attached transitive gadget. Here we take one
# realisation of every distinct dyadic triple of order <= 8, apply every cone step TT_m (m <= MMAX, entry p >= 2),
# and look for a dyadic port ANYWHERE in the output N (any vertex w, in T or in the gadget), and for Q(N) a power
# of two (a dyadic source port above N). Any hit is a cone step whose output again carries a cone.
import sys, collections
from a_lib import *
MMAX = int(sys.argv[1]) if len(sys.argv) > 1 else 4
rows = [tuple(map(int, l.split())) for l in open('a03_triples_le8.txt')]
witness = {}
for n, c, v, h, e, R in rows:
    if ispow2(R) and (h, e, R) not in witness: witness[(h, e, R)] = (n, c, v)
hits_T = []; hits_G = []; hits_Q = []; total = 0
for b, (n, c, v) in sorted(witness.items()):
    T = from_code(c, n)
    jobs = []; keys = []
    for m in range(2, MMAX + 1):
        if n + m > 14: continue
        for p in range(2, m + 1):
            jobs.append(one_arc_join(T, v, TT(m), p - 1)); keys.append((m, p))
    if not jobs: continue
    for (m, p), r in zip(keys, run(jobs)):
        total += 1
        assert r['H'] == b[0] + 2 ** (p - 1) * b[2]
        for w in range(n + m):
            if ispow2(r['R'][w]):
                rec = (b, (m, p), w, r['H'], r['e'][w], r['R'][w])
                (hits_T if w < n else hits_G).append(rec)
        if ispow2(r['Q']): hits_Q.append((b, (m, p), r['H'], r['Q']))
print('cone steps examined: %d (bases: %d dyadic triples of order <= 8; TT_m, m <= %d, p >= 2)' % (total, len(witness), MMAX))
print('dyadic ports of the output inside the gadget: %d' % len(hits_G))
for x in hits_G[:10]: print('   ', x)
print('dyadic ports of the output inside the base T (old vertices): %d' % len(hits_T))
for x in hits_T[:15]: print('   ', x)
print('outputs N with Q(N) a power of two (new source above N is dyadic): %d' % len(hits_Q))
for x in hits_Q[:10]: print('   ', x)
vals = sorted(set(x[3] for x in hits_T) | set(x[2] for x in hits_Q))
print('values H(N) of cone outputs that again carry a dyadic port outside the gadget:', vals[:40])
