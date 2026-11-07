"""THM-4891 / README section 7: how many pairs 'flip' (depend on the order of the other two elements) in the
n = 4 hexagon-coherent orientations with exactly 7 sinks?  Auditor's own enumeration over the 12 hexagon classes
(from aud_perm.Perm: geometric classes of the Cayley graph), with the flip count read off each orientation."""
import itertools, sys, os
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from aud_perm import Perm, sinks_and_acyclic, describe_edge

P = Perm(4)
cid, cpar, nc, contra, nonpar = P.classes(['hex'])
assert nc == 12 and contra == 0
# per class: (pair, order of others)
cls_key = {}
for e in range(len(P.E)):
    pair, before, others = describe_edge(P, e)
    cls_key[cid[e]] = (pair, others)
dist = {}
session_tab = {"01|23": 0, "01|32": 1, "02|13": 2, "02|31": 0, "03|12": 0, "03|21": 0,
               "12|03": 1, "12|30": 2, "13|02": 3, "13|20": 3, "23|01": 2, "23|10": 3}
sess_flips = sum(1 for pr in ('01', '02', '03', '12', '13', '23')
                 if len({session_tab[k] for k in session_tab if k.startswith(pr + '|')}) == 2)
count7 = 0
for M in range(1 << nc):
    bits = [((M >> cid[e]) & 1) ^ cpar[e] for e in range(len(P.E))]
    sinks, ac, coh = sinks_and_acyclic(P, bits)
    if len(sinks) != 7:
        continue
    count7 += 1
    # kappa(pair, others) = element first at the head of the class edges
    kap = {}
    for e in range(len(P.E)):
        u, v, i = P.E[e]
        head = P.V[v] if bits[e] else P.V[u]
        pair, others = cls_key[cid[e]]
        first = pair[0] if head.index(pair[0]) < head.index(pair[1]) else pair[1]
        kap[(pair, others)] = first
    flips = sum(1 for pair in itertools.combinations(range(4), 2)
                if len({kap[k] for k in kap if k[0] == pair}) == 2)
    dist[flips] = dist.get(flips, 0) + 1
print('n=4 hexagon-coherent orientations with 7 sinks: %d; distribution of the number of flipping pairs: %s'
      % (count7, dict(sorted(dist.items()))))
print('session witness table (README section 7): flipping pairs = %d' % sess_flips)
