# Audit check of HYP-9390 Remark 1 (coefficient ring of the class counts).  For a tournament tau on [k], the number of
# induced copies of tau in P_q is (1/|Aut tau|) * 2^-C(k,2) * sum_(G subset K_k) sign_tau(G) T_(G, order), with
# sign_tau(G) = prod over edges of G of (+1 if tau orients i -> j for i < j, else -1).  The audit's symbolic engine
# gives the exact polynomial and wheel parts.  Example classes: the 3-cycle (k = 3), the regular 5-tournament R5
# (|Aut| = 5), and the transitive TT5 (|Aut| = 1); numerical check against direct counts in P_19 and P_23.
import itertools, sys
from fractions import Fraction
from math import comb
sys.path.insert(0, '.')
import audit_census_symbolic as acs
from audit_graph_sums import S_norm, W4
import numpy as np

def class_count_symbolic(k, tau_arcs):
    arcset = set(tau_arcs)
    tot = {}
    for m in range(1 << comb(k, 2)):
        s = 1
        for n in acs.bits(m):
            i, j = acs.EDGES[k][n]
            s *= 1 if (i, j) in arcset else -1
        tot = acs.combo_add(tot, {key: acs.pscale(v, s) for key, v in acs.T_oriented(k, m).items()})
    # automorphisms of tau
    aut = sum(1 for p in itertools.permutations(range(k)) if all((p[i], p[j]) in arcset for (i, j) in tau_arcs))
    return tot, aut

def direct_count(p, k, tau_arcs):
    sq = {(x * x) % p for x in range(1, p)}
    arcset = set(tau_arcs); cnt = 0
    aut = sum(1 for pm in itertools.permutations(range(k)) if all((pm[i], pm[j]) in arcset for (i, j) in tau_arcs))
    for S in itertools.combinations(range(p), k):
        for pm in itertools.permutations(S):
            if all(((pm[j] - pm[i]) % p) in sq for (i, j) in tau_arcs): cnt += 1
    return cnt // aut

out = []
C3 = [(0, 1), (1, 2), (2, 0)]
R5 = [(i, (i + d) % 5) for i in range(5) for d in (1, 2)]
TT5 = [(i, j) for i in range(5) for j in range(i + 1, 5)]
for name, k, arcs in [('C3', 3, C3), ('R5 (regular 5-tournament)', 5, R5), ('TT5', 5, TT5)]:
    tot, aut = class_count_symbolic(k, arcs)
    lab = acs.express(tot, k)
    den = Fraction(1, aut * (1 << comb(k, 2)))
    poly = [Fraction(c) * den for c in lab.get('1', (0,))]
    w4 = [Fraction(c) * den for c in lab.get('T_W4', (0,))]
    dens = sorted({f.denominator for f in poly + w4})
    out.append(f'{name}: |Aut| = {aut}; count = sum_i c_i q^i + d T_W4 with c = {[str(f) for f in poly]}, d = {[str(f) for f in w4]}; '
               f'denominators {dens}')
    for p in (19, 23):
        val = sum(f * p ** i for i, f in enumerate(poly)) + sum(f * p ** i for i, f in enumerate(w4)) * p * (p - 1) * S_norm(p, 5, W4)
        if k <= 5 and p <= 23:
            out.append(f'    q = {p}: formula {val}, direct count {direct_count(p, k, arcs)}')
print('\n'.join(out))
open('audit_class_counts.out', 'w').write('\n'.join(out) + '\n')
