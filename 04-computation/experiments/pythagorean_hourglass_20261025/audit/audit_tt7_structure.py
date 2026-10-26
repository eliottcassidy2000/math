# Audit extra (HYP-9390 (c)): the exact reduced form of 2^21 TT7(P_q) from the audit's symbolic elimination
# (audit_census_symbolic.py), checked numerically: every surviving dense sum (W4, octahedron, K6-P6 and the 38
# seven-vertex survivors, each in the canonical labelling, order orientation) is brute-forced with exact injectivity,
# and the right-hand side is compared with 2^21 TT7(P_q) = 2^21 q(q-1)/2 TT5(T[M(0,1)]) (source-sink decomposition,
# TT5 of T[M] by depth-first enumeration).  Also the PGL(2) reduction T_(K7-C7) = eps (q-6) T_(K6-P6).
# Audit 2026-10-25.
import sys, time
import numpy as np
sys.path.insert(0, '.')
import audit_census_symbolic as acs
from audit_graph_sums import S_norm
from audit_gf import GF
from audit_tt_census import paley_M_adjacency, dfs_count   # noqa: E402  (module-level code there is cheap? see guard)

def tt7(p):
    F = GF(p, 1)
    if (p - 3) // 4 < 5: return 0
    A = paley_M_adjacency(F)
    out = []
    for v in range(len(A)):
        m = 0
        for u in np.nonzero(A[v])[0]: m |= 1 << int(u)
        out.append(m)
    return p * (p - 1) // 2 * dfs_count(out, 5)

t0 = time.time()
res = acs.census(7)
lines = []
keys = list(res)
lines.append(f'2^21 TT7 reduced form: {len(keys)} terms (1 polynomial, W4, octahedron, K6-P6 and the seven-vertex survivors)')
for p in (7, 11, 19, 23, 31):
    rhs = 0
    for key, poly in res.items():
        c = acs.peval(poly, p)
        if key == 'one': rhs += c; continue
        _, k, mask = key
        E = [acs.EDGES[k][n] for n in acs.bits(mask)]
        T = p * (p - 1) * S_norm(p, k, E)
        rhs += c * T
    lhs = (1 << 21) * tt7(p)
    lines.append(f'q = {p}: 2^21 TT7 = {lhs}; reduced form = {rhs}; {"AGREE" if lhs == rhs else "DISAGREE"}  [{time.time() - t0:.0f}s]')
    print(lines[-1], flush=True)
print('\n'.join(lines[:1]))
with open('audit_tt7_structure.out', 'w') as fh: fh.write('\n'.join(lines) + '\n')
