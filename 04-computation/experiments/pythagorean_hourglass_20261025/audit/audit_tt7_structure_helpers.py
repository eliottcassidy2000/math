# TT7(P_p) = p(p-1)/2 * TT5(T[M(0,1)]) by depth-first enumeration (audit helper).
import numpy as np
from audit_gf import GF
from audit_tt_census import paley_M_adjacency, dfs_count
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
