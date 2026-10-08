#!/usr/bin/env python3
"""audit H: results note 6b, "Why 5I - J keeps appearing": (a) in the S_5-symmetric metric (Q = 5I - J) every translation is
strictly balanced and every non-translation coupling is at exact equality, for a rank-4 Z_5 map; (b) count of lifts of the
constants r mod 25 (r_i = -m_i i mod 5 fixed) for which 5I - J balances every nonzero 2-step block, for two maps."""
import itertools, time
import numpy as np
from hdp import Tables, balance_exact
from hcore import std_r, adj_int, det_int
Q = [[4, -1, -1, -1], [-1, 4, -1, -1], [-1, -1, 4, -1], [-1, -1, -1, 4]]
tb = Tables(5, [1, 2, 3, 7, 11])
adjQ, detQ = adj_int(Q), det_int(Q)
for a in (1, 2, 3, 4):
    for b in range(5):
        if (a, b) == (1, 0): continue
        D = tb.Dab[a, b].astype(np.int64)
        S = (np.einsum('xy,yx->', D, np.array(adjQ)) * np.array(Q) - 2 * detQ * D)
        ev = np.linalg.eigvalsh(S.astype(float))
        st = 'strict' if ev[0] > 1e-9 else ('EQUALITY (PSD singular)' if ev[0] > -1e-9 else 'FAILS')
        if b in (0, 1): print(f"  coupling ({a},{b}): {st}")
for m in ([1, 2, 3, 7, 11], [1, 6, 11, 7, 2]):
    t0 = time.time(); base = std_r(m); good = 0; std_ok = None
    for t in itertools.product(range(5), repeat=5):
        r = [base[i] + 5 * t[i] for i in range(5)]
        tb = Tables(5, m, r)
        T = tb.table(2).reshape(-1, 4, 4)
        nz = np.any(T.reshape(len(T), -1) != 0, axis=1)
        U = np.unique(T[nz].reshape(-1, 16), axis=0).reshape(-1, 4, 4)
        ok = bool(balance_exact(Q, U).all())
        good += ok
        if t == (0, 0, 0, 0, 0): std_ok = ok
    print(f"m = {m}: 5I - J certifies k = 2 for {good} of 3125 lifts r mod 25 (standard constants: {std_ok})  [{time.time() - t0:.1f}s]")
