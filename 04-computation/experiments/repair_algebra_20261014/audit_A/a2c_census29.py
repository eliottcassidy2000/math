#!/usr/bin/env python3
"""Audit A: N(28), N(29) from S_27 alone (memory-lean).  Membership: x in T_e iff (x - 2^e) % 3 == 0 and
(x - 2^e)/3 in S_(e-1).  S_28 = S_27 u T_28; T_29 = 3 S_28 + 2^29 is enumerated as 3z + 2^29 and 9z + 3*2^28 + 2^29
for z in S_27 (chunked)."""
import math, time
import numpy as np

t0 = time.time()
S = np.array([-1], dtype=np.int64)
P = np.array([], dtype=np.int64)
for e in range(1, 28):
    T = 3 * S + (1 << e)
    idx = np.searchsorted(S, T); idx[idx == len(S)] = len(S) - 1
    P = np.union1d(P, T[S[idx] == T]); del idx
    S = np.union1d(S, T); del T
# S = S_27, P = P_27
print(f"N(27) = {len(P)}, |S_27| = {len(S)}  [{time.time()-t0:.1f}s]", flush=True)


def in_S(x):
    j = np.searchsorted(S, x); j[j == len(S)] = len(S) - 1
    return S[j] == x


def in_T(x, e):
    y = x - (1 << e)
    ok = (y % 3 == 0)
    res = np.zeros(len(x), dtype=bool)
    if ok.any():
        res[ok] = in_S(y[ok] // 3)
    return res


CH = 2_000_000
h28 = []; h29 = []
for i in range(0, len(S), CH):
    z = S[i:i + CH]
    x28 = 3 * z + (1 << 28)                       # T_28
    h28.append(x28[in_S(x28)])                    # T_28 n S_27
    for y in (z, x28):                            # S_28 = S_27 u T_28
        x = 3 * y + (1 << 29)                     # T_29
        h29.append(x[in_S(x) | in_T(x, 28)])      # T_29 n S_28
P28 = np.union1d(P, np.unique(np.concatenate(h28)))
P29 = np.union1d(P28, np.unique(np.concatenate(h29)))
for E, N in ((27, len(P)), (28, len(P28)), (29, len(P29))):
    print(f"E={E} N(E)={N}  N sqrt(E)/2^E = {N*math.sqrt(E)/2**E:.6f}")
for (E, N1, N2) in ((28, len(P), len(P28)), (29, len(P28), len(P29))):
    print(f"beta_local({E}) = {-math.log(N2/(2*N1))/math.log(E/(E-1)):.3f}")
print(f"[{time.time()-t0:.1f}s]")
