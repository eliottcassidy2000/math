# Audit of HYP-9391 (B): numerical check of the auditor's proof of Y04 = -Y02 and Y06 = Y05 (q = 3 mod 4).
# With A(u) = sum_{s != 0,1} phi(s) phi(u-s) L(s)^2 and A2(u) = sum_{s != 0,1} phi(s) phi(s-1) phi(u-s) L(s)^2:
#   -Y02 = sum_{u != 0,1} phi(u) L(u) A(u),          Y04 = sum_{u != 0,1} phi(u) phi(u-1) L(u) A(u),
#    Y05 = sum_{u != 0,1} phi(u) L(u) A2(u),          Y06 = sum_{u != 0,1} phi(u) phi(u-1) L(u) A2(u),
# and on the set U- = {u : phi(u-1) = -1} the involution u -> u/(u-1) has no fixed point, preserves U-, fixes L(u),
# A(u), A2(u) and flips phi(u); so the U- parts vanish and Y04 + Y02 = Y06 - Y05 = 0.  Checked termwise here; the
# same steps at q = 1 mod 4 are shown to fail (phi(-1) = +1 there).  Writes y_symmetry_proof_check.out.
import os
import numpy as np
from sympy import primerange
HERE = os.path.dirname(os.path.abspath(__file__))
out = []
def tabs(q):
    phi = np.full(q, -1, dtype=np.int64); phi[0] = 0; phi[[(x * x) % q for x in range(1, q)]] = 1
    w = np.arange(q)
    L = np.array([int((phi[(w * (w - 1)) % q] * phi[(w - l) % q]).sum()) for l in range(q)], dtype=np.int64)
    inv = np.zeros(q, dtype=np.int64); inv[1:] = [pow(int(x), q - 2, q) for x in range(1, q)]
    return phi, L, inv
def closed(q, phi, L, inv):
    Y02 = Y04 = Y05 = Y06 = 0
    for s in range(2, q):
        for t in range(2, q):
            if s == t: continue
            r = t * inv[s] % q
            Y02 += -phi[s] * phi[t] * phi[t - 1] * L[s] ** 2 * L[r]
            Y05 += -phi[s - 1] * phi[t] * phi[t - 1] * L[s] ** 2 * L[r]
            Y04 += phi[s] * phi[t] * phi[t - 1] * phi[(t - s) % q] * L[s] ** 2 * L[t]
            Y06 += phi[s] * phi[s - 1] * phi[t] * phi[t - 1] * phi[(t - s) % q] * L[s] ** 2 * L[t]
    return int(Y02), int(Y04), int(Y05), int(Y06)
for q in [p for p in primerange(7, 140)]:
    phi, L, inv = tabs(q)
    s = np.arange(2, q)
    A = {u: int((phi[s] * phi[(u - s) % q] * L[s] ** 2).sum()) for u in range(q)}
    A2 = {u: int((phi[s] * phi[s - 1] * phi[(u - s) % q] * L[s] ** 2).sum()) for u in range(q)}
    us = range(2, q)
    rY02 = -sum(phi[u] * L[u] * A[u] for u in us); rY04 = sum(phi[u] * phi[u - 1] * L[u] * A[u] for u in us)
    rY05 = sum(phi[u] * L[u] * A2[u] for u in us); rY06 = sum(phi[u] * phi[u - 1] * L[u] * A2[u] for u in us)
    Y02, Y04, Y05, Y06 = closed(q, phi, L, inv)
    red_ok = (rY02, rY04, rY05, rY06) == (Y02, Y04, Y05, Y06)
    Um = [u for u in us if phi[u - 1] == -1]
    inv_ok = all(phi[u * inv[u - 1] % q] == -phi[u] and L[u * inv[u - 1] % q] == L[u] and A[u * inv[u - 1] % q] == A[u]
                 and A2[u * inv[u - 1] % q] == A2[u] and phi[(u * inv[u - 1] - 1) % q] == -1 and u * inv[u - 1] % q != u for u in Um)
    out.append(f"q={q} ({q % 4} mod 4): reductions match closed forms: {red_ok}; involution properties on U-: {inv_ok}; "
               f"Y04+Y02 = {Y04 + Y02}, Y06-Y05 = {Y06 - Y05}")
print('\n'.join(out))
open(os.path.join(HERE, 'y_symmetry_proof_check.out'), 'w').write('\n'.join(out) + '\n')
