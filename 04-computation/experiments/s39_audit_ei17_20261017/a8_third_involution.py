#!/usr/bin/env python3
"""Audit of THM-4665, Open section: the 'third involution'.

Claim checked (sympy, then on every faithful solution found at several primes):
  v1 + v14 = v2 + v11 identically (VERT), so with m the common midpoint the map
      iota: (U0, U4) -> (U9, U11) = (d1 - U0, d2 - U4),  d1 = v14 - v1,  d2 = v11 - v2,  other free directions fixed,
  re-places v4 -> 2m - v4 and v6 -> 2m - v6 (and v13 = v4 + U2, v16 = v6 + 1 accordingly); it fixes U7, U10, U12,
  swaps U0 <-> U9 and U4 <-> U11, and sends U8 -> -U8, so it maps solutions of the 13 unit conditions to solutions.
  It has no fixed point at all: iota(s) = s forces v4 = m = v6, i.e. U8 = 0, not a unit.
Then: on every faithful solution at the listed primes, iota(s) is the unique other faithful solution sharing the 13
vertices outside {4, 6, 13, 16}; and which element of {M, A, MA} has fixed points at p = 383 (twisted, count 4).
"""
import os
os.environ.setdefault('NUMBA_CACHE_DIR', os.path.join(os.environ.get('TEMP', os.environ.get('TMPDIR', '/tmp')), 'nbc_s39_audit_ei17'))  # short path outside the repo (Windows MAX_PATH)
import sys

import sympy as sp

HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, HERE)
import a2_embedding_search as A2  # noqa: E402
import a6_verify_embeddings as A6  # noqa: E402

BASIS = ['1', 'U0', 'U1', 'U2', 'U3', 'U4', 'U6']
U = {b: (sp.Integer(1) if b == '1' else sp.Symbol(b)) for b in BASIS}
vt = {v: sum(A2.VERT[v][t] * U[BASIS[t]] for t in range(7)) for v in A2.VERT}
print('v1 + v14 - v2 - v11 =', sp.expand(vt[1] + vt[14] - vt[2] - vt[11]))
d1 = sp.expand(vt[14] - vt[1]); d2 = sp.expand(vt[11] - vt[2])
U0, U1, U2, U3, U4, U6 = [U[b] for b in BASIS[1:]]
F = {'U7': U1 + U2 - U3 - 1, 'U8': -U0 - U2 + U4 + 1, 'U9': -U0 + U1 - U2 + U3 + U6 + 1, 'U10': 2 * U2 - U3 - U6 - 1,
     'U11': U1 + U2 + U3 - U4 + U6 - 1, 'U12': -U1 - U2 - U6 + 2}
sub = {U0: d1 - U0, U4: d2 - U4}
print('d1 = v14 - v1 =', d1, '  equals U0 + U9:', sp.expand(U0 + F['U9'] - d1) == 0)
print('d2 = v11 - v2 =', d2, '  equals U4 + U11:', sp.expand(U4 + F['U11'] - d2) == 0)
for k, e in F.items():
    print('  iota(%-3s) = %s' % (k, sp.factor(sp.expand(e.subs(sub, simultaneous=True)))))
for v in (4, 6, 13, 16):
    print('  iota(v%d) = %s ; 2m - v%d = %s' % (v, sp.expand(vt[v].subs(sub, simultaneous=True)), v,
                                              sp.expand(vt[1] + vt[14] - vt[v])))
print('  other vertices fixed:', all(sp.expand(vt[v].subs(sub, simultaneous=True) - vt[v]) == 0
                                     for v in vt if v not in (4, 6, 13, 16)))

ST, _ = A2.pick_order(); ARR = A2.order_arrays(ST)


def check(model, p):
    b, sols = A2.solutions(model, p, ARR)
    S = {tuple(tuple(int(c) for c in s[v - 1]) for v in range(1, 18)) for s in sols}
    n_ok = 0
    for T in S:
        P = {v: T[v - 1] for v in range(1, 18)}
        Pi = dict(P)
        for v in (4, 6):
            Pi[v] = ((P[1][0] + P[14][0] - P[v][0]) % p, (P[1][1] + P[14][1] - P[v][1]) % p)
        U2v = ((P[17][0] - P[14][0]) % p, (P[17][1] - P[14][1]) % p)
        Pi[13] = ((Pi[4][0] + U2v[0]) % p, (Pi[4][1] + U2v[1]) % p)
        one = (1, 1) if model == 0 else (1, 0)
        Pi[16] = ((Pi[6][0] + one[0]) % p, (Pi[6][1] + one[1]) % p)
        Ti = tuple(Pi[v] for v in range(1, 18))
        n_ok += (Ti in S and Ti != T)
    return len(S), n_ok


for model, plist in ((0, (103, 139, 281, 313, 401, 677, 2861, 3191)), (1, (127, 149, 227, 383, 2441, 2551))):
    for p in plist:
        n, k = check(model, p)
        print('%s p = %5d: %2d faithful solutions, iota(s) faithful, distinct and in the set for %2d of them'
              % ('split  ' if model == 0 else 'twisted', p, n, k))

# p = 383 twisted: which element of the Klein group fixes solutions
p = 383
b, sols = A2.solutions(1, p, ARR)
F2 = A6.Fq2(p, b)
S = [tuple(tuple(int(c) for c in s[v - 1]) for v in range(1, 18)) for s in sols]
conj = lambda z: ((z[0] + z[1]) % p, (-z[1]) % p)
mirror = lambda T: tuple(conj(z) for z in T)


def autom(T):
    P = {v: T[v - 1] for v in range(1, 18)}
    r = F2.mul((p - 1, 0), F2.pow(P[14], p * p - 2))
    return tuple(F2.mul(r, P[A6.AUT[v]]) for v in range(1, 18))


for T in S:
    print('p = 383 solution: fixed by M %s, by A %s, by MA %s' % (mirror(T) == T, autom(T) == T, mirror(autom(T)) == T))
