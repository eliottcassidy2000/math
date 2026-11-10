#!/usr/bin/env python3
# Checks intentionally survive python -O.
import sys
if "--fail" in sys.argv:
    raise RuntimeError("deliberate checker failure")
"""Literal permutations, independent of probe.py DP and transfer evaluation."""
import itertools, json, math
from pathlib import Path

def paths(a, V):
    return [p for p in itertools.permutations(V) if all((a[x] >> y & 1 for x, y in zip(p, p[1:])))]

def H(a, V):
    return len(paths(a, V))

def z(a, u):
    V = set(range(len(a)))
    R = 0
    for k in range(1, len(a) + 1):
        for tail in itertools.permutations(sorted(V - {u}), k - 1):
            p = (u,) + tail
            if all((a[x] >> y & 1 for x, y in zip(p, p[1:]))):
                R += H(a, sorted(V - set(p)))
    p = paths(a, sorted(V))
    return (len(p), R, sum((w[-1] == u for w in p)), H(a, sorted(V - {u})))

def aut(a):
    return [p for p in itertools.permutations(range(len(a))) if all((a[i] >> j & 1 == a[p[i]] >> p[j] & 1 for i in range(len(a)) for j in range(len(a))))]

def strong(a):
    for i in range(len(a)):
        seen = 1 << i
        for _ in a:
            seen |= sum((1 << j for j in range(len(a)) if any((seen >> k & 1 and a[k] >> j & 1 for k in range(len(a))))))
        if seen != (1 << len(a)) - 1:
            return False
    return True
T = [24, 1, 19, 22, 2]
U = [58, 60, 57, 16, 32, 8]
if not (len(aut(T)) == 1 and strong(T)):
    raise RuntimeError('live mathematical check failed')
if not z(T, 1) == z(U, 3) == (9, 9, 3, 3):
    raise RuntimeError('live mathematical check failed')
if not z(T, 0) == (9, 11, 1, 1):
    raise RuntimeError('live mathematical check failed')
if not all((math.gcd(*z(U, v)) == 3 for v in range(6))):
    raise RuntimeError('live mathematical check failed')
E = [r | 32 for r in T] + [1]
E[0] ^= 32
if not H(E, list(range(6))) == 19:
    raise RuntimeError('live mathematical check failed')
R = [sum((1 << j for j in range(5) if T[j] >> i & 1)) for i in range(5)]
if not z(R, 1) == (9, 21, 1, 3):
    raise RuntimeError('live mathematical check failed')
Urev = [sum((1 << j for j in range(6) if U[j] >> i & 1)) for i in range(6)]
if not z(Urev, 3) == (9, 36, 0, 3):
    raise RuntimeError('live mathematical check failed')
X = [62] + [r << 1 for r in T]
Xrev = [sum((1 << j for j in range(6) if X[j] >> i & 1)) for i in range(6)]
if not (z(X, 2) == z(U, 3) and z(X, 1) == (9, 11, 1, 1)):
    raise RuntimeError('live mathematical check failed')
if not (z(Xrev, 2) == (9, 42, 0, 3) and len(aut(X)) == 1):
    raise RuntimeError('live mathematical check failed')
checks = 0
for l in range(3):
    for word in itertools.product([0, 1], repeat=l):
        values = []
        for a, u in [(T, 1), (U, 3)]:
            a = list(a)
            for b in word:
                n = len(a)
                a = [r | 1 << n for r in a] + [0]
                if b:
                    a[u] ^= 1 << n
                    a[n] |= 1 << u
            values.append(H(a, list(range(len(a)))))
        if not (values[0] == values[1] and values[0] % 3 == 0):
            raise RuntimeError('live mathematical check failed')
        checks += 1
out = {'T_rows': T, 'T_full_automorphism_order': len(aut(T)), 'T_strong': strong(T), 'T_root_states': [z(T, v) for v in range(5)], 'same_order_X_rows': X, 'X_root2': z(X, 2), 'X_converse_root2': z(Xrev, 2), 'U_rows': U, 'U_full_automorphism_order': len(aut(U)), 'U_root_states': [z(U, v) for v in range(6)], 'converse_T_rows': R, 'converse_root1': z(R, 1), 'U_converse_root3': z(Urev, 3), 'escape_H': 19, 'literal_equal_tail_pairs': checks}
print(json.dumps(out, indent=2))
