#!/usr/bin/env python3
"""Independent named one-auxiliary native compiler; no repository imports.
Finite universe: every labelled tournament and arc at orders 2 through 5.
"""
from itertools import combinations, permutations
from collections import Counter
from functools import lru_cache
import json


def need(p, message):
    if not p:
        raise RuntimeError(message)


def decode(n, mask):
    rows = [0] * n
    for k, (i, j) in enumerate(combinations(range(n), 2)):
        u, v = (i, j) if mask >> k & 1 else (j, i)
        rows[u] |= 1 << v
    return tuple(rows)


def flip(T, a, b):
    U = list(T)
    U[a] ^= 1 << b
    U[b] ^= 1 << a
    return tuple(U)


def legal(T, p):
    return sorted(p) == list(range(len(T))) and all(T[a] >> b & 1 for a, b in zip(p, p[1:]))


@lru_cache(None)
def paths(T):
    return tuple(p for p in permutations(range(len(T))) if legal(T, p))


def native(T, p):
    need(legal(T, p), 'native guard')
    kept = {frozenset((a, b)) for a, b in zip(p, p[1:])}
    U = T
    for a, b in combinations(range(len(T)), 2):
        if frozenset((a, b)) not in kept:
            U = flip(U, a, b)
    need(legal(U, p), 'guard survives native move')
    return U


def extend(T, out):
    """Temporary vertex z beats exactly old named vertices in out."""
    n = len(T)
    return tuple(r | (0 if i in out else 1 << n) for i, r in enumerate(T)) + (sum(1 << i for i in out),)


def restrict(T, n):
    return tuple(T[i] & ((1 << n) - 1) for i in range(n))


def traced_gadget(T, u, v, p):
    """Compile a traced arc, returning the complete named intermediate states."""
    n = len(T)
    need((u, v) in tuple(zip(p, p[1:])), 'chosen traced arc')
    k = p.index(u)
    if n == 2:
        U = extend(T, set())
        p1, q, case = p + (n,), (v, n, u), 'n2'
    elif v != p[-1]:
        U = extend(T, {u})
        p1, q, case = p + (n,), p[:k+1] + (n,) + p[k+1:], 'append_single_out'
    else:
        U = extend(T, set(range(n)) - {v})
        p1, q, case = (n,) + p, p[:k+1] + (n,) + p[k+1:], 'prepend_single_in'
    V = native(U, p1)
    W = native(V, q)
    need(restrict(W, n) == flip(T, u, v), 'only the requested old arc changes')
    if n > 2:
        if case == 'append_single_out':
            need(W[n] == (1 << v) | (1 << p[-1]), 'final two out-ports')
        else:
            need(W[n] == ((1 << n)-1) ^ (1 << u) ^ (1 << p[0]), 'final two in-ports')
    return U, V, W, p1, q, case


def compile_arc(T, a, b):
    u, v = (a, b) if T[a] >> b & 1 else (b, a)
    p = next((p for p in paths(T) if (u, v) in tuple(zip(p, p[1:]))), None)
    if p is not None:
        return traced_gadget(T, u, v, p) + (False,)
    S = flip(T, u, v)
    p = next((p for p in paths(S) if (v, u) in tuple(zip(p, p[1:]))), None)
    need(p is not None, 'inherited one-guard separation')
    U, V, W, p1, q, case = traced_gadget(S, v, u, p)
    need(restrict(W, len(T)) == T, 'inverse gadget source')
    need(native(W, q) == V and native(V, p1) == U, 'inverse legal native moves')
    return W, V, U, q, p1, case, True


def main():
    checks = Counter()
    examples = {}
    layers = []
    for n in range(2, 6):
        arcs = 0
        for mask in range(1 << (n*(n-1)//2)):
            T = decode(n, mask)
            for a, b in combinations(range(n), 2):
                U, V, W, p, q, case, reverse = compile_arc(T, a, b)
                need(restrict(U, n) == T, 'source retains every old arc and name')
                need(native(U, p) == V and native(V, q) == W, 'actual two-step replay')
                need(restrict(W, n) == flip(T, a, b), 'exact target')
                need(len(U) == n+1, 'peak clock')
                # At a fixed order, two native moves change an even number of arcs.
                differences = sum(((U[i] >> j) ^ (W[i] >> j)) & 1 for i, j in combinations(range(n+1), 2))
                need(differences % 2 == 0, 'two-native parity')
                key = case + ('_inverse' if reverse else '_forward')
                checks[key] += 1
                examples.setdefault(key, {'n':n, 'old_arc':[a,b], 'source_rows':T, 'birth_rows':U,
                    'intermediate_rows':V, 'final_rows':W, 'paths':[p,q]})
                arcs += 1
        layers.append({'n':n, 'labelled_sources':1 << (n*(n-1)//2), 'arc_compilations':arcs})
    # Literal adaptive-family exit: R4 has the spine as a marked path, and its
    # native partner is transitive. Names and the exact path are retained.
    TT4 = tuple(sum(1 << j for j in range(i+1,4)) for i in range(4))
    spine = tuple(range(4))
    R4 = native(TT4, spine)
    need(len(paths(TT4)) == 1 and len(paths(R4)) == 5, 'adaptive exit H controls')
    need(native(R4, spine) == TT4, 'native exit replay')
    return {'status':'FINITE-EXACT', 'universe':layers, 'branches':dict(checks), 'examples':examples,
        'adaptive_native_exit':{'adaptive_rows':R4, 'target_rows':TT4, 'path':spine, 'H':[5,1]},
        'scope':'One birth, two native moves, one named deletion; not fixed-order native connectivity.'}


if __name__ == '__main__':
    print(json.dumps(main(), indent=2, sort_keys=True))
