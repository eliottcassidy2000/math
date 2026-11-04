#!/usr/bin/env python3
"""Exact all-strong three-helper compiler controls, standalone.

Every strong labelled source/target pair through4 and every scaffold root;
every labelled old core through5 and old pair inside the scaffold; setup and
cleanup at every root of every strong old core through5. No producer import.
"""
from itertools import combinations, permutations
from pathlib import Path
import argparse
import hashlib
import json

GATES = 0
GADGETS = 0
INVERSE = 0


def need(ok, label):
    global GATES
    GATES += 1
    if not ok:
        raise RuntimeError(label)


def graph(n, code):
    T = [0] * n
    for i, (u, v) in enumerate(combinations(range(n), 2)):
        if code >> i & 1:
            T[u] |= 1 << v
        else:
            T[v] |= 1 << u
    return tuple(T)


def flip(T, u, v):
    out = list(T)
    out[u] ^= 1 << v
    out[v] ^= 1 << u
    return tuple(out)


def strong(T):
    R = [r | (1 << i) for i, r in enumerate(T)]
    for k in range(len(T)):
        for i in range(len(T)):
            if R[i] >> k & 1:
                R[i] |= R[k]
    return all(r == (1 << len(T)) - 1 for r in R)


def birth(T, out):
    n = len(T)
    return tuple(row | (0 if out >> i & 1 else 1 << n) for i, row in enumerate(T)) + (out,)


def delete_last(T):
    mask = (1 << (len(T) - 1)) - 1
    return tuple(row & mask for row in T[:-1])


def directed(T, P):
    return len(P) == len(T) and len(set(P)) == len(T) and set(P) == set(range(len(T))) and all(T[u] >> v & 1 for u, v in zip(P, P[1:]))


def native(T, P):
    need(directed(T, P), 'actual full native guard')
    protected = {frozenset((u, v)) for u, v in zip(P, P[1:])}
    out = list(T)
    for u, v in combinations(range(len(T)), 2):
        if frozenset((u, v)) not in protected:
            out[u] ^= 1 << v
            out[v] ^= 1 << u
    return tuple(out)


def native_sets(T, P):
    """Independent adjacency-set complement interpreter."""
    edges = {(u, v) for u, row in enumerate(T) for v in range(len(T)) if row >> v & 1}
    if len(P) != len(T) or len(set(P)) != len(T) or set(P) != set(range(len(T))) or any((u, v) not in edges for u, v in zip(P, P[1:])):
        raise ValueError('illegal native guard')
    keep = {frozenset((u, v)) for u, v in zip(P, P[1:])}
    result = {(u, v) if frozenset((u, v)) in keep else (v, u) for u, v in edges}
    return tuple(sum(1 << v for u0, v in result if u0 == u) for u in range(len(T)))


def pair_word(T, u, v):
    """Actual direct pair-block algorithm, applied only to the old core."""
    P = []
    for a in range(len(T)):
        if a in (u, v):
            continue
        at = next((i for i, b in enumerate(P) if T[a] >> b & 1), len(P))
        P.insert(at, a)
    at = next((i for i, b in enumerate(P) if (T[u] | T[v]) >> b & 1), len(P))
    block = (u, v) if at == len(P) or T[v] >> P[at] & 1 else (v, u)
    P[at:at] = block
    oriented = T if T[block[0]] >> block[1] & 1 else flip(T, u, v)
    need(directed(oriented, tuple(P)), 'old-core pair block in one orientation')
    return tuple(P), oriented, block


def hp_start(T, w):
    for tail in permutations(v for v in range(len(T)) if v != w):
        P = (w,) + tail
        if directed(T, P):
            return P
    raise RuntimeError('strong old core missing prescribed starting path')


def scaffold(T):
    n = len(T)
    return tuple(row | (1 << (n + 1)) for row in T) + ((1 << n) - 1, 1 << n)


def push(log, kind, T, guard=None, out=None):
    if log is not None:
        rec = {'command': kind, 'rows': T}
        if guard is not None:
            rec['guard'] = guard
        if out is not None:
            rec['helper_out'] = out
        log.append(rec)


def gadget(current, base, P, u, v, log=None):
    """Append gadget, possibly reversed; both endpoint restrictions strong."""
    global GADGETS, INVERSE
    N = len(current)
    at = P.index(u)
    need(at > 0 and at + 1 < N - 1 and P[at + 1] == v, 'selected pair strictly internal')
    need(directed(base, P) and strong(base) and strong(current), 'strong endpoint bases and oriented path')
    need(current in (base, flip(base, u, v)), 'actual forward or inverse base')
    A = birth(base, 1 << u)
    first = P + (N,)
    second = P[:at + 1] + (N,) + P[at + 1:]
    middle = native(A, first)
    Z = native(middle, second)
    need(delete_last(Z) == flip(base, u, v), 'single-pair endpoint identity')
    need(Z[-1] == (1 << v) | (1 << P[-1]), 'exact reverse helper cut')
    need(directed(middle, first) and middle[N] >> P[0] & 1, 'middle has literal Hamilton cycle Pz')
    inverse = current != base
    if inverse:
        A, Z = Z, A
        first, second = second, first
    need(delete_last(A) == current, 'helper birth preserves actual old relation')
    need(0 < A[-1] < (1 << N) - 1 and A[-1].bit_count() <= 2, 'actual helper birth proper sparse cut')
    need(native_sets(A, first) == middle and native_sets(middle, second) == Z, 'independent actual native replay')
    after = delete_last(Z)
    need(after == flip(current, u, v), 'actual forward/inverse named endpoint')
    need(all(strong(S) for S in (A, middle, Z, after)), 'every helper/native/deletion state strong')
    push(log, 'helper_birth', A, out=A[-1])
    push(log, 'native_complement', middle, guard=first)
    push(log, 'native_complement', Z, guard=second)
    push(log, 'helper_delete', after)
    GADGETS += 1
    INVERSE += inverse
    return after


def setup(T, w, log=None):
    n = len(T)
    A = birth(T, ((1 << n) - 1) ^ (1 << w))
    B = birth(A, 1 << n)
    need(strong(T) and strong(A) and strong(B), 'two scaffold setup births strong')
    push(log, 'scaffold_a_birth', A, out=A[-1])
    push(log, 'scaffold_b_birth', B, out=B[-1])
    S = scaffold(T)
    P = (n + 1, n) + hp_start(T, w)
    out = gadget(B, S, P, n, w, log)
    need(out == S, 'complete scaffold installed')
    return out


def cleanup(S, T, w, log=None):
    n = len(T)
    need(S == scaffold(T) and strong(T), 'cleanup has actual strong target core')
    P = (n + 1, n) + hp_start(T, w)
    out = gadget(S, S, P, n, w, log)
    out = delete_last(out)
    need(strong(out), 'delete b leaves target plus a strong')
    push(log, 'scaffold_b_delete', out)
    out = delete_last(out)
    need(out == T and strong(out), 'delete a leaves exact strong target')
    push(log, 'scaffold_a_delete', out)
    return out


def old_pair(S, u, v, log=None):
    n = len(S) - 2
    old = tuple(row & ((1 << n) - 1) for row in S[:n])
    R, oriented, (a, b) = pair_word(old, u, v)
    base = scaffold(oriented)
    return gadget(S, base, (n,) + R + (n + 1,), a, b, log)


def compile(T, U, w, trace=False):
    log = [] if trace else None
    push(log, 'initial', T)
    before = GADGETS
    S = setup(T, w, log)
    d = 0
    for a, b in combinations(range(len(T)), 2):
        if bool(T[a] >> b & 1) != bool(U[a] >> b & 1):
            S = old_pair(S, a, b, log)
            d += 1
    need(S == scaffold(U), 'all old discrepancies corrected inside scaffold')
    result = cleanup(S, U, w, log)
    need(result == U and GADGETS - before == d + 2, 'complete named compiler and gadget budget')
    if log is not None:
        births = sum(r['command'].endswith('_birth') for r in log)
        deletes = sum(r['command'].endswith('_delete') for r in log)
        complements = sum(r['command'] == 'native_complement' for r in log)
        need((births, deletes, complements) == (d + 4, d + 4, 2*d + 4), 'literal exact event budget')
        need(max(len(r['rows']) for r in log) == len(T) + 3, 'literal peak order')
    return d, log


def observer(T, w):
    n, size = len(T), 1 << len(T)
    D = [[0] * n for _ in range(size)]
    H = [0] * size
    H[0] = 1
    for v in range(n):
        D[1 << v][v] = 1
    for mask in range(1, size):
        for v in range(n):
            if mask >> v & 1 and mask != 1 << v:
                old = mask ^ (1 << v)
                D[mask][v] = sum(D[old][u] for u in range(n) if old >> u & 1 and T[u] >> v & 1)
        H[mask] = sum(D[mask])
    full = size - 1
    F = sum(D[s][w] * H[full ^ s] for s in range(size))
    return {'H':H[full], 'F_root':F, 'defect':F-H[full]}


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument('--output', type=Path, default=Path(__file__).with_suffix('.json'))
    args = ap.parse_args()
    banks = {n:[graph(n,c) for c in range(1 << (n*(n-1)//2))] for n in (3,4,5)}
    strong_banks = {n:[T for T in bank if strong(T)] for n,bank in banks.items()}
    need([len(strong_banks[n]) for n in (3,4,5)] == [2,24,544], 'complete labelled strong census')
    pair_cases = setup_cases = full_cases = 0
    for n, bank in banks.items():
        for T in bank:
            S = scaffold(T)
            need(strong(S), 'scaffold strong for arbitrary weak or strong old core')
            for a,b in combinations(range(n),2):
                need(old_pair(S,a,b) == scaffold(flip(T,a,b)), 'complete arbitrary-core single-pair scaffold')
                pair_cases += 1
        for T in strong_banks[n]:
            for w in range(n):
                need(cleanup(setup(T,w),T,w) == T, 'every strong core/root setup cleanup')
                setup_cases += 1
    for n in (3,4):
        for T in strong_banks[n]:
            for U in strong_banks[n]:
                for w in range(n):
                    compile(T,U,w)
                    full_cases += 1
    need((pair_cases,setup_cases,full_cases)==(10648,2822,2316), 'complete declared universes')
    T,U,w=(8,1,3,6),(2,4,9,3),0
    d,trace=compile(T,U,w,True)
    for r in trace:
        r['observer'] = observer(r['rows'],w)
    need(trace[0]['observer']['H']==trace[-1]['observer']['H']==5, 'H-neutral endpoint control')
    need(any(r['observer']['H']!=5 for r in trace), 'intermediate H is not promised invariant')
    # Cleanup in the opposite deletion order is a genuine weak-state hostile.
    A=birth(U,((1<<len(U))-1)^(1<<w));B=birth(A,1<<len(U))
    keep=tuple(range(len(U)))+(len(U)+1,)
    wrong=tuple(sum(((B[a]>>b)&1)<<j for j,b in enumerate(keep)) for a in keep)
    need(not strong(wrong), 'must delete b before a at cleanup')
    premature=delete_last(scaffold(U))
    need(not strong(premature) and premature[-1] == (1 << len(U)) - 1,
         'deleting b before restoring w-to-a leaves universal source a')
    result={'scope':'finite controls for universal all-strong compiler; every old label survives; no count invariance or optimal helper claim',
            'gates':GATES,'gadgets':GADGETS,'inverse_gadgets':INVERSE,
            'arbitrary_core_pair_cases':pair_cases,'strong_setup_cleanup_cases':setup_cases,'all_strong_endpoint_root_cases':full_cases,
            'example':{'source':T,'target':U,'retained_root':w,'d':d,'trace':trace},
            'cleanup_wrong_order_hostile':wrong,'cleanup_premature_b_hostile':premature}
    raw=(json.dumps(result,indent=2,sort_keys=True)+'\n').encode();args.output.write_bytes(raw)
    print(json.dumps({'gates':GATES,'bytes':len(raw),'sha256':hashlib.sha256(raw).hexdigest()},sort_keys=True))


if __name__=='__main__':
    main()
