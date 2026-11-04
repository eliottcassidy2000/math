import json, subprocess, itertools, functools, pathlib
D = pathlib.Path(__file__).resolve().parent
import tempfile
_tmp = tempfile.TemporaryDirectory(prefix='native-independent-')
EXE = pathlib.Path(_tmp.name) / 'count'
subprocess.run(['g++', '-O3', '-std=c++20', str(D / 'count.cpp'), '-o', str(EXE)], check=True)
source = json.load(open(D.parent / 'returns/source_interface.json'))
z = json.load(open(D.parent / 'returns/candidate.json'))
names = source['target_names']
S = tuple(source['target'])
N = len(S)

def strong(T):
    for i in range(len(T)):
        seen = 1 << i
        while True:
            q = seen
            for j in range(len(T)):
                if seen >> j & 1:
                    q |= T[j]
            if q == seen:
                break
            seen = q
        if seen != (1 << len(T)) - 1:
            return False
    return True

def flip(T, a, b):
    U = list(T)
    U[a] ^= 1 << b
    U[b] ^= 1 << a
    return tuple(U)

def birth(T, cut):
    n = len(T)
    return tuple((x + ((cut >> i & 1) << n) for i, x in enumerate(T))) + ((1 << n) - 1 ^ cut,)

def native(T, P):
    if not sorted(P) == list(range(len(T))):
        raise RuntimeError('Audit condition failed: sorted(P) == list(range(len(T)))')
    if not all((T[a] >> b & 1 for a, b in zip(P, P[1:]))):
        raise RuntimeError('Audit condition failed: all((T[a] >> b & 1 for a, b in zip(P, P[1:])))')
    U = [(1 << len(T)) - 1 ^ x ^ 1 << i for i, x in enumerate(T)]
    for a, b in zip(P, P[1:]):
        U[a] ^= 1 << b
        U[b] ^= 1 << a
    return tuple(U)

def guard(T, a, b):
    n = len(T)
    full = (1 << n) - 1

    @functools.lru_cache(None)
    def dfs(mask, last):
        if mask == full:
            return () if last != b else None
        avail = T[last] & ~mask if last >= 0 else full
        if last == a:
            avail &= 1 << b
        else:
            avail &= ~(1 << b)
        while avail:
            bit = avail & -avail
            avail -= bit
            v = bit.bit_length() - 1
            tail = dfs(mask | bit, v)
            if tail is not None:
                return (v,) + tail
        return None
    return dfs(0, -1)
changes = [tuple((names.index(a) for a in e)) for e in z['old_edge_flips']]
for order in itertools.permutations(changes):
    T = S
    trace = []
    for k, (a, b) in enumerate(order):
        U = flip(T, a, b)
        if not strong(U):
            break
        if not T[a] >> b & 1:
            a, b = (b, a)
        P = guard(T, a, b)
        if P is None:
            break
        i = P.index(a)
        cut = (1 << N) - 1 ^ 1 << a
        X = birth(T, cut)
        P1 = (*P, N)
        P2 = (*P[:i + 1], N, *P[i + 1:])
        V = native(native(X, P1), P2)
        result = tuple((x & (1 << N) - 1 for x in V[:N]))
        if not result == U:
            raise RuntimeError('Audit condition failed: result == U')
        local = names + [100 + k]
        trace.append(dict(arc=[names[a], names[b]], helper=100 + k, cut=cut, paths=[[local[v] for v in p] for p in (P1, P2)], rows=list(U), strong=True))
        T = U
    else:
        break
else:
    raise RuntimeError('No guarded strong schedule')
births = []
for i in range(N, len(z['rows'])):
    cut = sum((1 << j for j in range(i) if z['rows'][j] >> i & 1))
    T = birth(T, cut)
    if not strong(T):
        raise RuntimeError('Audit condition failed: strong(T)')
    births.append(dict(name=z['names'][i], cut=cut, rows=T))
if not list(T) == z['rows']:
    raise RuntimeError("Audit condition failed: list(T) == z['rows']")
if not len(set(z['names'])) == len(z['names']):
    raise RuntimeError("Audit condition failed: len(set(z['names'])) == len(z['names'])")
if not set(names) <= set(z['names']):
    raise RuntimeError("Audit condition failed: set(names) <= set(z['names'])")
counts = list(map(int, subprocess.check_output([str(EXE)], input=''.join((str(len(T)) + ' ' + ' '.join(map(str, T)) + '\n' for T in [S, *[t['rows'] for t in trace], *[b['rows'] for b in births]])), text=True).split()))
if counts[0] != z['source_H'] or counts[-1] != z['target_H']:
    raise RuntimeError('Count target mismatch')
result = dict(status='INDEPENDENT_EXACT_COUNT_AND_NATIVE_REPLAY', source_names=names, source=S, trace=trace, births=births, target_names=z['names'], target=T, counts=counts, native_complements=2 * len(trace), temporary_helpers=len(trace), deleted_old_names=[], source_H=counts[0], target_H=counts[-1], all_completed_states_strong=True, scope='Literal named program and independent subset DP. Root39 and spare6 stars change; no preserved-frame claim.')
(D / 'native_audit.json').write_text(json.dumps(result, indent=2) + '\n')
print('PASS', counts, 'arcs', [t['arc'] for t in trace], flush=True)
