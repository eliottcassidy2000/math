"""Exact boundary-valued vertex-insertion update, independently checked by path cutting."""
def need(condition, message):
    if not condition:
        raise ArithmeticError(message)


from itertools import permutations, combinations
import json
from pathlib import Path
from native_kernel import observer, direct_observer, det, sectors, rref_kernel

def mk(A, k):
    n = len(A)
    D = {}
    if k > n:
        return D
    for p in permutations(range(n)):
        for cuts in combinations(range(1, n), k - 1):
            cuts = (0,) + cuts + (n,)
            if all((A[p[j]][p[j + 1]] for i in range(k) for j in range(cuts[i], cuts[i + 1] - 1))):
                s = tuple((p[cuts[i]] for i in range(k)))
                t = tuple((p[cuts[i + 1] - 1] for i in range(k)))
                D[s, t] = D.get((s, t), 0) + 1
    return D

def insertion_update(A, incoming):
    n = len(A)
    x = n
    old = {k: mk(A, k) for k in (1, 2, 3)}
    ports = list(permutations(range(n + 1), 2))

    def get(k, s, t):
        return old[k].get((tuple(s), tuple(t)), 0)

    def entry(s, t):
        if x in s and x in t:
            j = s.index(x)
            h = t.index(x)
            return get(1, [s[1 - j]], [t[1 - j]]) if j == h else 0
        if x in s:
            j = s.index(x)
            out = 0
            for v in range(n):
                ss = list(s)
                ss[j] = v
                out += (1 - incoming[v]) * get(2, ss, t)
            return out
        if x in t:
            j = t.index(x)
            out = 0
            for u in range(n):
                tt = list(t)
                tt[j] = u
                out += incoming[u] * get(2, s, tt)
            return out
        out = 0
        for u in range(n):
            if not incoming[u]:
                continue
            for v in range(n):
                if incoming[v]:
                    continue
                out += get(3, [s[0], v, s[1]], [u, t[0], t[1]])
                out += get(3, [s[0], s[1], v], [t[0], u, t[1]])
        return out
    return (ports, [[entry(s, t) for t in ports] for s in ports])

def tournament(n, bits):
    A = [[0] * n for _ in range(n)]
    for q, (i, j) in enumerate(combinations(range(n), 2)):
        A[i][j] = bits >> q & 1
        A[j][i] = 1 - A[i][j]
    return A

def main():
    counts = []
    for n in range(2, 5):
        cases = 0
        for bits in range(1 << n * (n - 1) // 2):
            A = tournament(n, bits)
            need(observer(A) == direct_observer(A), 'exact check failed')
            for mask in range(1 << n):
                inc = [mask >> u & 1 for u in range(n)]
                B = [row + [inc[u]] for u, row in enumerate(A)] + [[1 - v for v in inc] + [0]]
                need(insertion_update(A, inc) == direct_observer(B), (n, bits, mask))
                cases += 1
        counts.append({'old_order': n, 'all_labelled_extensions': cases})
    for bits in range(64):
        ports, M = observer(tournament(4, bits))
        P, S = sectors(ports, M)
        need(det(M) == det(P) * det(S), 'exact check failed')
        need(len(rref_kernel(M)) == len(rref_kernel(P)) + len(rref_kernel(S)), 'exact check failed')
    result = {'status': 'FINITE-EXACT verification of proved universal boundary update', 'extension_checks': counts, 'total_extensions': sum((c['all_labelled_extensions'] for c in counts)), 'order4_sector_and_kernel_checks': 64, 'insertion_word_criterion': 'OPEN: exact b/e semantics unavailable; no criterion assumed'}
    Path(__file__).with_name('receipt.json').write_text(json.dumps(result, indent=2) + '\n')
    print(json.dumps(result, indent=2))
if __name__ == '__main__':
    main()
