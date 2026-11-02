"""Exact referee; Python stdlib only. Run normally and with -O.

The path kernel and child HP DP are algorithmically distinct. Assertions are
deliberately not used, so optimized mode executes every verification gate.
"""
from collections import deque
import hashlib
import json
from probe import adj, h_all, kernel, cut, strong


def require(condition, message):
    if not condition:
        raise RuntimeError(message)


def child(a, s):
    n = len(a)
    return [row | (0 if s >> i & 1 else 1 << n)
            for i, row in enumerate(a)] + [s]


def minflow(capacity, source, sink):
    n = len(capacity)
    residual = [row[:] for row in capacity]
    value = 0
    while True:
        prev = [-1] * n
        prev[source] = source
        q = deque([source])
        while q and prev[sink] < 0:
            i = q.popleft()
            for j in range(n):
                if residual[i][j] and prev[j] < 0:
                    prev[j] = i
                    q.append(j)
        if prev[sink] < 0:
            return value
        delta = None
        j = sink
        while j != source:
            i = prev[j]
            delta = residual[i][j] if delta is None else min(delta, residual[i][j])
            j = i
        value += delta
        j = sink
        while j != source:
            i = prev[j]
            residual[i][j] -= delta
            residual[j][i] += delta
            j = i


def compatible(a, b, full):
    return not (a & b and a & (full ^ b) and b & (full ^ a)
                and full ^ (a | b))


def check_parent(a, direct=False, flows=False, terminal=True):
    n = len(a)
    full = (1 << n) - 1
    hs = h_all(a)
    H = hs[full]
    r = kernel(a, hs)
    f = [H + 2 * cut(r, s) for s in range(full + 1)]
    rho = min(hs[full ^ (1 << i) ^ (1 << j)]
              for i in range(n) for j in range(i + 1, n))
    require(rho % 2 == 1 and rho > 0, 'rho parity')
    for i in range(n):
        for j in range(i + 1, n):
            w = r[i][j] + r[j][i]
            require(w >= hs[full ^ (1 << i) ^ (1 << j)], 'pair floor')
            require(f[1 << i] + f[1 << j] - f[(1 << i) | (1 << j)] - H == 2*w,
                    'strict pair curvature')
    if direct:
        for s in range(full + 1):
            require(h_all(child(a, s))[-1] == f[s], ('literal child DP', a, s))
    m = min(f[1:full])
    M = max(f[1:full])
    lows = [s for s in range(1, full) if f[s] <= m + rho - 1]
    require(len(lows) <= 4*n - 6, ('low-state bound', a))
    for i, s in enumerate(lows):
        for t in lows[i+1:]:
            require(compatible(s, t, full), ('compatible low splits', a, s, t))
    if n >= 4:
        require(M - m >= rho + 1, ('response range', a))
    if flows:
        global_flow = min(minflow(r, 0, t) for t in range(1, n))
        global_flow = min(global_flow, *(minflow(r, t, 0) for t in range(1, n)))
        require(H + 2*global_flow == m, 'global minflow')
    if terminal:
        for i in range(n):
            for j in range(n):
                if i == j:
                    continue
                feasible = [s for s in range(1, full) if s >> i & 1 and not s >> j & 1]
                minimum = min(f[s] for s in feasible)
                chain = [s for s in feasible if f[s] <= minimum + rho - 1]
                require(len(chain) <= n-1, ('terminal low-chain size', a, i, j))
                chain.sort(key=int.bit_count)
                for s, t in zip(chain, chain[1:]):
                    require(s & t == s, ('terminal low chain', a, i, j))
    minima = [s for s in range(1, full) if f[s] == m]
    band_values = set(f[s] for s in lows)
    missing = [v for v in range(m, m + rho, 2) if v not in band_values]
    forced = max(0, (rho+1)//2 - (4*n-6))
    require(len(missing) >= forced, ('forced hole count', a))
    return dict(H=H, rho=rho, minimum=m, maximum=M, minima=minima,
                singleton_minimum=any(s.bit_count() in (1, n-1) for s in minima),
                low_states=len(lows), low_values=len(band_values),
                low_band_missing=len(missing), forced_missing=forced), r, f


def main():
    digest = hashlib.sha256()
    rows = []
    total_parents = total_direct_cuts = 0
    for n in range(2, 6):
        count = strong_count = missing_single = 0
        max_low = 0
        for code in range(1 << (n*(n-1)//2)):
            a = adj(n, code)
            result, r, f = check_parent(a, direct=True, flows=True)
            count += 1
            total_direct_cuts += 1 << n
            max_low = max(max_low, result['low_states'])
            is_strong = strong(a)
            strong_count += is_strong
            missing_single += is_strong and not result['singleton_minimum']
            digest.update(json.dumps([n, code, result, r, f], sort_keys=True,
                                     separators=(',', ':')).encode())
        total_parents += count
        row = dict(n=n, all_labelled=count, strong=strong_count,
                   no_singleton_minimum=missing_single, max_low_states=max_low)
        rows.append(row)
        print(json.dumps(row, sort_keys=True), flush=True)
    n = 6
    count = missing_single = 0
    first_hostile = None
    max_low = 0
    for code in range(1 << 15):
        a = adj(n, code)
        if not strong(a):
            continue
        result, r, f = check_parent(a, direct=False, flows=False, terminal=True)
        count += 1
        max_low = max(max_low, result['low_states'])
        missing_single += not result['singleton_minimum']
        if not result['singleton_minimum'] and first_hostile is None:
            # The first hostile gets all literal child counts and independent flows.
            checked, _, _ = check_parent(a, direct=True, flows=True)
            first_hostile = dict(code=code, adjacency=a, result=checked, R=r,
                                 singleton_values=[f[s] for s in range(1, 63)
                                                   if s.bit_count() in (1, 5)])
        digest.update(json.dumps([n, code, result], sort_keys=True,
                                 separators=(',', ':')).encode())
    row = dict(n=6, strong=count, no_singleton_minimum=missing_single,
               max_low_states=max_low)
    rows.append(row)
    print(json.dumps(row, sort_keys=True), flush=True)
    print(json.dumps({'first_hostile': first_hostile}, sort_keys=True), flush=True)
    for n in (7, 9, 11):
        a = [sum(1 << ((i+d) % n) for d in range(1, (n+1)//2)) for i in range(n)]
        result, r, f = check_parent(a, direct=(n <= 9), flows=True, terminal=True)
        # Literal independent DP replays the extremal children also at n=11.
        for value in (result['minimum'], result['maximum']):
            s = f.index(value)
            require(h_all(child(a, s))[-1] == value, ('cyclic extremal child', n))
        print(json.dumps({'cyclic_order': n, 'result': result}, sort_keys=True), flush=True)
        digest.update(json.dumps([n, result, r, f], sort_keys=True,
                                 separators=(',', ':')).encode())
    print(json.dumps({'status': 'PASS', 'all_labelled_parents_through5': total_parents,
                      'literal_child_DP_checks_through5': total_direct_cuts,
                      'strong_labelled_order6': count,
                      'semantic_sha256': digest.hexdigest()}, sort_keys=True))


if __name__ == '__main__':
    main()
