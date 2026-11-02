"""Clean-room symmetry referee: imports no primary code and no path kernel.

Direct subset DP reconstructs the ear quadratic from singleton/pair child
values. Exhaustive literal children additionally check that reconstruction
for cyclic orders 5,7,9; order11 checks each extremum and all low-band states.
"""
import hashlib
import json


def check(ok, data):
    if not ok:
        raise ValueError(data)


def path_count(a):
    n = len(a)
    if n == 0:
        return 1
    table = [dict() for _ in range(1 << n)]
    for v in range(n):
        table[1 << v][v] = 1
    for mask, endpoints in enumerate(table):
        for end, multiplicity in endpoints.items():
            available = a[end] & ~mask
            while available:
                bit = available & -available
                available ^= bit
                nxt = bit.bit_length() - 1
                target = table[mask | bit]
                target[nxt] = target.get(nxt, 0) + multiplicity
    return sum(table[-1].values())


def restrict(a, vertices):
    return [sum(1 << j for j, v in enumerate(vertices) if a[u] >> v & 1)
            for u in vertices]


def ear(a, s):
    n = len(a)
    return [out | ((1 << n) if not (s >> i & 1) else 0)
            for i, out in enumerate(a)] + [s]


def rotate(s, n):
    return ((s << 1) & ((1 << n)-1)) | (s >> (n-1))


def floor_log3(n):
    k, power = 0, 1
    while power*3 <= n:
        power *= 3
        k += 1
    return k


def prime(n):
    return all(n % d for d in range(2, n) if d*d <= n)


def main():
    semantic = []
    for n in (5, 7, 9, 11):
        full = (1 << n)-1
        a = [sum(1 << ((i+d) % n) for d in range(1, (n+1)//2))
             for i in range(n)]
        H = path_count(a)
        singleton = [path_count(ear(a, 1 << i)) for i in range(n)]
        A = [v-H for v in singleton]
        K = [[0]*n for _ in range(n)]
        deletions = []
        for i in range(n):
            for j in range(i+1, n):
                f_pair = path_count(ear(a, (1 << i) | (1 << j)))
                K[i][j] = K[j][i] = singleton[i]+singleton[j]-f_pair-H
                deletion = path_count(restrict(a, [v for v in range(n) if v not in (i,j)]))
                deletions.append(deletion)
                check(K[i][j] >= 2*deletion, ('pair floor', n, i, j))
        rho = min(deletions)
        f = []
        for s in range(full+1):
            members = [i for i in range(n) if s >> i & 1]
            value = H + sum(A[i] for i in members)
            value -= sum(K[i][j] for i in members for j in members if i < j)
            f.append(value)
            if n <= 9:
                check(path_count(ear(a, s)) == value, ('literal child', n, s))
        m, M = min(f[1:full]), max(f[1:full])
        lows = {s for s in range(1, full) if f[s] <= m+rho-1}
        small = {s if s.bit_count() <= n//2 else full ^ s for s in lows}
        for s in lows:
            check(f[s] == f[full ^ s], ('complement equality', n, s))
            check(rotate(s, n) in lows, ('rotation invariance', n, s))
        for s in small:
            t = s
            orbit = set()
            for _ in range(n):
                orbit.add(t)
                t = rotate(t, n)
            check(sum(t.bit_count() for t in orbit) == n, ('block partition size', n, s))
            covered = 0
            for t in orbit:
                check(not (covered & t), ('disjoint orbit', n, s))
                covered |= t
            check(covered == full, ('orbit covers', n, s))
        sizes = sorted({s.bit_count() for s in small})
        check(all(n % d == 0 for d in sizes), ('block divisor', n))
        check(all(e % d == 0 and e//d >= 3 for d, e in zip(sizes, sizes[1:])),
              ('nested size divisibility', n))
        values = {f[s] for s in lows}
        check(len(values) <= len(sizes) <= floor_log3(n), ('logarithmic value bound', n))
        if prime(n):
            check(sizes == [1] and len(lows) == 2*n and len(values) == 1,
                  ('primitive corollary', n))
            twice_W = sum(K[i][j] for i in range(n) for j in range(i+1, n))
            check(n*(m-H) == twice_W, ('primitive global singleton mean', n))
        divisor_masks = {d: sum(1 << (j*(n//d)) for j in range(d))
                         for d in range(1, n) if n % d == 0}
        divisor_values = {d: f[s] for d, s in divisor_masks.items()}
        check(min(divisor_values.values()) == m, ('divisor-cut minimum', n))
        cosets_and_complements = set()
        for subgroup in divisor_masks.values():
            s = subgroup
            for _ in range(n):
                cosets_and_complements.update((s, full ^ s))
                s = rotate(s, n)
        check(all(s in cosets_and_complements for s in range(1, full) if f[s] == m),
              ('minimum subgroup cosets', n))
        if n == 11:
            for s in lows | {f.index(m), f.index(M)}:
                check(path_count(ear(a, s)) == f[s], ('order11 direct low/extremum', s))
        missing = (rho+1)//2-len(values)
        row = dict(order=n, H=H, rho=rho, minimum=m, maximum=M,
                   small_side_sizes=sizes, low_states=len(lows), low_values=len(values),
                   missing=missing, vertex_transitive_missing_floor=max(0, (rho+1)//2-floor_log3(n)),
                   primitive_control=prime(n), divisor_cut_values=divisor_values)
        print(json.dumps(row, sort_keys=True), flush=True)
        semantic.append(row)
    print(json.dumps(dict(status='PASS', semantic_sha256=hashlib.sha256(
        json.dumps(semantic, sort_keys=True, separators=(',', ':')).encode()).hexdigest()), sort_keys=True))


if __name__ == '__main__':
    main()
