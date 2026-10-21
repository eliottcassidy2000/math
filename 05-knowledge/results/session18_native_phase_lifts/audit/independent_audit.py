#!/usr/bin/env python3
"""Independent integer audit; no producer import, repository import or dependency."""
import argparse
import hashlib
import json
from collections import defaultdict
from pathlib import Path

COUNT = 0
MANIFEST = 'f20cb439a11829cd1ba0529d06f5e598f5fa496faf19ac930697b9ad88850ffa'

def require(value, message):
    global COUNT
    COUNT += 1
    if not value:
        raise ValueError(message)

def digest(raw):
    return hashlib.sha256(raw).hexdigest()

def encoded(value):
    return (json.dumps(value, sort_keys=True, separators=(',', ':'))+'\n').encode()

def primitive(s):
    return s > 0 and s % 8 in (1, 3, 7)

def val2(n):
    return (n & -n).bit_length()-1

def native(s):
    # Unlike the producer's strip loop, recover the strip count from valuation.
    a = val2(3*s+1)
    require(a in (1, 2), 'native source branch')
    x = (3*s+1) >> a
    r = (val2(3*x+1)-1)//2
    u = ((3*x+1)//(4**r)-1)//3
    require(primitive(u), 'native primitive target')
    require((4**r*(3*u+1)-1)//3 == x, 'exact strip reconstruction')
    return u, a+2*r, r

def lifted(s):
    return 64*s + (49 if s % 8 == 1 else 35)

def rotate(v, c):
    # Use the other two real readouts, rather than iterating the producer R.
    a, b = v
    return ((a, b), (-b, a-b), (b-a, -a))[c % 3]

def norm(v):
    a, b = v
    return a*a-a*b+b*b

def plus(v, w):
    return v[0]+w[0], v[1]+w[1]

def minus(v, w):
    return v[0]-w[0], v[1]-w[1]

def beta(v, w):
    a, b = v
    c, d = w
    return 2*a*c-a*d-b*c+2*b*d

def row_to_charge(rows):
    require(sum(rows) == 0, 'three readout sum')
    require((rows[0]-rows[1]) % 3 == 0, 'readout integrality')
    a = (rows[0]-rows[1])//3
    b = (rows[2]-rows[1])//3
    require(rows == [2*a-b, -a-b, -a+2*b], 'three readout reconstruction')
    return a, b

def weight(residue, lower):
    row = (0, -3, 3) if lower else (-1, -1, 2)
    return row_to_charge([row[(residue+k) % 3] for k in range(3)])

def charge(N, s):
    rows = [0, 0, 0]
    while s < 4*N:
        row = (0, -3, 3) if s < 2*N else (-1, -1, 2)
        for k in range(3):
            rows[k] += row[(s+k) % 3]
        s = 4*s+1
    return row_to_charge(rows)

def events(N):
    # Solve each branch's rational crossing interval with two integer bounds.
    # This does not use the producer's greatest-odd-source candidate procedure.
    answer = []
    p = 1
    j = 0
    while (4*p-1)//3 < 4*N:
        b = (p-1)//3
        for K in (2, 4):
            for c in (35, 49):
                numerator = 64*K*N - (c-21)*p - 21 - 64*b
                lower = -((-numerator)//(64*p))
                upper = (K*N-b-1)//p
                require(upper-lower <= 0, 'interval width less than one')
                for s in range(lower, upper+1):
                    if not primitive(s) or s >= 4*N or lifted(s)-64*s != c:
                        continue
                    m = p*s+b
                    ml = (p*(3*lifted(s)+1)-1)//3
                    require(m < K*N and ml >= 64*K*N, 'interval crossing')
                    require(j > 0 and N < m < 4*N, 'single B-scale localization')
                    answer.append([s, j, K, c, m, ml])
        p *= 4
        j += 1
    require(len({r[0] for r in answer}) == len(answer), 'unit-defect uniqueness')
    require(len(answer) <= 2*(j-1), 'sparse global bound')
    return sorted(answer)

def expected_defect(event):
    s, j, K, c, m, ml = event
    w = weight(m % 3, False)
    if K == 2:
        w = minus(weight(m % 3, True), w)
    w = rotate(w, c)
    return -w[0], -w[1]

def route(s, V):
    t = d = 0
    seen = set()
    while s not in V:
        require(s not in seen and d < 10000, 'actual finite route')
        seen.add(s)
        s, clock, _ = native(s)
        t += clock
        d += 1
    return s, t, d

def main():
    parser = argparse.ArgumentParser()
    parser.add_argument('packet', type=Path)
    args = parser.parse_args()
    packet = args.packet
    require(digest((packet/'MANIFEST.json').read_bytes()) == MANIFEST, 'frozen manifest')
    manifest = json.loads((packet/'MANIFEST.json').read_bytes())
    for item in manifest['files']:
        raw = (packet/item['path']).read_bytes()
        require(len(raw) == item['bytes'] and digest(raw) == item['sha256'], 'frozen '+item['path'])
    data = json.loads((packet/'frontier.json').read_bytes())
    saved = json.loads((packet/'phase_receipt.json').read_bytes())
    V = set(data['frontier'])
    ranks = dict(data['frontier_ranks'])
    require(len(V) == len(data['frontier']) == len(ranks) == 627, 'frontier unique')
    require(set(ranks) == V and ranks[1] == 0, 'frontier exact rank support/root')
    for s in V:
        require(primitive(s) and type(ranks[s]) is int and ranks[s] >= 0, 'frontier source/rank')
        if s != 1:
            u, _, _ = native(s)
            require(u in V and ranks[u] < ranks[s], 'frontier descent')
    first_outside = next(s for s in range(1, max(V)+1) if primitive(s) and s not in V)
    require(first_outside == 1009 and max(V) == 63665, 'literal frontier extrema')
    require(lifted(first_outside) > max(V) and 64*first_outside+35 > max(V), 'outside guard')

    source_digest = hashlib.sha256()
    all_phase_digest = hashlib.sha256()
    phase_rows = 0
    for N in range(1, 513):
        ev = events(N)
        by_source = {r[0]: r for r in ev}
        found = []
        for s in range(1, 4*N, 2):
            if not primitive(s):
                continue
            ls = lifted(s)
            c = ls-64*s
            u, t, r = native(s)
            require(native(ls) == (u, t+6, r+3), 'integer lift clock/target')
            require(primitive(ls) and val2(3*ls+1) == val2(3*s+1), 'stable native branch')
            require(ls < 256*N, 'admission forward')
            v, w = charge(N, s), charge(64*N, ls)
            defect = minus(w, rotate(v, c))
            expect = expected_defect(by_source[s]) if s in by_source else (0, 0)
            require(defect == expect, 'complete phase law')
            require(norm(defect) in (0, 1), 'unit defect')
            require(charge(64*N, s) == v, 'inherited fixed-source transport')
            if defect != (0, 0):
                found.append(s)
            row = [N, s, ls, c, list(v), list(w), list(defect)]
            all_phase_digest.update(encoded(row))
            if N <= 256:
                source_digest.update(encoded(row))
            phase_rows += 1
        require(found == sorted(by_source), 'complete sparse exception equality')
        for s in range(max(1, 4*N-9), 4*N+10, 2):
            if primitive(s):
                require((s < 4*N) == (lifted(s) < 256*N), 'both admission directions')
    require(source_digest.hexdigest() == saved['phase_source_rows_sha256'], 'independent saved source digest')
    for control in saved['sparse_large_controls']:
        N = control['N']
        ev = events(N)
        require(ev == sorted(control['events']), 'large sparse control')
        for event in ev:
            s, _, _, c, _, _ = event
            require(minus(charge(64*N, lifted(s)), rotate(charge(N, s), c)) == expected_defect(event), 'large exact defect')

    inverse_rows = 0
    inverse_digest = hashlib.sha256()
    for s in range(1, 262144, 2):
        if not primitive(s):
            continue
        a = val2(3*s+1)
        u, t, r = native(s)
        r0 = ((2 if a == 1 else 1)-u) % 3
        x = (4**r0*(3*u+1)-1)//3
        base = (2**a*x-1)//3
        require(primitive(base) and val2(3*base+1) == a, 'inverse branch admissible')
        require(r >= r0 and (r-r0) % 3 == 0, 'inverse residue')
        k = (r-r0)//3
        c = 49 if a == 2 else 35
        recovered = 64**k*base+c*(64**k-1)//63
        require(recovered == s, 'closed-form inverse reconstruction')
        if s < 65536:
            inverse_digest.update(encoded([s, u, a, r, base, r0, k]))
        inverse_rows += 1
    require(inverse_digest.hexdigest() == saved['inverse_address_sha256'], 'saved inverse digest')
    for u in range(1, 16384, 2):
        if not primitive(u):
            continue
        for a in (1, 2):
            r0 = ((2 if a == 1 else 1)-u) % 3
            for k in range(5):
                r = r0+3*k
                x = (4**r*(3*u+1)-1)//3
                require((2**a*x-1) % 3 == 0, 'target-side inverse integer')
                s = (2**a*x-1)//3
                require(primitive(s) and val2(3*s+1) == a, 'target-side inverse primitive')
                require(native(s) == (u, a+2*r, r), 'target-side exact inverse')

    for rec in saved['packets']:
        N, M = rec['base_N'], rec['scale_M']
        require(M == 4096*N, 'packet scale')
        rows = rec['rows']
        s0 = rows[0][0]
        require([r[0] for r in rows] == [s0, lifted(s0), lifted(lifted(s0))], 'packet labels')
        total = (0, 0)
        for row in rows:
            s, r15, r35, u, t, d, a, b, q = row
            require(s not in V and (r15, r35) == (s % 15, s % 35), 'source guard/registers')
            require(route(s, V) == (u, t, d), 'actual first arrival')
            require(charge(M, s) == (a, b) and norm((a, b)) == q, 'actual packet charge')
            require(native(s)[0] == rec['immediate_target'], 'common first merge')
            total = plus(total, ((1-2*(t % 2))*a, (1-2*(t % 2))*b))
        require(total == (0, 0), 'signed packet zero')
        require(len({r[5] for r in rows}) == len({r[3] for r in rows}) == 1, 'equal terminal/depth')
        require([r[4]-rows[0][4] for r in rows] == [0, 6, 12], 'actual integer clock shifts')
        require([native(r[0])[1] for r in rows] == rec['edge_clocks'], 'literal edge clocks')
        # At a clock threshold only the first nonzero packet member is visible.
        first_clock_sum = tuple((1-2*(rows[0][4] % 2))*x for x in rows[0][6:8])
        require(norm(first_clock_sum) > 0, 'actual integer-clock partial-packet hostile')
        # In contrast, every nonnegative P-depth arrival cut sees all or none.
        for cut in range(rows[0][5]+2):
            arrived = [row for row in rows if row[5] <= cut]
            require(len(arrived) in (0, 3), 'synchronous P-depth cut')

    hostile = saved['boundary_hostile']
    require(events(2019)[0] == [31, 4, 4, 35, 8021, 516949], 'second literal hostile event')
    require([r for r in events(2019) if r[0] == 1009] == [[1009, 1, 2, 49, 4037, 258501]], 'sharp native boundary')
    require(charge(2019, 1009) == (-1, -2) and charge(64*2019, 64625) == (1, 0), 'boundary charge witness')
    require(expected_defect([1009, 1, 2, 49, 4037, 258501]) == (-1, -1), 'boundary unit defect')
    for row in hostile['triple_rows'] + saved['frontier_hostile_rows']:
        s, _, _, u, t, d, a, b, q = row
        scale = 4096*2019 if row in hostile['triple_rows'] else 4096
        require(route(s, V) == (u, t, d), 'hostile actual terminal/clock/depth')
        require(charge(scale, s) == (a, b), 'hostile actual charge')
    require([route(s, V) for s in (1, 113, 7281)] == [(1, 0, 0), (113, 0, 0), (1, 14, 1)], 'stopped-frontier counterexample')
    require(tuple(map(sum, zip(*(charge(4096*2019, s) for s in (1009, 64625, 4136049))))) == (0, -1), 'unguarded triple residual')
    require(norm(plus(charge(4096, 1), charge(4096, 7281)))+norm(charge(4096, 113)) == 6, 'false cross-terminal packet has positive stopped energy')

    # Exercise the same-scale compiler, without importing the producer.
    # Isolated singletons cannot be packets, so only adjacency endpoints need storage.
    compiler = []
    for N in (1, 3, 7, 55, 256, 2019, 4096, 16384, 65536, 65543, 129251):
        M = 64*N
        bad_sources = {e[0] for e in events(N)}
        successor = {}
        predecessor = {}
        nodes = set()
        for s in range(1, 4*N, 2):
            if not primitive(s):
                continue
            z = lifted(s)
            require(z < 4*M and s < 4*N, 'same-scale pair admission')
            if s in V or z in V:
                continue
            nodes.update((s, z))
            if s in bad_sources:
                continue
            require(z not in predecessor, 'no inverse-edge reuse')
            successor[s] = z
            predecessor[z] = s
        used = set()
        packets = 0
        nonzero_packets = 0
        for first in sorted(nodes-set(predecessor)):
            segment = [first]
            while segment[-1] in successor:
                segment.append(successor[segment[-1]])
            require(len(set(segment)) == len(segment), 'acyclic lift segment')
            for start in range(0, len(segment)-2, 3):
                ss = segment[start:start+3]
                require(not (used & set(ss)), 'no source reused between packets')
                used.update(ss)
                vv = [charge(M, s) for s in ss]
                require(plus(plus(vv[0], vv[1]), vv[2]) == (0, 0), 'same-scale compiler zero')
                routes = [route(s, V) for s in ss]
                require(len({(u, d) for u, t, d in routes}) == 1, 'compiler equal arrival terminal/depth')
                require([t-routes[0][1] for u, t, d in routes] == [0, 6, 12], 'compiler full clocks')
                packets += 1
                nonzero_packets += norm(vv[0]) != 0
        if N == 129251:
            require(64625 in bad_sources and 1009 not in bad_sources, 'second-edge boundary hostile')
            require(not (used & {1009, 64625, 4136049}), 'bad edge retains the three residual sources')
            vv = [charge(M, s) for s in (1009, 64625, 4136049)]
            require(plus(plus(vv[0], vv[1]), vv[2]) == (1, 0), 'bad second-edge packet is not zero')
        compiler.append({'N': N, 'M': M, 'packets': packets, 'nonzero_packets': nonzero_packets, 'used_sources': len(used)})

    # Verify algebra beyond the twelve-source alphabet, including external coherence.
    for a in range(-4, 5):
        for b in range(-4, 5):
            v = (a, b)
            for c in (1, 2):
                w = rotate(v, c)
                w2 = rotate(v, 2*c)
                require(norm(v) == norm(w) and plus(plus(v, w), w2) == (0, 0), 'rotation algebra')
                require(norm(plus(v, w)) == norm(v) and beta(v, w) == -norm(v), 'pair diagonal payment')
                require(minus((0, 0), plus(v, w2)) == w, 'clock polynomial coefficient factorization')
                for x in range(-2, 3):
                    for y in range(-2, 3):
                        Z = (x, y)
                        require(norm(plus(plus(Z, v), w))-norm(plus(Z, v)) == beta(Z, w), 'external cross-term survives')
    require(norm(plus((1, 0), rotate((1, 0), 1))) == 1, 'clock-prefix nonzero hostile')
    result = {'status': 'PASS', 'scope': 'Scoped phase-lift and packet claims; no uniform annulus or entry inference',
              'packet_manifest_sha256': MANIFEST, 'packet_receipt_sha256': digest((packet/'phase_receipt.json').read_bytes()),
              'checks': COUNT, 'phase_N_range': [1, 512], 'phase_rows': phase_rows,
              'extended_phase_sha256': all_phase_digest.hexdigest(), 'inverse_source_limit_exclusive': 262144,
              'inverse_rows': inverse_rows, 'inverse_target_limit_exclusive': 16384, 'inverse_chain_k_range': [0, 4],
              'same_scale_compiler': compiler}
    print(json.dumps(result, sort_keys=True))

if __name__ == '__main__':
    main()
