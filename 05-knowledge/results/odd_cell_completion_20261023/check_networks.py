"""Exact small attachment-network orbits and the sharp nine-cell boundary.

No assertion is used for a load-bearing test, so -O runs the same checks.
Frozen output is compared by default; --write regenerates it intentionally.
"""
from collections import Counter, deque
from itertools import combinations, permutations, product
from math import comb
from pathlib import Path
import json
import sys

HERE = Path(__file__).resolve().parent
checks = 0


def require(condition, message):
    global checks
    checks += 1
    if not condition:
        raise RuntimeError(message)


def formula(p, vertices, degrees, counts):
    raw = 1
    fixed = 1
    for d in degrees:
        raw *= counts[d]
        fixed *= comb(p, d)
    numerator = raw + (p - 1) ** vertices * fixed
    sym_numerator = (p - 1) ** (vertices - 1) * fixed
    require(numerator % p ** vertices == 0, "orbit integrality")
    require(sym_numerator % p ** (vertices - 1) == 0, "symmetric integrality")
    total = numerator // p ** vertices
    symmetric = sym_numerator // p ** (vertices - 1)
    return dict(total=total, symmetric=symmetric, rigid=total-symmetric,
                labelled=raw)


def matching_orbits(p, vertices, edges):
    alphabet = list(permutations(range(p)))
    ids = {x: i for i, x in enumerate(alphabet)}
    size = len(alphabet)
    # Direct image under a rotation at the row/column cell.
    row = [ids[tuple(q[(x-1) % p] for x in range(p))] for q in alphabet]
    col = [ids[tuple((q[x]+1) % p for x in range(p))] for q in alphabet]
    edge_orbit = {}
    edge_slopes = {}
    for seed in range(size):
        if seed in edge_orbit: continue
        orbit = {seed}
        todo = [seed]
        for state in todo:
            for new in (row[state],col[state]):
                if new not in orbit: orbit.add(new); todo.append(new)
        for state in orbit: edge_orbit[state]=seed
        q=alphabet[seed]
        diffs={(q[(i+1)%p]-q[i])%p for i in range(p)}
        edge_slopes[seed]=next(iter(diffs)) if len(diffs)==1 else None
    powers = [size ** j for j in range(len(edges))]
    incident = [[(j, row if u == v else col)
                 for j, (u, w) in enumerate(edges) if u == v or w == v]
                for v in range(vertices)]
    total_states = size ** len(edges)
    seen = bytearray(total_states)
    histogram = Counter()
    local_fibres = Counter()
    for seed in range(total_states):
        if seen[seed]:
            continue
        seen[seed] = 1
        orbit = [seed]
        for state in orbit:
            coordinates = [(state // x) % size for x in powers]
            for at_vertex in incident:
                new = state + sum((action[coordinates[j]]-coordinates[j])*powers[j]
                                  for j, action in at_vertex)
                if not seen[new]:
                    seen[new] = 1
                    orbit.append(new)
        histogram[len(orbit)] += 1
        local_fibres[tuple(edge_orbit[(seed//x)%size] for x in powers)] += 1
    require(sum(k*v for k, v in histogram.items()) == total_states, "orbit mass")
    expected = formula(p, vertices, [1]*len(edges), {1:size})
    require(sum(histogram.values()) == expected['total'], "network orbit formula")
    require(histogram[p ** (vertices-1)] == expected['symmetric'], "symmetric orbits")
    require(histogram[p ** vertices] == expected['rigid'], "rigid orbits")
    require(set(histogram) <= {p ** (vertices-1), p ** vertices}, "stabilizer dichotomy")
    fibre_histogram=Counter()
    beta=len(edges)-vertices+1
    for signature,count in local_fibres.items():
        slopes=[edge_slopes[x] for x in signature]
        rigid=sum(s is None for s in slopes)
        balanced=False
        if rigid==0:
            gains={0:1};todo=[0];consistent=True
            for u in todo:
                for (a,b),s in zip(edges,slopes):
                    if u==a: v=b;candidate=gains[u]*s%p
                    elif u==b: v=a;candidate=gains[u]*pow(s,-1,p)%p
                    else: continue
                    if v not in gains: gains[v]=candidate;todo.append(v)
                    elif gains[v]!=candidate:consistent=False
            balanced=consistent
        dimension=beta+rigid-1+int(balanced)
        require(dimension>=0 and count==p**dimension,'fixed local orbit alignment dimension')
        fibre_histogram[(rigid,balanced,dimension,count)]+=1
    return dict(p=p, vertices=vertices, edges=edges,
                orbit_size_histogram=dict(sorted(histogram.items())),
                local_orbit_signatures=len(local_fibres),
                alignment_fibres=[dict(rigid_edges=k[0],balanced=k[1],dimension=k[2],
                                      classes_per_signature=k[3],signatures=v)
                                  for k,v in sorted(fibre_histogram.items())], **expected)


def regular_binary_counts(p):
    counts = {}
    for d in range(1, (p-1)//2+1):
        options = [tuple(int(j in s) for j in range(p)) for s in combinations(range(p), d)]
        count = 0
        # Independent row-word scan, not Burnside or a coefficient recurrence.
        for rows in product(options, repeat=p):
            if all(sum(row[j] for row in rows) == d for j in range(p)):
                count += 1
        counts[d] = count
        counts[p-d] = count
    return counts


def wreath_tournament():
    def beats(u, v):
        a, x = divmod(u, 3)
        b, y = divmod(v, 3)
        return int((y-x) % 3 == 1) if a == b else int((b-a) % 3 == 1)
    return [[0 if u == v else beats(u, v) for v in range(9)] for u in range(9)]


def autos(matrix):
    n = len(matrix)
    result = []
    def extend(partial, unused):
        u = len(partial)
        if u == n:
            result.append(tuple(partial))
            return
        for v in unused:
            if all(matrix[u][k] == matrix[v][partial[k]] for k in range(u)):
                extend(partial+[v], [x for x in unused if x != v])
    extend([], list(range(n)))
    return result


def refine(matrix, cells, marks=()):
    cells = [[v] for v in marks] + [[v for v in c if v not in marks] for c in cells]
    cells = [c for c in cells if c]
    while True:
        new = []
        for cell in cells:
            parts = {}
            for v in cell:
                key = tuple(sum(matrix[v][w] for w in c) for c in cells)
                parts.setdefault(key, []).append(v)
            new.extend(parts[key] for key in sorted(parts))
        if len(new) == len(cells):
            return new
        cells = new


def nine_boundary():
    r9 = wreath_tournament()
    group = autos(r9)
    require(len(group) == 81, 'R9 automorphism count')
    require(all(sum(r) == 4 for r in r9), 'R9 regular')
    block = [[int(u//3 == v//3) for v in range(9)] for u in range(9)]
    require(all(sum(r) == 3 for r in block), 'block row degree')
    require(all(sum(block[u][v] for u in range(9)) == 3 for v in range(9)), 'block column degree')
    compatible = [(g,h) for g in group for h in group
                  if all(block[u][v] == block[g[u]][h[v]] for u in range(9) for v in range(9))]
    require(len(compatible) == 2187, 'pair compatibility count')
    fibers = Counter(g for g,h in compatible)
    require(set(fibers.values()) == {27}, 'nonfunctional fibres')
    full = [[0]*18 for _ in range(18)]
    for u in range(9):
        for v in range(9):
            full[u][v] = full[9+u][9+v] = r9[u][v]
            full[u][9+v] = block[u][v]
            full[9+v][u] = 1-block[u][v]
    cells = [list(range(9)), list(range(9,18))]
    source_named = refine(full, cells, tuple(range(9)))
    require(sorted(map(len,source_named)) == [1]*9+[3]*3, 'naming entire source is insufficient')
    # Six independent internal C3 rotations force one owner per triple.
    base = (0,3,6,9,12,15)
    require(len(refine(full,cells,base)) == 18, 'six-owner refinement certificate')
    all_maps = [g+tuple(9+x for x in h) for g,h in compatible]
    require(sum(all(g[v] == v for v in base) for g in all_maps) == 1, 'six-owner base')
    for omitted in base:
        require(sum(all(g[v] == v for v in base if v != omitted) for g in all_maps) >= 3,
                'each internal triple still needs an owner')
    return dict(internal_auts=81, pair_auts=2187, targets_per_fixed_source_chart=27,
                source_named_partition=sorted(map(len,source_named)),
                exact_minimum_owners=6, base=base)


def prime_eleven_boundary():
    # Strong four-tournament: path 0->1->2->3->0 plus 0->2,1->3.
    r4 = [[0]*4 for _ in range(4)]
    for u,v in [(0,1),(1,2),(2,3),(3,0),(0,2),(1,3)]:
        r4[u][v] = 1
    low = [u for u in range(4) if sum(r4[u]) == 1]
    high = [u for u in range(4) if sum(r4[u]) == 2]
    r11 = [[0]*11 for _ in range(11)]
    for u in range(3):
        r11[u][(u+1)%3] = 1
    for u in range(4):
        for v in range(4):
            r11[3+u][3+v] = r11[7+u][7+v] = r4[u][v]
            bit = int((u,v) in set(zip(low,high)))
            r11[3+u][7+v] = bit
            r11[7+v][3+u] = 1-bit
        for v in range(3):
            r11[3+u][v] = r11[v][7+u] = 1
    require(all(sum(r) == 5 for r in r11), 'prime11 regularity')
    require(all(r11[u][v]+r11[v][u] == 1 for u in range(11) for v in range(u+1,11)),
            'prime11 tournament')
    rotation = (1,2,0)+tuple(range(3,11))
    require(all(r11[u][v] == r11[rotation[u]][rotation[v]] for u in range(11) for v in range(11)),
            'prime11 hidden C3 module')
    block = [[int(u<3 and v<3) if u<3 or v<3 else int((v-u)%8 in (0,1,2))
              for v in range(11)] for u in range(11)]
    require(all(sum(r) == 3 for r in block), 'prime11 row degree')
    require(all(sum(block[u][v] for u in range(11)) == 3 for v in range(11)), 'prime11 column degree')
    require(all(block[u][v] == block[u][rotation[v]] for u in range(11) for v in range(11)),
            'prime11 nontrivial one-sided kernel')
    return dict(internal_rows=[''.join(map(str,r)) for r in r11],
                block_rows=[''.join(map(str,r)) for r in block], kernel_generator=rotation,
                scope='Prime cell size alone does not imply functional chart attachments.')


def main():
    counts = regular_binary_counts(5)
    require(counts == {1:120,4:120,2:2040,3:2040}, 'independent regular-five census')
    networks = [matching_orbits(3,3,[(0,1),(1,2)]),
                matching_orbits(3,3,[(0,1),(1,2),(0,2)]),
                matching_orbits(3,4,[(0,1),(1,2),(2,3),(0,3)]),
                matching_orbits(5,3,[(0,1),(1,2)]),
                matching_orbits(5,3,[(0,1),(1,2),(0,2)])]
    mixed = formula(5,3,[1,2,1],counts)
    boundary = nine_boundary()
    prime_boundary = prime_eleven_boundary()
    output = dict(scope='Named cells, actual shared rotations, fixed nonuniform support; no quotient by cell permutations.',
                  regular_five_counts=counts, complete_network_censuses=networks,
                  mixed_degree_triangle_formula=mixed, nine_boundary=boundary,
                  prime_eleven_boundary=prime_boundary, checks=checks)
    target = HERE/'networks.json'
    text = json.dumps(output,indent=2,sort_keys=True)+'\n'
    if '--write' in sys.argv:
        target.write_text(text,encoding='utf-8')
    else:
        require(json.loads(target.read_text()) == json.loads(text),'frozen receipt')
    print(json.dumps(dict(checks=checks,network_orbits=[x['total'] for x in networks],
                          nine_pair_auts=boundary['pair_auts']),sort_keys=True))


if __name__ == '__main__':
    main()
