"""Full n=8 port replay and small incidence/electrical certificates.

No removable assertions; Python -O produces the same receipt. The C++ producer
uses Hamiltonian expansion; this consumer checks independent orbit accounting,
the complement ledger, profile decoding, and binary linear algebra.
"""
from pathlib import Path
from itertools import combinations, permutations, product
from collections import Counter, defaultdict, deque
import argparse, csv, gzip, hashlib, json, math, struct, subprocess, tempfile

HERE = Path(__file__).resolve().parent
ROOT = HERE.parents[2]
CENSUS = ROOT / '04-computation/experiments/tiling_line_metagraph_census_20261020_classes_n8.tsv'
CHECKS = 0

def need(ok, message):
    global CHECKS
    CHECKS += 1
    if not ok:
        raise ValueError(message)

def digest(data):
    return hashlib.sha256(data).hexdigest()

def adjacency(code):
    A = [[0] * 8 for _ in range(8)]
    for k, (i, j) in enumerate(combinations(range(8), 2)):
        A[i][j] = (code >> (27-k)) & 1
        A[j][i] = 1-A[i][j]
    return A

def profiles(A):
    n = len(A)
    s = list(map(sum, A)); c = [0]*n; f = [0]*n
    for i, j, k in combinations(range(n), 3):
        if A[i][j] == A[j][k] == A[k][i]:
            for u in (i, j, k): c[u] += 1
    for B in combinations(range(n), 4):
        if sorted(sum(A[u][v] for v in B) for u in B) == [1, 1, 2, 2]:
            for u in B: f[u] += 1
    tags = [(s[u], c[u], f[u]) for u in range(n)]
    rho = [(tags[u], tuple(sorted(tags[v] for v in range(n) if A[u][v]))) for u in range(n)]
    types = sorted(set(rho)); rank = [types.index(x) for x in rho]
    cells = [[u for u in range(n) if rank[u] == j] for j in range(len(types))]
    E = tuple(tuple(sum(A[u][v] for u in C for v in D) for D in cells) for C in cells)
    return tuple(sorted(rho)), rho, rank, cells, E

def overlap(A, D, v, w):
    return sum(A[u][v]*A[u][w] for u in D)

def binary_kernel(rows, n):
    """Reduced row echelon form and a nullspace basis over F2."""
    piv = {}
    for row in rows:
        while row:
            p = (row & -row).bit_length()-1
            if p in piv: row ^= piv[p]
            else:
                piv[p] = row
                break
    for p in sorted(piv, reverse=True):
        for q in piv:
            if q < p and (piv[q] >> p) & 1: piv[q] ^= piv[p]
    basis = []
    for f in range(n):
        if f not in piv:
            v = 1 << f
            for p, row in piv.items():
                if (row >> f) & 1: v |= 1 << p
            basis.append(v)
    need(all(not ((r & v).bit_count() & 1) for r in rows for v in basis), 'kernel verification')
    return basis, len(piv)

def binary_rank(rows):
    piv = {}
    for row in rows:
        while row:
            p = row.bit_length()-1
            if p in piv: row ^= piv[p]
            else:
                piv[p] = row
                break
    return len(piv)

def binary_solve(rows, rhs, n):
    piv = {}
    for row, b in zip(rows, rhs):
        while row:
            p = (row & -row).bit_length()-1
            if p in piv:
                r, s = piv[p]; row ^= r; b ^= s
            else:
                piv[p] = (row, b); break
        else: need(not b, 'inconsistent binary lift')
    x = 0
    for p in sorted(piv, reverse=True):
        row, b = piv[p]
        if ((row & x).bit_count() & 1) != b: x |= 1 << p
    need(all(((r & x).bit_count() & 1) == b for r, b in zip(rows, rhs)), 'binary solve residual')
    return x

def bits(x):
    while x:
        p = (x & -x).bit_length()-1
        yield p
        x &= x-1

def components(nodes, edges):
    parent = {v: v for v in nodes}
    def root(v):
        while parent[v] != v:
            parent[v] = parent[parent[v]]; v = parent[v]
        return v
    for a, b, w in edges:
        if w: parent[root(a)] = root(b)
    groups = defaultdict(list)
    for v in nodes: groups[root(v)].append(v)
    return sorted((sorted(g) for g in groups.values()), key=lambda g: (len(g), g))

def decode_doubleton(A, rho, cells, E):
    """Reconstruct with singleton edges, margins and a star of overlaps."""
    double = [j for j, C in enumerate(cells) if len(C) == 2]
    need(len(double) == 1 and all(len(C) <= 2 for C in cells), 'doubleton scope')
    d = double[0]; D = cells[d]
    single = [j for j, C in enumerate(cells) if len(C) == 1]
    mixed = [j for j in single if E[d][j] == 1]
    need(bool(mixed), 'equal-score doubleton needs mixed observers')
    root = mixed[0]
    # Name one row by its arrow to the root observer: this fixes the swap gauge.
    assigned = {root: 1}
    for j in mixed[1:]: assigned[j] = overlap(A, D, cells[root][0], cells[j][0])
    B = [[0]*len(A) for _ in A]
    u, v = D
    for i, j in combinations(single, 2):
        a, b = cells[i][0], cells[j][0]
        B[a][b] = E[i][j]; B[b][a] = 1-B[a][b]
    for j in single:
        w = cells[j][0]
        if E[d][j] == 1:
            B[u][w] = assigned[j]; B[v][w] = 1-assigned[j]
        else: B[u][w] = B[v][w] = E[d][j]//2
        B[w][u] = 1-B[u][w]; B[w][v] = 1-B[v][w]
    score = rho[u][0][0]
    B[u][v] = score-sum(B[u]); B[v][u] = 1-B[u][v]
    need(B[u][v] in (0, 1) and sum(B[v]) == score, 'internal degree feasibility')
    swap = list(range(len(A))); swap[u], swap[v] = swap[v], swap[u]
    need(B == A or [[B[swap[i]][swap[j]] for j in range(len(A))] for i in range(len(A))] == A, 'doubleton decoder lost incidence')
    for i, j in combinations(mixed, 2):
        need(overlap(A, D, cells[i][0], cells[j][0]) == int(assigned[i] == assigned[j]), 'overlap parity equation')
    return {'doubleton_type': d, 'mixed_observer_types': mixed, 'tree_bits': len(mixed)-1}

def main():
    parser = argparse.ArgumentParser()
    parser.add_argument('--prefix', type=Path, help='reuse outputs of the C++ producer')
    parser.add_argument('--write', action='store_true')
    args = parser.parse_args()
    rows = list(csv.DictReader(CENSUS.open(), delimiter='\t'))
    fields = ('class', 'canon', 'H', 'Aut', 'tc', 'bm', 'sigmaClass')
    data = [{k: int(r[k]) for k in fields} for r in rows]
    out = {'status': 'FINITE-EXACT', 'scope': 'all 6880 ordinary n=8 tournaments and all 2^21 ascending fixed-path tilings',
           'sources': {str(CENSUS.relative_to(ROOT)): digest(CENSUS.read_bytes())}}
    matrices = []; sigs = []; prof = []; groups = defaultdict(list); canons = set(); mass = 0; perms = 0
    for r in data:
        A = adjacency(r['canon']); matrices.append(A)
        s = list(map(sum, A)); blocks = [[u for u in range(8) if s[u] == a] for a in sorted(set(s))]
        best = None; aut = 0
        for bp in product(*(permutations(B) for B in blocks)):
            p = sum(bp, ()); x = 0
            for i, j in combinations(range(8), 2): x = 2*x + A[p[i]][p[j]]
            perms += 1
            if best is None or x < best: best, aut = x, 1
            elif x == best: aut += 1
        need(best == r['canon'] and aut == r['Aut'], 'independent score canonicalization')
        need(best not in canons, 'duplicate class'); canons.add(best); mass += math.factorial(8)//aut
        P = profiles(A); prof.append(P); sigs.append(P[0]); groups[P[0]].append(r['class'])
    need(mass == 1 << 28 and len(data) == 6880, 'labelled tournament orbit mass')
    out['orbit_audit'] = {'labelled_mass': mass, 'score_preserving_permutations': perms}
    collisions = []
    for sig, ids in groups.items():
        if len(ids) == 1: continue
        need(len(ids) == 2, 'unexpected larger profile collision')
        p, q = ids; P, Q = prof[p], prof[q]
        need(P[4] == Q[4], 'full profile quotient should still collide')
        d = next(j for j, C in enumerate(P[3]) if len(C) == 2)
        candidates = []
        singles = [j for j, C in enumerate(P[3]) if len(C) == 1]
        for i, j in combinations(singles, 2):
            values = [overlap(matrices[c], prof[c][3][d], prof[c][3][i][0], prof[c][3][j][0]) for c in ids]
            if values[0] != values[1]: candidates.append({'observer_types': [i, j], 'values': values})
        need(bool(candidates) and all(set(v['values']) == {0, 1} for v in candidates), 'one-bit separator')
        dec = [decode_doubleton(matrices[c], prof[c][1], prof[c][3], prof[c][4]) for c in ids]
        collisions.append({'classes': ids, 'canons': [data[c]['canon'] for c in ids],
                           'H': [data[c]['H'] for c in ids], 'doubleton_type': d,
                           'equal_complete_type_incidence': True, 'separators': candidates, 'decoders': dec})
    need(len(groups) == 6878 and len(collisions) == 2, 'known J4 boundary')
    # Recheck the intrinsic selector under relabelling and converse.
    for entry in collisions:
        for c in entry['classes']:
            A = matrices[c]; original = prof[c]
            for shift in range(8):
                p = [(i+shift)%8 for i in range(8)]
                B = [[A[p[i]][p[j]] for j in range(8)] for i in range(8)]
                P = profiles(B)
                need(P[0] == original[0] and P[4] == original[4], 'profile relabelling')
                decode_doubleton(B, P[1], P[3], P[4])
            B = [list(x) for x in zip(*A)]
            need(profiles(B)[0] == sigs[data[c]['sigmaClass']], 'converse profile transport')
    out['incidence'] = {'J4_signatures': len(groups), 'repaired_signatures': len(data), 'collisions': collisions}
    print('profiles and orbit audit complete', flush=True)

    with tempfile.TemporaryDirectory(prefix='incidence-n8-') as td:
        tmp = Path(td)
        if args.prefix is None:
            exe = tmp/'wiring'; inp = tmp/'input.txt'; prefix = tmp/'n8'
            inp.write_text(''.join(' '.join(str(r[k]) for k in fields)+'\n' for r in data))
            subprocess.run(['clang++', '-O2', '-std=c++17', str(HERE/'wiring.cpp'), '-o', str(exe)], check=True)
            producer = json.loads(subprocess.check_output([str(exe), str(inp), str(prefix)]))
        else:
            prefix = args.prefix
            producer = {'classes':6880, 'tilings':2097152, 'literal_paths':2143304, 'closed_ports':1048576,
                        'blue_tilings':4096, 'aggregated_port_entries':956339, 'merged_colored_entries':438374}
        owners_raw = prefix.with_suffix('.owners').read_bytes()
        ports_raw = prefix.with_suffix('.ports').read_bytes()
        edges_raw = prefix.with_suffix('.edges').read_bytes()
    owners = [v[0] for v in struct.iter_unpack('<H', owners_raw)]
    ports = Counter()
    for key, count in struct.iter_unpack('<QI', ports_raw):
        c = key & 8191; t = (key >> 13) & 8191; a = (key >> 26) & 7; b = (key >> 29) & 7; blue = (key >> 32) & 1
        ports[c, t, a, b, blue] += count
    merged = Counter()
    for key, count in struct.iter_unpack('<QI', edges_raw):
        merged[key & 8191, (key >> 13) & 8191, (key >> 26) & 1] += count
    need(len(owners) == 1 << 21 and len(ports) == 956339 and len(merged) == 438374, 'full output sizes')
    need(sum(ports.values()) == 1 << 20 and sum(merged.values()) == 1 << 20, 'full complement pairing mass')
    need(sum(r['H'] for r in data) == producer['literal_paths'], 'producer HP mass')
    need(sum(r['bm'] for r in data) == producer['blue_tilings'], 'producer blue mass')
    fib = Counter(owners)
    need(all(fib[r['class']] == r['tc'] for r in data), 'independent owner fibres')
    merged_from_ports = Counter(); colcap = Counter(); sourcecap = Counter()
    def node(c): return min(c, data[c]['sigmaClass'])
    for (c, t, a, b, blue), count in ports.items():
        u, v = sorted((node(c), node(t))); merged_from_ports[u, v, blue] += count
        colcap[t, 8-a, 6-b, blue] += count; sourcecap[c, a, b, blue] += count
        conv = (data[c]['sigmaClass'], data[t]['sigmaClass'], 7-b, 7-a, blue)
        need(ports[conv] == count, 'converse-reflected exact port count')
    need(merged_from_ports == merged, 'independent closed-to-open line ledger')
    need(all(count % 2 == 0 for (u, v, blue), count in merged.items() if not blue), 'black twin parity')
    # Independently recover endpoint/open capacities from the owners, not the port output.
    cells = [(i, j) for i in range(8) for j in range(i+2, 8)]
    first = sum(1 << k for k, (i, j) in enumerate(cells) if i == 0)
    last = sum(1 << k for k, (i, j) in enumerate(cells) if j == 7)
    closing_bit = 1 << cells.index((0, 7)); F = (1 << 21)-1
    independent_cols = Counter(); independent_source = Counter(); independent_ports = Counter()
    reflection = {k: cells.index((7-j, 7-i)) for k, (i, j) in enumerate(cells)}
    orbits = sorted({tuple(sorted({k, reflection[k]})) for k in range(21)})
    blue_masks = {sum(sum(1 << k for k in C) for j, C in enumerate(orbits) if word >> j & 1)
                  for word in range(1 << len(orbits))}
    need(len(blue_masks) == 4096, 'independent blue reflection masks')
    for m, c in enumerate(owners):
        a = 1+(m & first).bit_count(); b = ((~m) & last).bit_count()
        if not (m & closing_bit):
            independent_source[c, a, b] += 1
            independent_ports[c, owners[m ^ F], a, b, int(m in blue_masks)] += 1
        else: independent_cols[c, a, b] += 1
    need(independent_cols == Counter({k: sum(colcap[k+(b,)] for b in (0, 1)) for k in independent_cols}), 'open endpoint capacities')
    need(independent_source == Counter({k: sum(sourcecap[k+(b,)] for b in (0, 1)) for k in independent_source}), 'closed endpoint capacities')
    need(independent_ports == ports, 'independent literal mask-to-port replay')
    del independent_ports
    # Only the four collided targets need a nontrivial transportation test.
    ambiguity = []; counterfeits = []
    for entry in collisions:
        ids = entry['classes']; records = defaultdict(dict)
        for (c, t, a, b, blue), count in ports.items():
            if t in ids: records[c, a, b, blue][t] = count
        bytype = defaultdict(list)
        for row, values in records.items(): bytype[row[1:]].append((row, values))
        unresolved_mass = 0; residual_rows = 0; eligible_rows = 0
        for (a, b, blue), rr in bytype.items():
            available = [t for t in ids if colcap[t, 8-a, 6-b, blue]]
            if len(available) == 2:
                eligible_rows += len(rr)
                if len(rr) > 1:
                    unresolved_mass += sum(sum(v.values()) for _, v in rr); residual_rows += len(rr)
                    # An actual two-row swap is a certificate of nonuniqueness.
                    p, q = ids
                    left = next(((r, v) for r, v in rr if v.get(p, 0)), None)
                    right = next(((r, v) for r, v in rr if v.get(q, 0) and r != left[0]), None) if left else None
                    if right and not blue and not counterfeits:
                        x, y = left[0], right[0]
                        edits = [(x, p, -1), (x, q, 1), (y, q, -1), (y, p, 1)]
                        both = list(edits)
                        for (c, a, b, blue), t, v in edits:
                            both.append(((data[c]['sigmaClass'], 7-b, 7-a, blue), data[t]['sigmaClass'], v))
                        delta = Counter()
                        for (c, a, b, blue), t, v in both: delta[c, t, a, b, blue] += v
                        need(all(ports[e]+v >= 0 for e, v in delta.items()), 'negative reflected counterfeit')
                        dr = Counter(); dc = Counter(); dm = Counter()
                        for (c, t, a, b, blue), v in delta.items():
                            dr[c, a, b, blue, sigs[t]] += v; dc[t, 8-a, 6-b, blue] += v
                            u, w = sorted((node(c), node(t))); dm[u, w, blue] += v
                        need(not any(dr.values()) and not any(dc.values()), 'counterfeit changes capacities')
                        need(any(dm.values()) and all(v % 2 == 0 for v in dm.values()), 'counterfeit must change graph and preserve twins')
                        counterfeits.append({'port_edits': [[*e, v] for e, v in sorted(delta.items()) if v],
                                             'merged_edits': [[*e, v] for e, v in sorted(dm.items()) if v]})
        ambiguity.append({'target_classes': ids, 'candidate_rows': len(records), 'two_candidate_rows': eligible_rows,
                          'unresolved_rows': residual_rows, 'unresolved_mass': unresolved_mass})
    need(bool(counterfeits), 'need genuine n8 ambiguity witness')
    out['ports'] = {'producer': producer, 'binary_sha256': {'owners': digest(owners_raw), 'ports': digest(ports_raw), 'edges': digest(edges_raw)},
                    'J4_ambiguity': ambiguity, 'J4_unresolved_mass': sum(x['unresolved_mass'] for x in ambiguity),
                    'one_bit_repair_unresolved_mass': 0, 'reflected_counterfeit': counterfeits[0]}
    print('full port and counterfeit checks complete', flush=True)

    nodes = sorted({node(c) for c in range(6880)})
    blue_edges = [(u, v, w) for (u, v, b), w in merged.items() if b]
    black_edges = [(u, v, w//2) for (u, v, b), w in merged.items() if not b]
    full_edges = [(u, v, w) for (u, v, b), w in merged.items()]
    cc = components(nodes, full_edges); bc = components(nodes, black_edges); lc = components(nodes, blue_edges)
    blue_nodes = sorted({u for e in blue_edges for u in e[:2]})
    need(len(blue_nodes) == 176 and len(cc) == 1, 'n8 connected and selfconverse blue inventory')
    isolated = sorted(set(nodes)-set(blue_nodes)); bi = {u: i for i, u in enumerate(blue_nodes)}
    need(len([C for C in lc if len(C)>1]) == 1, 'blue nontrivial component connected')
    L = [[0]*176 for _ in range(176)]
    for u, v, w in blue_edges:
        if u != v:
            i, j = bi[u], bi[v]; L[i][i] += w; L[j][j] += w; L[i][j] -= w; L[j][i] -= w
    # Fix one blue root; the other blue null vectors give the bicycle coordinates.
    reduced = [sum((L[i][j] & 1) << (j-1) for j in range(1, 176)) for i in range(1, 176)]
    blue_kernel, _ = binary_kernel(reduced, 175)
    blue_kernel = [v << 1 for v in blue_kernel]
    need(len(blue_kernel) == 11, 'blue bicycle dimension')
    dim = len(isolated)+len(blue_kernel)
    coords = {u: 1 << i for i, u in enumerate(isolated)}
    for u in blue_nodes:
        coords[u] = sum(((v >> bi[u]) & 1) << (len(isolated)+j) for j, v in enumerate(blue_kernel))
    Q = [0]*dim
    for u, v, w in black_edges:
        if w & 1:
            z = coords[u]^coords[v]
            for i in bits(z): Q[i] ^= z
    blue_half = [[0]*11 for _ in range(11)]
    for i, a in enumerate(blue_kernel):
        for j, b in enumerate(blue_kernel):
            total = sum(L[x][y] for x in bits(a) for y in bits(b))
            need(total % 2 == 0, 'half energy integrality')
            blue_half[i][j] = (total//2) & 1
            if blue_half[i][j]: Q[len(isolated)+i] ^= 1 << (len(isolated)+j)
    need(all(((Q[i] >> j)&1) == ((Q[j] >> i)&1) for i in range(dim) for j in range(i)), 'full half form symmetry')
    rankQ = binary_rank(Q)
    # Lift the five radical directions to actual integer potentials with
    # Lap(G)v divisible by four; compute the next form without floating point.
    radical, rankQ2 = binary_kernel(Q, dim)
    need(rankQ2 == rankQ, 'independent pivot order half-form rank')
    ni = {u: i for i, u in enumerate(nodes)}
    indexed_edges = [(ni[u], ni[v], w) for u, v, w in full_edges if u != v]
    def laplacian(v):
        y = [0]*len(nodes)
        for a, b, w in indexed_edges:
            z = w*(v[a]-v[b]); y[a] += z; y[b] -= z
        return y
    lifts = []
    laplifts = []
    for r in radical:
        x = [(coords[u] & r).bit_count() & 1 for u in nodes]
        Ax = laplacian(x)
        need(all(v % 2 == 0 for v in Ax), 'radical integral half image')
        need(all((Ax[ni[u]]//2) % 2 == 0 for u in isolated), 'radical image at blue isolates')
        rhs = [(Ax[ni[u]]//2) & 1 for u in blue_nodes[1:]]
        y = binary_solve(reduced, rhs, 175) << 1
        v = list(x)
        for u in blue_nodes: v[ni[u]] += 2*((y >> bi[u]) & 1)
        Av = laplacian(v)
        need(all(a % 4 == 0 for a in Av), 'quarter-energy lift')
        lifts.append(v); laplifts.append(Av)
    quarter = [[(sum(a*b for a, b in zip(v, Aw))//4) & 1 for Aw in laplifts] for v in lifts]
    quarter_rows = [sum(b << j for j, b in enumerate(row)) for row in quarter]
    quarter_rank = binary_rank(quarter_rows)
    quarter_radical, _ = binary_kernel(quarter_rows, len(radical))
    eighth_lifts = []; eighth_images = []
    for r in quarter_radical:
        v = [sum(lifts[j][i] for j in bits(r)) for i in range(len(nodes))]
        Av = laplacian(v)
        # First cancel the obstruction against all of V using q2. Vanishing
        # against rad(q2) is exactly the q4-radical condition.
        obstruction = 0
        for u in nodes:
            if (Av[ni[u]]//4)&1: obstruction ^= coords[u]
        coeff = binary_solve(Q, [(obstruction >> i)&1 for i in range(dim)], dim)
        w = [(coords[u] & coeff).bit_count() & 1 for u in nodes]
        v = [a+2*b for a, b in zip(v, w)]; Av = laplacian(v)
        need(all(a % 4 == 0 for a in Av), 'eighth lift intermediate integrality')
        need(all((Av[ni[u]]//4)%2 == 0 for u in isolated), 'eighth lift isolate solvability')
        y = binary_solve(reduced, [(Av[ni[u]]//4)&1 for u in blue_nodes[1:]], 175) << 1
        for u in blue_nodes: v[ni[u]] += 4*((y >> bi[u])&1)
        Av = laplacian(v)
        need(all(a % 8 == 0 for a in Av), 'eighth-energy lift')
        eighth_lifts.append(v); eighth_images.append(Av)
    eighth = [[(sum(a*b for a, b in zip(v, Aw))//8)&1 for Aw in eighth_images] for v in eighth_lifts]
    eighth_rank = binary_rank([sum(b << j for j, b in enumerate(row)) for row in eighth])
    stages = {2: (lifts, quarter_rows), 3: (eighth_lifts, [sum(b << j for j, b in enumerate(row)) for row in eighth])}
    higher = []
    for k in range(3, 13):
        current_lifts, current_form = stages[k]
        surviving, _ = binary_kernel(current_form, len(current_lifts))
        if not surviving: break
        next_lifts = []; next_images = []
        for r in surviving:
            v = [sum(current_lifts[j][i] for j in bits(r)) for i in range(len(nodes))]
            Av = laplacian(v)
            # Cancel obstructions against the nested surviving spaces from
            # deepest to shallowest, keeping A*v divisible by 2^k throughout.
            for j in range(k-1, 1, -1):
                old_lifts, old_form = stages[j]
                obstruction = [(sum(a*b for a, b in zip(w, Av))//(1 << k)) & 1 for w in old_lifts]
                coeff = binary_solve(old_form, obstruction, len(old_lifts))
                for i in range(len(nodes)):
                    v[i] += (1 << (k-j))*sum(old_lifts[t][i] for t in bits(coeff))
                Av = laplacian(v)
                need(all(a % (1 << k) == 0 for a in Av), 'nested lift divisibility')
            obstruction = 0
            for u in nodes:
                if (Av[ni[u]]//(1 << k)) & 1: obstruction ^= coords[u]
            coeff = binary_solve(Q, [(obstruction >> i)&1 for i in range(dim)], dim)
            w = [(coords[u] & coeff).bit_count() & 1 for u in nodes]
            v = [a+(1 << (k-1))*b for a, b in zip(v, w)]; Av = laplacian(v)
            need(all(a % (1 << k) == 0 for a in Av), 'nested half-form correction')
            need(all((Av[ni[u]]//(1 << k))%2 == 0 for u in isolated), 'nested isolate lift')
            y = binary_solve(reduced, [(Av[ni[u]]//(1 << k))&1 for u in blue_nodes[1:]], 175) << 1
            for u in blue_nodes: v[ni[u]] += (1 << k)*((y >> bi[u])&1)
            Av = laplacian(v)
            need(all(a % (1 << (k+1)) == 0 for a in Av), 'nested full lift equation')
            next_lifts.append(v); next_images.append(Av)
        M = [[(sum(a*b for a, b in zip(v, Aw))//(1 << (k+1)))&1 for Aw in next_images] for v in next_lifts]
        mrows = [sum(b << j for j, b in enumerate(row)) for row in M]
        stages[k+1] = next_lifts, mrows
        higher.append({'valuation': k+1, 'form': M, 'rank': binary_rank(mrows), 'dimension': len(M)})
    # Direct integer energy on a deterministic collection of basis and mixed vectors.
    sample = [(0, 0), (0, 1), (97, 97), (97, 194)]
    sample += [(len(isolated)+i, len(isolated)+j) for i, j in ((0, 0), (0, 1), (3, 7), (10, 10))]
    sample += [(97*i, len(isolated)+i) for i in (1, 3, 7, 10)]
    for i, j in sample:
        total = sum(w*(((coords[u]>>i)&1)-((coords[v]>>i)&1))*(((coords[u]>>j)&1)-((coords[v]>>j)&1)) for u, v, w in full_edges)
        need(total % 2 == 0 and ((total//2)&1) == ((Q[i]>>j)&1), 'direct integer half energy')
    fake = counterfeits[0]['merged_edits']
    deltaQ = [0]*dim
    for u, v, blue, w in fake:
        need(not blue and w % 2 == 0, 'black counterfeit only')
        if (w//2) & 1:
            z = coords[u]^coords[v]
            for i in bits(z): deltaQ[i] ^= z
    doubled_ok = all(ports[tuple(e[:-1])]+2*e[-1] >= 0 for e in counterfeits[0]['port_edits'])
    need(doubled_ok, 'doubled native-capacity counterfeit feasibility')
    quarter_delta = [[sum((w//2)*(a[ni[u]]-a[ni[v]])*(b[ni[u]]-b[ni[v]]) for u, v, blue, w in fake) & 1
                      for b in lifts] for a in lifts]
    out['electrical'] = {'merged_nodes': len(nodes), 'full_components': len(cc), 'paired_black_components': len(bc),
                         'paired_black_component_sizes': list(map(len, bc)), 'blue_components': len(lc), 'blue_nodes': len(blue_nodes),
                         'blue_bicycle_dimension': len(blue_kernel), 'critical_group_2_rank': dim,
                         'half_energy_rank': rankQ, 'smith_factors_divisible_by_4': dim-rankQ,
                         'quarter_energy_matrix': quarter, 'quarter_energy_rank': quarter_rank,
                         'smith_factors_divisible_by_8': dim-rankQ-quarter_rank,
                         'eighth_energy_matrix': eighth, 'eighth_energy_rank': eighth_rank,
                         'smith_factors_divisible_by_16': dim-rankQ-quarter_rank-eighth_rank,
                         'higher_dyadic_layers': higher,
                         'blue_half_rank': binary_rank([sum(v<<j for j, v in enumerate(row)) for row in blue_half]),
                         'nonzero_half_diagonal_entries': sum((r>>i)&1 for i, r in enumerate(Q)),
                         'counterfeit_half_difference_rank': binary_rank(deltaQ),
                         'doubled_counterfeit_feasible': doubled_ok, 'doubled_counterfeit_half_difference_rank': 0,
                         'doubled_counterfeit_quarter_difference_rank': binary_rank([sum(v << j for j, v in enumerate(row)) for row in quarter_delta]),
                         'literal_blue_lines': sum(w for u, v, w in blue_edges),
                         'literal_black_lines': sum(2*w for u, v, w in black_edges),
                         'blue_loops': sum(w for u, v, w in blue_edges if u == v),
                         'literal_black_loops': sum(2*w for u, v, w in black_edges if u == v)}
    exact_orders = {1: rankQ, 2: quarter_rank, 3: eighth_rank}
    exact_orders.update({r['valuation']: r['rank'] for r in higher})
    last_dimension = higher[-1]['dimension'] if higher else len(eighth)
    last_rank = higher[-1]['rank'] if higher else eighth_rank
    out['electrical']['filtration_terminated'] = last_rank == last_dimension
    out['electrical']['exact_2_primary_exponents_counts'] = [[j, r] for j, r in sorted(exact_orders.items()) if r]
    out['electrical']['tree_count_2_adic_valuation'] = sum(j*r for j, r in exact_orders.items()) if last_rank == last_dimension else None
    need(binary_rank(deltaQ) == 2, 'new counterfeit must be detected by half energy')
    out['checks'] = CHECKS
    print('electrical:', json.dumps(out['electrical']), flush=True)
    dest = HERE/'receipt.json'
    if args.write:
        dest.write_text(json.dumps(out, indent=2)+'\n')
        # Inspectable complete merged wiring; deterministic gzip metadata.
        table = 'source_class_min\ttarget_class_min\tblue\tliteral_lines\n'+''.join(f'{u}\t{v}\t{b}\t{w}\n' for (u, v, b), w in sorted(merged.items()))
        (HERE/'n8_merged_wiring.tsv.gz').write_bytes(gzip.compress(table.encode(), mtime=0))
    else: need(json.loads(dest.read_text()) == out, 'frozen receipt differs')
    print('receipt', 'written' if args.write else 'matched', flush=True)

if __name__ == '__main__': main()
