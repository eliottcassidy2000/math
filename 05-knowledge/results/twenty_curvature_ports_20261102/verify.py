"""Exact, finite audits of the separately proved carrier and geometry lemmas.

Run from any directory; normal and -O executions must agree. No asserts.
"""
from collections import Counter
from fractions import Fraction as F
from functools import lru_cache
from itertools import combinations, product
from math import comb, gcd, isqrt
from pathlib import Path
import argparse
import hashlib
import importlib.util
import json

HERE = Path(__file__).resolve().parent


def need(ok, why):
    if not ok:
        raise RuntimeError(why)


def module(name, relative):
    spec = importlib.util.spec_from_file_location(name, HERE.parent / relative)
    obj = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(obj)
    return obj


K = module('inherited_path_kernel', 'joint_incidence_clocked_frames_20261102/path_kernel.py')


def subsets(mask):
    sub = mask
    while True:
        yield sub
        if sub == 0:
            return
        sub = (sub-1) & mask


def quotient(states, observations, successors):
    """Coarsest deterministic observation congruence, with named partial actions."""
    def classify(keys):
        table = {}
        return [table.setdefault(key, len(table)) for key in keys]
    labels = classify(observations)
    rounds = [len(set(labels))]
    while True:
        new = classify((observations[i], tuple((a, labels[j]) for a, j in successors[i]))
                       for i in range(len(states)))
        rounds.append(len(set(new)))
        if new == labels:
            break
        labels = new
    # Stable equivalence explicitly commutes with every legal action.
    class_signature = {}
    for i, label in enumerate(labels):
        key = (observations[i], tuple((a, labels[j]) for a, j in successors[i]))
        need(label not in class_signature or class_signature[label] == key, 'closed quotient')
        class_signature[label] = key
    return labels, rounds


def return_recurrence(labels, successors, maximum_degree=24):
    """Prove a return-count recurrence by annihilating the whole Krylov vector."""
    size = max(labels)+1
    weights = Counter(labels)
    matrix = [None]*size
    for i, row in enumerate(successors):
        value = Counter(labels[j] for _,j in row)
        if matrix[labels[i]] is None:
            matrix[labels[i]] = value
        need(matrix[labels[i]] == value, 'command matrix well-defined')
    for i,row in enumerate(matrix):
        for j,m in row.items():
            need(matrix[j][i] == m and weights[i] == weights[j], 'symmetric involutive command matrix')
    def multiply(v):return [sum(m*v[j] for j,m in row.items()) for row in matrix]
    basis=[]
    v=[1]*size
    sequence=[]
    polynomial=None
    for degree in range(min(size,maximum_degree)+1):
        sequence.append(sum(weights[i]*v[i] for i in range(size)))
        w=list(map(F,v)); poly=[F(0)]*degree+[F(1)]
        for pivot,vector,previous in basis:
            factor=w[pivot]
            if factor:
                w=[a-factor*b for a,b in zip(w,vector)]
                for j,p in enumerate(previous):poly[j]-=factor*p
        pivot=next((i for i,x in enumerate(w) if x),None)
        if pivot is None:
            polynomial=poly
            break
        scale=w[pivot]
        basis.append((pivot,[x/scale for x in w],[x/scale for x in poly]))
        v=multiply(v)
    need(polynomial is not None, 'bounded Krylov search found an exact annihilator')
    # Independent literal state-walk counts through six commands.
    literal=[1]*len(labels)
    for t in range(7):
        if t < len(sequence):need(sum(literal)==sequence[t], 'literal return-word total')
        literal=[sum(literal[j] for _,j in row) for row in successors]
    while len(sequence)<16:
        v=multiply(v)
        sequence.append(sum(weights[i]*v[i] for i in range(size)))
    d=len(polynomial)-1
    need(all(sum(polynomial[j]*sequence[t+j] for j in range(d+1))==0 for t in range(len(sequence)-d)), 'exact scalar recurrence')
    return {'active_word_counts':sequence,'minimal_recurrence_order':d,
            'annihilator_low_to_high':[str(x) for x in polynomial]}


def window_model(sizes):
    offsets = [0]
    for size in sizes:
        offsets.append(offsets[-1]+size)
    blocks = [((1 << size)-1) << start for size, start in zip(sizes, offsets)]
    A, B, C, D = blocks
    windows = [A | B, B | C, C | D]
    states = sorted(set().union(*(set(subsets(w)) for w in windows)))
    index = {s: i for i, s in enumerate(states)}
    obs = [tuple(i for i, w in enumerate(windows) if s & ~w == 0) for s in states]
    successors = [[(i, index[windows[i] ^ s]) for i in obs[j]] for j, s in enumerate(states)]
    labels, rounds = quotient(states, obs, successors)
    a, b, c, d = [int(size >= 2) for size in sizes]
    expected = 8+a+d+3*b+3*c+b*c
    need(len(set(labels)) == expected, 'three-window minimum formula')
    weights = Counter(labels)
    A_size,B_size,C_size,D_size = sizes
    b_mass,c_mass = (1 << B_size)-2,(1 << C_size)-2
    constant = ((1 << A_size)-2)*(1 << B_size)+((1 << D_size)-2)*(1 << C_size)+b_mass*c_mass
    fib=[0,1]
    for _ in range(20):fib.append(fib[-1]+fib[-2])
    walk=[1]*len(states)
    counts=[]
    for t in range(16):
        total=sum(walk)
        formula=constant+(b_mass+c_mass)*(2**((t+1)//2)+2*2**(t//2))+3**((t+1)//2)+3*3**(t//2)+2*fib[t+3]
        need(total==formula,'all-size three-window weighted walk formula')
        counts.append(total)
        walk=[sum(walk[j] for _,j in row) for row in successors]
    return {'block_sizes': sizes, 'active_sources': len(states), 'minimum': expected,
            'refinement_rounds': rounds, 'class_weight_histogram': dict(sorted(Counter(weights.values()).items())),
            'active_word_counts':counts,'minimal_count_recurrence_order':4+2*int(b_mass+c_mass>0)+int(constant>0)}


def count_paths_dp(n):
    """Independent endpoint/subset recurrence, without the inherited leaf filter."""
    adj = [tuple(j for j in range(n) if i != j and isqrt(i+j+2)**2 == i+j+2) for i in range(n)]
    full = (1 << n)-1
    @lru_cache(None)
    def count(last, seen):
        if seen == full:
            return 1
        return sum(count(j, seen | (1 << j)) for j in adj[last] if not (seen >> j) & 1)
    return sum(count(i, 1 << i) for i in range(n)), count.cache_info().misses


def atom_component(pats, edge_count):
    """Independent signed-incidence compiler; masks are not enumerated here."""
    groups = {}
    for e in range(edge_count):
        column = tuple(1 if p >> e & 1 else 0 if q >> e & 1 else -1 for p, q in pats)
        if all(v == -1 for v in column):
            continue
        first = next(v for v in column if v != -1)
        normalized = tuple(v ^ first if v != -1 else -1 for v in column)
        groups.setdefault(normalized, []).append(e)
    profiles = sorted(groups)
    sizes = [len(groups[p]) for p in profiles]
    guards = [tuple(p[i] for p in profiles) for i in range(len(pats))]
    active = set()
    for guard in guards:
        choices = [(v,) if v != -1 else tuple(range(2+(size >= 2))) for v, size in zip(guard, sizes)]
        active.update(product(*choices))
    states = sorted(active)
    index = {s:i for i,s in enumerate(states)}
    obs = [tuple(i for i,g in enumerate(guards) if all(v == -1 or v == s[j] for j,v in enumerate(g))) for s in states]
    successors = []
    for s, legal in zip(states, obs):
        row = []
        for i in legal:
            target = tuple((1-v if v != 2 else 2) if g == -1 else v for v,g in zip(s,guards[i]))
            row.append((i,index[target]))
        successors.append(row)
    labels, rounds = quotient(states, obs, successors)
    outside = edge_count-sum(sizes)
    weights = Counter()
    for s,label in zip(states,labels):
        mass = 1 << outside
        for v,size in zip(s,sizes):
            if v == 2:
                mass *= (1 << size)-2
        weights[label] += mass
    return {'named_paths': len(pats), 'atoms': sizes, 'compressed_active':len(states),
            'minimum_active':len(weights), 'source_mass':sum(weights.values()),
            'class_weight_histogram':dict(sorted(Counter(weights.values()).items())), 'refinement_rounds':rounds}


def path_carrier(n, edges, paths, pats):
    full = (1 << len(edges))-1
    states = sorted(set().union(*(set(p | s for s in subsets(full ^ (p | q))) for p, q in pats)))
    index = {s: i for i, s in enumerate(states)}
    obs = [tuple(j for j, (p, q) in enumerate(pats) if s & p == p and s & q == 0) for s in states]
    successors = [[(j, index[s ^ (full ^ (pats[j][0] | pats[j][1]))]) for j in obs[i]]
                  for i, s in enumerate(states)]
    labels, rounds = quotient(states, obs, successors)
    weights = Counter(labels)
    pairs = [(i, j) for i, j in combinations(range(len(pats)), 2)
             if not (pats[i][0] & pats[j][1] or pats[i][1] & pats[j][0])]
    unseen = set(range(len(paths)))
    components = []
    while unseen:
        part = {min(unseen)}
        while True:
            grown = part | {j for i, j in pairs if i in part} | {i for i, j in pairs if j in part}
            if grown == part:
                break
            part = grown
        unseen -= part
        components.append(sorted(part))
    triangles = []
    predicted_active = 0
    atom_models = []
    atom_weights = Counter()
    for part in components:
        atom = atom_component([pats[i] for i in part], len(edges))
        atom_models.append(atom)
        atom_weights.update(atom['class_weight_histogram'])
        if len(part) == 1:
            predicted_active += 1
        elif len(part) == 2:
            i, j = part
            d = ((pats[i][0] | pats[i][1]) ^ (pats[j][0] | pats[j][1])).bit_count() // 2
            predicted_active += 3+2*(d >= 2)
        elif n == 26 and len(part) == 3 and all((i,j) in pairs for i,j in combinations(part,2)):
            need(len(part) == 3 and all((i, j) in pairs for i, j in combinations(part, 2)), 'triangle component')
            edge_sets = [p | q for p, q in (pats[i] for i in part)]
            U = edge_sets[0] | edge_sets[1] | edge_sets[2]
            windows = [U ^ e for e in edge_sets]
            P, Q, R = windows
            need(P & R == 0 and Q == (P & Q) | (Q & R), 'three-window arrangement')
            sizes = [(P & ~Q).bit_count(), (P & Q).bit_count(), (Q & R).bit_count(), (R & ~Q).bit_count()]
            need(sizes == [2]*4, 'order26 block sizes')
            predicted_active += window_model(sizes)['minimum']
            triangles.append({'path_ids': part, 'paths': [[v+1 for v in paths[i]] for i in part],
                              'window_masks': windows, 'union_edges': U.bit_count(), 'block_sizes': sizes})
        else:
            predicted_active += atom['minimum_active']
    need(predicted_active == len(weights), 'independent component minimum')
    need(sum(a['minimum_active'] for a in atom_models) == len(weights), 'atom minimum agreement')
    need(atom_weights == Counter(weights.values()), 'atom source multiplicities agree')
    need(len(states) < 1 << len(edges), 'dead source exists')
    if n <= 26:
        need(len(weights)+1 == {25: 33, 26: 71}[n], 'global carrier size')
    digest = hashlib.sha256(json.dumps(list(zip(states, labels)), separators=(',', ':')).encode()).hexdigest()
    result = {'n': n, 'active_sources': len(states), 'minimum_including_dead': len(weights)+1,
            'refinement_rounds_active': rounds, 'component_size_histogram': dict(sorted(Counter(map(len, components)).items())),
            'class_weight_histogram_active': dict(sorted(Counter(weights.values()).items())),
            'source_partition_sha256': digest, 'triangles': triangles, 'atom_models':atom_models}
    if n <= 26:
        result['return_words'] = return_recurrence(labels,successors)
        fib=[0,1]
        for _ in range(20):fib.append(fib[-1]+fib[-2])
        for t,total in enumerate(result['return_words']['active_word_counts']):
            A=2**((t+1)//2)+2*2**(t//2)
            formula=4096+256*A if n==25 else 8320+800*A+64*(3**((t+1)//2)+3*3**(t//2))+128*fib[t+3]
            need(total==formula,'actual guard all-length return formula')
    return result


def path_tests():
    rows, carriers = [], []
    for n in range(25, 31):
        edges, paths = K.graph_paths(n, True)
        pats = K.patterns(edges, paths)
        moments, visits = K.moment_spectrum(pats, len(edges))
        shannon, calls = K.shannon_spectrum(pats, len(edges))
        need(moments == shannon, 'independent weighted spectra')
        need(sum(shannon.values()) == 1 << len(edges), 'all guard orientations counted')
        need(sum(k*v for k, v in shannon.items()) == len(paths)*2**(len(edges)-n+1), 'first moment')
        row = {'n': n, 'edges': len(edges), 'directed_paths': len(paths), 'spectrum': shannon,
               'outside_free_bits': comb(n, 2)-len(edges), 'clique_nodes': visits, 'shannon_states': calls}
        content=0
        for value in shannon.values():content=gcd(content,value)
        need(content & (content-1) == 0,'binary coefficient content')
        row['full_content_exponent']=row['outside_free_bits']+content.bit_length()-1
        row['primitive_coefficients']={k:v//content for k,v in shannon.items()}
        if n in (25,26):need(row['full_content_exponent']=={25:276,26:296}[n],'incoming primitive normalization')
        if n <= 26:
            count, dp_states = count_paths_dp(n)
            need(count == len(paths), 'independent path enumeration')
            row['independent_dp_states'] = dp_states
        if n <= 27:
            carriers.append(path_carrier(n, edges, paths, pats))
        rows.append(row)
    windows = [window_model(list(sizes)) for sizes in product((1, 2), repeat=4)]
    windows.append(window_model([2, 3, 4, 2]))
    return {'spectra': rows, 'carriers': carriers, 'three_window_models': windows}


def points(base, depth):
    out = {(): (F(1, 2), F(1))}
    for d in range(depth):
        for address, (x, y) in list(out.items()):
            if len(address) == d:
                for bit in (0, 1):
                    out[address+(bit,)] = (x+(2*bit-1)*y*(base-1)/(2*base), y/base)
    return out


def cosh_distance(p, q):
    x, y = p
    X, Y = q
    return 1+((x-X)**2+(y-Y)**2)/(2*y*Y)


def lca(a, b):
    d = 0
    while d < min(len(a), len(b)) and a[d] == b[d]:
        d += 1
    return a[:d]


def cosh_log(q):
    return (q+1/q)/2


def geometry_tests():
    rows = []
    for base in (2, 3, 4):
        pts = points(base, 6)
        edge = F(5*base*base-2*base+5, 8*base)
        chebyshev = [F(1), edge]
        for d in range(2, 13):
            chebyshev.append(2*edge*chebyshev[-1]-chebyshev[-2])
        pairs = 0
        for a, b in combinations(pts, 2):
            q = lca(a, b)
            D = len(a)+len(b)-2*len(q)
            cd = cosh_distance(pts[a], pts[b])
            need(cd >= edge and (cd == edge) == (D == 1), 'induced fixed-distance tree')
            need(cd <= chebyshev[D], 'triangle upper bound')
            need(cd <= cosh_log(F(25,16)*base**D), 'additive LCA upper bound')
            if base > 2:
                # exp(D log b + 2 log g), positive lower bound when >1.
                exp_lower = F(base**D)*(1-F(2, base))**2
                if exp_lower >= 1:
                    need(cd >= cosh_log(exp_lower), 'gapped metric lower bound')
            pairs += 1
        for a in pts:
            for cut in range(len(a)+1):
                q = a[:cut]
                D = len(a)-cut
                cd = cosh_distance(pts[a], pts[q])
                need(cosh_log(F(base**D)) <= cd <= cosh_log(F(5, 4)*base**D), 'ancestor bounds')
        central = []
        for h in range(1, 21):
            cd = 1+(F(base-2, base)*base**h+1)**2/2
            if h <= 6:
                need(cd == cosh_distance(pts[(0,)+(1,)*(h-1)], pts[(1,)+(0,)*(h-1)]), 'central pair formula')
            central.append({'depth': h, 'tree_distance': 2*h, 'cosh_distance': str(cd)})
        rows.append({'base': base, 'vertices': len(pts), 'pairs': pairs, 'cosh_edge_length': str(edge), 'central_pairs': central})
    need(11**2+8**2 == 13**2+4**2 == 185, 'dyadic geodesic circle')
    need(7**2+6**2 == 9**2+2**2 == 85, 'ternary geodesic circle')
    need(F(5,4)**2+F(3,4)**2 == F(17,8) and 2*F(5,4)*F(3,4) == F(15,8), 'boost composition')
    return rows


def collatz_tests():
    def C(n):
        return n//2 if n % 2 == 0 else 3*n+1
    def children(n):
        exponents = (2,4,6) if n % 3 == 1 else (3,5,7)
        result = [((2**a*n-1)//3, a) for a in exponents if (2**a*n-1) % 9]
        need(len(result) == 2, 'two guarded children')
        return result
    nodes = {(): 5}
    clocks = {}
    for h in range(6):
        for address in [a for a in nodes if len(a) == h]:
            for bit, (n, a) in enumerate(children(nodes[address])):
                nodes[address+(bit,)] = n
                clocks[address+(bit,)] = a+1
    need(len(set(nodes.values())) == len(nodes), 'injective selected source addresses')
    paths = {}
    for address, n in nodes.items():
        path = [n]
        while path[-1] != 5:
            need(len(path) <= 8*len(address), 'stored finite return bound')
            path.append(C(path[-1]))
        paths[address] = path
    pts = points(3, 6)
    edge = F(11,6)
    cheb = [F(1), edge]
    for i in range(40):
        cheb.append(2*edge*cheb[-1]-cheb[-2])
    pairs = 0
    for a, b in combinations(nodes, 2):
        q = lca(a, b)
        D = len(a)+len(b)-2*len(q)
        pb = {n: i for i, n in enumerate(paths[b])}
        ia, meet = next((i,n) for i,n in enumerate(paths[a]) if n in pb)
        native = ia+pb[meet]
        need(3*D-2 <= native <= 8*D, 'native versus selected metric')
        wa = sum(clocks[a[:i]] for i in range(len(q)+1, len(a)+1))
        wb = sum(clocks[b[:i]] for i in range(len(q)+1, len(b)+1))
        shift = 0 if q in (a,b) else min(clocks[a[:len(q)+1]], clocks[b[:len(q)+1]])-1
        need(native == wa+wb-2*shift and meet == 2**shift*nodes[q], 'literal earlier even meeting')
        pairs += 1
    seen, frontier, counts = {5}, {5}, [1]
    for depth in range(1, 27):
        new = {2*n for n in frontier}
        new |= {(n-1)//3 for n in frontier if n > 1 and n % 6 == 4}
        need(not (seen & new), 'acyclic root5 inverse levels')
        need(all(C(n) in frontier for n in new), 'native inverse edge')
        seen |= new
        frontier = new
        counts.append(len(seen))
    need(counts[21] == 704 and 704-(1+3*21*22//2) == 10, 'correct honeycomb deficit')
    need(counts[26] == 2263 and counts[26] > 1+3*26*27, 'triangular obstruction')
    def G(n, sign):
        residues = (4,7) if sign == -1 else (2,5)
        k = next(k for k in range(6) if 2**k*n % 9 in residues)
        return (2**k*n+sign)//3, k
    need(G(4,1) == (11,3) and G(11,1) == (4,0), 'signed greedy cycle with clocks')
    need(G(13,-1) == (4,0) and G(4,-1) == (1,0) and C(4) == 2, 'distinct maps')
    return {'selected_vertices': len(nodes), 'native_pairs': pairs, 'inverse_root5_counts': counts,
            'honeycomb_deficit': 10, 'signed_greedy_cycle': [[4,11,3],[11,4,0]],
            'positive_greedy_edges': [[13,4,0],[4,1,0]]}


def arithmetic_tests():
    bound = 10000
    sums = [0]*(bound+1)
    tau = [1]*(bound+1)
    for d in range(1, bound//2+1):
        for n in range(2*d, bound+1, d):
            sums[n] += d
            tau[n] += 1
    five_divisor_overlap = []
    for n in range(1, bound+1):
        pairs = [(a,n//a) for a in range(1,isqrt(n)+1) if n % a == 0]
        actual = any(sums[n] == 1+a+b for a,b in pairs)
        need(actual == (tau[n] == 4 or sums[n] == n+2), 'inherited all-size criterion finite hostile audit')
        if actual and tau[n] == 6:
            five_divisor_overlap.append(n)
    need(five_divisor_overlap == [20], 'unique five-proper-divisor overlap')
    sigma42=[n for n in range(1,42) if sums[n]+n==42]
    need(sigma42==[20,26,41],'complete sigma42 fibre; n>1 implies sigma(n)>=n+1')
    need(sums[16] == 15 and not any(1+a+16//a == 15 for a in (1,2,4)), '16 excluded')
    rows = [{'n': n, 's': sums[n], 'proper_divisors': [d for d in range(1,n) if n%d == 0]}
            for n in (8,12,16,20,24,28,44,52)]
    factor_rows = []
    for n in (20,44,52):
        pairs = [(a,n//a) for a in range(1,isqrt(n)+1) if n%a == 0]
        for a,b in pairs:
            need((a+b)**2-(b-a)**2 == 4*n, 'factor-pair Lorentz norm')
        factor_rows.append({'n': n, 'factor_pairs': pairs, 'lorentz_coordinates': [(a+b,b-a) for a,b in pairs]})
    need(pow(2,3,7) == 1 and {pow(2,i,7) for i in range(3)} == {1,2,4}, 'Frobenius exponent orbit')
    trace_rows = []
    for a, modulus in ((5,13),(2,11),(5,29)):
        def mul(x,y):
            r,s=x; t,u=y
            return ((r*t+s*u)%modulus,(r*u+s*t+a*s*u)%modulus)
        power=(1,0)
        for _ in range(12):power=mul(power,(0,1))
        def trace(x):return (2*x[0]+a*x[1])%modulus
        domain=list(product(range(modulus),repeat=2))
        primitive={(trace(x),trace(mul((0,1),x))) for x in domain}
        if modulus in (11,13):
            need(power==((1 if modulus==13 else modulus-1),0), 'twelfth clock scalar')
            sampled={(trace(x),trace(mul(power,x))) for x in domain}
            need(len(sampled)==modulus and len(primitive)==modulus**2, 'primitive clock repair')
        else:
            sampled=set()
            need(len(primitive)==29, 'ramified trace radical hostile')
        trace_rows.append({'a':a,'modulus':modulus,'u_power12':power,'primitive_pairs':len(primitive),
                           'twelfth_pairs':len(sampled) if sampled else None})
    return {'bound': bound, 'five_divisor_overlap': five_divisor_overlap, 'rows': rows,
            'factor_pairs': factor_rows, 'recovered_trace_clocks':trace_rows,'sigma42_fibre':sigma42}


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument('--output', default='verification.json')
    args = parser.parse_args()
    hostile = module('inherited_symmetry', 'regular_parent_symmetry_20261029/verify.py')
    result = {'paths': path_tests(), 'geometry': geometry_tests(), 'collatz': collatz_tests(),
              'arithmetic': arithmetic_tests(), 'hostile_81_27': hostile.hostile_checks()}
    target = HERE / args.output
    target.write_text(json.dumps(result, indent=2, sort_keys=True)+'\n', encoding='utf-8')
    print(json.dumps({'output': target.name, 'carriers': [r['minimum_including_dead'] for r in result['paths']['carriers']],
                      'geometry_pairs': sum(r['pairs'] for r in result['geometry']),
                      'native_pairs': result['collatz']['native_pairs'], 'hostile': result['hostile_81_27']}, sort_keys=True))


if __name__ == '__main__':
    main()
