"""Complete native invariant-form census and certified clock-isometry reduction."""
from collections import Counter
from itertools import combinations, product
from pathlib import Path
import hashlib
import json

from primary_native_clocks import (alternating_matrix, dual_map, identity,
                                  invariant_form_space, modular_kernel, mul, row_reduce)
from three_run_glue import rref_mod, hostile_singleton_slices
from paley_pencil import diagonal, local_layers

HERE = Path(__file__).resolve().parent
GATES = 0


def need(ok, why):
    global GATES
    GATES += 1
    if not ok:
        raise ValueError(why)


def valuation(a, p):
    need(a != 0, 'finite radical Smith factor')
    v = 0
    while a % p == 0:
        a //= p
        v += 1
    return v


def plane_actions():
    """All native-clock actions on the necessary extra symplectic F3 plane."""
    def mm(a,b):
        return [[sum(x*y for x,y in zip(row,col)) % 3 for col in zip(*b)] for row in a]
    I = identity(2)
    mats = [[[a,b],[c,d]] for a,b,c,d in product(range(3),repeat=4) if (a*d-b*c) % 3 == 1]
    need(len(mats) == 24, 'complete SL2(F3) enumeration')
    cubes = [A for A in mats if mm(mm(A,A),A) == I]
    involutions = [A for A in mats if mm(A,A) == I]
    need(len(cubes) == 9 and involutions == [I,[[2,0],[0,2]]], 'only scalar symplectic involutions')
    def key(F,R,S):
        return tuple(x for A in (F,R,S) for row in A for x in row)
    valid = {}
    for F,R,S in product(cubes,cubes,involutions):
        if (mm(F,R) == mm(R,F) and mm(S,F) == mm(F,S)
                and mm(mm(S,R),S) == mm(R,R)):
            need(R == I, 'run cycle is trivial on every symplectic defect plane')
            valid[key(F,R,S)] = (F,R,S)
    need(len(valid) == 18, 'all literal native-clock plane actions')
    remaining = set(valid)
    orbits = []
    while remaining:
        rep = min(remaining)
        orbit = set()
        for A in mats:
            Ai = [[A[1][1],-A[0][1] % 3],[-A[1][0] % 3,A[0][0]]]
            orbit.add(key(*(mm(mm(A,C),Ai) for C in valid[rep])))
        need(orbit <= remaining, 'disjoint complete symplectic plane conjugacy classes')
        remaining -= orbit
        orbits.append({'representative':valid[rep],'size':len(orbit)})
    need(len(orbits) == 6 and sorted(o['size'] for o in orbits) == [1,1,4,4,4,4], 'six complete plane-action types')
    J = [[0,1],[2,0]]
    reflection = [[1,0],[0,2]]
    shear = [[1,1],[0,1]]
    need(mm(mm(reflection,shear),reflection) == mm(shear,shear), 'hostile has the correct group relation')
    need(mm(mm(reflection,J),reflection) != J, 'hostile reflection does not preserve the alternating form')
    return {'all_symplectic_matrices':24,'valid_actions':18,'conjugacy_classes':orbits,
            'hostile':{'run_cycle':shear,'run_swap':reflection,'failure':'run swap is anti-symplectic'}}


def run():
    src = json.loads((HERE/'primary_native_clocks.json').read_text())
    case = src['primary_cases'][1]
    es, ds, clocks = case['exponents'], case['orders'], case['clocks']
    n, q, E = 9, 81, 4
    slots = list(combinations(range(n), 2))
    basis, info = invariant_form_space(es, clocks, 3, True)
    need(basis == src['three_primary_invariant_radical']['complete_lift_basis'], 'unchanged complete source form module')
    # Uniform coordinates have a precisely known kernel of representation.
    kernel_log = sum(info['kernel_pivot_valuations'])+(len(slots)-len(info['kernel_pivot_valuations']))*E
    redundant_log = sum(E-min(es[i], es[j]) for i, j in slots)
    need(kernel_log-redundant_log == 6, 'exact invariant-form module has order3^6')
    flattened = []
    for v in basis:
        J = alternating_matrix(v, es, 3)
        flattened.append(tuple(J[i][j] for i, j in slots))
    zero = (0,)*len(slots)
    space, generators, orders = {zero}, [], []
    for g in flattened:
        if g in space:
            continue
        order = next(k for k in (3, 9, 27, 81) if not any(k*x % q for x in g))
        new = {tuple((x+t*y) % q for x, y in zip(v, g)) for v in space for t in range(order)}
        need(len(new) == len(space)*order, 'chosen form generators independent')
        space = new
        generators.append(g)
        orders.append(order)
    need(orders == [3, 3, 3, 3, 9] and len(space) == 729, 'full C3^4 plus C9 parameterization')

    def matrix(v):
        J = [[0]*n for _ in range(n)]
        for x, (i, j) in zip(v, slots):
            J[i][j], J[j][i] = x, -x % q
        return J

    forms, parameters, records = [], {}, []
    distribution = Counter()
    D = [[ds[i]*int(i == j) for j in range(n)] for i in range(n)]
    for co in product(*(range(d) for d in orders)):
        v = tuple(sum(a*g[i] for a, g in zip(co, generators)) % q for i in range(len(slots)))
        need(v not in parameters, 'unique marked form parameters')
        parameters[v] = list(co)
        forms.append(v)
        B = dual_map(matrix(v), es, 3)
        presentation = [B[i]+D[i] for i in range(n)]
        layers = local_layers(presentation, 3, 5)
        smith = diagonal(presentation)
        valuations = [valuation(a, 3) for a in smith]
        need(Counter(valuations) == Counter({i:c for i,c in enumerate(layers) if c}),
             'independent integer Smith and local radical calculation')
        r = sum(valuations)
        distribution[r] += 1
        records.append({'parameters': list(co), 'radical_log_order': r, 'radical_exponents': [v for v in valuations if v]})
    need(dict(distribution) == {18:1,14:70,16:10,12:324,8:324}, 'complete native radical strata')
    need(10 not in distribution, 'no source form with radical order3^10')

    # Generate a certified subgroup of the full clock centralizer.
    scale = [3**max(0, es[i]-es[j]) for i in range(n) for j in range(n)]
    equations = []
    for C in clocks.values():
        for i, j in product(range(n), repeat=2):
            row = [0]*(n*n)
            for b in range(n):
                row[i*n+b] += scale[i*n+b]*C[b][j]
            for a in range(n):
                row[a*n+j] -= C[i][a]*scale[a*n+j]
            equations.append([x*(q//ds[i]) % q for x in row])
    ends, _ = modular_kernel(equations, 3, E)
    units, actions, signatures = [], [], set()
    candidates, residue_basis, residue_flat = [], [], []
    for v in ends:
        M = row_reduce([[scale[i*n+j]*v[i*n+j] for j in range(n)] for i in range(n)], ds)
        flat = [x % 3 for row in M for x in row]
        if len(rref_mod(residue_flat+[flat],3)[1]) > len(residue_flat):
            residue_flat.append(flat)
            residue_basis.append(M)
        for t in (1, 2):
            candidates.append(row_reduce([[int(i == j)+t*M[i][j] for j in range(n)] for i in range(n)], ds))
    need(len(residue_basis) == 6, 'complete six-dimensional residue endomorphism image')
    residue_units = 0
    for co in product(range(3),repeat=len(residue_basis)):
        A = row_reduce([[sum(a*M[i][j] for a,M in zip(co,residue_basis)) for j in range(n)] for i in range(n)],ds)
        if len(rref_mod(A,3)[1]) == n:
            residue_units += 1
            candidates.append(A)
    for A in candidates:
        if len(rref_mod(A, 3)[1]) != n:
            continue
        need(all(ds[j]*A[i][j] % ds[i] == 0 for i, j in product(range(n), repeat=2)), 'heterogeneous source homomorphism')
        need(all(mul(A,C,ds) == mul(C,A,ds) for C in clocks.values()), 'isometry change intertwines all source clocks')
        action = [[(A[a][i]*A[b][j]-A[b][i]*A[a][j]) % q for a,b in slots] for i,j in slots]
        def apply(v):
            return tuple(sum(x*y for x,y in zip(row,v)) % q for row in action)
        signature = tuple(apply(g) for g in generators)
        if signature == tuple(generators) or signature in signatures:
            continue
        signatures.add(signature)
        units.append(A)
        actions.append(action)
    need(bool(units), 'nontrivial source-isometry subgroup')
    record_by_form = dict(zip(forms, records))
    unseen = set(forms)
    orbits = []
    while unseen:
        root = min(unseen)
        queue, seen = [root], {root}
        edges = []
        for v in queue:
            for g, action in enumerate(actions):
                w = tuple(sum(x*y for x,y in zip(row,v)) % q for row in action)
                need(w in parameters, 'automorphism preserves complete invariant-form space')
                need(record_by_form[w]['radical_log_order'] == record_by_form[v]['radical_log_order'], 'radical order constant under source isometry')
                if w not in seen:
                    seen.add(w)
                    queue.append(w)
                    edges.append({'parent': parameters[v], 'generator': g, 'child': parameters[w]})
        need(seen <= unseen, 'isometry orbits disjoint')
        unseen -= seen
        orbits.append({'representative': parameters[root], 'radical_log_order': record_by_form[root]['radical_log_order'],
                       'size': len(seen), 'spanning_tree': edges})
    optimal = [o for o in orbits if o['radical_log_order'] == 8]
    need(sum(o['size'] for o in optimal) == 324, 'every optimal source form covered by certified source isometries')
    planes = plane_actions()
    hostile = hostile_singleton_slices()
    return {'status': 'FINITE-EXACT complete form census / sufficient isometry representatives',
            'source_sha256_canonical': hashlib.sha256(json.dumps(src,sort_keys=True,separators=(',',':')).encode()).hexdigest(),
            'uniform_kernel_log_order': kernel_log, 'coordinate_redundancy_log_order': redundant_log,
            'parameter_orders': orders, 'form_generators_upper_triangle': generators,
            'radical_distribution': dict(sorted(distribution.items())), 'forms': records,
            'clock_centralizer_generators': units, 'certified_subgroup_orbits': orbits,
            'residue_endomorphism_dimension': len(residue_basis), 'all_residue_units_lifted': residue_units,
            'optimal_representative_count': len(optimal), 'optimal_orbit_sizes': sorted(o['size'] for o in optimal),
            'necessary_extra_plane_actions': planes,
            'full_centralizer_orbit_minimality': 'NOT CLAIMED', 'hostile_81_27': hostile, 'gates': GATES}


if __name__ == '__main__':
    receipt = run()
    (HERE/'ternary_form_strata.json').write_text(json.dumps(receipt,indent=2)+'\n')
    print(json.dumps({k:receipt[k] for k in ('gates','radical_distribution','optimal_representative_count','optimal_orbit_sizes')}))
