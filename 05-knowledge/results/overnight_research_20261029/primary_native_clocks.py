"""Induce all native clocks on saved primary coordinates; test actual completions."""
from itertools import combinations, product
from pathlib import Path
from random import Random
import hashlib
import json

from primary_clock_layers import Ring, identity
from three_run_glue import three_runs, mm_mod, rref_mod, hostile_singleton_slices
from paley_pencil import local_layers, diagonal

HERE = Path(__file__).resolve().parent
GATES = 0


def need(ok, why):
    global GATES
    GATES += 1
    if not ok:
        raise ValueError(why)


def transpose(a):
    return list(map(list, zip(*a)))


def row_reduce(a, orders):
    return [[x % q for x in row] for q, row in zip(orders, a)]


def mul(a, b, orders):
    # Every order divides their maximum, so uniform reduction is safe first.
    return row_reduce(mm_mod(a, b, max(orders)), orders)


def power(a, n, orders):
    b = identity(len(a))
    for _ in range(n):
        b = mul(b, a, orders)
    return b


def native_actions(ports):
    ix = {a:i for i, a in enumerate(ports)}
    return {
        'translation': [ix[tuple((x+1) % 7 for x in a)] for a in ports],
        'vertex_frobenius': [ix[tuple(2*x % 7 for x in a)] for a in ports],
        'run_swap': [ix[(a[1], a[0], a[2])] for a in ports],
        'run_cycle': [ix[(a[1], a[2], a[0])] for a in ports]}


def check_relations(clocks, orders):
    T, F, S, R = [clocks[x] for x in ('translation','vertex_frobenius','run_swap','run_cycle')]
    I = identity(len(T))
    for C, n in ((T,7),(F,3),(S,2),(R,3)):
        need(power(C,n,orders) == I, 'actual generator order')
        need(all(orders[j]*C[i][j] % orders[i] == 0
                 for i,j in product(range(len(orders)), repeat=2)), 'well-defined heterogeneous group map')
    need(mul(mul(F,T,orders),power(F,2,orders),orders) == power(T,2,orders), 'affine clock relation')
    need(mul(mul(S,R,orders),S,orders) == power(R,2,orders), 'joint S3 clock relation')
    for A, B in product((T,F),(S,R)):
        need(mul(A,B,orders) == mul(B,A,orders), 'vertex and run clocks commute')


def induce(case, data, M, actions):
    p, k, q = case['prime'], case['precision'], case['modulus']
    orbits = data['translation_orbits']
    F, Fi = case['source_orbit_chart'], case['source_orbit_chart_inverse']
    B, columns, names, exponents = [], [], [], []
    offset = 0
    for sector, record in enumerate(case['factors']):
        ring = Ring(p,k,record['factor_low_coefficients'])
        U, Ui = record['source_row_map_U'], record['source_row_map_inverse']
        for i,e in enumerate(record.get('pivot_valuations',record.get('all_30_pivot_valuations'))):
            if not e:
                continue
            for coord in range(ring.d):
                row, col = [0]*210, [0]*210
                for j, orbit in enumerate(orbits):
                    block, inverse = ring.times(tuple(U[i][j])), ring.times(tuple(Ui[j][i]))
                    for phase, v in enumerate(orbit):
                        row[v] = sum(block[coord][b]*F[offset+b][phase] for b in range(ring.d)) % (p**e)
                        col[v] = sum(Fi[phase][offset+b]*inverse[b][coord] for b in range(ring.d)) % q
                B.append(row)
                columns.append(col)
                names.append({'sector':sector,'smith_row':i,'coefficient':coord})
                exponents.append(e)
        offset += ring.d
    section = transpose(columns)
    orders = [p**e for e in exponents]
    # Integer coordinate lifts, not a homomorphic splitting into the free lattice.
    need(mul(B,section,orders) == identity(len(B)), 'coordinate lifts recover every quotient generator')
    need(mul(B,M,orders) == [[0]*210 for _ in B], 'every native relation is zero in the quotient')
    clocks = {}
    for name, perm in actions.items():
        moved = [[row[j] for j in perm] for row in B]
        C = mul(moved,section,orders)
        need(mul(C,B,orders) == moved, 'clock intertwines every named port at its exact exponent')
        clocks[name] = C
    check_relations(clocks,orders)
    return {'prime':p,'coordinate_names':names,'exponents':exponents,'orders':orders,
            'projection':B,'coordinate_lifts':section,'clocks':clocks}


def modular_kernel(matrix, p, k):
    """Full kernel over Z/(p^k), via reversible rows and retained column map."""
    q = p**k
    a = [[x % q for x in row] for row in matrix]
    nr, nc = len(a), len(a[0])
    V, vals = identity(nc), []
    def val(x):
        for h in range(k):
            if x % p**(h+1):
                return h
        return k
    for r in range(min(nr,nc)):
        v,i,j = min((val(a[i][j]),i,j) for i in range(r,nr) for j in range(r,nc))
        if v == k:
            break
        vals.append(v)
        a[r],a[i] = a[i],a[r]
        for row in a:
            row[r],row[j] = row[j],row[r]
        for row in V:
            row[r],row[j] = row[j],row[r]
        pv = p**v
        inv = pow(a[r][r]//pv,-1,q)
        a[r] = [x*inv % q for x in a[r]]
        for i in range(r+1,nr):
            need(a[i][r] % pv == 0, 'local kernel row division')
            c = a[i][r]//pv
            a[i] = [(x-c*y) % q for x,y in zip(a[i],a[r])]
        for j in range(r+1,nc):
            need(a[r][j] % pv == 0, 'local kernel column division')
            c = a[r][j]//pv
            for row in a:
                row[j] = (row[j]-c*row[r]) % q
            for row in V:
                row[j] = (row[j]-c*row[r]) % q
    need(all(a[i][j] == (p**vals[i] if i == j and i < len(vals) else 0)
             for i in range(nr) for j in range(nc)), 'entire local relation matrix diagonalized')
    need(len(rref_mod(V,p)[1]) == nc, 'kernel column chart invertible')
    basis = []
    for j in range(nc):
        multiplier = p**(k-vals[j]) if j < len(vals) else 1
        col = [row[j]*multiplier % q for row in V]
        if any(col):
            need(all(sum(x*y for x,y in zip(row,col)) % q == 0 for row in matrix), 'literal kernel generator')
            basis.append(col)
    return basis, vals


def alternating_matrix(vector, exponents, p):
    n, E = len(exponents), max(exponents)
    q = p**E
    J = [[0]*n for _ in range(n)]
    for x,(i,j) in zip(vector,combinations(range(n),2)):
        y = x*p**(E-min(exponents[i],exponents[j])) % q
        J[i][j],J[j][i] = y,-y % q
    return J


def dual_map(J, exponents, p):
    E = max(exponents)
    out = []
    for e,row in zip(exponents,J):
        divisor = p**(E-e)
        need(all(x % divisor == 0 for x in row), 'pairing defines a homomorphism to the actual dual')
        out.append([x//divisor for x in row])
    return out


def invariant_form_space(exponents,clocks,p,all_generators=False):
    n,E = len(exponents),max(exponents)
    q,slots = p**E,list(combinations(range(n),2))
    factors = [p**(E-min(exponents[i],exponents[j])) for i,j in slots]
    equations = []
    for C in clocks.values():
        for i,j in slots:
            equations.append([s*(C[a][i]*C[b][j]-C[b][i]*C[a][j]-int((a,b)==(i,j))) % q
                              for s,(a,b) in zip(factors,slots)])
    basis,vals = modular_kernel(equations,p,E)
    # Only the mod-p dual maps affect perfection (finite-module Nakayama).
    projected, chosen, rank = [], [], 0
    for v in basis:
        J = alternating_matrix(v,exponents,p)
        flat = [x % p for row in dual_map(J,exponents,p) for x in row]
        newrank = len(rref_mod(projected+[flat],p)[1])
        if newrank > rank:
            projected.append(flat)
            chosen.append(v)
            rank = newrank
    return (basis if all_generators else chosen), {'kernel_pivot_valuations':vals,'kernel_generator_count':len(basis),
                    'mod_p_dual_map_space_dimension':rank}


def radical_profile(case):
    p,es,clocks = case['prime'],case['exponents'],case['clocks']
    n,E = len(es),max(es)
    basis,info = invariant_form_space(es,clocks,p,True)
    forms = [alternating_matrix(v,es,p) for v in basis]
    maps = [dual_map(J,es,p) for J in forms]
    D = [[p**es[i]*int(i==j) for j in range(n)] for i in range(n)]
    def kernel_log(ms):
        joined = [sum((B[i] for B in ms),[])+D[i] for i in range(n)]
        layers = local_layers(joined,p,E+1)
        need(sum(layers) == n,'complete cokernel image index')
        return sum(e*r for e,r in enumerate(layers)),layers
    common,common_layers = kernel_log(maps)
    joined = [sum((B[i] for B in maps),[])+D[i] for i in range(n)]
    expected = [p**e for e,count in enumerate(common_layers) for _ in range(count)]
    need(diagonal(joined)==expected,'independent integer Smith check of the universal radical')
    rng,best,chosen = Random(113),sum(es),None
    for trial in range(2000):
        cs = [rng.randrange(p**E) for _ in basis]
        v = [sum(c*b[j] for c,b in zip(cs,basis)) % (p**E) for j in range(n*(n-1)//2)]
        J = alternating_matrix(v,es,p)
        value,layers = kernel_log([dual_map(J,es,p)])
        if value < best:
            best,chosen = value,{'coefficients':cs,'form':J,'radical_layers':layers}
        if best == common:
            break
    if chosen:
        for C in clocks.values():
            need(mm_mod(mm_mod(transpose(C),chosen['form'],p**E),C,p**E) == chosen['form'],
                 'literal best restriction form is jointly invariant')
        best_map=dual_map(chosen['form'],es,p)
        need(diagonal([best_map[i]+D[i] for i in range(n)])==
             [p**e for e,count in enumerate(chosen['radical_layers']) for _ in range(count)],
             'independent integer Smith check of the attaining radical')
    return {**info,'complete_lift_basis':basis,'common_radical_log_order':common,
            'common_radical_layers':common_layers,'best_radical_log_order':best,
            'best_form':chosen,'trials':trial+1,'attains_common_radical':best==common}


def find_perfect(exponents,clocks,p,seed=29):
    basis,info = invariant_form_space(exponents,clocks,p)
    rng,n,q = Random(seed),len(exponents),p**max(exponents)
    if len(basis) <= 7:
        candidates = product(range(p),repeat=len(basis))
        exhaustive = True
    else:
        candidates = ([rng.randrange(p) for _ in basis] for _ in range(2000))
        exhaustive = False
    tested,best = 0,0
    for coefficients in candidates:
        tested += 1
        v = [sum(c*b[j] for c,b in zip(coefficients,basis)) % q for j in range(n*(n-1)//2)]
        J = alternating_matrix(v,exponents,p)
        D = dual_map(J,exponents,p)
        rank = len(rref_mod(D,p)[1])
        best = max(rank,best)
        if rank == n:
            for C in clocks.values():
                need(mm_mod(mm_mod(transpose(C),J,q),C,q) == J,'perfect form preserves every native clock')
            need(all(J[i][i] == 0 for i in range(n)),'alternating includes zero diagonal')
            info.update({'tested':tested,'perfect':True,'form':J,'dual_map':D,
                         'chosen_coefficients':list(coefficients),'lift_basis':basis})
            return info
    info.update({'tested':tested,'perfect':False,'best_mod_p_rank':best,
                 'exhaustive_mod_p_form_space':exhaustive,'lift_basis':basis})
    return info


def three_primary_completion(case):
    p,es = case['prime'],case['exponents']
    need(p == 3 and es == [1,1,2,2,2,2,2,2,4],'actual native3-primary ordered group')
    results = []
    for frob,sgn in product((1,28,55),(1,80)):
        characters = dict(translation=1,vertex_frobenius=frob,run_swap=sgn,run_cycle=1)
        clocks = {name:[row+[0] for row in C]+[[0]*len(es)+[characters[name]]]
                  for name,C in case['clocks'].items()}
        check_relations(clocks,[3**e for e in es]+[81])
        form = find_perfect(es+[4],clocks,3)
        results.append({'extra_cyclic_order':81,'character':characters,'clocks':clocks,**form})
        print(json.dumps({'prime':3,'character':characters,'perfect':form['perfect'],
                          'space':form['mod_p_dual_map_space_dimension'],
                          'tested':form['tested'],'best_rank':form.get('best_mod_p_rank')}),flush=True)
    return results


def heisenberg_control():
    out = []
    for p in (2,3,5,11,13):
        J = [[0,1,0,0],[-1,0,0,0],[0,0,0,1],[0,0,-1,0]]
        commute = 0
        for a,b in product(range(p),repeat=2):
            A = [[1,0,0,0],[0,1,0,-1],[1,0,1,a],[0,0,0,1]]
            B = [[1,0,0,1],[0,1,0,0],[0,1,1,b],[0,0,0,1]]
            for C in (A,B):
                need(mm_mod(mm_mod(transpose(C),J,p),C,p) == [[x%p for x in row] for row in J],
                     'each individual hostile clock extends isometrically')
                need(power(C,p,[p]*4) == identity(4),'each hostile clock has the required order')
            AB,BA = mm_mod(A,B,p),mm_mod(B,A,p)
            expected = [[(2 if (i,j)==(2,3) else 0) % p for j in range(4)] for i in range(4)]
            need([[(x-y)%p for x,y in zip(r,s)] for r,s in zip(AB,BA)] == expected,
                 'all possible lifts have the same joint central error')
            commute += AB == BA
        need(commute == (p*p if p == 2 else 0),'sharp odd-prime versus dyadic boundary')
        out.append({'prime':p,'lift_pairs':p*p,'commuting_lift_pairs':commute,
                    'commutator_displacement':2%p})
    return out


def main():
    data = json.loads((HERE/'primary_clock_layers.json').read_text())
    rows = [sum(1 << ((v+s)%7) for s in (1,2,4)) for v in range(7)]
    _,ports,M = three_runs(rows)
    need([list(x) for x in ports] == data['ports'],'saved source labels match native ports')
    actions = native_actions(ports)
    for perm in actions.values():
        need(all(M[perm[i]][perm[j]] == M[i][j] for i,j in product(range(210),repeat=2)),
             'every extra clock preserves the actual native operator')
    cases = [induce(case,data,M,actions) for case in data['primary_cases']]
    print(json.dumps({'induced_dimensions':[len(c['orders']) for c in cases]}),flush=True)
    three = three_primary_completion(cases[1])
    radical = radical_profile(cases[1])
    print(json.dumps({'prime':3,'common_radical':radical['common_radical_log_order'],
                      'best_radical':radical['best_radical_log_order'],
                      'trials':radical['trials']}),flush=True)
    hostile = heisenberg_control()
    twins = hostile_singleton_slices()
    out = {'status':'FINITE-EXACT actual native primary clocks and explicit completion tests',
           'primary_source_canonical_sha256':hashlib.sha256(
               json.dumps(data,sort_keys=True,separators=(',',':')).encode()).hexdigest(),
           'primary_cases':cases,'three_primary_completion_tests':three,
           'three_primary_invariant_radical':radical,'joint_clock_hostile':hostile,
           'hostile_81_27':twins,'gates':GATES}
    (HERE/'primary_native_clocks.json').write_text(json.dumps(out,indent=2)+'\n')
    print(json.dumps({'gates':GATES,'joint_clock_hostile':hostile}),flush=True)


if __name__ == '__main__':
    main()
