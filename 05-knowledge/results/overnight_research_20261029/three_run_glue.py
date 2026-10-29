"""Native three-run responses and their exact integral S3-sector lattice."""
from collections import Counter
from itertools import combinations, permutations, product
from math import gcd, lcm, prod, isqrt
from pathlib import Path
import hashlib
import argparse
import json
import sys

HERE = Path(__file__).resolve().parent
sys.path.insert(0, str(HERE.parent/'overnight_rank_return_20261029'))
from two_run import subset_paths, two_runs, walk_pair_inclusion_exclusion
from atom_precision import hostile, automorphisms, determinant
from paley_pencil import diagonal, local_layers
from port_queries import mm
from sympy import Matrix

GATES = 0
GROUP = list(permutations(range(3)))


def need(ok, why):
    global GATES
    GATES += 1
    if not ok:
        raise RuntimeError(why)


def compose(g,h):
    return tuple(g[h[i]] for i in range(3))


def inverse(g):
    return tuple(g.index(i) for i in range(3))


def representation(g):
    return [[int(g[j] == i)-int(g[2] == i) for j in range(2)] for i in range(2)]


SIGNS = [(-1)**sum(g[i] > g[j] for i in range(3) for j in range(i+1,3)) for g in GROUP]
REPS = [representation(g) for g in GROUP]
FOURIER = list(map(list, zip(*[[1,sgn,*sum(R,[])] for sgn,R in zip(SIGNS,REPS)])))
INVERSE6 = [[int(x) for x in row] for row in (6*Matrix(FOURIER).inv()).tolist()]


def lattice_checks():
    need(abs(determinant(FOURIER)) == 54, 'integral Fourier lattice index')
    need(diagonal(FOURIER) == [1,1,1,3,3,6], 'all six observer invariant factors')
    need(mm(INVERSE6,FOURIER) == [[6*int(i == j) for j in range(6)] for i in range(6)],
         'literal inverse and complete modulo6 realization guard')
    for g,h in product(GROUP,repeat=2):
        need(representation(compose(g,h)) == mm(representation(g),representation(h)),
             'the integral standard representation is a homomorphism')
    for c in product(range(-1,2),repeat=6):
        y = mm(FOURIER,[[x] for x in c])
        need(mm(INVERSE6,y) == [[6*x] for x in c], 'source coefficient recovery')
    # Complete image enumeration at the bad primes, not a generic-rank claim.
    modular = []
    for m in (2,3,4,6):
        images = Counter(tuple(sum(F*x for F,x in zip(row,c)) % m for row in FOURIER)
                         for c in product(range(m),repeat=6))
        expected = m**6//(gcd(m,3)**2*gcd(m,6))
        need(len(images) == expected, 'complete modular Fourier image count')
        need(set(images.values()) == {m**6//expected}, 'all source fibers accounted for')
        modular.append({'modulus':m,'states':m**6,'image':expected,'fiber':m**6//expected})
    for u,v,a,b,c,d in product(range(6),repeat=6):
        explicit = ((u-v) % 2 == 0 and (a+b-u) % 3 == 0 and
                    (c+d-u) % 3 == 0 and (a-c-v) % 3 == 0)
        y = [u,v,a,b,c,d]
        inverse_guard = all(sum(x*z for x,z in zip(row,y)) % 6 == 0 for row in INVERSE6)
        need(explicit == inverse_guard, 'four explicit congruences give the entire integral image')
    return {'matrix':FOURIER,'inverse_times_six':INVERSE6,'smith':[1,1,1,3,3,6],
            'index':54,'modular_censuses':modular}


def three_runs(rows):
    n = len(rows)
    subsets = list(combinations(range(n),3))
    ports = [tuple(s[i] for i in g) for s in subsets for g in GROUP]
    index = {p:i for i,p in enumerate(ports)}
    out = [[0]*len(ports) for _ in ports]
    table = subset_paths(rows)
    entries = [[(s,t,x) for s,row in enumerate(block) for t,x in enumerate(row) if x]
               for block in table]
    full = (1 << n)-1
    for mask1 in range(1,full):
        rest = full ^ mask1
        mask2 = (rest-1) & rest
        while mask2:
            mask3 = rest ^ mask2
            for a,b,x in entries[mask1]:
                for c,d,y in entries[mask2]:
                    xy = x*y
                    for e,f,z in entries[mask3]:
                        out[index[a,c,e]][index[b,d,f]] += xy*z
            mask2 = (mask2-1) & rest
    # Independent complete permutations with two cuts, including singleton runs.
    brute = [[0]*len(ports) for _ in ports]
    for order in permutations(range(n)):
        mandatory = {i+1 for i in range(n-1) if not(rows[order[i]] >> order[i+1] & 1)}
        if len(mandatory) > 2:
            continue
        for u,v in combinations(range(1,n),2):
            if mandatory <= {u,v}:
                brute[index[order[0],order[u],order[v]]][index[order[u-1],order[v-1],order[-1]]] += 1
    need(out == brute, 'complete independent permutation-and-two-cut census')
    return subsets, ports, out


def split_operator(M, subsets):
    b = len(subsets)
    plus = [[0]*b for _ in range(b)]
    minus = [[0]*b for _ in range(b)]
    standard = [[0]*(2*b) for _ in range(2*b)]
    for i,j in product(range(b),repeat=2):
        kernel = M[6*i][6*j:6*j+6]
        plus[i][j] = sum(kernel)
        minus[i][j] = sum(s*x for s,x in zip(SIGNS,kernel))
        L = [[sum(x*representation(inverse(g))[r][s] for x,g in zip(kernel,GROUP))
              for s in range(2)] for r in range(2)]
        for r,s in product(range(2),repeat=2):
            standard[2*i+r][2*j+s] = L[r][s]
        block = [row[6*j:6*j+6] for row in M[6*i:6*i+6]]
        for r,s in product(range(6),repeat=2):
            need(block[r][s] == kernel[GROUP.index(compose(GROUP[s],inverse(GROUP[r])))],
                 'literal run permutation covariance')
        predicted = [[0]*6 for _ in range(6)]
        predicted[0][0], predicted[1][1] = plus[i][j],minus[i][j]
        for col in range(2):
            for r,s in product(range(2),repeat=2):
                predicted[2+2*r+col][2+2*s+col] = L[r][s]
        need(mm(mm(FOURIER,block),INVERSE6) == [[6*x for x in row] for row in predicted],
             'exact rational splitting with two standard-column copies')
    return plus,minus,standard


def fast_sector_smith(s):
    # Declared certificate prime bank: exact determinant division proves
    # completeness, rather than presuming that no other prime occurs.
    det = abs(determinant(s))
    need(det != 0, 'the large-sector certificate requires nonsingularity')
    remaining, factors = det, {}
    for p in (2,3,5,11,17,79,5861,392809):
        need(all(p % d for d in range(2,isqrt(p)+1)), 'prime certificate by trial division')
        e = 0
        while remaining % p == 0:
            remaining //= p
            e += 1
        if e:
            factors[p] = e
    need(remaining == 1, 'the entire determinant factors over the certified bank')
    ds = [1]*len(s)
    for p,e in factors.items():
        layers = local_layers(s,p,e+1)
        vals = [v for v,multiplicity in enumerate(layers) for _ in range(multiplicity)]
        need(len(vals) == len(s) and sum(vals) == e, 'complete local factor certificate')
        ds = [x*p**v for x,v in zip(ds,vals)]
    need(prod(ds) == det and all(y % x == 0 for x,y in zip(ds,ds[1:])),
         'certified full sector Smith chain')
    return ds


def rref_mod(matrix,p):
    a = [[x % p for x in row] for row in matrix]
    pivots, row = [],0
    for col in range(len(a[0])):
        target = next((i for i in range(row,len(a)) if a[i][col]),None)
        if target is None:
            continue
        a[row],a[target] = a[target],a[row]
        iv = pow(a[row][col],-1,p)
        a[row] = [x*iv % p for x in a[row]]
        for i in range(len(a)):
            if i != row and a[i][col]:
                c = a[i][col]
                a[i] = [(x-c*y) % p for x,y in zip(a[i],a[row])]
        pivots.append(col)
        row += 1
        if row == len(a):
            break
    return a,pivots


def mm_mod(a,b,p):
    return [[x % p for x in row] for row in mm(a,b)]


def clock_pairing_modules(M,ports,ds):
    index = {port:i for i,port in enumerate(ports)}
    move = [index[tuple((x+1) % 7 for x in port)] for port in ports]
    out = []
    for p in (5,11,17,79,5861,392809):
        dim = sum(d % p == 0 for d in ds)
        need(all(d % (p*p) != 0 for d in ds), 'the selected primary group has exponent p')
        reduced,pivots = rref_mod(list(map(list,zip(*M))),p)
        free = [i for i in range(len(M)) if i not in pivots]
        need(len(free) == dim, 'independent cokernel dimension')
        projection = []
        for f in free:
            v = [0]*len(M)
            v[f] = 1
            for i,j in enumerate(pivots):
                v[j] = -reduced[i][f] % p
            projection.append(v)
        need(all(x == 0 for row in mm_mod(projection,M,p) for x in row), 'literal named cokernel projection')
        clock = [[row[move[f]] for f in free] for row in projection]
        need(mm_mod(clock,projection,p) == [[row[j] for j in move] for row in projection],
             'translation clock intertwines every named source port')
        identity = [[int(i == j) for j in range(dim)] for i in range(dim)]
        power = identity
        for _ in range(7):
            power = mm_mod(power,clock,p)
        need(power == identity, 'literal induced clock has order dividing7')
        nontrivial_rank = len(rref_mod([[clock[i][j]-identity[i][j] for j in range(dim)]
                                      for i in range(dim)],p)[1])
        trivial = dim-nontrivial_rank
        degree = next(k for k in range(1,7) if pow(p,k,7) == 1)
        if degree == 6:
            factors = [[1]*7]
        else:
            need(degree == 3, 'declared prime bank has cubic reciprocal factors')
            root = next(x for x in range(p) if x*x % p == (-7) % p)
            eta,bar = ((-1+root)*pow(2,-1,p)) % p,((-1-root)*pow(2,-1,p)) % p
            factors = sorted([[1,(-eta) % p,bar,p-1],[1,(-bar) % p,eta,p-1]])
            convolution = [0]*7
            for i,x in enumerate(factors[0]):
                for j,y in enumerate(factors[1]):
                    convolution[i+j] = (convolution[i+j]+x*y) % p
            need(convolution == [1]*7, 'literal quadratic-period cubic factorization')
        blocks = []
        for coefficients in factors:
            need(len(coefficients) == degree+1, 'cyclotomic factors have Frobenius-orbit degree')
            value = [[0]*dim for _ in range(dim)]
            for coeff in coefficients:
                value = mm_mod(value,clock,p)
                value = [[(value[i][j]+coeff*identity[i][j]) % p for j in range(dim)]
                         for i in range(dim)]
            nullity = dim-len(rref_mod(value,p)[1])
            need(nullity % degree == 0, 'complete isotypic block dimension')
            reciprocal = [(x*pow(coefficients[-1],-1,p)) % p for x in coefficients[::-1]]
            blocks.append({'polynomial':coefficients,'multiplicity':nullity//degree,
                           'reciprocal_polynomial':reciprocal})
        need(trivial+sum(degree*b['multiplicity'] for b in blocks) == dim,
             'every source character accounted for')
        missing = trivial % 2
        for b in blocks:
            partner = next(c for c in blocks if c['polynomial'] == b['reciprocal_polynomial'])
            if b['polynomial'] < partner['polynomial']:
                missing += degree*abs(b['multiplicity']-partner['multiplicity'])
        bilinear_dimension = None
        if missing == dim:
            equations = [[(clock[k][i]*clock[l][j]-int(i == k and j == l)) % p
                          for k,l in product(range(dim),repeat=2)]
                         for i,j in product(range(dim),repeat=2)]
            bilinear_dimension = dim*dim-len(rref_mod(equations,p)[1])
            need(bilinear_dimension == 0, 'independent check: every clock-invariant bilinear form vanishes')
            invclock = identity
            for _ in range(6):
                invclock = mm_mod(invclock,clock,p)
            doubled = [[0]*(2*dim) for _ in range(2*dim)]
            J = [[0]*(2*dim) for _ in range(2*dim)]
            for i,j in product(range(dim),repeat=2):
                doubled[i][j] = clock[i][j]
                doubled[i+dim][j+dim] = invclock[j][i]
                J[i][j+dim],J[i+dim][j] = identity[i][j],-identity[i][j]
            need(mm_mod(list(map(list,zip(*doubled))),mm_mod(J,doubled,p),p) ==
                 [[x % p for x in row] for row in J], 'hyperbolic dual clock preserves the named pairing')
        out.append({'prime':p,'source_dimension':dim,'trivial_multiplicity':trivial,
                    'Frobenius_orbit_degree':degree,'cyclotomic_blocks':blocks,
                    'minimum_C7_equivariant_extra_dimension':missing,
                    'abstract_extra_dimension':dim % 2,'clock_matrix':clock,
                    'invariant_bilinear_dimension_when_totally_unpaired':bilinear_dimension,
                    'projection_free_columns':free,'projection':projection})
    return out


def actual_case(label, rows, direct_smith=False):
    subsets,ports,M = three_runs(rows)
    sectors = split_operator(M,subsets)
    ds = [fast_sector_smith(s) if len(s)>30 else diagonal(s) for s in sectors]
    reference = json.loads((HERE/'three_run_direct_smith_reference.json').read_text())
    need(ds == reference['sectors'][label], 'match the independently computed direct Smith reference')
    if direct_smith:
        need(ds == [diagonal(s) for s in sectors], 'fresh direct SymPy Smith replay')
    sector_product = prod(ds[0])*prod(ds[1])*prod(ds[2])**2
    split_factors = ds[0]+ds[1]+ds[2]*2
    layers = {}
    for p in (2,3,5):
        actual = local_layers(M,p,12)
        sector_layers = [local_layers(s,p,12) for s in sectors]
        predicted = [sum(vals[k]*(2 if i == 2 else 1)
                         for i,vals in enumerate(sector_layers)) for k in range(12)]
        layers[str(p)] = {'actual':actual,'unglued':predicted,
                         'actual_unresolved':len(M)-sum(actual),
                         'unglued_unresolved':len(M)-sum(predicted)}
        if p == 5:
            need(actual == predicted, 'outside the projector primes the split is integral')
    if len(rows) <= 5:
        direct = diagonal(M)
        need(prod(direct) == sector_product, 'literal full Smith determinant equals sector product')
        for p in (2,3,5):
            vals = []
            for x in direct:
                if not x:
                    vals.append(None)
                    continue
                v = 0
                while x % p == 0:
                    x //= p
                    v += 1
                vals.append(v)
            need(layers[str(p)]['actual'] == [vals.count(k) for k in range(12)],
                 'direct full Smith versus separate modular elimination')
    else:
        direct = None
    reconstructed = None
    if sector_product:
        away = []
        for x in split_factors:
            for p in (2,3):
                while x % p == 0:
                    x //= p
            away.append(x)
        for i in range(len(away)):
            for j in range(i+1,len(away)):
                x,y = away[i],away[j]
                away[i],away[j] = gcd(x,y),lcm(x,y)
        exponents = []
        for p in (2,3):
            vals = [e for e,n in enumerate(layers[str(p)]['actual']) for _ in range(n)]
            need(len(vals) == len(M), 'every exceptional-prime Smith pivot is certified')
            exponents.append(vals)
        reconstructed = [x*2**a*3**b for x,a,b in zip(away,*exponents)]
        need(all(y % x == 0 for x,y in zip(reconstructed,reconstructed[1:])), 'assembled divisibility chain')
        need(prod(reconstructed) == sector_product, 'full Smith product from exact rational splitting')
        if direct is not None:
            need(reconstructed == direct, 'primary reconstruction versus direct full Smith')
    group = automorphisms(rows)
    idx = {x:i for i,x in enumerate(ports)}
    for g in group:
        moved = [idx[tuple(g[x] for x in p)] for p in ports]
        need(all(M[moved[i]][moved[j]] == M[i][j] for i,j in product(range(len(M)),repeat=2)),
             'actual tournament transporter covariance')
    clock_modules = clock_pairing_modules(M,ports,reconstructed) if label == 'Paley7' else None
    return {'label':label,'rows':rows,'matrix_order':len(M),'port_orbits':len(subsets),
            'total_three_run_covers':sum(map(sum,M)),'sector_smith':ds,'full_smith_small':direct,
            'full_smith_reconstructed':reconstructed,
            'odd_primary_translation_modules':clock_modules,
            'sector_determinant_product':sector_product,'local_layers':layers,
            'matrix_sha256':hashlib.sha256(json.dumps(M,separators=(',',':')).encode()).hexdigest(),
            'automorphism_order':len(group)}


def delete(rows,v):
    keep = [i for i in range(len(rows)) if i != v]
    return [sum(((rows[i] >> j) & 1) << k for k,j in enumerate(keep)) for i in keep]


def hostile_singleton_slices():
    out = []
    triples = []
    baseline = json.loads((HERE.parent/'overnight_rank_return_20261029/two_run.json').read_text())
    for flip in (False,True):
        rows = hostile(flip)
        group = automorphisms(rows)
        need(len(group) == (27 if flip else 81), 'original automorphism-order hostile')
        triples.append([s for s in combinations(range(12),3)
                        if all(sum((rows[v] >> w) & 1 for w in s) == 1 for v in s)])
        inherited = baseline['objects'][-1 if flip else -2]
        need(inherited['rows'] == rows, 'same labelled source as the two-run certificate')
        old_ranks = {str(p):local_layers(inherited['response'],p,1)[0] for p in (2,3,5)}
        orbits = []
        unused = set(range(len(rows)))
        while unused:
            v = min(unused)
            orbit = sorted({g[v] for g in group})
            unused.difference_update(orbit)
            child = delete(rows,v)
            ports,M,_ = two_runs(child)
            probes = [(0,0),(0,len(M)-1)]
            nz = next((i,j) for i,row in enumerate(M) for j,x in enumerate(row) if x)
            probes.append(nz)
            for i,j in probes:
                need(M[i][j] == walk_pair_inclusion_exclusion(child,ports[i],ports[j]),
                     'independent singleton-slice endpoint probe')
            orbits.append({'vertices':orbit,'representative':v,'slice_order':len(M),
                           'total_slice_covers':sum(map(sum,M)),
                           'ranks':{str(p):local_layers(M,p,1)[0] for p in (2,3,5)},
                           'matrix_sha256':hashlib.sha256(json.dumps(M,separators=(',',':')).encode()).hexdigest()})
        out.append({'control':81 if not flip else 27,'automorphism_order':len(group),
                    'two_run_baseline_ranks':old_ranks,
                    'singleton_run_slices':orbits})
    need(triples[0] == triples[1] and len(triples[0]) == 67, 'same 67 labelled cyclic triples')
    return out


def earlier_orders():
    pairs = list(combinations(range(4),2))
    reps = {}
    for code in range(64):
        rows = [0]*4
        for k,(i,j) in enumerate(pairs):
            u,v = (i,j) if code >> k & 1 else (j,i)
            rows[u] |= 1 << v
        canonical = min(sum(((rows[g[i]] >> g[j]) & 1) << k for k,(i,j) in enumerate(pairs))
                        for g in permutations(range(4)))
        reps.setdefault(canonical,{'rows':rows,'labelled_multiplicity':0})['labelled_multiplicity'] += 1
    need(len(reps) == 4, 'all 64 labelled order4 tournaments reach four classes')
    for code,item in sorted(reps.items()):
        subsets,_,M = three_runs(item['rows'])
        sectors = split_operator(M,subsets)
        ds = diagonal(M)
        ss = [diagonal(s) for s in sectors]
        split = ss[0]+ss[1]+ss[2]*2
        for i in range(len(split)):
            for j in range(i+1,len(split)):
                x,y = split[i],split[j]
                split[i],split[j] = gcd(x,y),lcm(x,y)
        need(ds == split, 'no Smith discrepancy in any order4 tournament')
        item['canonical_code'] = code
        item['full_smith'] = ds
    return list(reps.values())


def pairing_obstruction(actual):
    ds = actual[-1]['full_smith_reconstructed']
    bank = (2,3,5,11,17,79,5861,392809)
    profile, repairs, order_factors = {}, [], {}
    for p in bank:
        exponents = []
        for d in ds:
            e = 0
            while d % p == 0:
                d //= p
                e += 1
            if e:
                exponents.append(e)
        multiplicities = Counter(exponents)
        profile[str(p)] = dict(sorted(multiplicities.items()))
        order_factors[str(p)] = sum(exponents)
        for e,m in sorted(multiplicities.items()):
            if m % 2:
                repairs.append(p**e)
            need((m + m%2) % 2 == 0, 'minimal summand repair makes each cyclic multiplicity even')
    squarefree = prod(p for p in bank if order_factors[str(p)] % 2)
    need(squarefree == 79*392809 == 31031911, 'odd-prime obstruction to square or twice square')
    need(repairs == [2,8,16,81,79,392809], 'complete minimal abstract alternating repair')
    need(prod(int(p)**e for p,e in order_factors.items()) == prod(ds), 'entire finite group order certified')
    # General hyperbolic doubling on a complete small integer matrix box.
    cases = 0
    for a,b,c,d in product(range(-2,3),repeat=4):
        M = [[a,b],[c,d]]
        if a*d-b*c == 0:
            continue
        K = [[0,0,a,b],[0,0,c,d],[-a,-c,0,0],[-b,-d,0,0]]
        expected = sorted(diagonal(M)*2)
        need(diagonal(K) == expected, 'hyperbolic alternating presentation duplicates all Smith factors')
        cases += 1
    return {'primary_cyclic_multiplicities':profile,'order_factorization':order_factors,
            'squarefree_part_of_order':squarefree,'minimum_abstract_alternating_summands':repairs,
            'minimum_abstract_square_or_C2_paired_summands':repairs[1:],
            'hyperbolic_small_matrix_cases':cases,
            'scope':'abstract pairing obstruction and repair, not arithmetic realization or clock-minimal repair'}


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument('--direct-smith',action='store_true',help='also repeat the slower general integer Smith calculation')
    args = parser.parse_args()
    lattice = lattice_checks()
    small = earlier_orders()
    cases = [('C3',[2,4,1]),('strong4',[6,4,8,3])]
    for n,steps,label in ((5,{1,2},'T5'),(7,{1,2,4},'Paley7')):
        cases.append((label,[sum(1 << ((u+s) % n) for s in steps) for u in range(n)]))
    actual = []
    for label,rows in cases:
        actual.append(actual_case(label,rows,args.direct_smith))
        print(json.dumps({'progress':label,'mass':actual[-1]['total_three_run_covers'],
                          'local_layers':actual[-1]['local_layers']}),flush=True)
    controls = hostile_singleton_slices()
    obstruction = pairing_obstruction(actual)
    out = {'gates':GATES,'lattice':lattice,'order4_complete_census':small,
           'actual_tournaments':actual,'hostile_controls':controls,
           'pairing_obstruction':obstruction,
           'scope':'three-run endpoint response; complement relation still separate'}
    (HERE/'three_run_glue.json').write_text(json.dumps(out,indent=2)+'\n')
    print(json.dumps({'gates':GATES,'hostile_controls':controls}),flush=True)


if __name__ == '__main__':
    main()
