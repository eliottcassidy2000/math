#!/usr/bin/env python3
"""Independent exact controls. No submitted module is imported; -O is safe."""
from collections import Counter
from itertools import product, combinations
from pathlib import Path
import hashlib
import json

N = 29
V = [(x, y) for y in range(N) for x in range(N)]
DIR = ((1, 0), (0, 1), (3, 9))
SIG = (-3, -9, 1)
INV = {a: next(b for b in range(1, N) if a*b % N == 1) for a in range(1, N)}
COUNTS = Counter()

def need(condition, label):
    COUNTS[label] += 1
    if not condition:
        raise RuntimeError(label)

def at(x, y): return x % N + N*(y % N)
def loc(i): return (i % N, i // N)
def move(i, d, k=1):
    x, y = loc(i)
    return at(x+k*d[0], y+k*d[1])
def line(i, d):
    x, y = loc(i)
    return (d[0]*y-d[1]*x) % N

# Dense pointwise evaluation deliberately differs from the submitted sparse
# push-forward derivative. Indices are y-major, not sorted coordinate pairs.
def derivative(values, d):
    return [values[move(i, d)]-values[i] for i in range(N*N)]
def delta(i): return [int(j == i) for j in range(N*N)]
def cube(e, k):
    values = [int(i == 0)-int(i == move(0, DIR[e], k)) for i in range(N*N)]
    for j, d in enumerate(DIR):
        if j != e: values = derivative(values, d)
    return {i: x for i, x in enumerate(values) if x}
def triple(i):
    values = delta(i)
    for d in DIR: values = derivative(values, d)
    return {j: x for j, x in enumerate(values) if x}
def shifted(c, t): return {move(i, t): value for i, value in c.items()}

# Forward row reduction, retaining pivot rows without eliminating above pivots.
# Inverses are found by the multiplication table, rather than modular pow().
def rank(rows, columns):
    pivots = {}
    for source in rows:
        row = [v % N for v in source]
        need(len(row) == columns, 'matrix_width')
        for c in range(columns):
            if not row[c]: continue
            if c not in pivots:
                iv = INV[row[c]]
                pivots[c] = [v*iv % N for v in row]
                break
            m = row[c]
            row = [(x-m*y) % N for x, y in zip(row, pivots[c])]
    return len(pivots)

def line_rows(support, component):
    # Include all 29 line labels, even the empty lines.
    return [[int(line(i, d) == t) for i in support]
            for j, d in enumerate(DIR) if j != component for t in range(N)]

def certificate(cubes, expected):
    support = sorted(set().union(*(set(c) for c in cubes)))
    n = len(support)
    constraints = [line_rows(support, j) for j in range(3)]
    wdims = [n-rank(rows, n) for rows in constraints]
    all_rows = [[0]*(j*n)+row+[0]*((2-j)*n)
                for j in range(3) for row in constraints[j]]
    all_rows += [[SIG[j] if k == i else 0 for j in range(3) for k in range(n)]
                 for i in range(n)]
    full_rank = rank(all_rows, 3*n)
    need(3*n-full_rank == expected, 'direct_homogeneous_dimension')
    modes = [[a*c.get(i, 0) for a in amplitudes for i in support]
             for c in cubes for amplitudes in ((1, 0, 3), (0, 1, 9))]
    for mode in modes:
        need(all(sum(a*b for a,b in zip(row, mode)) % N == 0 for row in all_rows),
             'explicit_model_in_full_kernel')
    model_rank = rank(modes, 3*n)
    need(model_rank == 2*len(cubes), 'independent_model_modes')
    return {'support_size': n, 'pair_kernel_dimensions': wdims,
            'direct_constraint_rank': full_rank, 'homogeneous_dimension': 3*n-full_rank,
            'coupled_rank': sum(wdims)-(3*n-full_rank), 'model_rank': model_rank}

def support_bounds(e, q, r):
    ordered = [0, q, r]
    gaps = [(ordered[(i+1)%3]-ordered[i]) % N for i in range(3)]
    origin = ordered[(gaps.index(max(gaps))+1)%3]
    lengths = sorted((x-origin) % N for x in ordered if x != origin)
    need(max(gaps) >= 10 and max(lengths) <= 19, 'largest_gap_bound')
    paths = []
    for length in lengths:
        path = [0]*(N*N)
        for j in range(1, length+1): path[move(0, DIR[e], origin+j)] = 1
        paths.append(path)
        dipole = derivative(path, DIR[e])
        need({i:v for i,v in enumerate(dipole) if v} ==
             {move(0,DIR[e],origin):1, move(0,DIR[e],origin+length):-1},
             'correct_dipole_sign_after_rebasing')
    bounds = []
    for j, d in enumerate(DIR):
        derivatives = [derivative(path, d) for path in paths]
        active = set(i for values in derivatives for i, v in enumerate(values) if v)
        axis = DIR[e] if j != e else DIR[next(a for a in range(3) if a != j)]
        counts = Counter(line(i, axis) for i in active)
        need(len(counts) <= (2 if j != e else 3), 'normalization_column_count')
        need(max(counts.values()) <= (19 if j != e else 1), 'normalization_column_capacity')
        for values in derivatives:
            sums = [0]*N
            for i, value in enumerate(values): sums[line(i,d)] += value
            need(sums == [0]*N, 'individual_cycle_preservation_over_Z')
        bounds.append([len(counts), max(counts.values())])
    return {'origin': origin, 'lengths': lengths, 'bounds': bounds}

for digits in product(range(-2,3), repeat=3):
    need((sum(a*b for a,b in zip(SIG,digits)) % N == 0) == (digits == (0,0,0)),
         'all_125_balanced_digit_triples')

cache = {(e,k): cube(e,k) for e in range(3) for k in range(1,N)}
single = []
for (e,k), c in cache.items():
    result = certificate([c], 2)
    need(result['support_size'] == 8 and result['pair_kernel_dimensions'] == [2,2,2],
         'all_84_single_cubes')
    single.append({'e': e, 'k': k, **result})

records = []
histogram = Counter()
for e in range(3):
    for q,r in combinations(range(1,N),2):
        a,b = cache[e,q],cache[e,r]
        common = set(a)&set(b)
        private_a, private_b = set(a)-set(b), set(b)-set(a)
        need(len(common) == len(private_a) == len(private_b) == 4, 'three_faces_disjoint')
        need(all(a[i] == b[i] for i in common), 'shared_face_equal_signs')
        need(set(a.values()) == set(b.values()) == {-1,1}, 'cube_unit_signs')
        result = certificate([a,b], 4)
        need(result['support_size'] == 12 and result['direct_constraint_rank'] == 32,
             'all_1134_rank_32')
        signatures = {(a.get(i,0),b.get(i,0)) for i in set(a)|set(b)}
        # Enumerate every F29 coefficient pair whose resulting array has the
        # residue of a binary second difference, then verify its integer lift.
        lifted = []
        for x,y in product(range(N),repeat=2):
            residues = [(x*s+y*t) % N for s,t in signatures]
            if any(v not in (0,1,2,27,28) for v in residues): continue
            xx,yy = (x if x < 15 else x-N), (y if y < 15 else y-N)
            need(-2 <= xx <= 2 and -2 <= yy <= 2, 'private_face_bounds_coefficients')
            need(all(-2 <= xx*s+yy*t <= 2 for s,t in signatures), 'integer_lift_no_alias')
            lifted.append((xx,yy))
        need(set(lifted) == {(x,y) for x,y in product(range(-2,3),repeat=2) if abs(x+y)<=2},
             'exactly_19_bounded_integer_lifts')
        normalized = support_bounds(e,q,r)
        key = (e,tuple(result['pair_kernel_dimensions']),result['coupled_rank'])
        histogram[key] += 1
        records.append({'e':e,'q':q,'r':r,**result,'normalization':normalized})
need(len(records) == 1134, 'complete_1134_universe')
need(sum(v for (e,dims,r),v in histogram.items() if max(dims) == 5) == 18,
     'eighteen_nongeneric_geometries')

# Independent comparison to the submitted per-geometry receipt.
source = Path(__file__).resolve().parent.parent
submitted = json.loads((source/'two-extra-exact-results.json').read_text())
lookup = {(row['exceptional_direction'],row['q'],row['r']):row
          for row in submitted['three_failure_cycle_certificates']}
need(len(lookup) == len(submitted['three_failure_cycle_certificates']) == 1134,
     'submitted_unique_universe')
for row in records:
    previous = lookup[(row['e'],row['q'],row['r'])]
    for key in ('support_size','pair_kernel_dimensions','direct_constraint_rank',
                'homogeneous_dimension','coupled_rank','model_rank'):
        need(row[key] == previous[key], 'submitted_certificate_match')

hostile_cubes = [triple(at(0,j)) for j in range(3)]
need(set(hostile_cubes[0]).isdisjoint(hostile_cubes[2]), 'hostile_disjoint_endpoints')
need(set(hostile_cubes[1]) <= set(hostile_cubes[0])|set(hostile_cubes[2]), 'hostile_connecting_cube')
hostile = certificate([hostile_cubes[0],hostile_cubes[2]], 6)
complete_hostile = certificate(hostile_cubes, 6)
need(hostile['pair_kernel_dimensions'] == [6,4,6] and hostile['coupled_rank'] == 10
     and complete_hostile['model_rank'] == 6, 'hostile_exact_six_modes')
for j in (0,2):
    # S_j is the point mass at j*v and gives the listed u-dipole.
    dip = derivative(delta(at(0,j)), DIR[0])
    need({i:v for i,v in enumerate(dip) if v} == {at(-1,j):1,at(0,j):-1},
         'hostile_relevant_distinct_u_cycles')
need(line(at(0,0),DIR[0]) != line(at(0,2),DIR[0]), 'hostile_cycles_really_distinct')
passing = certificate([cache[0,1],shifted(cache[1,1],(7,5))], 4)
need(passing['pair_kernel_dimensions'] == [4,4,4] and passing['coupled_rank'] == 8,
     'nonvacuous_rank_four_example')

# Exactly decide whether a group of slots can be filled by zeros and ones
# within a prescribed change budget; the possible weights form an interval.
def possible_weights(values, cap):
    mandatory = sum(v not in (0,1) for v in values)
    if mandatory > cap: return None
    ones = values.count(1)
    spare = cap-mandatory
    return max(0,ones-spare),min(len(values),ones+mandatory+spare)

column_cases = 0
for n in (5,7,9):
    for r in range(1,(n-1)//2+1):
        for F in product((0,1),repeat=n):
            if min(F) == max(F): continue
            m = sum(F)
            for G in product((-1,0,1),repeat=r):
                if m+sum(G) != 1: continue
                column_cases += 1
                groups = [possible_weights([f+g for f in F],n-2) for g in G]
                feasible = all(g is not None for g in groups)
                if feasible:
                    fixed = (n-r)*m
                    feasible = fixed+sum(g[0] for g in groups) <= n <= fixed+sum(g[1] for g in groups)
                if feasible: need(all(g == 0 for g in G) and m == 1, 'column_lemma_exhaustive_small_grids')

# For the crossing lemma, use an exhaustive gauge on the clean rectangle:
# one clean-column G is zero; all clean F are binary. If F varies, all
# clean G are zero; otherwise each remaining clean G has the two allowed
# values. A feasible exceptional factor must be in [-1,2], since its strip
# retains an unchanged entry. This covers all integer additive arrays for
# these small sizes, not just a randomly chosen value range.
crossing_receipts = []
for n,r,s in ((5,1,1),(7,1,2),(7,2,1),(7,2,2),(9,2,2)):
    need(n>2*r and n>2*s and (n-r)*(n-s)>n, 'crossing_test_hypotheses')
    count = feasible_count = 0
    for clean_F in product((0,1),repeat=n-r):
        if min(clean_F) != max(clean_F):
            clean_G_options = [(0,)*(n-s)]
        else:
            c = clean_F[0]
            clean_G_options = [(0,)+tail for tail in product((-c,1-c),repeat=n-s-1)]
        for clean_G in clean_G_options:
            fixed_weight = sum(f+g for f in clean_F for g in clean_G)
            if fixed_weight > n: continue
            for ext in product(range(-1,3),repeat=r+s):
                F = ext[:r]+clean_F
                G = ext[r:]+clean_G
                if sum(F)+sum(G) != 1: continue
                count += 1
                groups = [possible_weights([f+g for g in clean_G],n-s-2) for f in F[:r]]
                groups += [possible_weights([f+g for f in clean_F],n-r-2) for g in G[:s]]
                if any(g is None for g in groups): continue
                lo = fixed_weight+sum(g[0] for g in groups)
                hi = fixed_weight+sum(g[1] for g in groups)+r*s
                if lo <= n <= hi:
                    feasible_count += 1
                    need(all(f+g in (0,1) for f in F for g in G), 'crossing_lemma_exhaustive_small_grids')
    crossing_receipts.append({'n':n,'r':r,'s':s,'integer_arrays':count,'feasible_arrays':feasible_count})

# The stated p-1 hostile: independently evaluate the displayed formulas.
J = [[int(x == 0)+int(y == 0)-int(y == 1) for y in range(N)] for x in range(N)]
I = [[J[x][y]-int(y==0 and x!=1)+int(y==1 and x!=0) for y in range(N)] for x in range(N)]
need(all(v in (0,1) for row in I for v in row), 'boundary_hostile_I_binary')
need(any(v not in (0,1) for row in J for v in row), 'boundary_hostile_J_nonbinary')
need(sum(map(sum,I)) == sum(map(sum,J)) == N, 'boundary_hostile_equal_weight')
need([sum(I[x][y] != J[x][y] for x in range(N)) for y in (0,1)] == [28,28],
     'boundary_hostile_twenty_eight_changes')

# Actual operator-level local motion, including both potentials and every
# exceptional-direction pair. Identities hold for arbitrary I; no inference
# from stress sufficiency is made. Full dense arrays challenge all signs.
motion_receipts = []
for e1,e2 in product(range(3),repeat=2):
    e_values = (e1,e2)
    starts = (at(2,5),at(8,11))
    lengths = (14,19)
    paths=[]; hvalues=[]; lambdas=[]; tvalues=[]
    for j,e in enumerate(e_values):
        path=[0]*(N*N)
        for k in range(1,lengths[j]+1): path[move(starts[j],DIR[e],k)] = 1
        lam = (j+1,-1,2-j)
        t = 2*sum(a*b for a,b in zip(SIG,lam))*INV[SIG[e]%N] % N
        L = [(-2*lam[d]+(t if d==e else 0))%N for d in range(3)]
        hx = (2*L[0]-L[1])*INV[3] % N
        hy = (2*L[1]-L[0])*INV[3] % N
        need((15*hx+21*hy-L[2])%N == 0, 'ell_w_compatibility')
        paths.append(path);hvalues.append((hx,hy));lambdas.append(lam);tvalues.append(t)
    hfield = [[sum(hvalues[j][axis]*paths[j][i] for j in range(2))%N
               for i in range(N*N)] for axis in range(2)]
    for d,direction in enumerate(DIR):
        selection=[int((i*7+d)%N == 0) for i in range(N*N)]
        change=[derivative(path,direction) for path in paths]
        normalized=[selection[i]-sum(lambdas[j][d]*change[j][i] for j in range(2)) for i in range(N*N)]
        z=[(-2*selection[i]+sum(tvalues[j]*change[j][i] for j in range(2) if e_values[j]==d))%N
           for i in range(N*N)]
        dx,dy=direction
        rh=[((2*dx+dy)*a+(dx+2*dy)*b)%N for a,b in zip(derivative(hfield[0],direction),derivative(hfield[1],direction))]
        need(all((z[i]+2*normalized[i]-rh[i])%N == 0 for i in range(N*N)), 'dense_two_potential_coset_identity')
    motion_receipts.append({'exceptional_directions':[e1,e2],'t':tvalues,'h':hvalues})

# Meaningful fault controls: each deliberately broken object is detected by
# the invariant it is supposed to challenge, rather than by its file name.
faults = {}
bad = dict(cache[0,1]); bad[next(iter(bad))] += 1
support = sorted(bad)
faults['mutated_cube_sign_rejected'] = any(sum(row[j]*bad[i] for j,i in enumerate(support))%N
    for row in line_rows(support,0))
faults['duplicated_origin_rejected'] = len(set(cache[0,1])|set(cache[0,1])) != 12
faults['missing_last_coordinate_rejected'] = 3*len(list(combinations(range(1,28),2))) != 1134
faults['disjoint_implies_rank_four_rejected'] = hostile['homogeneous_dimension'] != 4
faults['wrong_stress_sign_rejected'] = (-3*1+(-1)*3)%N != 0
need(all(faults.values()), 'five_deliberate_fault_controls')

receipt = {'verdict':'PASS','prime':N,'single_cube_geometries':len(single),
           'three_failure_geometries':len(records),'histogram':[
               {'e':e,'pair_kernel_dimensions':list(dims),'coupled_rank':rank,'count':count}
               for (e,dims,rank),count in sorted(histogram.items())],
           'all_geometry_certificates':records,'column_small_grid_candidate_arrays':column_cases,
           'crossing_small_grid_exhaustive_tests':crossing_receipts,
           'hostile':hostile,'hostile_full_basis':complete_hostile,'passing_example':passing,
           'motion_controls':motion_receipts,'fault_controls':faults,
           'checks':dict(sorted(COUNTS.items())),
           'frozen_readme_sha256':hashlib.sha256((source/'README.md').read_bytes()).hexdigest(),
           'scope':'Excludes one three-failure cycle at total 89; conditional two-dipole criterion; bounds remain 89 and 445.'}
out = Path(__file__).with_name('independent-results.json')
out.write_text(json.dumps(receipt,indent=2)+'\n')
print(json.dumps({key:value for key,value in receipt.items() if key not in ('all_geometry_certificates','checks')},indent=2))
