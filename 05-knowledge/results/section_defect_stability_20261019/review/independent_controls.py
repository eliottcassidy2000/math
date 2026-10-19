#!/usr/bin/env python3
"""Independent exact audit; imports no submitted code; checks remain active under -O."""
from itertools import product
from fractions import Fraction
from collections import Counter
from pathlib import Path
import hashlib, json

p = 29
vectors = ((1, 0), (0, 1), (3, 9))
stress = (-3, -9, 1)
points = [(x, y) for y in range(p) for x in range(p)]
checks = Counter()

def require(test, label):
    checks[label] += 1
    if not test:
        raise RuntimeError('FAILED: ' + label)

def plus(a, b): return ((a[0]+b[0]) % p, (a[1]+b[1]) % p)
def times(n, a): return ((n*a[0]) % p, (n*a[1]) % p)
def det(a, b): return (a[0]*b[1]-a[1]*b[0]) % p

def diff(f, a):
    out = {}
    for x, value in f.items():
        for loc, sign in ((plus(x, times(-1, a)), 1), (x, -1)):
            out[loc] = out.get(loc, 0) + sign*value
    return {x: v for x, v in out.items() if v}

def dense_diff(f, a): return {x: f[plus(x, a)]-f[x] for x in points}

def elimination(matrix, n, modulus=None):
    A = [[(v % modulus if modulus else Fraction(v)) for v in row] for row in matrix]
    r = 0; pivots = []
    for col in range(n):
        pivot = next((j for j in range(r, len(A)) if A[j][col]), None)
        if pivot is None: continue
        A[r], A[pivot] = A[pivot], A[r]
        inverse = pow(A[r][col], -1, modulus) if modulus else 1/A[r][col]
        A[r] = [(v*inverse % modulus if modulus else v*inverse) for v in A[r]]
        for j in range(len(A)):
            if j == r or not A[j][col]: continue
            mult = A[j][col]
            A[j] = [((a-mult*b) % modulus if modulus else a-mult*b) for a, b in zip(A[j], A[r])]
        pivots.append(col); r += 1
        if r == len(A): break
    return A[:r], pivots

def kernel(matrix, n, modulus=None):
    reduced, pivots = elimination(matrix, n, modulus)
    result = []
    for free in range(n):
        if free in pivots: continue
        v = [0]*n; v[free] = 1
        for row, pivot in zip(reduced, pivots):
            v[pivot] = -row[free] % modulus if modulus else -row[free]
        result.append(v)
    return result

def constraints(support, omitted):
    rows = []
    for i, direction in enumerate(vectors):
        if i == omitted: continue
        for level in sorted({det(direction, x) for x in support}):
            rows.append([int(det(direction, x) == level) for x in support])
    return rows

for a,b,c in product(range(-2,3),repeat=3):
    require(((a-3*b-9*c) % p == 0) == (a == b == c == 0), 'all_125_balanced_digit_triples')

geometry_rows = []
short_data = []
for exceptional in range(3):
    e = vectors[exceptional]
    for separation in range(1, p):
        q = times(separation, e)
        segment = {times(j, e): 1 for j in range(1, separation+1)}
        dipole = diff(segment, e)
        require(dipole == {(0, 0): 1, q: -1}, 'dipole_orientation')
        cube = segment
        for d in vectors: cube = diff(cube, d)
        support = sorted(cube)
        require(len(support) == 8 and set(cube.values()) == {-1, 1}, 'eight_distinct_cube_vertices')
        axes = [times(separation if i == exceptional else 1, d) for i, d in enumerate(vectors)]
        shared_counts = []
        rational_bases = []
        for omitted in range(3):
            others = [i for i in range(3) if i != omitted]
            faces = []; face_nodes = []
            for b in (0, 1):
                face = {plus(q, times(-1, plus(times(b, axes[omitted]), plus(times(s, axes[others[0]]), times(t, axes[others[1]]))))) for s, t in product((0, 1), repeat=2)}
                faces.append(face)
                nodes = {(i, det(vectors[i], x)) for i in others for x in face}
                face_nodes.append(nodes)
                require(len(face) == 4 and len(nodes) == 4, 'face_is_incidence_four_cycle')
            shared = len(face_nodes[0] & face_nodes[1]); shared_counts.append(shared)
            require(not (faces[0] & faces[1]) and faces[0] | faces[1] == set(support) and shared <= 1, 'two_faces_share_at_most_one_vertex')
            A = constraints(support, omitted)
            basis = kernel(A, 8)
            require(len(basis) == len(kernel(A, 8, p)) == 2, 'integer_and_field_face_nullity_two')
            require(all(value.denominator == 1 for v in basis for value in map(Fraction, v)), 'saturated_integer_face_basis')
            basis = [[int(x) for x in v] for v in basis]
            rational_bases.append(basis)
            require(all(sum(row[j]*cube[x] for j, x in enumerate(support)) == 0 for row in A), 'integer_cube_line_sums')
            # Confirm actual integer kernel is precisely the two alternating face cycles.
            face_vectors = [[cube[x] if x in face else 0 for x in support] for face in faces]
            require(len(elimination(face_vectors, 8)[1]) == 2, 'two_independent_integer_face_cycles')
            require(all(sum(row[j]*fv[j] for j in range(8)) == 0 for row in A for fv in face_vectors), 'integer_face_cycle_kernel')
        coupled = [[stress[i]*v[j] for i in range(3) for v in rational_bases[i]] for j in range(8)]
        require(len(elimination(coupled, 6)[1]) == len(elimination(coupled, 6, p)[1]) == 4, 'coupled_rank_four_over_Q_and_F29')
        # An exact enumeration of every bounded integer face amplitude, not numerical rank alone.
        feasible = 0; modular_only = 0
        for coeffs in product(range(-2, 3), repeat=6):
            fields = [[sum(coeffs[2*i+b]*rational_bases[i][b][j] for b in range(2)) for j in range(8)] for i in range(3)]
            if any(abs(value) > 2 for field in fields for value in field): continue
            signed = [sum(stress[i]*fields[i][j] for i in range(3))*cube[x] for j, x in enumerate(support)]
            if any((v-signed[0]) % p for v in signed[1:]): continue
            feasible += 1
            if any(v != signed[0] for v in signed[1:]): modular_only += 1
            require(all(len({fields[i][j]*cube[x] for j, x in enumerate(support)}) == 1 for i in range(3)), 'exhaustive_bounded_coupling_forces_integer_scalars')
        require(feasible == 125 and modular_only == 0, 'bounded_coupling_exactly_125_scalar_triples')
        geometry_rows.append({'exceptional': exceptional, 'separation': separation, 'shared_vertices': shared_counts, 'bounded_face_solutions': feasible, 'modular_only_solutions': modular_only})
        if separation <= 14:
            # Check column support in actual coordinates of the derivative grid, independent of rank calculations.
            for d_index, d in enumerate(vectors):
                derivative = diff(segment, d)
                other = [j for j in range(3) if j != d_index]
                if d_index != exceptional:
                    other = [exceptional] + [j for j in other if j != exceptional]
                a, b = [vectors[j] for j in other]
                denom = pow(det(a, b), -1, p)
                columns = Counter((det(a, x)*denom) % p for x in derivative)
                require(len(columns) <= 2 and max(columns.values(), default=0) <= 14, 'actual_derivative_grid_short_columns')
                require(sum(derivative.values()) == 0, 'derivative_total_zero')
                cycle_sums = Counter()
                for x, value in derivative.items(): cycle_sums[det(d, x)] += value
                require(all(value == 0 for value in cycle_sums.values()), 'derivative_each_own_cycle_zero')
            short_data.append((exceptional, separation, segment, q))

# Exhaustive dipole normal forms over all ordered errors on all short cycle separations.
for exceptional, separation, segment, q in short_data:
    for beta in range(1, p):
        alpha = (-2-beta) % p
        if alpha == 0: continue
        t = (-beta) % p
        require(t != 0 and ((-2+t) % p, -t % p) == (alpha, beta), 'all_nonzero_exceptional_error_pairs')

# Binary normalization: exhaust all nonconstant F masks and all potentially admissible exceptional shifts.
# With an unmodified column fixing F as binary, any exceptional shift with an unmodified entry lies in {-1,0,1}.
normalization_rows = []
for size in (5, 6, 7, 8, 9):
    solutions = 0
    masks = list(product((0, 1), repeat=size))
    for F in masks:
        m = sum(F)
        if m in (0, size): continue  # Constant F is directly binary column-by-column at an unmodified entry.
        for g, h in product((-1, 0, 1), repeat=2):
            Jcols = [tuple(f+shift for f in F) for shift in (g, h)]
            Jtotal = size*m + size*(g+h)
            if Jtotal != size: continue
            options = []
            for Jcol in Jcols:
                weights = Counter(sum(Icol) for Icol in masks if sum(a != b for a, b in zip(Icol, Jcol)) <= size-2)
                options.append(weights)
            target = size-(size-2)*m
            for first_weight, multiplicity in options[0].items():
                count = multiplicity*options[1].get(target-first_weight, 0)
                if not count: continue
                solutions += count
                require(all(value in (0, 1) for col in Jcols for value in col), 'exhaustive_normalization_no_nonbinary_solution')
    normalization_rows.append({'p_or_grid_size': size, 'admissible_configurations': solutions, 'nonbinary_counterexamples': 0})

# Three explicit boundary hostiles: p-1 support, missing weight-p, and p=4.
def examine_hostile(size, F, G, changes, expected_weight):
    J = {(s, t): F[s]+G[t] for s in range(size) for t in range(size)}
    I = dict(J)
    for x, value in changes.items(): I[x] = value
    require(set(I.values()) <= {0, 1}, 'hostile_I_binary')
    require(not set(J.values()) <= {0, 1}, 'hostile_J_nonbinary')
    require(sum(I.values()) == sum(J.values()) == expected_weight, 'hostile_equal_weight')
    changed = [sum(I[s,t] != J[s,t] for s in range(size)) for t in range(size)]
    return {'size': size, 'I_and_J_weight': expected_weight, 'changed_positions_by_column': changed}

boundary_size = 5
boundary_F = [1,0,0,0,0]; boundary_G = [1,-1,0,0,0]
boundary_changes = {(s,0): boundary_F[s] for s in range(5) if s != 1}
boundary_changes.update({(s,1): boundary_F[s] for s in range(5) if s != 0})
hostile_boundary = examine_hostile(5, boundary_F, boundary_G, boundary_changes, 5)
weight_F = [1,1,0,0,0]; weight_G = [1,-1,0,0,0]
weight_changes = {(0,0):1, (1,0):1, (2,0):0, (2,1):0, (3,1):0, (4,1):0}
hostile_weight = examine_hostile(5, weight_F, weight_G, weight_changes, 10)
small_F = [1,1,0,0]; small_G = [0,-1,0,0]
small_changes = {(0,0):0, (1,0):0, (2,1):0, (3,1):0}
hostile_small = examine_hostile(4, small_F, small_G, small_changes, 4)

# Independently computed carry and period functional, with all rows retained.
def norm(x): return x[0]*x[0]+x[0]*x[1]+x[1]*x[1]
def linear(direction, value): return (2*direction[0]+direction[1])*value[0]+(direction[0]+2*direction[1])*value[1]
carry = []
for d in vectors:
    c = {}
    for x in points:
        y = plus(x,d); difference = (y[0]-x[0], y[1]-x[1])
        numerator = norm(difference)-1
        require(numerator % p == 0, 'literal_carry_divisibility')
        c[x] = numerator//p % p
    carry.append(c)
    sums=Counter()
    for x,value in c.items(): sums[det(d,x)]+=value
    require(len(sums)==p and all(value % p == p-2 for value in sums.values()),'all_direction_carry_cycle_sums')

def periods(z):
    r = [{x:(z[i].get(x,0)-carry[i][x]) % p for x in points} for i in range(3)]
    a = {}; b = {}
    for y in range(p):
        value = 0
        for x in range(p): a[x,y] = value; value = (value+r[0][x,y]) % p
        require(value == 0, 'period_u_integrability')
    for x in range(p):
        value = 0
        for y in range(p): b[x,y] = value; value = (value+r[1][x,y]) % p
        require(value == 0, 'period_v_integrability')
    return tuple(sum(r[2][x,y]-3*(a[plus((x,y),vectors[2])]-a[x,y])-9*(b[plus((x,y),vectors[2])]-b[x,y]) for x in range(p)) % p for y in range(p))

seams = []
for i in range(3):
    choices = []
    for parallel in range(3):
        if parallel == i: continue
        for offset in range(p):
            support = frozenset(x for x in points if det(vectors[parallel], x) == offset)
            choices.append((parallel, offset, support))
            require(len(support) == p and len({det(vectors[i], x) for x in support}) == p, 'all_seam_offsets_one_per_cycle')
    seams.append(choices)
base = [{x:-2 for x in seams[i][0][2]} for i in range(3)]
base_period = periods(base)
period_changes = []
for i in range(3):
    changes = []
    for _, _, support in seams[i]:
        z = [dict(f) for f in base]; z[i] = {x:-2 for x in support}
        result = periods(z)
        changes.append(tuple((a-b) % p for a,b in zip(result,base_period)))
    period_changes.append(changes)
seam_period_values = Counter()
for a,b,c in product(range(58), repeat=3):
    result = tuple((base_period[y]+period_changes[0][a][y]+period_changes[1][b][y]+period_changes[2][c][y]) % p for y in range(p))
    require(len(set(result)) == 1 and result[0] != 0, 'all_offsets_direct_integrated_period_nonzero')
    seam_period_values[result[0]] += 1
require(set(seam_period_values) == {2,4,13,14,15,16,25,27}, 'all_eight_expected_period_values')

# Positive controls: actual sections and arbitrary local rigidity motions have zero period.
for seed in range(3):
    h = {x: ((x[0]**2 + (seed+2)*x[1] + x[0]*x[1]) % p, (x[1]**2 + (seed+3)*x[0] + seed*x[0]*x[1]) % p) for x in points}
    z = []
    for i,d in enumerate(vectors):
        f = {}
        for x in points:
            y = plus(x,d); motion=(h[y][0]-h[x][0],h[y][1]-h[x][1])
            f[x]=(carry[i][x]+linear(d,motion)) % p
            X=(x[0]+p*h[x][0],x[1]+p*h[x][1]); Y=(y[0]+p*h[y][0],y[1]+p*h[y][1])
            direct=(norm((Y[0]-X[0],Y[1]-X[1]))-1)//p % p
            require(f[x] == direct, 'direct_precision_two_section_error')
        z.append(f)
    require(periods(z) == (0,)*p, 'actual_section_period_zero')
    for i,d in enumerate(vectors):
        sums=Counter()
        for x,value in z[i].items(): sums[det(d,x)]+=value
        require(len(sums)==p and all(value % p == p-2 for value in sums.values()),'actual_section_all_cycle_sums')
    mixed=[diff(diff(z[i],vectors[(i+1)%3]),vectors[(i+2)%3]) for i in range(3)]
    require(all(sum(stress[i]*mixed[i].get(x,0) for i in range(3)) % p == 0 for x in points),'actual_section_full_affine_stress')

# Exhaust every normalized seam/dipole candidate with all lambda choices and all seam offsets.
# This is exhaustive within the proved normal form, not a search of arbitrary sections.
normalized_rows = []
normalized_total = 0
amplitudes = Counter()
for exceptional, separation, segment, q in short_data:
    derivatives = [diff(segment,d) for d in vectors]
    options = []
    for i in range(3):
        family = []
        for seam_index, (parallel, offset, support) in enumerate(seams[i]):
            for lam in range(-2,3):
                I = set(support)
                good = True
                for x, value in derivatives[i].items():
                    value = int(x in support)+lam*value
                    if value not in (0,1): good = False; break
                    if value: I.add(x)
                    else: I.discard(x)
                if not good: continue
                if i == exceptional and ((0,0) not in I or q in I): continue
                require(len(I) == p and len({det(vectors[i],x) for x in I}) == p, 'normalized_candidate_still_one_per_cycle')
                # The mixed derivatives have exactly lambda*C as integer arrays.
                field = {x:1 for x in I}
                for j in range(3):
                    if j != i: field=diff(field,vectors[j])
                cube=segment
                for d in vectors: cube=diff(cube,d)
                require(field == {x:lam*v for x,v in cube.items() if lam*v}, 'normalized_candidate_integer_mixed_derivative')
                family.append((seam_index,lam,frozenset(I)))
        options.append(family)
    count=0; sample=None
    for combination in product(*options):
        lam=[item[1] for item in combination]
        amplitude=2*sum(stress[i]*lam[i] for i in range(3))*pow(stress[exceptional] % p,-1,p) % p
        if amplitude in (0,2): continue
        count+=1; amplitudes[amplitude]+=1
        L=[(-2*lam[i]+(amplitude if i == exceptional else 0)) % p for i in range(3)]
        H=((2*L[0]-L[1])*pow(3,-1,p) % p, (2*L[1]-L[0])*pow(3,-1,p) % p)
        require(all(linear(vectors[i],H) % p == L[i] for i in range(3)), 'all_normalized_candidates_actual_local_motion')
        period=tuple((base_period[y]+sum(period_changes[i][combination[i][0]][y] for i in range(3))) % p for y in range(p))
        require(len(set(period)) == 1 and period[0] != 0, 'all_normalized_candidates_global_period_nonzero')
        if sample is None: sample=(combination,amplitude,L,H,period)
    require(count > 0, 'normalized_hostiles_exist_every_short_geometry')
    # Independently expand one nontrivial sample per geometry to direct z and run the full period calculation.
    combination,amplitude,L,H,expected_period=sample
    z=[]
    for i in range(3):
        values={x:(-2*int(x in combination[i][2]) + (amplitude*diff(segment,vectors[exceptional]).get(x,0) if i == exceptional else 0)) % p for x in points}
        support=seams[i][combination[i][0]][2]
        require(all((values[x]+2*int(x in support)-L[i]*derivatives[i].get(x,0)) % p == 0 for x in points), 'expanded_hostile_actual_local_motion')
        z.append(values)
    require(sum(bool(v) for f in z for v in f.values()) == 88, 'expanded_hostile_has_88_errors')
    require(periods(z) == expected_period, 'expanded_hostile_direct_period_preserved')
    mixed=[diff(diff({x:z[i][x] for x in points if z[i][x]},vectors[(i+1)%3]),vectors[(i+2)%3]) for i in range(3)]
    require(all(sum(stress[i]*mixed[i].get(x,0) for i in range(3)) % p == 0 for x in points), 'expanded_hostile_compact_stress_zero')
    normalized_total+=count
    normalized_rows.append({'exceptional':exceptional,'separation':separation,'family_option_counts':list(map(len,options)),'88_error_candidates':count})

# Count unit pairs and rotate triples using the actual quadratic-ring multiplication.
def multiply(a,b): return ((a[0]*b[0]-a[1]*b[1]) % p,(a[0]*b[1]+a[1]*b[0]+a[1]*b[1]) % p)
def pair(a): return min(a,times(-1,a))
units=[x for x in points if norm(x) % p == 1]
pairs=sorted({pair(x) for x in units})
require(len(units)==30 and len(pairs)==15,'norm_circle_has_fifteen_pairs')
for g in units:
    k=(norm(g)-1)//p % p
    multiplier=(1-p*k*pow(2,-1,p)) % (p*p)
    lift=tuple(coordinate*multiplier % (p*p) for coordinate in g)
    require(tuple(coordinate % p for coordinate in lift)==g and norm(lift) % (p*p)==1,'all_thirty_norm_one_lifts')
occurrences=Counter(pair(multiply(g,d)) for g in pairs for d in vectors)
require(set(occurrences)==set(pairs) and set(occurrences.values())=={3},'rotated_triples_count_every_family_three_times')
require(15*89==3*445 and p*p*15==12615,'final_445_count')

out={'verdict':'PASS','modulus':p,'geometry_cases':len(geometry_rows),'geometries':geometry_rows,'normalization_exhaustion':normalization_rows,'normalization_boundary_hostiles':{'relaxed_support_p_minus_1':hostile_boundary,'removed_weight_p':hostile_weight,'removed_p_at_least_5':hostile_small},'seam_offset_combinations_directly_integrated_by_additivity':58**3,'seam_period_counts':dict(sorted(seam_period_values.items())),'normalized_88_error_candidates':normalized_total,'normalized_geometry_cases':normalized_rows,'dipole_amplitude_counts':dict(sorted(amplitudes.items())),'checks':dict(sorted(checks.items())),'proof_scope':'Exact controls supplement a universal proof; no enumeration of arbitrary sections, sharpness claim, or chromatic implication.'}
path=Path(__file__).with_name('independent-results.json')
path.write_text(json.dumps(out,indent=2)+'\n')
print(json.dumps({key:value for key,value in out.items() if key not in ('geometries','normalized_geometry_cases')},indent=2))
