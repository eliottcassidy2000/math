"""Source-fixed odd-primary root recurrence with two commuting clocks.

Every check uses explicit exceptions. Counts are marked presentations,
not unlabelled extension classes. The native search keeps a specified form.
"""
from itertools import product
from pathlib import Path
import hashlib
import json

from sympy import Matrix, ZZ
from sympy.polys.matrices import DomainMatrix
from sympy.polys.matrices.normalforms import smith_normal_decomp
from primary_native_clocks import modular_kernel, identity, mul, power, dual_map
from paley_pencil import local_layers
from three_run_glue import rref_mod

HERE = Path(__file__).resolve().parent
GATES = 0


def need(ok, why):
    global GATES
    GATES += 1
    if not ok:
        raise ValueError(why)


def mv(a, v):
    return [sum(x*y for x, y in zip(row, v)) for row in a]


def tr(a):
    return list(map(list, zip(*a)))


def affine_solve(matrix, rhs, p, k):
    """Particular solution and COMPLETE homogeneous kernel over Z/p^k."""
    q = p**k
    a = [[x % q for x in row] for row in matrix]
    b = [x % q for x in rhs]
    nr, nc = len(a), len(a[0])
    V, vals = identity(nc), []

    def val(x):
        for h in range(k):
            if x % p**(h+1):
                return h
        return k

    for r in range(min(nr, nc)):
        v, i, j = min((val(a[i][j]), i, j)
                      for i in range(r, nr) for j in range(r, nc))
        if v == k:
            break
        vals.append(v)
        a[r], a[i], b[r], b[i] = a[i], a[r], b[i], b[r]
        for row in a:
            row[r], row[j] = row[j], row[r]
        for row in V:
            row[r], row[j] = row[j], row[r]
        pv = p**v
        inv = pow(a[r][r]//pv, -1, q)
        a[r] = [x*inv % q for x in a[r]]
        b[r] = b[r]*inv % q
        for i in range(r+1, nr):
            need(a[i][r] % pv == 0, 'exact affine row division')
            c = a[i][r]//pv
            a[i] = [(x-c*y) % q for x, y in zip(a[i], a[r])]
            b[i] = (b[i]-c*b[r]) % q
        for j in range(r+1, nc):
            need(a[r][j] % pv == 0, 'exact affine column division')
            c = a[r][j]//pv
            for row in a:
                row[j] = (row[j]-c*row[r]) % q
            for row in V:
                row[j] = (row[j]-c*row[r]) % q
    need(all(a[i][j] == (p**vals[i] if i == j and i < len(vals) else 0)
             for i in range(nr) for j in range(nc)), 'affine diagonalization')
    if any(b[i] % p**v for i, v in enumerate(vals)) or any(b[len(vals):]):
        return None, [], vals
    y = [b[i]//p**vals[i] if i < len(vals) else 0 for i in range(nc)]
    particular = [x % q for x in mv(V, y)]
    basis = []
    for j in range(nc):
        m = p**(k-vals[j]) if j < len(vals) else 1
        col = [row[j]*m % q for row in V]
        if any(col):
            basis.append(col)
    need([x % q for x in mv(matrix, particular)] == [x % q for x in rhs],
         'literal affine particular solution')
    for col in basis:
        need(all(x % q == 0 for x in mv(matrix, col)), 'literal affine kernel')
    return particular, basis, vals


def radical(case):
    es, J = case['exponents'], case['form']
    p, n, E = case['prime'], len(es), max(es)
    D = [[p**es[i]*int(i == j) for j in range(n)] for i in range(n)]
    B = dual_map(J, es, p)
    layers = local_layers([B[i]+D[i] for i in range(n)], p, E+1)
    gens, _ = modular_kernel(J, p, E)
    orders = [p**e for e in es]
    gens = [[x % d for x, d in zip(g, orders)] for g in gens]
    gens = [g for g in gens if any(g)]
    return sum(e*c for e, c in enumerate(layers)), layers, gens


def root_classes(case, eps):
    p, es, clocks = case['prime'], case['exponents'], case['clocks']
    n = len(es)
    F, R, S = [clocks[x] for x in ('vertex_frobenius', 'run_cycle', 'run_swap')]
    eq = [[C[i][j]-(a if i == j else 0) for j in range(n)]
          for C, a in ((F, 1), (R, 1), (S, eps)) for i in range(n)]
    basis, _ = modular_kernel(eq, p, 1)
    orders = [p**e for e in es]
    for coeffs in product(range(p), repeat=len(basis)):
        v = [sum(a*g[i] for a, g in zip(coeffs, basis)) % p for i in range(n)]
        sv = mv(S, v)
        v = [(x+eps*y)*pow(2, -1, d) % d for x, y, d in zip(v, sv, orders)]
        need(all((x-eps*y) % d == 0 for x, y, d in zip(mv(S, v), v, orders)),
             'chosen root is an exact involution eigenvector')
        yield list(coeffs), v


def extension_equations(case, eps, v):
    p, es, J = case['prime'], case['exponents'], case['form']
    q, n = p**max(es), len(es)
    ds = [p**e for e in es]
    scale = [q//d for d in ds]
    F, R, S = [case['clocks'][x] for x in ('vertex_frobenius', 'run_cycle', 'run_swap')]
    FF, RR, RS = mul(F, F, ds), mul(R, R, ds), mul(R, S, ds)

    def residual(u):
        wf, wr, z0 = u[:n], u[n:2*n], u[2*n:]
        z = [s*t for s, t in zip(scale, z0)]
        out = []
        def group_eq(xs):
            out.extend(x*s for x, s in zip(xs, scale))
        for C, w in ((F, wf), (R, wr)):
            cv = mv(C, v)
            group_eq([p*a-b+c for a, b, c in zip(w, cv, v)])
        for C, CC, w in ((F, FF, wf), (R, RR, wr)):
            group_eq([a+b+c for a, b, c in zip(w, mv(C, w), mv(CC, w))])
        group_eq([x-eps*y for x, y in zip(mv(S, wf), wf)])
        group_eq([x+eps*y for x, y in zip(wr, mv(RS, wr))])
        group_eq([a-b-c+d for a, b, c, d in zip(mv(F, wr), wr, mv(R, wf), wf)])
        out.extend(p*x-y for x, y in zip(z, mv(J, v)))
        out.extend(eps*x-y for x, y in zip(mv(tr(S), z), z))
        for C, w in ((F, wf), (R, wr)):
            out.extend(a+b-c for a, b, c in zip(mv(tr(C), z), mv(tr(C), mv(J, w)), z))
        out.append(sum(x*y for x, y in zip(v, z)))
        return [x % q for x in out]

    zero = residual([0]*(3*n))
    columns = []
    for j in range(3*n):
        unit = [int(i == j) for i in range(3*n)]
        columns.append([(a-b) % q for a, b in zip(residual(unit), zero)])
    return tr(columns), [(-x) % q for x in zero]


def extend(case, eps, v, solution):
    """Retain the source injection, root, pairing, and ALL three clocks."""
    p, es, J = case['prime'], case['exponents'], case['form']
    n, q = len(es), p**max(es)
    ds = [p**e for e in es]
    wf = [a % d for a, d in zip(solution[:n], ds)]
    wr = [a % d for a, d in zip(solution[n:2*n], ds)]
    z = [(q//d)*a % q for a, d in zip(solution[2*n:], ds)]
    presentation = [[ds[i]*int(i == j) for j in range(n)]+[-v[i]] for i in range(n)]
    presentation.append([0]*n+[p])
    dm = DomainMatrix.from_Matrix(Matrix(presentation)).convert_to(ZZ)
    diag, U, V = smith_normal_decomp(dm)
    need(diag == U*dm*V, 'literal integer Smith identity')
    U = U.to_Matrix()
    Ui = U.inv()
    need(all(x.q == 1 for x in Ui), 'unimodular integer source chart')
    allorders = [abs(int(diag.to_Matrix()[i, i])) for i in range(n+1)]
    keep = [i for i, d in enumerate(allorders) if d > 1]
    orders = [allorders[i] for i in keep]
    newes = []
    for d in orders:
        e, t = 0, d
        while t % p == 0:
            e, t = e+1, t//p
        need(t == 1, 'root group stays primary')
        newes.append(e)
    newq = p**max(newes)
    K = Matrix([J[i]+[z[i]] for i in range(n)]+[[-x for x in z]+[0]])
    need(all(int(x) % q == 0 for x in Matrix(presentation).T*K), 'pairing descends to presentation')
    K = Ui.T*K*Ui
    newJ = [[int(K[i, j])*newq//q % newq for j in keep] for i in keep]
    clocks = {}
    for name, C in case['clocks'].items():
        w, t = ([0]*n, 1)
        if name == 'vertex_frobenius':
            w = wf
        elif name == 'run_cycle':
            w = wr
        elif name == 'run_swap':
            t = eps
        lifted = Matrix([C[i]+[w[i]] for i in range(n)]+[[0]*n+[t]])
        transformed = U*lifted*Ui
        clocks[name] = [[int(transformed[i, j]) % allorders[i] for j in keep] for i in keep]
    from primary_native_clocks import check_relations
    check_relations(clocks, orders)
    for C in clocks.values():
        mat = Matrix(C)
        need(all(int(x) % newq == 0 for x in mat.T*Matrix(newJ)*mat-Matrix(newJ)),
             'new full-clock invariant alternating form')
    injection = [[int(U[i, j]) % allorders[i] for j in range(n)] for i in keep]
    im = Matrix(injection)
    need(all(int(x) % newq == 0 for x in im.T*Matrix(newJ)*im-Matrix(J)*(newq//q)),
         'old pairing recovered on every source pair')
    for name, C in case['clocks'].items():
        diff = Matrix(clocks[name])*im-im*Matrix(C)
        need(all(int(diff[i, j]) % orders[i] == 0 for i in range(len(keep)) for j in range(n)),
             'source injection intertwines actual clocks')
    out = {'prime': p, 'exponents': newes, 'clocks': clocks, 'form': newJ}
    receipt = {'sign': eps, 'root': v, 'frobenius_carry': wf, 'run_carry': wr,
               'pairing_column': z, 'old_orders': ds, 'new_orders': orders,
               'source_injection': injection, 'root_coordinate': [int(U[i, n]) % allorders[i] for i in keep],
               'presentation': presentation, 'smith_left': [list(map(int, row)) for row in U.tolist()]}
    return out, receipt


def first_layer(case, take_first=False):
    """Exact root-class scan; affine kernels are complete, not sampled."""
    p, es = case['prime'], case['exponents']
    q, n = p**max(es), len(es)
    radlog, _, radgens = radical(case)
    rows, successes = [], []
    for eps in (1, -1):
        for coeffs, v in root_classes(case, eps):
            A, b = extension_equations(case, eps, v)
            sol, kernel, pivots = affine_solve(A, b, p, max(es))
            record = {'sign': eps, 'root_coefficients': coeffs, 'root': v,
                      'solvable': sol is not None, 'radical_reducing': False}
            if sol is not None:
                def nonzero_restriction(u):
                    z = [(q//p**e)*a for e, a in zip(es, u[2*n:])]
                    return any(sum(a*b for a, b in zip(z, g)) % q for g in radgens)
                witness = sol if nonzero_restriction(sol) else next(
                    ([a+b for a, b in zip(sol, g)] for g in kernel if nonzero_restriction(g)), None)
                record['marked_solution_log_order'] = (sum(pivots)+(3*n-len(pivots))*max(es)
                                                       -3*sum(max(es)-e for e in es))
                if witness is not None:
                    nxt, receipt = extend(case, eps, v, witness)
                    nextlog, _, _ = radical(nxt)
                    need(nextlog == radlog-1, 'successful root lowers radical order by exactly p')
                    need(sum(nxt['exponents']) == sum(es)+1, 'source index exactly p')
                    record['radical_reducing'] = True
                    receipt['radical_log_before'] = radlog
                    receipt['radical_log_after'] = nextlog
                    successes.append((nxt, receipt))
                    if take_first:
                        return rows+[record], successes
            rows.append(record)
    return rows, successes


def span_basis(vectors, p):
    out = []
    for v in vectors:
        if len(rref_mod(out+[v], p)[1]) > len(out):
            out.append([x % p for x in v])
    return out


def combinations(basis, p, size):
    for co in product(range(p), repeat=len(basis)):
        yield list(co), [sum(a*g[i] for a, g in zip(co, basis)) % p for i in range(size)]


def universal_obstruction(case, form_vectors):
    """Exclude coisotropic first steps for EVERY invariant source form.

    The first census is over all root classes. The second is over the
    81 cocycle/gauge classes per sign, independently checked on all
    729+2187 literal torsion carry pairs. All finite linear witnesses
    are saved; no randomized search or fixed-form inference is used.
    """
    from primary_native_clocks import alternating_matrix
    p, es = case['prime'], case['exponents']
    E, n = max(es), len(es)
    q, ds = p**E, [p**e for e in es]
    F, R, S = [case['clocks'][x] for x in ('vertex_frobenius', 'run_cycle', 'run_swap')]
    radlog, _, radgens = radical(case)
    need(radlog == 8, 'native common radical order')
    radmod = span_basis(radgens, p)
    need(len(radmod) == 1 and radmod[0][5] != 0 and
         all(x == 0 for i, x in enumerate(radmod[0]) if i != 5),
         'only the sixth source direction survives from R in V/3V')
    forms, flattened = [], []
    allforms = [alternating_matrix(v, es, p) for v in form_vectors]
    for J in allforms:
        need(all(x % q == 0 for g in radgens for x in mv(J, g)),
             'displayed radical lies in EVERY invariant form radical')
        entries = [J[i][j]*(ds[j]//p) for i in range(n) for j in range(n)]
        need(all(x % (q//p) == 0 for x in entries), 'torsion pairings have order at most p')
        flat = [x//(q//p) % p for x in entries]
        if len(rref_mod(flattened+[flat], p)[1]) > len(flattened):
            flattened.append(flat)
            forms.append(J)
    need(len(forms) == 3, 'complete torsion pairing image dimension')
    torsion = []
    for C in (F, R, S):
        need(all(C[i][j]*ds[j] % ds[i] == 0 for i in range(n) for j in range(n)),
             'actual torsion action integral')
        torsion.append([[C[i][j]*ds[j]//ds[i] % p for j in range(n)] for i in range(n)])
    ft, rt, st = torsion
    result = {'radical_generators': radgens, 'radical_mod_p_basis': radmod,
              'torsion_form_basis': forms, 'signs': []}
    for eps in (1, -1):
        zeroA, _ = extension_equations(case, eps, [0]*n)
        groupA = [row[:2*n] for row in zeroA[:7*n]]
        leftnull, _ = modular_kernel(tr(groupA), p, E)
        root_receipts = []
        for co, v in root_classes(case, eps):
            a, b = extension_equations(case, eps, v)
            rhs = b[:7*n]
            need([row[:2*n] for row in a[:7*n]] == groupA, 'root-independent group equations')
            if any(co):
                y = next((y for y in leftnull if sum(x*z for x, z in zip(y, rhs)) % q), None)
                need(y is not None, 'every nonzero root class has an exact incompatibility witness')
                need(all(x % q == 0 for x in mv(tr(groupA), y)), 'left annihilator certificate')
                root_receipts.append({'coefficients': co, 'root': v, 'left_annihilator': y,
                                      'rhs_pairing': sum(x*z for x, z in zip(y, rhs)) % q})
            else:
                need(not any(rhs), 'zero root has zero group obstruction')
        entries = [[row[j]*(ds[j % n]//p) for j in range(2*n)] for row in groupA]
        need(all(x % (q//p) == 0 for row in entries for x in row), 'torsion equations descend to Fp')
        torsion_eq = [[x//(q//p) % p for x in row] for row in entries]
        cocycles, _ = modular_kernel(torsion_eq, p, 1)
        eigen, _ = modular_kernel([[st[i][j]-eps*int(i == j) for j in range(n)] for i in range(n)], p, 1)
        coboundaries = span_basis([[(a-b) % p for C in (ft, rt) for a, b in zip(mv(C, v), v)]
                                 for v in eigen], p)
        need(all(x % p == 0 for v in coboundaries for x in mv(torsion_eq, v)),
             'root gauge changes are literal cocycles')
        joined, representatives = list(coboundaries), []
        for v in cocycles:
            if len(rref_mod(joined+[v], p)[1]) > len(joined):
                joined.append(v)
                representatives.append(v)
        need(len(joined) == len(cocycles) and len(representatives) == 4,
             'four-dimensional complete quotient of carries by source-fixed gauge')

        def pairing_rows(w):
            wf = [a*d//p for a, d in zip(w[:n], ds)]
            wr = [a*d//p for a, d in zip(w[n:], ds)]
            rows = [[eps*S[j][i]-int(i == j) for j in range(n)]+[0]*len(forms) for i in range(n)]
            for C, ww in ((F, wf), (R, wr)):
                cols = [mv(tr(C), mv(J, ww)) for J in forms]
                need(all(x % (q//p) == 0 for col in cols for x in col), 'clock/form term has order p')
                rows.extend([[C[j][i]-int(i == j) for j in range(n)]+
                             [a[i]//(q//p) % p for a in cols] for i in range(n)])
            return [[x % p for x in row] for row in rows]

        target = [int(i == 5) for i in range(n+len(forms))]
        certs = []
        for co, w in combinations(representatives, p, 2*n):
            rows = pairing_rows(w)
            y, _, _ = affine_solve(tr(rows), target, p, 1)
            need(y is not None, 'every source form and carry class forces z(e_6)=0')
            need([x % p for x in mv(tr(rows), y)] == target, 'literal field row certificate')
            certs.append({'carry_coefficients': co, 'row_combination': y})
        # Independent Gaussian row-rank path on every literal carry pair,
        # including all gauge representatives omitted from the saved census.
        literal_count = 0
        for _, w in combinations(cocycles, p, 2*n):
            rows = pairing_rows(w)
            need(len(rref_mod(rows, p)[1]) == len(rref_mod(rows+[target], p)[1]),
                 'independent exhaustive carry-pair row-space audit')
            literal_count += 1
        result['signs'].append({'sign': eps, 'excluded_nonzero_roots': root_receipts,
                                'cocycle_dimension': len(cocycles), 'coboundary_dimension': len(coboundaries),
                                'carry_class_basis': representatives, 'row_certificates': certs,
                                'literal_carry_pairs_checked': literal_count})
    result['conclusion'] = 'No invariant perfect completion has index 3^8; square order raises lower bound to 3^10.'
    return result


def affine_audit():
    totals = []
    for p, k in ((3, 1), (2, 2)):
        q = p**k
        vectors = list(product(range(q), repeat=2))
        checked = 0
        for entries in product(range(q), repeat=4):
            a = [list(entries[:2]), list(entries[2:])]
            fibres = {}
            for v in vectors:
                image = tuple(x % q for x in mv(a, v))
                fibres.setdefault(image, set()).add(v)
            for b in vectors:
                sol, ker, _ = affine_solve(a, b, p, k)
                literal = fibres.get(b, set())
                if sol is None:
                    need(not literal, 'exhaustive affine nonexistence audit')
                else:
                    generated = {tuple(sol)}
                    for g in ker:
                        generated = {tuple((x+t*y) % q for x, y in zip(v, g))
                                     for v in generated for t in range(q)}
                    need(generated == literal, 'exhaustive complete affine fibre audit')
                checked += 1
        totals.append({'modulus': q, 'all_2_by_2_matrices': q**4, 'all_rhs_checked': checked})
    return totals


def controls():
    # Two genuine non-split roots: C3^2 -> C3+C9 -> C9^2.
    # This checks the code does not accidentally forbid root attachments.
    case = {'prime': 3, 'exponents': [1, 1], 'form': [[0, 0], [0, 0]],
            'clocks': {name: identity(2) for name in
                       ('translation', 'vertex_frobenius', 'run_swap', 'run_cycle')}}
    first, r1 = extend(case, 1, [1, 0], [0, 0, 0, 0, 0, 1])
    need(radical(first)[0] == 1, 'first non-split root reduces positive-control radical')
    rows, candidates = first_layer(first, True)
    need(bool(candidates), 'second non-split root is found by recurrence')
    final, r2 = candidates[0]
    need(radical(final)[0] == 0 and final['exponents'] == [2, 2], 'positive control reaches perfect C9 squared')
    # Source-free coefficient mistakes are caught by the old joint-clock
    # hostile; this remains a different group from the native S3 example.
    from primary_native_clocks import heisenberg_control
    return {'nonsplit_root_chain': [r1, r2], 'final': final,
            'heisenberg': heisenberg_control()}


def main():
    data = json.loads((HERE/'primary_native_clocks.json').read_text())
    case = dict(data['primary_cases'][1])
    case['form'] = data['three_primary_invariant_radical']['best_form']['form']
    rows, successes = first_layer(case)
    need(not successes, 'saved optimal form has no radical-reducing first step')
    universal = universal_obstruction(case, data['three_primary_invariant_radical']['complete_lift_basis'])
    audit = affine_audit()
    positive = controls()
    out = {'status': 'PROVED finite-certificate exclusion / all-size marked root recurrence',
           'source_sha256_canonical': hashlib.sha256(json.dumps(data, sort_keys=True, separators=(',', ':')).encode()).hexdigest(),
           'specified_form_first_layer': rows, 'universal_first_layer_obstruction': universal,
           'affine_solver_audit': audit, 'controls': positive, 'gates': GATES}
    (HERE/'ternary_root_extensions.json').write_text(json.dumps(out, indent=2)+'\n')
    print(json.dumps({'gates': GATES, 'native_root_classes': len(rows),
                      'all_carry_pairs': sum(s['literal_carry_pairs_checked'] for s in universal['signs']),
                      'ternary_index_lower_exponent': 10, 'affine_audit': audit}))


if __name__ == '__main__':
    main()
