"""Read-only certificate verifier, independent of the extension solver.

Uses field Gaussian elimination and literal matrix products only.
It imports neither ternary_root_extensions nor its affine solver.
"""
from copy import deepcopy
from itertools import product
from pathlib import Path
import hashlib
import json

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


def kernel(a, p):
    a = [[x % p for x in row] for row in a]
    nc, pivots, r = len(a[0]), [], 0
    for j in range(nc):
        at = next((i for i in range(r, len(a)) if a[i][j]), None)
        if at is None:
            continue
        a[r], a[at] = a[at], a[r]
        inv = pow(a[r][j], -1, p)
        a[r] = [x*inv % p for x in a[r]]
        for i in range(len(a)):
            if i != r:
                c = a[i][j]
                a[i] = [(x-c*y) % p for x, y in zip(a[i], a[r])]
        pivots.append(j)
        r += 1
    basis = []
    for j in range(nc):
        if j in pivots:
            continue
        v = [int(i == j) for i in range(nc)]
        for i, col in enumerate(pivots):
            v[col] = -a[i][j] % p
        basis.append(v)
    return basis, r


def verify(out, source):
    digest = hashlib.sha256(json.dumps(source, sort_keys=True, separators=(',', ':')).encode()).hexdigest()
    need(digest == out['source_sha256_canonical'], 'unchanged canonical native source')
    c = source['primary_cases'][1]
    es, ds = c['exponents'], c['orders']
    F, R, S = [c['clocks'][k] for k in ('vertex_frobenius', 'run_cycle', 'run_swap')]
    u = out['universal_first_layer_obstruction']
    forms = u['torsion_form_basis']
    need(kernel(u['radical_generators'], 3)[1] == 1, 'radical image dimension')
    need(all(g[i] % 3 == 0 for g in u['radical_generators'] for i in range(9) if i != 5),
         'radical image supported only at e6')
    flat = [[J[i][j]*(ds[j]//3)//27 % 3 for i in range(9) for j in range(9)] for J in forms]
    need(kernel(flat, 3)[1] == 3, 'three independent torsion forms')
    for v in source['three_primary_invariant_radical']['complete_lift_basis']:
        J = [[0]*9 for _ in range(9)]
        at = 0
        for i in range(9):
            for j in range(i+1, 9):
                z = v[at]*3**(4-min(es[i], es[j])) % 81
                J[i][j], J[j][i] = z, -z % 81
                at += 1
        need(kernel(flat+[[J[i][j]*(ds[j]//3)//27 % 3 for i in range(9) for j in range(9)]], 3)[1] == 3,
             'all saved invariant forms represented on torsion')
        need(all(x % 81 == 0 for g in u['radical_generators'] for x in mv(J, g)), 'common radical retained')

    for rec in u['signs']:
        ep = rec['sign']
        rootset = set()
        for v in product(range(3), repeat=9):
            if any((x-y) % 3 for C in (F, R) for x, y in zip(mv(C, v), v)):
                continue
            if all((x-ep*y) % 3 == 0 for x, y in zip(mv(S, v), v)):
                rootset.add(v)
        need(len(rootset) == 9, 'independent exhaustive root-class census')

        def relations(w):
            f, r = w[:9], w[9:]
            blocks = [[3*x for x in f], [3*x for x in r],
                      [x+y+z for x, y, z in zip(f, mv(F, f), mv(F, mv(F, f)))],
                      [x+y+z for x, y, z in zip(r, mv(R, r), mv(R, mv(R, r)))],
                      [x-ep*y for x, y in zip(mv(S, f), f)],
                      [x+ep*y for x, y in zip(r, mv(R, mv(S, r)))],
                      [a-b-c+d for a, b, c, d in zip(mv(F, r), r, mv(R, f), f)]]
            return [x*(81//d) % 81 for block in blocks for x, d in zip(block, ds)]

        cols = [relations([int(i == j) for i in range(18)]) for j in range(18)]
        seen = set()
        for cert in rec['excluded_nonzero_roots']:
            v, y = cert['root'], cert['left_annihilator']
            residue = tuple(x % 3 for x in v)
            need(residue in rootset and any(residue) and residue not in seen, 'distinct actual root class')
            seen.add(residue)
            rhs = [(x-z)*(81//d) % 81 for C in (F, R) for x, z, d in zip(mv(C, v), v, ds)]+[0]*45
            need(all(sum(a*b for a, b in zip(y, col)) % 81 == 0 for col in cols), 'literal obstruction annihilator')
            value = sum(a*b for a, b in zip(y, rhs)) % 81
            need(value != 0 and value == cert['rhs_pairing'], 'literal contradictory RHS')
        need(seen == rootset-{(0,)*9}, 'every nonzero root excluded')

        torsioncols = [relations([int(i == j)*(ds[i % 9]//3) for i in range(18)]) for j in range(18)]
        need(all(x % 27 == 0 for col in torsioncols for x in col), 'torsion normalization exact')
        equations = tr([[x//27 % 3 for x in col] for col in torsioncols])
        cocycles, rank = kernel(equations, 3)
        need(len(cocycles) == rec['cocycle_dimension'], 'complete carry dimension')
        ft, rt, st = [[[C[i][j]*ds[j]//ds[i] % 3 for j in range(9)] for i in range(9)] for C in (F, R, S)]
        eigen, _ = kernel([[st[i][j]-ep*int(i == j) for j in range(9)] for i in range(9)], 3)
        cb = [[(x-y) % 3 for C in (ft, rt) for x, y in zip(mv(C, a), a)] for a in eigen]
        need(kernel(cb, 3)[1] == rec['coboundary_dimension'], 'actual source gauge rank')
        reps = rec['carry_class_basis']
        need(all(x % 3 == 0 for v in cb+reps for x in mv(equations, v)), 'all representatives are legal carries')
        need(kernel(cb+reps, 3)[1] == len(cocycles) and len(reps) == 4, 'complete carry quotient')
        seen = set()
        for cert in rec['row_certificates']:
            co = cert['carry_coefficients']
            need(tuple(co) not in seen and len(co) == 4 and all(x in range(3) for x in co), 'distinct carry class')
            seen.add(tuple(co))
            w = [sum(a*g[i] for a, g in zip(co, reps)) % 3 for i in range(18)]
            wf = [a*d//3 for a, d in zip(w[:9], ds)]
            wr = [a*d//3 for a, d in zip(w[9:], ds)]
            rows = [[ep*S[j][i]-int(i == j) for j in range(9)]+[0]*3 for i in range(9)]
            for C, ww in ((F, wf), (R, wr)):
                terms = [mv(tr(C), mv(J, ww)) for J in forms]
                need(all(x % 27 == 0 for term in terms for x in term), 'exact order-three pairing terms')
                rows.extend([[C[j][i]-int(i == j) for j in range(9)]+[t[i]//27 for t in terms] for i in range(9)])
            need([x % 3 for x in mv(tr(rows), cert['row_combination'])] == [int(i == 5) for i in range(12)],
                 'explicit row certificate forces radical character zero')
        need(len(seen) == 81, 'all carry classes certified')


def main():
    source = json.loads((HERE/'primary_native_clocks.json').read_text())
    out = json.loads((HERE/'ternary_root_extensions.json').read_text())
    verify(out, source)
    mutations = []
    for name in ('root_annihilator', 'character_row'):
        bad = deepcopy(out)
        sign = bad['universal_first_layer_obstruction']['signs'][0]
        if name == 'root_annihilator':
            sign['excluded_nonzero_roots'][0]['left_annihilator'] = [0]*63
        else:
            sign['row_certificates'][0]['row_combination'] = [0]*27
        try:
            verify(bad, source)
        except ValueError:
            mutations.append(name)
        else:
            raise ValueError('mutated certificate accepted: '+name)
    from three_run_glue import hostile_singleton_slices
    control = hostile_singleton_slices()
    receipt = {'status': 'VERIFIED read-only certificate replay and hostile mutations',
               'gates': GATES, 'rejected_mutations': mutations, 'hostile_81_27': control}
    (HERE/'ternary_obstruction_verify.json').write_text(json.dumps(receipt, indent=2)+'\n')
    print(json.dumps({'gates': GATES, 'rejected_mutations': mutations, 'hostile_81_27': 'passed'}))


if __name__ == '__main__':
    main()
