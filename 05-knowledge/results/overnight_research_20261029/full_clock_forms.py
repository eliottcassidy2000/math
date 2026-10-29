"""Actual invariant alternating forms for the remaining elementary primes5/17."""
from itertools import combinations, permutations, product
from pathlib import Path
import hashlib
import json

from three_run_glue import three_runs, rref_mod, mm_mod

HERE = Path(__file__).resolve().parent
GATES = 0


def need(ok, why):
    global GATES
    GATES += 1
    if not ok:
        raise ValueError(why)


def transpose(a):
    return list(map(list, zip(*a)))


def nullspace(a, p, n):
    reduced, pivots = rref_mod(a, p)
    free = [i for i in range(n) if i not in pivots]
    basis = []
    for f in free:
        row = [int(i == f) for i in range(n)]
        for i, j in enumerate(pivots):
            row[j] = -reduced[i][f] % p
        basis.append(row)
    return basis


def form(a, n, p):
    J = [[0]*n for _ in range(n)]
    for x, (i, j) in zip(a, combinations(range(n), 2)):
        J[i][j], J[j][i] = x % p, -x % p
    return J


def main():
    source_path = HERE/'three_run_glue.json'
    source_data = json.loads(source_path.read_text())
    source = source_data['actual_tournaments'][-1]
    _, ports, M = three_runs(source['rows'])
    ix = {x:i for i, x in enumerate(ports)}
    actions = {
        'translation': [ix[tuple((x+1) % 7 for x in a)] for a in ports],
        'multiplication_two': [ix[tuple(2*x % 7 for x in a)] for a in ports],
        'run_swap_01': [ix[(a[1], a[0], a[2])] for a in ports],
        'run_cycle_012': [ix[(a[1], a[2], a[0])] for a in ports]}
    group = {tuple(ix[tuple((pow(2, b, 7)*a[g[j]]+c) % 7 for j in range(3))] for a in ports)
             for b, c, g in product(range(3), range(7), permutations(range(3)))}
    need(source['automorphism_order'] == 21 and len(group) == 126,
         'generators cover the independently counted full vertex group times S3')
    for perm in actions.values():
        need(sorted(perm) == list(range(210)), 'literal native source permutation')
        need(all(M[perm[i]][perm[j]] == M[i][j] for i, j in product(range(210), repeat=2)),
             'native source and target covariance, including run clocks')
    records = []
    for data in source['odd_primary_translation_modules']:
        p = data['prime']
        if p not in (5, 17):
            continue
        B, free = data['projection'], data['projection_free_columns']
        n = len(B)
        need(mm_mod(B, M, p) == [[0]*210 for _ in range(n)], 'actual named quotient relations')
        need([[row[j] for j in free] for row in B] ==
             [[int(i == j) for j in range(n)] for i in range(n)], 'surjective named projection')
        clocks = {}
        for name, perm in actions.items():
            C = [[row[perm[j]] for j in free] for row in B]
            need(mm_mod(C, B, p) == [[row[j] for j in perm] for row in B],
                 'every native clock intertwines all210 port columns')
            clocks[name] = C
        slots = list(combinations(range(n), 2))
        equations = []
        for C in clocks.values():
            for i, j in slots:
                equations.append([(C[k][i]*C[l][j]-C[l][i]*C[k][j]-int((k, l) == (i, j))) % p
                                  for k, l in slots])
        basis = nullspace(equations, p, len(slots))
        need(len(basis) <= 6, 'bounded exhaustive form search')
        chosen, tested = None, 0
        for coefficients in product(range(p), repeat=len(basis)):
            tested += 1
            a = [sum(c*b[i] for c, b in zip(coefficients, basis)) % p for i in range(len(slots))]
            J = form(a, n, p)
            if len(rref_mod(J, p)[1]) == n:
                chosen = (coefficients, J)
                break
        need(chosen is not None, 'perfect form for all actual vertex and run clock generators')
        coefficients, J = chosen
        need(all(J[i][i] == 0 and all((J[i][j]+J[j][i]) % p == 0 for j in range(n))
                 for i in range(n)), 'literal alternating form')
        for C in clocks.values():
            need(mm_mod(mm_mod(transpose(C), J, p), C, p) == J, 'literal invariant form')
        records.append({'prime':p, 'dimension':n, 'projection':B, 'free_columns':free,
                        'clocks':clocks, 'invariant_alternating_space_dimension':len(basis),
                        'space_basis_upper_triangles':basis, 'chosen_coefficients':coefficients,
                        'perfect_alternating_form':J, 'candidates_tested':tested,
                        'minimum_extra_dimension':0})
        print(json.dumps({'prime':p, 'dimension':n, 'form_space_dimension':len(basis),
                          'candidates_tested':tested}), flush=True)
    out = {'status':'FINITE-EXACT full native vertex/run clock forms', 'gates':GATES,
           'source_canonical_json_sha256':hashlib.sha256(
               json.dumps(source_data, sort_keys=True, separators=(',', ':')).encode()).hexdigest(),
           'native_clock_group_order':len(group),
           'ports':ports, 'native_clock_permutations':actions, 'cases':records,
           'scope':'Aut(Paley7) x S3, elementary primary groups5/17; native complement separate'}
    (HERE/'full_clock_forms.json').write_text(json.dumps(out, indent=2)+'\n')
    print(json.dumps({'gates':GATES}), flush=True)


if __name__ == '__main__':
    main()
