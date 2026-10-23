"""Exact restricted-observer controls; standard library only, no asserts.

Run: python -B check.py [--write]. The supplied observer is (0,...,n-1).
The native target has n vertices; the auxiliary target has n-1 vertices.
"""
from collections import Counter
from fractions import Fraction as F
import importlib.util
import json
from pathlib import Path
import sys

HERE = Path(__file__).resolve().parent
PRIOR = HERE.parents[1] / 'attachment_lattices_square_sum_20261023' / 'complement' / 'check.py'
spec = importlib.util.spec_from_file_location('observer_canon', PRIOR)
old = importlib.util.module_from_spec(spec)
spec.loader.exec_module(old)
CHECKS = 0


def need(test, message):
    global CHECKS
    CHECKS += 1
    if not test:
        raise RuntimeError(message)


def psd(matrix):
    """Exact rational LDL, including the zero-pivot obstruction."""
    a = [list(map(F, row)) for row in matrix]
    for k in range(len(a)):
        p = a[k][k]
        if p < 0:
            return False
        if not p:
            if any(a[k][j] for j in range(k+1, len(a))):
                return False
            continue
        for i in range(k+1, len(a)):
            for j in range(i, len(a)):
                a[i][j] -= a[i][k] * a[k][j] / p
                a[j][i] = a[i][j]
    return True


def coupling(cells):
    gs = sorted({g for g, t in cells})
    ts = sorted({t for g, t in cells})
    w = [[F(cells[g, t]) for t in ts] for g in gs]
    r = [sum(row) for row in w]
    c = [sum(row[j] for row in w) for j in range(len(ts))]
    gram = [[sum(w[i][j]*w[k][j]/c[j] for j in range(len(ts)))
             for k in range(len(gs))] for i in range(len(gs))]
    k = [[gram[i][j]/r[i] for j in range(len(gs))] for i in range(len(gs))]
    need(all(sum(row) == 1 for row in k), 'stochastic')
    need(all(r[i]*k[i][j] == r[j]*k[j][i]
             for i in range(len(gs)) for j in range(len(gs))), 'reversible')
    need(psd(gram), 'PSD kernel')
    unseen = set(range(len(gs)))
    components = []
    while unseen:
        seen = {min(unseen)}
        todo = list(seen)
        while todo:
            i = todo.pop()
            for j, p in enumerate(k[i]):
                if p and j not in seen:
                    seen.add(j)
                    todo.append(j)
        unseen -= seen
        components.append(sorted(seen))
    return gs, ts, w, r, c, gram, k, components


def gap_certificate(data, gamma, component=None):
    gs, ts, w, r, c, gram, k, components = data
    inds = list(range(len(gs))) if component is None else component
    z = sum(r[i] for i in inds)
    return psd([[(1-gamma)*r[i]*(i == j) + gamma*r[i]*r[j]/z - gram[i][j]
                 for j in inds] for i in inds])


def finite_summary(data):
    gs, ts, w, r, c, gram, k, components = data
    return {'source_classes': len(gs), 'target_classes': len(ts),
            'component_source_sizes': [len(s) for s in components],
            'global_gap_at_least_8_over_9': gap_certificate(data, F(8, 9)),
            'component_gap_at_least_8_over_9': [gap_certificate(data, F(8, 9), s)
                                               for s in components]}


def run():
    census = []
    exact_four = None
    for n in range(3, 7):
        m = n-1
        mpairs, midx, mperms, mactions = old.setup(m)
        npairs, nidx, nperms, nactions = old.setup(n)
        size = 1 << len(mpairs)
        gc, greps, _ = old.orbit_table(size, mactions, False)
        uc, ureps, _ = old.orbit_table(size, mactions, True)
        tc, treps, _ = old.orbit_table(1 << len(npairs), nactions, True)
        free_to_native = [nidx[i, j+1] for i, j in mpairs]
        native = [old.relabel(x, free_to_native) for x in range(size)]
        fulln = (1 << len(npairs))-1
        need(all(not (t >> nidx[i, i+1] & 1) for t in native for i in range(n-1)),
             'all native states contain the supplied Hamiltonian path')
        need([x for x, t in enumerate(native) if tc[t] == tc[0]] == [0],
             'unique transitive presentation')
        targets = {
            'native': [tc[t] for t in native],
            'native_converse': [min(tc[t], tc[t ^ fulln]) for t in native],
            'auxiliary': uc,
        }
        report = {'n': n, 'compressed_vertices': m, 'path_cube_size': size}
        data_by_name = {}
        for name, target in targets.items():
            cells = Counter((gc[x], target[x]) for x in range(size))
            data = coupling(cells)
            data_by_name[name] = data
            report[name] = finite_summary(data)
            if name.startswith('native'):
                need(len(data[-1]) >= 2, 'native coupling disconnected')
                need(data[-1][0] == [0], 'empty graph isolated')
                need(not gap_certificate(data, F(8, 9)), 'native global bound rejected')
            else:
                need(len(data[-1]) == 1, 'auxiliary connected')
                need(gap_certificate(data, F(1) if m == 2 else F(8, 9)),
                     'whole-cube auxiliary gap')
                if m >= 3:
                    need(not gap_certificate(data, F(8, 9)+F(1, 1000)),
                         'auxiliary sharp gap')

        # Coarsening both observations on the SAME measure retains the gap.
        fullm = size-1
        coarse = coupling(Counter((min(gc[x], gc[x ^ fullm]),
                                   min(uc[x], uc[x ^ fullm])) for x in range(size)))
        need(gap_certificate(coarse, F(8, 9)), 'safe coarsening')
        report['auxiliary_complement_converse'] = finite_summary(coarse)

        # Soft, strictly positive, non-invariant reweighting; w ranges in [1,3].
        softcells = Counter()
        for x in range(size):
            softcells[gc[x], uc[x]] += 1 + x.bit_count() % 3
        soft = coupling(softcells)
        need(gap_certificate(soft, F(8, 27)), 'positive weight comparison')
        report['soft_weight_gap_at_least'] = '8/27'

        if n == 4:
            data = data_by_name['native']
            need([len(s) for s in data[-1]] == [1, 3], 'order-four components')
            s = data[-1][1]
            need(gap_certificate(data, F(7, 15), s), 'order-four lower bound')
            need(not gap_certificate(data, F(7, 15)+F(1, 1000), s),
                 'order-four bound sharp')
            # Exact eigenvector, in increasing graph-class canonical mask order.
            f = [F(4), F(-3), F(-3)]
            need(sum(data[3][i]*f[a] for a, i in enumerate(s)) == 0,
                 'centered exact native hostile')
            need(all(sum(data[6][i][j]*f[b] for b, j in enumerate(s)) == F(8, 15)*f[a]
                     for a, i in enumerate(s)), 'native eigenvalue 8/15')
            exact_four = {
                'source_masks': [data[0][i] for i in s],
                'source_masses': [str(data[3][i]) for i in s],
                'target_masks': data[1],
                'incidence': [[str(v) for v in row] for row in data[2]],
                'nontrivial_component_kernel': [[str(data[6][i][j]) for j in s] for i in s],
                'nonconstant_eigenvector': list(map(str, f)),
                'nonconstant_eigenvalue': '8/15', 'component_gap': '7/15',
            }
            paired_native = coupling(Counter((min(gc[x], gc[x ^ fullm]),
                                              min(tc[native[x]], tc[native[x] ^ fulln]))
                                             for x in range(size)))
            need(len(paired_native[-1]) == 1, 'paired native connected')
            need(gap_certificate(paired_native, F(8, 15)), 'paired native gap lower bound')
            need(not gap_certificate(paired_native, F(8, 15)+F(1, 1000)),
                 'paired native exact gap')
            exact_four['paired_source_complement_native_converse_gap'] = '8/15'
        census.append(report)

    # Actual first-vertex marking creates a shared degree observable on the full cube.
    marked = []
    for n in range(3, 6):
        pairs, idx, perms, actions = old.setup(n)
        fixed = [a for p, a in zip(perms, actions) if p[0] == 0]
        size = 1 << len(pairs)
        gc, _, _ = old.orbit_table(size, fixed, False)
        tc, _, _ = old.orbit_table(size, fixed, True)
        degrees = [sum(x >> idx[0, j] & 1 for j in range(1, n)) for x in range(size)]
        for classes in (gc, tc):
            seen = {}
            for x, g in enumerate(classes):
                seen.setdefault(g, degrees[x])
                need(seen[g] == degrees[x], 'first-root degree is jointly observable')
        data = coupling(Counter((gc[x], tc[x]) for x in range(size)))
        need(not gap_certificate(data, F(8, 9)), 'root marking breaks full-cube bound')
        need(len(data[-1]) >= n, 'at least one marked component per degree')
        boundary = []
        for degree in (0, n-1):
            local = coupling(Counter((gc[x], tc[x]) for x in range(size) if degrees[x] == degree))
            need(gap_certificate(local, F(8, 9)), 'source/sink root stratum retains gap')
            boundary.append(finite_summary(local))
        marked.append({'n': n, 'components': len(data[-1]), 'root_source_sink_strata': boundary})
    return {'status': 'FINITE-EXACT', 'checks': CHECKS,
            'universe': 'every supplied-path free-bit word at native n=3..6; every first-root-marked full word at n=3..5',
            'native_and_auxiliary': census, 'exact_order_four_hostile': exact_four,
            'root_marking': marked,
            'scope': 'native electrical mod-8 receipts are not replaced by a toy residue; these are exact supplied-model controls'}


if __name__ == '__main__':
    result = run()
    receipt = HERE / 'receipt.json'
    if '--write' in sys.argv:
        receipt.write_text(json.dumps(result, indent=2, sort_keys=True)+'\n', encoding='utf-8')
    else:
        need(result == json.loads(receipt.read_text(encoding='utf-8')), 'saved exact receipt')
    print(json.dumps({'status': 'PASS', 'checks': CHECKS,
                      'exact_order_four_hostile': result['exact_order_four_hostile'],
                      'native_components': [(x['n'], x['native']['component_source_sizes'])
                                            for x in result['native_and_auxiliary']]}, indent=2))
