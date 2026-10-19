#!/usr/bin/env python3
"""Exact support controls for the two-path complement repair at p=29.
Only the 348 potentially crowded geometries are enumerated. No solver.
"""
from collections import Counter
from pathlib import Path
import json

P = 29
DS = ((1, 0), (0, 1), (3, 9))


def add(x, y):
    return ((x[0] + y[0]) % P, (x[1] + y[1]) % P)


def scale(t, x):
    return ((t * x[0]) % P, (t * x[1]) % P)


def cross(x, y):
    return (x[0] * y[1] - x[1] * y[0]) % P


def derivative(f, d):
    out = {}
    for x, value in f.items():
        for y, sign in ((add(x, scale(-1, d)), 1), (x, -1)):
            out[y] = out.get(y, 0) + sign * value
    return {x: value for x, value in out.items() if value}


def triple(f):
    for d in DS:
        f = derivative(f, d)
    return f


def columns(paths, e_index, component):
    # For d!=e use the derivative-grid columns parallel to e. For d=e,
    # either other derivative direction gives a valid column coordinate.
    d = DS[component]
    axis = DS[e_index] if component != e_index else DS[(e_index+1) % 3]
    supports = [set(derivative(path, d)) for path in paths]
    counts = Counter(cross(axis, x) for x in set().union(*supports))
    return {'columns': len(counts), 'maximum_positions': max(counts.values()),
            'shared_columns': len({cross(axis, x) for x in supports[0]} &
                                  {cross(axis, x) for x in supports[1]})}


# Quotient projections are nonzero and unequal up to sign.
quotients = []
for e in DS:
    projections = [cross(e, d) for d in DS if d != e]
    assert all(projections)
    assert projections[0] not in (projections[1], -projections[1] % P)
    quotients.append(projections)

records = []
repaired_count = 0
for ei, e in enumerate(DS):
    for di, d in enumerate(DS):
        if ei == di:
            continue
        for sign in (-1, 1):
            for phase in range(P):
                # Normalize P_1=0. Potential crowding requires k_1=k_2=14
                # and the two e-lines to differ by +/-d modulo <e>.
                origin2 = add(scale(sign, d), scale(phase, e))
                first = {scale(j, e): 1 for j in range(1, 15)}
                second = {add(origin2, scale(j, e)): 1 for j in range(1, 15)}
                before = [columns([first, second], ei, i) for i in range(3)]
                shared_components = [i for i in range(3) if i != ei and before[i]['shared_columns']]
                assert shared_components == [di]
                crowded = before[di]['maximum_positions'] == 28
                if crowded:
                    repaired_count += 1
                    assert phase in (14, 15)
                    # first - indicator(full e-line): minus the 15-point complement.
                    repaired = {scale(j, e): -1 for j in range(P) if j not in range(1, 15)}
                    assert len(repaired) == 15
                    assert derivative(first, e) == derivative(repaired, e)
                    assert triple(first) == triple(repaired)
                    after_paths = [repaired, second]
                else:
                    after_paths = [first, second]
                after = [columns(after_paths, ei, i) for i in range(3)]
                for i in range(3):
                    assert after[i]['columns'] <= 4
                    assert after[i]['maximum_positions'] <= P-2
                    assert P > 2*after[i]['columns']
                    if crowded and i != ei:
                        assert after[i]['maximum_positions'] <= 15
                    for path in after_paths:
                        delta = derivative(path, DS[i])
                        assert sum(delta.values()) == 0
                        assert all(sum(value for x, value in delta.items() if cross(DS[i], x) == t) == 0
                                   for t in range(P))
                records.append({'exceptional_direction': ei, 'possibly_shared_component': di,
                                'line_shift_sign': sign, 'relative_phase': phase,
                                'repair_needed': crowded, 'before': before, 'after': after})
assert len(records) == 3*2*2*29 == 348
assert repaired_count == 3*2*2*2 == 24
out = {
    'p': P, 'directions': DS, 'quotient_projections': quotients,
    'potentially_crowded_geometries': len(records), 'complement_repairs': repaired_count,
    'claim': 'Two same-direction dipoles on distinct cycles always admit paths satisfying column normalization.',
    'consequence': 'For two disjoint cubes at p=29, homogeneous incidence dimension four suffices for exclusion.',
    'numerical_lower_bounds_unchanged': {'triple': 89, 'full_graph': 445},
    'certificates': records
}
Path(__file__).with_name('path-complement-exact-results.json').write_text(json.dumps(out, indent=2)+'\n')
print('348 potentially crowded geometries checked; exactly 24 require repair.')
print('At most one component shares correction columns, by unequal quotient projections.')
print('Every repair preserves its endpoint dipole and mixed cube exactly over the integers.')
print('Every repaired nonexceptional component has at most 15 positions per column.')
print('All 348 final support patterns satisfy the r-column normalization bounds.')
print('Two disjoint cubes with homogeneous incidence dimension four are therefore excluded without an extra sparsity hypothesis.')
print('Numerical bounds remain 89 and 445; interacting or overlapping two-cycle cases remain open.')
