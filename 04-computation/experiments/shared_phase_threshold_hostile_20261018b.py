#!/usr/bin/env python3
"""Exact metric hostile: thirteen quotient constraints, no retraction, full support.

The translations are engineered mixed radii, not the fixed b/c search bank.
All graph edges are authenticated; no drawing or numerical tolerance is used.
"""
import argparse
import hashlib
import itertools
import json
from pathlib import Path

from cubic_boundary_search_20261018a import Context
from exceptional_factor_support_20261018a import BASE_SET, partitions


def build():
    ctx = Context()
    f, C, j = ctx.f, ctx.ports, ctx.geo.j
    terminal = (2, 4)
    W = [0, 2, 3, 4, 5, 6, 7]
    added = sorted(set(itertools.combinations(W, 2)) - BASE_SET - {terminal})
    assert len(added) == 13 and (2, 6) in added
    pairs = [(6, 2) if e == (2, 6) else e for e in added]
    shifts = [f.zero]
    ms, spokes = [], []
    for layer, (a, b) in enumerate(pairs, 1):
        m = 1000 + 37 * layer
        # Rational parametrization in Q(sqrt(-11)); complex norm is exactly one.
        unit = f.mul(f.add(f.scalar(m), j), f.inv(f.sub(f.scalar(m), j)))
        assert ctx.norm(unit) == f.one
        shifts.append(f.sub(unit, f.sub(C[b], C[a])))
        ms.append(m)
        spokes.append((a, 8 * layer + b))
    points = [f.add(c, t) for t in shifts for c in C]
    assert len(points) == len(set(points)) == 112
    # The residue rejection is exact only at places where every denominator is a unit.
    for p, _, _ in ctx.geo.pl:
        assert all(co.denominator % p for z in points for co in z)
    edges = set(ctx.geo.unit_edges(points))
    expected = {(8 * k + a, 8 * k + b) for k in range(14) for a, b in BASE_SET}
    expected.update(spokes)
    if edges != expected:
        raise AssertionError({'extra': sorted(edges - expected),
                              'missing': sorted(expected - edges)})
    assert len(edges) == 139
    projected = {tuple(sorted((a % 8, b % 8))) for a, b in edges}
    assert projected - BASE_SET == set(added)
    return ctx, points, shifts, pairs, ms, sorted(edges), terminal


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument('--out', type=Path, required=True)
    args = ap.parse_args()
    ctx, points, shifts, pairs, ms, edges, terminal = build()
    print('exact geometry: 112 points, 139 edges, 13 new projected pairs', flush=True)
    pats = partitions()
    quotient = BASE_SET | {tuple(sorted(e)) for e in pairs}
    quotient_support = [list(p) for p in pats if all(p[a] != p[b] for a, b in quotient)]
    assert quotient_support and all(p[2] == p[4] for p in quotient_support)
    witnesses = []
    for p in pats:
        # Three pair labels remain fixed; independently flip each copy's shade.
        bits = [0] + [int(p[a] == p[b]) for a, b in pairs]
        colors = [c ^ bit for bit in bits for c in p]
        assert colors[:8] == list(p)
        assert all(colors[a] != colors[b] for a, b in edges)
        witnesses.append({'source': list(p), 'copy_shade_flips': bits,
                          'coloring': colors})
    assert len(witnesses) == 570
    # Independent finite check behind the nonretraction argument.
    adjacency = [set() for _ in range(8)]
    for a, b in BASE_SET:
        adjacency[a].add(b); adjacency[b].add(a)
    cycle_maps = [list(w) for w in itertools.product(range(8), repeat=5)
                  if all(w[(i+1) % 5] in adjacency[w[i]] for i in range(5))]
    assert len(cycle_maps) == 10
    assert all(set(w) == set(range(5)) for w in cycle_maps)
    assert adjacency[6] == {5, 7}
    assert not adjacency[6] & set(range(5))
    # Positive/hostile support controls are separate from metric realizability.
    separating = sum(p[2] != p[4] for p in pats)
    assert separating > 0 and separating < len(pats)
    # Record the complete exact basis coefficients, including zero entries.
    encoded = [[str(co) for co in z] for z in points]
    geometry_bytes = json.dumps(encoded, separators=(',', ':')).encode()
    out = {
        'status': 'FINITE-EXACT',
        'scope': 'engineered mixed-radius shared phase q=1; not the fixed b/c bank',
        'basis_dependency': 'golden_forcing_search_20261017.Geo / make(sqrt3)',
        'point_count': len(points), 'exterior_to_original_eight': len(points)-8,
        'edge_count': len(edges), 'edges': edges,
        'new_projected_pairs': sorted(quotient-BASE_SET), 'terminal': terminal,
        'unit_parameters_m': ms, 'oriented_spokes': pairs,
        'exact_coordinates': encoded,
        'coordinates_sha256': hashlib.sha256(geometry_bytes).hexdigest(),
        'good_places': [p for p, _, _ in ctx.geo.pl],
        'quotient_support': quotient_support,
        'physical_support_count': len(witnesses),
        'equal_terminal_support': len(pats)-separating,
        'separating_terminal_support': separating,
        'source_fixed_retraction': False,
        'pentagon_closed_walks': cycle_maps,
        'witnesses': witnesses,
    }
    args.out.parent.mkdir(parents=True, exist_ok=True)
    args.out.write_text(json.dumps(out, indent=2)+'\n', encoding='utf-8')
    print(json.dumps({k: out[k] for k in ('point_count','edge_count',
        'physical_support_count','equal_terminal_support','separating_terminal_support',
        'source_fixed_retraction')}, indent=2), flush=True)


if __name__ == '__main__':
    main()
