#!/usr/bin/env python3
"""All-phase event sieve for source-owned cubic rotation layers.

Finite residue roots are ALL quadratic-extension roots, not a finite norm-one
circle.  Residue rejection proves a bound on actual characteristic-zero
constraints; a residue hit is never promoted to a geometric edge.
"""
from collections import Counter
from dataclasses import dataclass
from fractions import Fraction
from itertools import combinations
from math import isqrt
import hashlib
import json

from cubic_boundary_search_20261018a import Context, encode_element
from golden_boundary_20261017c import PORT_EDGES, patterns


def need(test, message):
    if not test:
        raise ValueError(message)


def integer(x, low=0):
    need(type(x) is int and x >= low, 'exact integer in declared domain')
    return x


def prime(p):
    integer(p, 3)
    need(p % 2 and all(p % d for d in range(3, isqrt(p)+1, 2)), 'odd prime')
    return p


@dataclass(frozen=True)
class Place:
    p: int
    z: int
    j: int
    w: int


PLACES = (Place(10321, 3678, 1363, 394),
          Place(10531, 846, 3442, 9367),
          Place(100621, 2204, 9259, 47934))


class QuadraticRoots:
    def __init__(self, p):
        self.p = prime(p)
        # Fixed canonical square roots. No root-finding heuristic is used.
        self.squares = {}
        for x in range(p):
            self.squares.setdefault(x*x % p, x)
        self.nu = next(x for x in range(2, p) if x not in self.squares)
        self.inv_nu = pow(self.nu, -1, p)

    def mul(self, x, y):
        a, b = x; c, d = y
        p = self.p
        return ((a*c+self.nu*b*d) % p, (a*d+b*c) % p)

    def evaluate(self, coeff, z):
        out = (0, 0)
        for c in reversed(coeff):
            out = self.mul(out, z)
            out = ((out[0]+c) % self.p, out[1])
        return out

    def roots(self, coeff):
        need(type(coeff) is tuple and len(coeff) == 3
             and all(type(x) is int for x in coeff), 'typed quadratic coefficients')
        p = self.p
        A, B, C = (x % p for x in coeff)
        need(C != 0, 'unit leading coefficient')
        inv = pow(2*C, -1, p)
        center = -B*inv % p
        disc = (B*B-4*A*C) % p
        if disc in self.squares:
            root = self.squares[disc]*inv % p
            out = ((center+root, 0), (center-root, 0))
        else:
            root = self.squares[disc*self.inv_nu % p]*inv % p
            out = ((center, root), (center, -root))
        return tuple(sorted({(x % p, y % p) for x, y in out}))

    def field_norm(self, root):
        a, b = root
        return (a*a-self.nu*b*b) % self.p


class Engine:
    def __init__(self):
        self.ctx = Context()
        self.f = self.ctx.f
        self.base = self.ctx.ports
        self.edges = tuple(self.ctx.geo.unit_edges(self.base))
        self.missing = tuple(e for e in combinations(range(8), 2) if e not in self.edges)
        self.pair_index = {e: i for i, e in enumerate(self.missing)}
        # Exact source difference catalogue, retaining every oriented address.
        cat = {}
        for a in range(8):
            for b in range(8):
                e = tuple(sorted((a, b)))
                if a != b and e in self.pair_index:
                    d = self.f.sub(self.base[b], self.base[a])
                    cat.setdefault(d, []).append((a, b))
        self.catalogue = tuple((d, tuple(addresses),
                                sum(1 << self.pair_index[tuple(sorted(ab))]
                                    for ab in addresses))
                               for d, addresses in sorted(cat.items()))
        # In this source each oriented difference has a single provenance.
        need(len(self.catalogue) == 38 and all(len(x[1]) == 1 for x in self.catalogue),
             'marked original difference catalogue')

    def place_basis(self, place):
        need(type(place) is Place, 'marked place receipt')
        p = prime(place.p)
        roots = (place.z, place.j, place.w)
        for x in roots:
            integer(x)
            need(x < p, 'normalized place generator')
        for poly, x in zip(self.f.polys, roots):
            need(sum(c*pow(x, i, p) for i, c in enumerate(poly)) % p == 0,
                 'all source field relations at the place')
        need(pow(place.z, 15, p) == 1 and pow(place.z, 3, p) != 1
             and pow(place.z, 5, p) != 1, 'primitive cyclotomic marking')
        conjugates = (pow(place.z, -1, p), -place.j % p, place.w)
        def basis(rs):
            return tuple(pow(rs[0], a, p)*pow(rs[1], b, p)*pow(rs[2], c, p) % p
                         for a, b, c in self.f.basis)
        return basis(roots), basis(conjugates)

    def reduce(self, value, p, basis):
        self.ctx.element(value)
        answer = 0
        for q, e in zip(value, basis):
            need(q.denominator % p != 0, 'coefficient integral at chosen place')
            answer += q.numerator*pow(q.denominator, -1, p)*e
        return answer % p

    def residue_source(self, place, scales, length):
        need(type(scales) is tuple and scales in (('b',), ('c',), ('b', 'c')),
             'named ordered radius bank')
        integer(length, 1)
        need(length <= 64, 'declared finite bank cap64')
        basis, conjugates = self.place_basis(place)
        p = place.p
        def red(x):
            return self.reduce(x, p, basis), self.reduce(x, p, conjugates)
        U, Ub = red(self.ctx.cubic_unit)
        need(U*Ub % p == 1, 'source unit and its actual conjugate')
        base = tuple(red(x) for x in self.base)
        layers = []
        labels = []
        for name in scales:
            k, kb = red(self.ctx.coefficients[name])
            need(k and kb, 'scale is a place unit')
            for i in range(length):
                layers.append((k*pow(U, i, p) % p, kb*pow(Ub, i, p) % p))
                labels.append((name, i))
        differences = []
        for d, provenance, mask in self.catalogue:
            x, xb = red(d)
            need(x and xb, 'all oriented nonedge differences are place units')
            differences.append((x, xb, provenance, mask))
        return base, tuple(layers), tuple(labels), tuple(differences)

    def compile(self, place, scales, length):
        base, layers, labels, catalogue = self.residue_source(place, scales, length)
        p = place.p
        roots = QuadraticRoots(p)
        buckets = {}
        loop_hits, collision_hits, bad_leads = [], [], []
        count = 0
        first_outside = None
        anchor_norms = tuple(((a, b), (base[b][0]-base[a][0]) *
                              (base[b][1]-base[a][1]) % p)
                             for a, b in combinations(range(8), 2))
        for i, j in combinations(range(len(layers)), 2):
            v = (layers[j][0]-layers[i][0]) % p
            vb = (layers[j][1]-layers[i][1]) % p
            if not v or not vb:
                bad_leads.append((i, j))
                continue
            nv = v*vb % p
            if nv == 1:
                loop_hits.append((i, j))
            collision_hits.extend((i, j, *ab) for ab, nd in anchor_norms if nd == nv)
            for d, db, provenance, mask in catalogue:
                coeff = (d*vb % p, (d*db+nv-1) % p, v*db % p)
                need(coeff[0] and coeff[2], 'constant and leading unit guards')
                rr = roots.roots(coeff)
                count += 1
                for beta in rr:
                    buckets[beta] = buckets.get(beta, 0) | mask
                    if first_outside is None and roots.field_norm(beta) != 1:
                        first_outside = {'root': beta, 'field_norm': roots.field_norm(beta),
                                         'coefficients': coeff, 'layers': (i, j),
                                         'anchor_addresses': provenance}
        histogram = Counter(mask.bit_count() for mask in buckets.values())
        maximum = max(histogram, default=0)
        max_roots = [(beta, mask) for beta, mask in sorted(buckets.items())
                     if mask.bit_count() == maximum]
        digest = hashlib.sha256()
        for (a, b), mask in sorted(buckets.items()):
            digest.update(f'{a},{b}:{mask}\n'.encode('ascii'))
        certified = not bad_leads and not loop_hits and not collision_hits
        # A convenient maximum bucket is only a residue obstruction, not geometry.
        return {'place': {'p': p, 'z': place.z, 'j': place.j, 'w': place.w},
                'scales': scales, 'length': length, 'layer_labels': labels,
                'vertex_count': 8*len(layers), 'layer_count': len(layers),
                'event_count': count, 'distinct_residue_roots': len(buckets),
                'root_bucket_histogram': dict(sorted(histogram.items())),
                'maximum_new_pairs': maximum, 'maximum_buckets': max_roots,
                'residue_root_extension_nonsquare': roots.nu,
                'all_bucket_sha256': digest.hexdigest(),
                'outside_finite_norm_one_circle': first_outside,
                'bad_leading_layers': bad_leads, 'loop_candidates': loop_hits,
                'collision_candidates': collision_hits,
                'all_phase_projection_certified': certified,
                'separating_support_certified': certified and maximum <= 12}


def fractional_packet(ctx):
    """One exact radius-adjusted phase; no assertion about other phases."""
    f = ctx.f
    k = f.scalar(2)
    delta = f.sub(ctx.geo.eta, f.one)
    root5 = f.add(f.one, f.scale(f.add(ctx.geo.alpha, f.power(ctx.geo.alpha, 4)), 2))
    root_minus15 = f.mul(root5, f.add(f.one, f.scale(delta, 2)))
    q = f.scale(f.mul(delta, f.sub(root_minus15, f.one)), Fraction(1, 4))
    need(ctx.norm(q) == f.one, 'exact unit phase')
    T = ctx.reflection_link(f.one, k)
    S = ctx.reflection_link(ctx.geo.eta, k)
    R = ctx.rotation(ctx.cubic_unit)
    maps = (ctx.identity, T, S, ctx.compose(T, S), ctx.compose(S, T),
            R, ctx.compose(T, R), ctx.compose(S, R))
    values = tuple(ctx.map_value(m, q) for m in maps)
    shifts = tuple(f.mul(k, f.sub(v, values[0])) for v in values)
    need(shifts[1] == delta and ctx.norm(delta) == f.one, 'fractional unit-translation identity')
    points, aliases, lookup = [], [], {}
    for shift in shifts:
        for c in ctx.ports:
            point = f.add(c, shift)
            if point not in lookup:
                lookup[point] = len(points)
                points.append(point)
            aliases.append(lookup[point])
    need(tuple(points[:8]) == ctx.ports, 'literal original ports retained after coalescence')
    edges = tuple(ctx.geo.unit_edges(points))
    adjacent = [set() for _ in points]
    for a, b in edges:
        adjacent[a].add(b); adjacent[b].add(a)
    # A relative degeneracy certificate proves every source colouring extends.
    live = set(range(len(points)))
    order = []
    while len(live) > 8:
        choices = [v for v in sorted(live) if v >= 8 and len(adjacent[v] & live) <= 5]
        need(choices, 'this packet needs a stronger boundary solver')
        v = choices[0]
        order.append(v)
        live.remove(v)
    witnesses = []
    for source in patterns(8, PORT_EDGES):
        colours = list(source) + [None]*(len(points)-8)
        for v in reversed(order):
            used = {colours[w] for w in adjacent[v] if colours[w] is not None}
            colours[v] = next(c for c in range(6) if c not in used)
        need(tuple(colours[:8]) == source and all(colours[a] != colours[b] for a, b in edges),
             'literal source-bound greedy six-colour lift')
        witnesses.append(''.join(map(str, colours)))
    return {'radius': 2, 'phase': encode_element(q),
            'map_words': ('id', 'T1', 'Teta', 'T1*Teta', 'Teta*T1', 'Ru', 'T1*Ru', 'Teta*Ru'),
            'map_coefficients': [[encode_element(x) for x in m] for m in maps],
            'raw_vertex_count': len(aliases), 'vertex_count': len(points),
            'raw_labels': [[i, c] for i in range(8) for c in range(8)],
            'raw_to_vertex': aliases, 'coordinates': [encode_element(x) for x in points],
            'edges': edges, 'relative_five_degeneracy_order': order,
            'six_port_witnesses': witnesses}, points, maps, q


def verify():
    engine = Engine()
    checks = 0
    def check(ok, message):
        nonlocal checks
        checks += 1
        if not ok:
            raise RuntimeError(message)
    def rejects(call):
        try:
            call()
        except ValueError:
            check(True, 'invalid packet rejected')
        else:
            check(False, 'invalid packet accepted')
    check(len(engine.missing) == 19 and len(engine.edges) == 9, 'source graph')
    # Independent complete Fp2 enumeration for every small quadratic.
    for p in (3, 5, 7):
        R = QuadraticRoots(p)
        points = tuple((a, b) for a in range(p) for b in range(p))
        for A in range(p):
            for B in range(p):
                for C in range(1, p):
                    coeff = (A, B, C)
                    brute = tuple(x for x in points if R.evaluate(coeff, x) == (0, 0))
                    check(R.roots(coeff) == brute, 'all quadratic roots independently enumerated')
    # A common-root plant with13 distinct pair labels must survive the threshold.
    R = QuadraticRoots(10321)
    plant = {}
    for label in range(13):
        a = label+2
        coeff = (a, -a-1, 1)  # (q-1)(q-a)
        for beta in R.roots(coeff):
            check(R.evaluate(coeff, beta) == (0, 0), 'planted root authenticity')
            plant[beta] = plant.get(beta, 0) | (1 << label)
    check(plant[(1, 0)].bit_count() == 13, 'positive common-factor threshold control')
    # Opposite guard hostile: a relevant root need not have finite-field norm1.
    check((2, 0) in R.roots((2, -3, 1)) and R.field_norm((2, 0)) != 1,
          'finite-circle restriction loses a legitimate event root')
    # Dropped leading units can hide a nonintegral root at infinity.
    p = 3
    nonintegral = (1, -1, 3)  # one root of negative3-adic valuation.
    check(nonintegral[2] % p == 0 and nonintegral[0] % p != 0,
          'leading-unit hypothesis is material')
    rejects(lambda: QuadraticRoots(3).roots(nonintegral))

    rows = []
    for place in PLACES[:2]:
        for scales in (('b',), ('c',), ('b', 'c')):
            result = engine.compile(place, scales, 16)
            check(result['separating_support_certified'], 'independent small-bank source certificate')
            check(result['event_count'] == 38*result['layer_count']*(result['layer_count']-1)//2,
                  'complete oriented difference-event universe')
            row = result['outside_finite_norm_one_circle']
            R = QuadraticRoots(place.p)
            check(row is not None and R.evaluate(tuple(row['coefficients']), tuple(row['root'])) == (0, 0)
                  and row['field_norm'] != 1, 'actual source-event finite-circle hostile')
            rows.append(result)
    for length in (32, 64):
        result = engine.compile(PLACES[2], ('b', 'c'), length)
        check(result['separating_support_certified'], 'large-bank all-phase source certificate')
        check(result['maximum_new_pairs'] == (6 if length == 32 else 10), 'declared exact residue bound')
        check(result['event_count'] == 38*result['layer_count']*(result['layer_count']-1)//2,
              'large-bank full event universe')
        rows.append(result)

    # Source field/place authentication before any finite-field computation.
    rejects(lambda: engine.compile(PLACES[0], ('b',), True))
    rejects(lambda: engine.compile(PLACES[0], ('b',), 65))
    rejects(lambda: engine.compile(PLACES[0], ('c', 'b'), 16))
    rejects(lambda: engine.compile(Place(10321.0, 3678, 1363, 394), ('b',), 16))
    rejects(lambda: engine.compile(Place(10321, 3678, 1363, 395), ('b',), 16))
    rejects(lambda: QuadraticRoots(9))
    rejects(lambda: QuadraticRoots(3).roots((True, 1, 1)))
    ctx = engine.ctx
    fractional, points, maps, phase = fractional_packet(ctx)
    check(fractional['raw_vertex_count'] == 64 and len(points) == 56
          and len(fractional['edges']) == 74, 'complete fractional candidate census')
    check(fractional['raw_to_vertex'][32:40] == list(range(8)),
          'phase-specific Teta*T1 aliases retained')
    check(not ctx.same_map(maps[4], ctx.identity), 'phase fixed point is not a rational map identity')
    check(tuple(e for e in fractional['edges'] if e[1] < 16) == ctx.exact_edges(points[:16]),
          'independent complete two-layer incidence')
    check(len([e for e in fractional['edges'] if e[1] < 16]) == 27,
          'two-layer loop-escape count')
    for i in range(8):
        check((i, 8+i) in fractional['edges'], 'actual same-anchor unit edge')
    check(len(fractional['relative_five_degeneracy_order']) == 48, 'complete extension certificate')
    for source, text in zip(patterns(8, PORT_EDGES), fractional['six_port_witnesses']):
        colours = tuple(map(int, text))
        check(colours[:8] == source, 'all labelled source patterns')
        for a, b in fractional['edges']:
            check(colours[a] != colours[b], 'explicit fractional edge witness')
    return {'status': 'FINITE-EXACT residue certificate proving bounded-family all-phase exclusions',
            'checks': checks,
            'source': {'field_polynomials': ctx.f.polys, 'basis': ctx.f.basis,
                       'anchors': [encode_element(x) for x in engine.base],
                       'edges': engine.edges, 'nonedges': engine.missing,
                       'cubic_unit': encode_element(ctx.cubic_unit),
                       'radii': {name: encode_element(k) for name, k in ctx.coefficients.items()},
                       'difference_catalogue': [{'value': encode_element(d), 'provenance': ab, 'pair_mask': mask}
                                                for d, ab, mask in engine.catalogue]},
            'synthetic_common_root_control': {'distinct_pairs': 13, 'phase': (1, 0),
                                              'geometric_claim': False},
            'banks': rows, 'fractional_candidate': fractional}


if __name__ == '__main__':
    print(json.dumps(verify(), indent=2, sort_keys=True))
