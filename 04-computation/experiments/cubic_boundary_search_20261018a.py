#!/usr/bin/env python3
"""Exact generic-phase obstruction and exceptional cubic-layer controls.

No search result is promoted to a six-colour equality gadget.  The field and
eight literal ports are inherited from golden_forcing_search_20261017.
"""
from fractions import Fraction as Q
from functools import lru_cache
from itertools import combinations
from pathlib import Path
import hashlib
import json

from golden_forcing_search_20261017 import Geo
from golden_boundary_20261017c import PORT_EDGES, patterns


def require(ok, message):
    if not ok:
        raise ValueError(message)


class Context:
    """All coefficients belong to the fixed, marked 48-dimensional field F."""
    def __init__(self):
        self.geo = Geo()
        self.f = f = self.geo.f
        self.zero, self.one = f.zero, f.one
        self.identity = (f.one, f.zero, f.zero, f.one)
        self.ports = tuple(self.geo.ports)
        a = self.geo.bar(self.geo.d)
        m2 = f.scale(f.add(f.scalar(20), f.scale(self.geo.j, 11)), Q(-1, 19))
        m1 = f.scale(f.sub(f.scalar(20), f.scale(self.geo.j, 11)), Q(1, 19))
        c = f.mul(f.add(f.sub(f.add(f.mul(a, m2), f.one),
                                  f.mul(f.mul(a, a), m1)), f.power(a, 3)),
                  f.inv(f.sub(f.mul(self.geo.d, a), f.one)))
        b = f.add(f.mul(c, self.geo.d), f.sub(f.mul(a, m1), f.mul(a, a)))
        self.coefficients = {'b': b, 'c': c}
        w = f.gen(2)
        self.cubic_unit = f.mul(f.sub(f.add(w, f.scalar(3)), self.geo.j),
                               f.inv(f.add(f.add(w, f.scalar(3)), self.geo.j)))

    def element(self, x):
        require(type(x) is tuple and len(x) == self.f.n
                and all(type(v) is Q for v in x), 'exact marked field element')
        return x

    def bar(self, x):
        self.element(x)
        return self._bar(x)

    @lru_cache(None)
    def _bar(self, x):
        return self.geo.bar(x)

    def norm(self, x):
        self.element(x)
        return self._norm(x)

    @lru_cache(None)
    def _norm(self, x):
        return self.f.mul(x, self.bar(x))

    def trim(self, p):
        p = list(p)
        while p and p[-1] == self.zero:
            p.pop()
        return tuple(p)

    def add(self, p, q):
        return self.trim(self.f.add(p[i] if i < len(p) else self.zero,
                                    q[i] if i < len(q) else self.zero)
                         for i in range(max(len(p), len(q))))

    def scale(self, p, c):
        return self.trim(self.f.mul(x, c) for x in p)

    def sub(self, p, q):
        return self.add(p, self.scale(q, self.f.scalar(-1)))

    def mul(self, p, q):
        if not p or not q:
            return ()
        out = [self.zero] * (len(p) + len(q) - 1)
        for i, a in enumerate(p):
            for j, b in enumerate(q):
                out[i+j] = self.f.add(out[i+j], self.f.mul(a, b))
        return self.trim(out)

    def star(self, p, grade):
        require(type(grade) is int and grade >= 0 and len(p) <= grade+1,
                'fixed polynomial grade')
        return self.trim(self.bar(p[grade-i]) if grade-i < len(p) else self.zero
                         for i in range(grade+1))

    def evaluate(self, p, z):
        self.element(z)
        out = self.zero
        for c in reversed(p):
            self.element(c)
            out = self.f.add(self.f.mul(out, z), c)
        return out

    def map(self, m):
        require(type(m) is tuple and len(m) == 4, 'four marked map coefficients')
        for x in m:
            self.element(x)
        a, b, c, d = m
        require(self.f.sub(self.f.mul(a, d), self.f.mul(b, c)) != self.zero,
                'nonconstant Mobius map')
        n, den = (b, a), (d, c)
        require(not self.sub(self.mul(n, self.star(n, 1)),
                             self.mul(den, self.star(den, 1))),
                'rational unit-circle identity')
        return m

    def compose(self, left, right):
        self.map(left); self.map(right)
        a, b, c, d = left
        e, h, i, j = right
        f = self.f
        return self.map((f.add(f.mul(a, e), f.mul(b, i)),
                         f.add(f.mul(a, h), f.mul(b, j)),
                         f.add(f.mul(c, e), f.mul(d, i)),
                         f.add(f.mul(c, h), f.mul(d, j))))

    def rotation(self, u):
        self.element(u)
        require(self.norm(u) == self.one, 'unit rotation')
        return self.map((u, self.zero, self.zero, self.one))

    def same_map(self, left, right):
        self.map(left); self.map(right)
        i = next(i for i, x in enumerate(left) if x != self.zero)
        return all(self.f.mul(x, right[i]) == self.f.mul(y, left[i])
                   for x, y in zip(left, right))

    def reflection_link(self, d, k):
        self.element(d); self.element(k)
        require(self.norm(d) == self.one and k != self.zero
                and self.norm(k) != self.one, 'nondegenerate reflected unit link')
        n = self.norm(k)
        f = self.f
        return self.map((n, f.scale(f.mul(d, self.bar(k)), -1),
                         f.scale(f.mul(self.bar(d), k), -1), n))

    def map_value(self, m, z):
        self.map(m); self.element(z)
        a, b, c, d = m
        den = self.f.add(self.f.mul(c, z), d)
        require(den != self.zero, 'phase is not a pole')
        return self.f.mul(self.f.add(self.f.mul(a, z), b), self.f.inv(den))

    @lru_cache(None)
    def _pair(self, left, right):
        # The public caller validates before this cache lookup.
        a, b, c, d = left
        e, h, i, j = right
        den = self.mul((d, c), (j, i))
        diff = self.sub(self.mul((h, e), (d, c)), self.mul((b, a), (j, i)))
        return diff, den

    def displacement_numerator(self, d, k, left, right):
        self.element(d); self.element(k); self.map(left); self.map(right)
        diff, den = self._pair(left, right)
        return self.add(self.scale(den, d), self.scale(diff, k)), den

    def event_polynomial(self, d, k, left, right):
        """Ascending coefficients of N(d+k(g_right-g_left))-1, cleared.

        This is a rational-identity numerator, not an existence certificate.
        At a proposed phase one must also check norm1, poles and coincidences.
        """
        n, den = self.displacement_numerator(d, k, left, right)
        return self.sub(self.mul(n, self.star(n, 2)),
                        self.mul(den, self.star(den, 2)))

    def specialize(self, base, k, maps, phase):
        require(type(base) is tuple and len(base) > 0, 'nonempty exact base')
        require(type(maps) is tuple and len(maps) > 0, 'nonempty exact maps')
        self.element(k); self.element(phase)
        require(k != self.zero and self.norm(phase) == self.one, 'nonzero scale/unit phase')
        for c in base:
            self.element(c)
        require(len(set(base)) == len(base), 'distinct source points')
        values = [self.map_value(m, phase) for m in maps]
        shifts = [self.f.mul(k, self.f.sub(x, values[0])) for x in values]
        pts = tuple(self.f.add(c, shift) for shift in shifts for c in base)
        require(len(set(pts)) == len(pts), 'specialized vertices must be distinct')
        require(pts[:len(base)] == base, 'literal source retained')
        edges = tuple(self.geo.unit_edges(pts))
        return pts, edges

    def exact_edges(self, pts):
        for p in pts:
            self.element(p)
        return tuple((i, j) for i, j in combinations(range(len(pts)), 2)
                     if self.norm(self.f.sub(pts[j], pts[i])) == self.one)


def encode_element(x):
    return [[i, c.numerator, c.denominator] for i, c in enumerate(x) if c]


def encode_poly(p):
    return [encode_element(c) for c in p]


def retraction(n, edges, source_count, source_edges, image):
    require(type(n) is int and type(source_count) is int and 0 < source_count <= n,
            'exact vertex counts')
    require(type(image) is tuple and len(image) == n
            and all(type(x) is int and 0 <= x < source_count for x in image),
            'source-owned vertex image')
    for collection in (edges, source_edges):
        require(type(collection) is tuple, 'canonical edge tuple')
        require(all(type(e) is tuple and len(e) == 2 and
                    all(type(x) is int for x in e) and 0 <= e[0] < e[1] < n
                    for e in collection), 'canonical typed edge')
        require(len(set(collection)) == len(collection), 'distinct edge records')
    require(all(v < source_count for e in source_edges for v in e), 'source edges')
    require(image[:source_count] == tuple(range(source_count)), 'fix all source vertices')
    E = set(source_edges)
    require(E == {e for e in edges if e[1] < source_count}, 'induced source edge list')
    require(all(tuple(sorted((image[a], image[b]))) in E for a, b in edges),
            'every actual edge projects to an actual source edge')
    return True


def designed_packet(ctx, coefficient, copies):
    require(coefficient in ('b', 'c') and type(coefficient) is str, 'named coefficient')
    require(type(copies) is int and copies in (2, 4), 'audited finite packet size')
    f = ctx.f
    k = ctx.coefficients[coefficient]
    powers = tuple(f.power(ctx.cubic_unit, i) for i in range(copies))
    shifts = tuple(f.mul(k, f.sub(u, f.one)) for u in powers)
    base = ctx.ports + tuple(f.sub(powers[i], shifts[i]) for i in range(1, copies))
    maps = tuple(ctx.rotation(u) for u in powers)
    points, edges = ctx.specialize(base, k, maps, f.one)
    return base, maps, points, edges


def verify():
    ctx = Context(); f = ctx.f
    checks = 0
    def check(ok, label):
        nonlocal checks
        checks += 1
        if not ok:
            raise RuntimeError(label)
    def rejects(fn, label):
        try:
            fn()
        except ValueError:
            check(True, label)
        else:
            check(False, label)
    check(f.n == 48, 'inherited field marking')
    check(ctx.norm(ctx.cubic_unit) == f.one, 'cubic unit')
    quintic = (f.scalar(-1), f.scale(f.add(f.one, ctx.geo.j), Q(-1, 2)),
               f.one, f.scalar(-1), f.scale(f.sub(f.one, ctx.geo.j), Q(1, 2)), f.one)
    check(ctx.mul(quintic, tuple(ctx.bar(x) for x in quintic)) == (f.one,)*11,
          'two quintic factors of Phi11')
    pats = tuple(patterns(8, PORT_EDGES))
    old_edges = tuple(ctx.geo.unit_edges(ctx.ports))
    check(len(pats) == 570 and old_edges == tuple(sorted(tuple(sorted(e)) for e in PORT_EDGES)),
          'literal old interface')
    packets = []
    for name in ('b', 'c'):
        k = ctx.coefficients[name]
        check(k != f.zero and ctx.norm(k) != f.one, 'scale hypotheses')
        T1 = ctx.reflection_link(f.one, k)
        Teta = ctx.reflection_link(ctx.geo.eta, k)
        R = ctx.rotation(ctx.cubic_unit)
        # Complete maps; generic theorem applies to all six, independently of samples.
        maps = (ctx.identity, R, T1, Teta, ctx.compose(T1, R), ctx.compose(Teta, R))
        for m in maps:
            check(ctx.map(m) == m, 'circle-map coefficient identity')
        check(all(not ctx.same_map(x, y) for x, y in combinations(maps, 2)),
              'six distinct rational phase maps')
        for d, m in ((f.one, T1), (ctx.geo.eta, Teta)):
            check(not ctx.event_polynomial(d, k, ctx.identity, m), 'reflected link identity')
        check(not ctx.event_polynomial(f.one, k, ctx.identity, ctx.identity), 'old edge identity')
        check(bool(ctx.event_polynomial(f.scalar(2), k, ctx.identity, T1)), 'nonunit projected displacement')
        # Same rational identity after a common right composition.
        check(not ctx.event_polynomial(f.one, k, R, maps[4]), 'right composition retained')
        for copies in (2, 4):
            base, ms, points, edges = designed_packet(ctx, name, copies)
            nB = len(base)
            source_edges = tuple(ctx.geo.unit_edges(base))
            check(source_edges == old_edges, 'all added source anchors isolated')
            check(edges == ctx.exact_edges(points), 'independent unsieved full incidence')
            expected = tuple(sorted([(i*nB+a, i*nB+b) for i in range(copies)
                                      for a, b in old_edges] +
                                     [(0, i*nB+7+i) for i in range(1, copies)]))
            check(edges == expected, 'complete translated copies plus designed spokes')
            check(len(points) == copies*(7+copies), 'distinct labelled source geometry')
            broken = tuple((i, j) for i, j in edges
                           if tuple(sorted((i % nB, j % nB))) not in set(source_edges))
            check(len(broken) == copies-1, 'canonical projection genuinely fails')
            image = tuple(i if i < nB else (i % nB if i % nB < 8 else 1)
                          for i in range(len(points)))
            check(retraction(len(points), edges, nB, source_edges, image), 'repaired full-source retraction')
            witnesses = []
            for p in pats:
                source = p + (0,)*(nB-8)
                colours = tuple(source[v] for v in image)
                check(colours[:8] == p, 'literal labelled port colours')
                for a, b in edges:
                    check(colours[a] != colours[b], 'explicit six-colour edge witness')
                witnesses.append(''.join(map(str, colours)))
            events = []
            for i in range(1, copies):
                d = base[7+i]
                poly = ctx.event_polynomial(d, k, ms[0], ms[i])
                check(bool(poly) and len(poly) <= 5, 'nonidentity low-degree event')
                check(ctx.evaluate(poly, f.one) == f.zero, 'exceptional phase is actual root')
                check(ctx.norm(d) != f.one, 'new edge has nonedge source labels')
                # Fixed-grade clearing adds a zero at 0 to this quadratic event.
                reduced = poly
                while reduced and reduced[0] == f.zero:
                    reduced = reduced[1:]
                check(len(reduced) == 3, 'actual rotation event quadratic')
                events.append({'phase_index': i, 'anchor_index': 7+i,
                               'cleared_event': encode_poly(poly),
                               'after_removing_zero_phase': encode_poly(reduced)})
            # Full graph is 3-colourable and contains the original C5.
            small = (0, 1, 0, 1, 2, 1, 0, 1) + (0,)*(nB-8)
            three = tuple(small[v] for v in image)
            check(all(three[a] != three[b] for a, b in edges), 'exact chromatic3 upper witness')
            check(all(tuple(sorted((i, (i+1) % 5))) in set(edges) for i in range(5)), 'C5 lower witness')
            packets.append({'coefficient': name, 'copies': copies,
                            'source_count': nB, 'vertex_count': len(points),
                            'exterior_to_full_source': len(points)-nB,
                            'beyond_eight_ports': len(points)-8,
                            'labels': [[i, j] for i in range(copies) for j in range(nB)],
                            'coordinates': [encode_element(x) for x in points],
                            'source_edges': source_edges, 'edges': edges,
                            'canonical_projection_failures': broken,
                            'repaired_retraction': image,
                            'three_colouring': three,
                            'six_port_witnesses': witnesses,
                            'events': events})
    second_roots = []
    for name in ('b', 'c'):
        k = ctx.coefficients[name]
        base, maps, _, _ = designed_packet(ctx, name, 4)
        nB = len(base)
        for h in range(1, 4):
            d = base[7+h]
            v = f.mul(k, f.sub(f.power(ctx.cubic_unit, h), f.one))
            A, C = f.mul(d, ctx.bar(v)), f.mul(v, ctx.bar(d))
            phase = f.mul(A, f.inv(C))
            check(ctx.norm(phase) == f.one and phase != f.one, 'distinct second unit root')
            poly = ctx.event_polynomial(d, k, maps[0], maps[h])
            check(ctx.evaluate(poly, phase) == f.zero, 'second event root')
            points, edges = ctx.specialize(base, k, maps, phase)
            check(all(ctx.norm(f.sub(points[b], points[a])) == f.one for a, b in edges),
                  'second-phase exact retained edges; nonedges have good-prime witnesses')
            expected = tuple(sorted([(i*nB+a, i*nB+b) for i in range(4)
                                      for a, b in old_edges] + [(0, h*nB+7+h)]))
            check(edges == expected, 'one and only one exceptional spoke survives')
            image = tuple(i if i < nB else (i % nB if i % nB < 8 else 1)
                          for i in range(len(points)))
            check(retraction(len(points), edges, nB, old_edges, image),
                  'second-phase full-source retraction')
            witnesses = []
            for p in pats:
                source = p + (0,)*(nB-8)
                colours = tuple(source[v] for v in image)
                check(colours[:8] == p and all(colours[a] != colours[b] for a, b in edges),
                      'complete second-phase labelled support')
                witnesses.append(''.join(map(str, colours)))
            second_roots.append({'coefficient': name, 'event_index': h,
                                 'phase': encode_element(phase),
                                 'source_count': nB, 'vertex_count': len(points),
                                 'labels': [[i, j] for i in range(4) for j in range(nB)],
                                 'coordinates': [encode_element(x) for x in points],
                                 'source_edges': old_edges, 'edges': edges,
                                 'repaired_retraction': image,
                                 'six_port_witnesses': witnesses})
    # The no-chord hypothesis is material: in this triangle-containing field,
    # k=1 and omega=eta give N(omega-1)=1, an edge projecting to a loop.
    omega = ctx.geo.eta
    check(ctx.norm(f.sub(omega, f.one)) == f.one, 'allowed translation chord hostile')
    check(not ctx.event_polynomial(f.zero, f.one, ctx.identity, ctx.rotation(omega)),
          'same-anchor unit edge is a rational identity without no-chord gate')
    rejects(lambda: ctx.map((f.one, f.zero, f.zero, f.zero)), 'constant map rejected')
    rejects(lambda: ctx.map((f.scalar(2), f.zero, f.zero, f.one)), 'nonunit circle map rejected')
    rejects(lambda: ctx.element((True,)+f.zero[1:]), 'boolean field alias rejected')
    rejects(lambda: ctx.element(f.one[:-1]), 'wrong field dimension rejected')
    ctx.norm(f.one); ctx.bar(f.one)
    rejects(lambda: ctx.norm(tuple(float(x) for x in f.one)), 'warm-cache norm float alias rejected')
    rejects(lambda: ctx.bar(tuple(float(x) for x in f.one)), 'warm-cache bar float alias rejected')
    rejects(lambda: ctx.norm((True,)+f.zero[1:]), 'warm-cache norm boolean alias rejected')
    rejects(lambda: ctx.bar((True,)+f.zero[1:]), 'warm-cache bar boolean alias rejected')
    rejects(lambda: ctx.specialize(ctx.ports, f.one, (ctx.identity, ctx.identity), f.one),
            'coincident specialized layers rejected')
    rejects(lambda: designed_packet(ctx, 'b', True), 'boolean packet count rejected')
    first = packets[0]
    badimage = tuple(False if i == 0 else x for i, x in enumerate(first['repaired_retraction']))
    rejects(lambda: retraction(first['vertex_count'], tuple(first['edges']), first['source_count'],
                               tuple(first['source_edges']), badimage), 'boolean vertex alias rejected')
    return {'status': 'PROVED generic theorem; FINITE-EXACT complete exceptional packets',
            'checks': checks,
            'field': {'polynomials': ctx.f.polys, 'dimension': ctx.f.n,
                      'basis_exponents': ctx.f.basis,
                      'cubic_unit': encode_element(ctx.cubic_unit),
                      'scales': {s: encode_element(k) for s, k in ctx.coefficients.items()}},
            'generic': {'phase_degree_threshold': 5, 'event_degree_bound': 4,
                        'maps_per_scale': 6, 'source_port_count': 8,
                        'formal_vertices': 48, 'old_patterns_retained': 570,
                        'no_chord_dependency': 'THM4700 for b and c'},
            'packets': packets, 'second_roots': second_roots}


if __name__ == '__main__':
    data = verify()
    print(json.dumps(data, indent=2, sort_keys=True))
