#!/usr/bin/env python3
"""Exact controls for atomic conductor debt, palette interfaces, and inverse clocks."""
import hashlib
import itertools as it
import json
from fractions import Fraction as F
from pathlib import Path
import sympy as s

ROOT = Path(__file__).resolve().parents[2]


def need(ok, message):
    if not ok:
        raise AssertionError(message)


def atomic_controls():
    packets = [{0: 1}, {0: 1, 1: -1}, {0: 1, 1: -2, 2: 1},
               {0: 1, 27: -1}, {2: 3, 11: -7, 76: 4}]
    for a in range(16, 32, 4):
        packet = {}
        for n in range(1, 129):
            c = max(0, min(a+3, n-1)-max(a, (n+3)//4)+1)
            if c:
                packet[n] = c*(2 if n % 3 == 2 else -1)
        packets.append(packet)
    count = 0
    for packet in packets:
        A, B = sum(abs(v) for v in packet.values()), sum(v*v for v in packet.values())
        need(A > 0 and B > 0, 'nonzero atomic packet')
        maxv = -1
        for n, m in it.combinations(packet, 2):
            d, v = abs(n-m), 0
            while d % 3 == 0:
                d //= 3
                v += 1
            maxv = max(maxv, v)
        q0 = max(1, maxv+2)
        for q in range(1, 10):
            N = 3**q
            def coarse_square(M):
                totals = {}
                for n, a in packet.items():
                    totals[n % M] = totals.get(n % M, 0)+a
                return M*sum(v*v for v in totals.values())
            energy = coarse_square(N)-coarse_square(N//3)
            def ramanujan(d):
                return 2*(N//3) if d % N == 0 else -(N//3) if d % (N//3) == 0 else 0
            pair_energy = sum(a*b*ramanujan(n-m) for n, a in packet.items()
                              for m, b in packet.items())
            need(energy == pair_energy and energy >= 0, 'primitive Fourier energy two paths')
            if q >= q0:
                need(energy == 2*3**(q-1)*B, 'eventual exact primitive energy')
                need(F(energy, A) > 0, 'nonvanishing L1 lower bound')
            count += 1
    return dict(packets=len(packets), exact_conductor_checks=count,
                hostile_vanishing_moments=[1, -2, 1])


def hazard_controls():
    count = 0
    delta = F(3, 10)
    for u in range(1, 101):
        for v in range(u+1):
            h = F(u-v, u)
            r = max(F(0), F(10*v-7*u, 3))
            integral = F(0) if h >= delta else 2*u*u/delta**2 * (delta-h)**2/2
            need(r*r == integral, 'exact lower-tail hazard identity')
            count += 1
    # Identical conductor and clock marginals do not specify collision fibres.
    need(abs(5-5) == 0 and abs(5)+abs(-5) == 10, 'same marginal mass, different joint fibre')
    return dict(exact_hazard_checks=count, marginal_collision_hostile=[0, 10])


def palette_controls():
    colours = tuple(range(6))
    pairs = list(it.combinations(colours, 2))
    results = []
    for k, ell in [(3, 3), (4, 2), (5, 1)]:
        options = set()
        for perm in it.permutations(colours):
            chunks = tuple(tuple(sorted(perm[2*i:2*i+2])) for i in range(ell))
            chunks += tuple((v,) for v in perm[2*ell:])
            options.add(chunks)
        expected = 720//2**ell
        need(len(options) == expected, 'all labelled disjoint palette assignments')
        blockers = {}
        for j in range(ell):
            for pair in pairs:
                bad = {opt for opt in options if opt[j] == pair}
                need(len(bad)*15 == len(options), 'one distinct class-pair blocker fraction')
                blockers[j, pair] = bad
        # Union-bound witness guarantee for all patterns follows from the exact
        # option count; these controls include duplicates and concentrated debt.
        for offset in range(15):
            marked = [(i % ell, pairs[(i+offset) % 15]) for i in range(9)]
            forbidden = set().union(*(blockers[z] for z in marked))
            need(len(options-forbidden) >= 2*len(options)//5, 'nine blockers retain >=2/5')
        need(set().union(*(blockers[0, p] for p in pairs)) == options,
             '15-blocker obstruction to this palette method only')
        results.append(dict(exterior_colours=k, interface_classes=ell,
                            options=len(options), blocked_per_vertex=len(options)//15,
                            guaranteed_survivors_with_nine=2*len(options)//5))
    # Complete finite set-family census, not a geometric chromatic assertion.
    intersection = []
    for a in pairs:
        intersection.append(sum(1 << j for j, b in enumerate(pairs) if set(a) & set(b)))
    best, maximizers = 0, 0
    allpairs = (1 << 15)-1
    for mask in range(1, 1 << 15):
        other = allpairs
        for j in range(15):
            if mask >> j & 1:
                other &= intersection[j]
        if other:
            total = mask.bit_count()+other.bit_count()
            if total > best:
                best, maximizers = total, 1
            elif total == best:
                maximizers += 1
    need(best == 10, 'two nonempty cross-intersecting pair families exact maximum')
    return dict(palette_census=results, nonempty_pair_family_subsets_checked=32767,
                cross_intersecting_maximum=best, maximizers=maximizers)


def jacobian_controls():
    x, y, z = s.symbols('x y z')
    checked = return_checks = 0
    maps = [(x, y), (x+y*y, y+(x+y*y)**2),
            (x+y**3, y+(x+y**3)**2)]
    for f, g in maps:
        fx, fy, gx, gy = [s.diff(p, v) for p, v in [(f,x), (f,y), (g,x), (g,y)]]
        need(s.expand(fx*gy-fy*gx) == 1, 'actual Keller source')
        def D(h):
            return s.expand(fx*s.diff(h, y)-fy*s.diff(h, x))
        def E(h):
            return s.expand(gy*s.diff(h, x)-gx*s.diff(h, y))
        def R(h):
            return s.expand(fx*s.integrate(h, (x, 0, x)))
        for r in [1, x, y, x*y, x*x+y*y]:
            r = s.sympify(r)
            P = s.integrate(r, (x, 0, x))
            need(s.expand(E(R(r))-r-D(-gx*P)) == 0, 'explicit cokernel right inverse')
            # If r=D(h), R(r) is D-exact; use a direct polynomial primitive of
            # the closed one-form produced in the proof, without quotient guesses.
            h = r
            b = R(D(h))
            a = -gx*s.integrate(D(h), (x, 0, x))+h
            # d(a df + b dg)=0; integrate in original Cartesian coordinates.
            px, py = s.expand(a*fx+b*gx), s.expand(a*fy+b*gy)
            primitive = s.integrate(px, (x, 0, x))
            primitive += s.integrate(s.expand(py-s.diff(primitive,y)), (y, 0, y))
            need(s.expand(D(primitive)-b) == 0, 'inverse respects actual image D(A)')
            power, invoice = r, s.Integer(0)
            for m in range(1, 5):
                term = -gx*s.integrate(power, (x,0,x))
                for _ in range(m-1):
                    term = E(term)
                invoice = s.expand(invoice+term)
                power = R(power)
                fs = f.subs(x,z)
                closed = fx*s.integrate((f-fs)**(m-1)*r.subs(x,z), (z,0,x))/s.factorial(m-1)
                need(s.expand(power-closed) == 0, 'iterated Volterra formula')
                if m <= 3:
                    returned = power
                    for _ in range(m):
                        returned = E(returned)
                    need(s.expand(returned-r-D(invoice)) == 0, 'literal clock-return invoice')
                    return_checks += 1
                checked += 1
    # A regular-singular, finite-monodromy Weyl-module control, not a Keller
    # cokernel realization: nabla z^n=(n+1/2)z^(n-1) is bijective on Laurent polys.
    alpha = F(1,2)
    inverse_coeff = F(1)
    for m in range(1, 21):
        inverse_coeff /= m+alpha
        forward = inverse_coeff
        for j in range(1, m+1):
            forward *= j+alpha
        need(forward == 1 and inverse_coeff > 0, 'nonzero small inverse-clock iterates')
    # A second abstract Weyl module retains a nonzero torsion arm while E is
    # bijective on the whole module: E=z multiplication, t=-d/dz.
    for n in range(-8,9):
        h = z**n
        need(s.expand(z*(-s.diff(h,z))+s.diff(z*h,z)-h) == 0, 'Laurent Weyl relation')
    for m in range(1,13):
        h = 1/z
        for _ in range(m-1):
            h = -s.diff(h,z)
        need(s.simplify(h/s.factorial(m-1)-z**(-m)) == 0, 'torsion head to free inverse tower')
    return dict(actual_source_maps=len(maps), volterra_checks=checked,
                exact_return_invoices=return_checks, laurent_weyl_checks=17,
                torsion_to_free_clock_checks=12,
                kummer_inverse_clock_checks=20,
                kummer_last_coefficient=[inverse_coeff.numerator, inverse_coeff.denominator])


if __name__ == '__main__':
    result = dict(atomic=atomic_controls(), hazard=hazard_controls(),
                  palettes=palette_controls(), jacobian=jacobian_controls())
    data = json.dumps(result, indent=2, sort_keys=True)+'\n'
    path = ROOT/'05-knowledge/results/marked_debt_controls_20261015f.json'
    path.write_text(data)
    print(data)
    print(hashlib.sha256(data.encode()).hexdigest())
