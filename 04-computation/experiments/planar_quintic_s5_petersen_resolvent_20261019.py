"""Exact coefficient and finite-graph certificates for packet 133.

The generic Galois-group proof is analytic/algebraic in the companion note;
this small checker verifies only its explicit formulas and graph labeling.
"""
import itertools
import json
import sympy as sp

checks = []


def check(name, condition):
    if not bool(condition):
        raise RuntimeError(name)
    checks.append(name)


v, y, p, q, e, z, lam, s, w = sp.symbols("v y p q E z lambda s w")
f = v**5 + p*v**2 + q
r = y**10 + p*y**7 - 11*q*y**5 - p**2*y**4 + 7*p*q*y**2 - p**3*y - q**2
diagonal = y**5 + 8*p*y**2 + 32*q
check("pair_sum_resultant", sp.expand(sp.resultant(f, f.subs(v, y-v), v) - diagonal*r**2) == 0)
check("diagonal_factor", sp.expand(32*f.subs(v, y/2) - diagonal) == 0)
check("monic_degree_ten", sp.Poly(r, y).degree() == 10 and sp.Poly(r, y).LC() == 1)
check("quintic_discriminant", sp.expand(sp.discriminant(f, v) - q*(108*p**5 + 3125*q**3)) == 0)

a = 3*e - 5*z**2 - 15*lam*z**5
quintic = a*v**5 + 5*z**5*v**2 - 3*e*z**5
kappa = 1 + 3*lam*z**3
check("zero_parameter_factorization", sp.expand(quintic.subs(e, 0) - 5*z**2*v**2*(z**3-kappa*v**3)) == 0)
lift = (3*s**2 - 5*z**2*kappa)*s**3*w**5 + 5*z**5*w**2 - 3*z**5
check("ramified_local_equation", sp.expand(quintic.subs({e: s**2, v: s*w}) - s**2*lift) == 0)
check("simple_lift_initial_quadratic", sp.expand(lift.subs(s, 0) - z**5*(5*w**2-3)) == 0)

vertices = tuple(frozenset(pair) for pair in itertools.combinations(range(5), 2))
edges = {frozenset((a, b)) for a, b in itertools.combinations(vertices, 2) if not a & b}
check("ten_vertices_fifteen_edges", len(vertices) == 10 and len(edges) == 15)
check("three_regular", all(sum(vtx in edge for edge in edges) == 3 for vtx in vertices))
stars = {frozenset(vtx for vtx in vertices if i in vtx) for i in range(5)}
independent_four = {frozenset(group) for group in itertools.combinations(vertices, 4)
                    if all(a & b for a, b in itertools.combinations(group, 2))}
check("five_maximum_stars", independent_four == stars and len(stars) == 5)
check("no_independent_five", not any(all(a & b for a, b in itertools.combinations(group, 2))
                                    for group in itertools.combinations(vertices, 5)))
outer = [frozenset((2*i % 5, (2*i+1) % 5)) for i in range(5)]
inner = [frozenset(((2*i+2) % 5, (2*i+4) % 5)) for i in range(5)]
check("standard_petersen_labels", set(outer+inner) == set(vertices))
standard_edges = {frozenset((outer[i], outer[(i+1) % 5])) for i in range(5)}
standard_edges |= {frozenset((outer[i], inner[i])) for i in range(5)}
standard_edges |= {frozenset((inner[i], inner[(i+2) % 5])) for i in range(5)}
check("standard_petersen_edges", standard_edges == edges)

h, phi, xi = sp.symbols("h phi xi")
p_h = 5*z**5*h/(3-(5*z**2+15*lam*z**5)*h)
q_h = -3*z**5/(3-(5*z**2+15*lam*z**5)*h)
check("normalized_limit_coefficients", p_h.subs(h, 0) == 0 and q_h.subs(h, 0) == -z**5)
r_limit = y**10 + 11*z**5*y**5 - z**10
check("pair_resolvent_at_infinity", sp.expand(r.subs({p: 0, q: -z**5}) - r_limit) == 0)
golden_poly = phi**2-phi-1
check("fifth_lucas_value", sp.rem(phi**5-(phi-1)**5-11, golden_poly, phi) == 0)
check("golden_factorization", sp.rem(sp.expand((y**5+phi**5*z**5)*(y**5-(phi-1)**5*z**5)-r_limit), golden_poly, phi) == 0)
cyclotomic = xi**4+xi**3+xi**2+xi+1
phi_xi = 1+xi+xi**4
check("side_fifth_power", sp.rem((1+xi)**5+phi_xi**5, cyclotomic, xi) == 0)
check("diagonal_fifth_power", sp.rem((1+xi**2)**5-(phi_xi-1)**5, cyclotomic, xi) == 0)
side_orbit = {frozenset((i, (i+1) % 5)) for i in range(5)}
diagonal_orbit = {frozenset((i, (i+2) % 5)) for i in range(5)}
check("cyclic_pair_orbits", len(side_orbit) == len(diagonal_orbit) == 5
      and side_orbit.isdisjoint(diagonal_orbit) and side_orbit | diagonal_orbit == set(vertices))
check("two_induced_five_cycles", all(sum(vtx in edge and edge <= orbit for edge in edges) == 2
      for orbit in (side_orbit, diagonal_orbit) for vtx in orbit))
cross_edges = {edge for edge in edges if len(edge & side_orbit) == 1}
check("five_cross_matching_edges", len(cross_edges) == 5
      and all(sum(vtx in edge for edge in cross_edges) == 1 for vtx in vertices))

print(json.dumps({"packet": 133, "status": "PASS", "checks": len(checks),
                  "verified": checks, "resolvent": str(r),
                  "discriminant": "q*(108*p**5+3125*q**3)",
                  "scope": "Explicit symbolic identities and finite graph only; universal Galois proof is in the note."},
                 ensure_ascii=False, indent=2, sort_keys=True))
