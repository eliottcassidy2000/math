#!/usr/bin/env python3
"""Focused replay of incoming THM-4662 after repairing its K0 square test.

No point search and no database: replay the local descent, then use the two
displayed points and a separate F2 calculation.  This is a focused audit of
the incoming mechanism, not an independent implementation of local descent.
"""
import importlib.util
import itertools
import json
from fractions import Fraction as Q
from pathlib import Path

spec = importlib.util.spec_from_file_location(
    "descent", Path(__file__).with_name("two_descent_q_sqrt5_fish_curve_20261016.py"))
d = importlib.util.module_from_spec(spec)
spec.loader.exec_module(d)
K = d.K
checks = 0


def require(ok, label):
    global checks
    checks += 1
    if not ok:
        raise RuntimeError(label)


def exact_square(x):
    r = d.is_square_K(x)
    if r is not None:
        require((r * r - x).iszero(), "returned square root")
    return r is not None


def xor_span(gens):
    ans = {(0, 0)}
    for a, b in gens:
        ans |= {(x ^ a, y ^ b) for x, y in list(ans)}
    return ans


def hub_audit():
    # Independent exact arithmetic in Q(sqrt7,i), basis 1,sqrt7,i,i sqrt7.
    def add(x, y):
        return tuple(a + b for a, b in zip(x, y))

    def neg(x):
        return tuple(-a for a in x)

    def mul(x, y):
        out = [Q(0)] * 4
        for a in range(4):
            for b in range(4):
                factor = (7 if a & b & 1 else 1) * (-1 if a & b & 2 else 1)
                out[a ^ b] += factor * x[a] * y[b]
        return tuple(out)

    def norm(x):
        return mul(x, (x[0], x[1], -x[2], -x[3]))

    z = (Q(0),) * 4
    one = (Q(1), Q(0), Q(0), Q(0))
    a = one
    b = (Q(0), Q(1, 4), Q(-3, 4), Q(0))
    c = (Q(0), Q(0), Q(-1), Q(0))
    d0 = (Q(3, 4), Q(0), Q(0), Q(1, 4))
    coords = [z, a, b, c, d0, add(a, b), add(a, c), add(b, d0), add(c, d0)]
    names = ['h', 'a', 'b', 'c', 'd', 'w_ab', 'w_ac', 'w_bd', 'w_cd']
    wanted = {(0, 1), (0, 2), (0, 3), (0, 4), (1, 5), (2, 5),
              (1, 6), (3, 6), (2, 7), (4, 7), (3, 8), (4, 8), (5, 8), (6, 7)}
    require(len(set(coords)) == 9, "injective hub coordinates")
    edges = {(i, j) for i in range(9) for j in range(i + 1, 9)
             if norm(add(coords[i], neg(coords[j]))) == one}
    require(wanted <= edges, "all fourteen required hub edges")
    require(not any(all(tuple(sorted(e)) in edges for e in itertools.combinations(t, 2))
                    for t in itertools.combinations(range(9), 3)), "no triangles")
    colours = next((c for c in itertools.product(range(3), repeat=9)
                    if c[0] == 0 and all(c[i] != c[j] for i, j in edges)), None)
    require(colours is not None, "three-colour positive control")
    two = any(all(c[i] != c[j] for i, j in edges)
              for c in itertools.product(range(2), repeat=9))
    require(not two, "nonbipartite hostile control")
    # Explicit point of C3 at D=1/2 over F4(sqrt7,sqrt11).
    require(11 * Q(1, 44) == Q(1, 2) * (1 - Q(1, 2)), "z0 squared")
    require(11 * Q(7, 44) == Q(1, 2) * (4 - Q(1, 2)), "z1 squared")
    require(11 * Q(7, 44) == (1 - Q(1, 2)) * (Q(1, 2) + 3), "z2 squared")
    # Hensel witnesses for the bad local places. At 7, 3 sqrt5 is a simple root.
    require((Q(-11, 5).numerator * pow(Q(-11, 5).denominator, -1, 8)) % 8 == 1,
            "2-adic square ratio")
    require((-11 - 1) % 3 == 0, "3-adic simple square")
    require((-11 - 45) % 7 == 0 and 6 % 7 != 0, "7-adic simple square in F49")
    return {"basis": ["1", "sqrt7", "i", "i*sqrt7"], "names": names,
            "coordinates": [[str(x) for x in v] for v in coords],
            "edges": sorted(edges), "additional_edges": sorted(edges - wanted),
            "three_colouring": colours, "chromatic_number": 3}


def replay(roots, gens, places, seed):
    d.random.seed(seed)
    images, gclasses, encoded = [], [], []
    for p in places:
        want = 2 if p.name.startswith("real") else 16 if p.name == "2(inert)" else 4
        im = d.local_image(p, roots, want)
        require(len(im) == want, "complete local image " + p.name)
        require(d.close_group(im) == im, "closed local image " + p.name)
        images.append(im)
        gc = [p.cls(g) for g in gens]
        allm = []
        for mask in range(1 << len(gens)):
            bits = [0] * p.dim
            for i in range(len(gens)):
                if mask >> i & 1:
                    bits = [a ^ b for a, b in zip(bits, gc[i])]
            allm.append(tuple(bits))
        gclasses.append(allm)
        encoded.append({"place": p.name, "size": want, "image": sorted(im)})
    survivors = [(a, b) for a in range(1 << len(gens))
                 for b in range(1 << len(gens))
                 if all((gc[a], gc[b]) in im for gc, im in zip(gclasses, images))]
    return {"seed": seed, "local_images": encoded, "classes": survivors}


def main():
    # Minimal omitted-root controls, and every square from this declared box.
    for q in (5, 20, 45, Q(5, 4), Q(5, 9)):
        require(exact_square(K(q)), "pure sqrt5 root")
    require(not exact_square(K(2)), "nonsquare 2")
    require(not exact_square(K(-1)), "negative nonsquare")
    for a, b, den in itertools.product(range(-8, 9), range(-8, 9), range(1, 5)):
        x = K(Q(a, den), Q(b, den))
        require(exact_square(x * x), "square box")
    # Pell near-cancellation: the old float comparison could erase this sign.
    a, b = 1, 1
    for _ in range(200):
        a, b = 9 * a + 20 * b, 4 * a + 9 * b
        # norm -4, a,b positive: a-b sqrt5 < 0 exactly.
        require(a * a - 5 * b * b == -4, "Pell norm")
        require(d.RealPlace(1).cls(K(a, -b)) == (1,), "exact real sign")
        require(d.RealPlace(-1).cls(K(a, -b)) == (0,), "conjugate real sign")

    delta = K(Q(55, 2), Q(11, 2))
    gens = [K(-1), d.PHI, K(2), K(3), K(7), d.SQRT5, K(4, 1), K(4, -1)]
    roots = [K(0), K(9) * delta, K(-7) * delta]
    small = replay([K(0), K(9), K(-7)], gens[:5],
                   d.places_K0([], [3, 7], False, True), 314159)
    require(len(small["classes"]) == 4, "rank-zero Selmer bound")
    runs = [replay(roots, gens, d.places_K0([(11, 7), (11, 4)], [3, 7], True, True), seed)
            for seed in (20261016, 314159)]
    require(runs[0]["classes"] == runs[1]["classes"], "two seeds same Selmer classes")
    require(runs[0]["local_images"] == runs[1]["local_images"], "two seeds same local images")
    require(len(runs[0]["classes"]) == 16, "twist rank upper bound two")
    points = [(K(484), K(1573, -3751), (0, 192)),
              (K(Q(-65, 2), Q(19, 2)), K(160, 400), (35, 129))]
    for x, y, masks in points:
        require((y * y - x * (x - roots[1]) * (x - roots[2])).iszero(), "exact twist point")
        for value, mask in zip((x, x - roots[1]), masks):
            require(exact_square(value / d.kprod(gens, mask)), "point descent class")
    torsion = [(17, 227), (226, 0)]
    raw_torsion = [((roots[0] - roots[1]) * (roots[0] - roots[2]), roots[0] - roots[1]),
                   (roots[1], roots[1] * (roots[1] - roots[2]))]
    for pair, masks in zip(raw_torsion, torsion):
        for value, mask in zip(pair, masks):
            require(exact_square(value / d.kprod(gens, mask)), "torsion class")
    span = xor_span(torsion + [p[2] for p in points])
    require(span == set(runs[0]["classes"]), "points fill all twist Selmer classes")

    def square_f4(x):
        return exact_square(x) or exact_square(x / delta)

    f4gens = [(17, 1), (9, 9), (0, 192), (35, 129)]
    f4span = sorted(xor_span(f4gens))
    reduced = {(min(a, a ^ 226), min(b, b ^ 226)) for a, b in f4span}
    require(len(reduced) == 16, "F2 independence modulo delta")
    for i, p in enumerate(f4span):
        for q in f4span[:i]:
            equal = all(square_f4(d.kprod(gens, a ^ b)) for a, b in zip(p, q))
            require(not equal, "independent exact field square test")
    hits = []
    for p in f4span:
        equal = all(square_f4(K(33) / d.kprod(gens, a)) for a in p)
        by_mask = all((a ^ 200) in (0, 226) for a in p)
        require(equal == by_mask, "field and F2 target tests agree")
        if equal:
            hits.append(p)
    require(not hits, "covering class excluded")
    hub = hub_audit()
    print(json.dumps({"status": "FINITE-EXACT focused replay; classical descent inputs inherited",
                      "square_box": [-8, 8, 4], "pell_sign_tests": 400,
                      "E_K0": small, "twist_runs": runs,
                      "twist_point_span": sorted(span), "F4_image": f4span,
                      "target_hits": hits, "even_extension_hub": hub,
                      "checks": checks}, sort_keys=True, separators=(",", ":")))


if __name__ == "__main__":
    main()
