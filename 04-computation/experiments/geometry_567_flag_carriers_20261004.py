"""Exact rotation-memory controls: Heawood maps and (2,3,5)/(2,3,7) quotients.

Standard library, exhaustive declared finite universes, checks survive -O.
"""
from collections import Counter
from fractions import Fraction
from itertools import permutations, product
import json


def need(ok, message):
    if not ok:
        raise ArithmeticError(message)


def compose(p, q):
    return tuple(p[q[i]] for i in range(len(p)))


def inverse(p):
    result = [None]*len(p)
    for i, j in enumerate(p):
        result[j] = i
    return tuple(result)


def cycles(p):
    unseen, result = set(range(len(p))), []
    while unseen:
        start = min(unseen)
        block, current = [], start
        while current in unseen:
            unseen.remove(current)
            block.append(current)
            current = p[current]
        need(current == start, 'permutation cycles close')
        result.append(tuple(block))
    return tuple(result)


D = (1, 2, 4)
LINES = tuple(frozenset((t+d) % 7 for d in D) for t in range(7))
NEIGHBORS = tuple(tuple(7+(p-d) % 7 for d in D) for p in range(7)) + tuple(
    tuple((t+d) % 7 for d in D) for t in range(7))
DARTS = tuple((v, w) for v in range(14) for w in NEIGHBORS[v])
DART_INDEX = {dart: i for i, dart in enumerate(DARTS)}
ALPHA = tuple(DART_INDEX[w, v] for v, w in DARTS)


def rotation(signs):
    need(len(signs) == 14 and set(signs) <= {-1, 1}, 'fourteen local orientations')
    return tuple(DART_INDEX[v, NEIGHBORS[v][(NEIGHBORS[v].index(w)+signs[v]) % 3]]
                 for v, w in DARTS)


def uniform_rotation(point_sign, line_sign):
    return rotation((point_sign,)*7+(line_sign,)*7)


def graph_automorphisms():
    # A connected bipartite graph automorphism globally preserves or swaps parts.
    line_lookup = {line: t for t, line in enumerate(LINES)}
    collineations = []
    for pi in permutations(range(7)):
        images = [frozenset(pi[p] for p in line) for line in LINES]
        if all(image in line_lookup for image in images):
            collineations.append(tuple(pi)+tuple(7+line_lookup[image] for image in images))
    duality = tuple(7+(-p) % 7 for p in range(7))+tuple((-t) % 7 for t in range(7))
    return tuple(collineations)+tuple(compose(duality, pi) for pi in collineations)


def dart_action(pi):
    return tuple(DART_INDEX[pi[v], pi[w]] for v, w in DARTS)


def map_profile(sigma, actions):
    faces = cycles(compose(sigma, ALPHA))
    face_of = {d: f for f, face in enumerate(faces) for d in face}
    edges = cycles(ALPHA)
    dual = Counter(tuple(sorted((face_of[d], face_of[ALPHA[d]]))) for d, _ in edges)
    positive, negative = [], []
    back = inverse(sigma)
    for action in actions:
        transported = compose(action, sigma)
        if transported == compose(sigma, action):
            positive.append(action)
        if transported == compose(back, action):
            negative.append(action)
    chi = 14-21+len(faces)
    need(chi % 2 == 0, 'orientable Euler characteristic')
    return dict(faces=faces, face_lengths=dict(sorted(Counter(map(len, faces)).items())),
                euler=chi, genus=1-chi//2, positive=positive, negative=negative,
                dual_edges=dual)


def orbit_sizes(objects, actions, action_fn):
    unseen, sizes = set(objects), []
    while unseen:
        obj = min(unseen)
        orbit = {action_fn(g, obj) for g in actions}
        need(orbit <= unseen, 'complete group orbits partition carrier')
        unseen -= orbit
        sizes.append(len(orbit))
    return sorted(sizes)


def matrix_product(a, b, p):
    return ((a[0]*b[0]+a[1]*b[2]) % p, (a[0]*b[1]+a[1]*b[3]) % p,
            (a[2]*b[0]+a[3]*b[2]) % p, (a[2]*b[1]+a[3]*b[3]) % p)


def psl_triangle(p):
    canon = lambda m: min(m, tuple(-x % p for x in m))
    group = sorted({canon(m) for m in product(range(p), repeat=4)
                    if (m[0]*m[3]-m[1]*m[2]) % p == 1})
    index = {g: i for i, g in enumerate(group)}
    s, t = (0, p-1, 1, 0), (0, p-1, 1, 1)
    alpha = tuple(index[canon(matrix_product(g, s, p))] for g in group)
    sigma = tuple(index[canon(matrix_product(g, t, p))] for g in group)
    face = compose(sigma, alpha)
    vc, ec, fc = cycles(sigma), cycles(alpha), cycles(face)
    n = len(group)
    need({len(c) for c in vc} == {3}, 'uniform ternary vertex rotations')
    need({len(c) for c in ec} == {2}, 'fixed-point-free edge involution')
    need({len(c) for c in fc} == {p}, 'uniform face order')
    reached, stack = {0}, [0]
    while stack:
        i = stack.pop()
        for j in (alpha[i], sigma[i]):
            if j not in reached:
                reached.add(j)
                stack.append(j)
    need(len(reached) == n, 'connected regular dart carrier')
    # Check the graph rather than assuming the group quotient is a simple map.
    vertex_of = {d: v for v, c in enumerate(vc) for d in c}
    undirected = [tuple(sorted((vertex_of[c[0]], vertex_of[c[1]]))) for c in ec]
    need(all(v != w for v, w in undirected) and len(set(undirected)) == len(ec), 'simple quotient graph')
    need(all(len({vertex_of[d] for d in c}) == p for c in fc), 'simple face boundaries')
    chi = len(vc)-len(ec)+len(fc)
    excess = Fraction(1, 2)+Fraction(1, 3)+Fraction(1, p)-1
    need(chi == n*excess, 'triangle excess equals Euler characteristic per dart')
    return dict(p=p, darts=n, flags=2*n, vertices=len(vc), edges=len(ec), faces=len(fc),
                genus=1-chi//2, excess=str(excess))


def main():
    need(len(DARTS) == 42 and len(set(DARTS)) == 42, 'Heawood darts')
    need(all(len(LINES[a] & LINES[b]) == 1 for a in range(7) for b in range(a)), 'Fano line intersection')
    aut = graph_automorphisms()
    need(len(aut) == len(set(aut)) == 336, 'all Heawood graph automorphisms')
    actions = tuple(dart_action(pi) for pi in aut)
    need(all(compose(g, ALPHA) == compose(ALPHA, g) for g in actions), 'edge reversal is intrinsic')
    paley = tuple(tuple((a*p+b) % 7 for p in range(7)) + tuple(7+(a*t+b) % 7 for t in range(7))
                  for a in D for b in range(7))
    paley_actions = tuple(dart_action(g) for g in paley)
    need(set(paley) <= set(aut), 'retained Paley group acts on incidences')
    need(all(((g[v]-g[u]) % 7 in D) == ((v-u) % 7 in D)
             for g in paley for u in range(7) for v in range(7) if u != v), 'Paley tournament orientation retained')
    canonical = uniform_rotation(-1, 1)
    for i, (_, neighbor) in enumerate(DARTS):
        next_neighbor = DARTS[canonical[i]][1]
        need((next_neighbor-neighbor) % 7 in D,
             'canonical star rotation follows each induced Paley three-cycle')
    reports = []
    for point_sign, line_sign in product((1, -1), repeat=2):
        sigma = uniform_rotation(point_sign, line_sign)
        profile = map_profile(sigma, actions)
        same = point_sign == line_sign
        need(profile['face_lengths'] == ({6: 7} if same else {14: 3}), 'relative rotation bit controls face lengths')
        need((len(profile['positive']), len(profile['negative'])) == ((42, 0) if same else (21, 21)), 'oriented versus reversing map symmetries')
        need(set(paley_actions) <= set(profile['positive']), 'same Paley21 group preserves every map')
        need(all(len({DARTS[d][0] for d in f}) == len(f) for f in profile['faces']), 'each face is a simple cycle')
        full = profile['positive']+profile['negative']
        face_sets = {frozenset(f) for f in profile['faces']}
        def face_action(g, f):
            # An orientation-reversing map sends the left face to the right face.
            if g in profile['positive']:
                return frozenset(g[d] for d in f)
            return frozenset(ALPHA[g[d]] for d in f)
        vertex_sets = {frozenset(DART_INDEX[v,w] for w in NEIGHBORS[v]) for v in range(14)}
        edge_sets = {frozenset(c) for c in cycles(ALPHA)}
        orbits = dict(vertices=orbit_sizes(vertex_sets, full, lambda g,c:frozenset(g[d] for d in c)),
                      edges=orbit_sizes(edge_sets, full, lambda g,c:frozenset(g[d] for d in c)),
                      faces=orbit_sizes(face_sets, full, face_action))
        need(orbits == dict(vertices=[14], edges=[21], faces=[7 if same else 3]), 'both maps are vertex/edge/face transitive')
        need(not any(a == b for a,b in profile['dual_edges']), 'dual has no loops')
        need(sorted(profile['dual_edges'].values()) == ([1]*21 if same else [7]*3), 'dual is K7 or seven-parallel K3')
        # 84 flags; map automorphisms act freely, but have only order42.
        need(84//len(full) == 2, 'cell-transitivity does not imply flag-transitivity')
        reports.append(dict(signs=[point_sign,line_sign], face_lengths=profile['face_lengths'],
                            genus=profile['genus'], orientation_preserving=len(profile['positive']),
                            orientation_reversing=len(profile['negative']), cell_orbits=orbits,
                            flag_orbits=2, dual='K7' if same else '7 parallel edges per pair of K3',
                            face_vertex_cycles=[[DARTS[d][0] for d in f] for f in profile['faces']]))
    census = Counter()
    invariant = []
    for mask in range(1 << 14):
        sigma = rotation(tuple(1 if (mask>>i)&1 else -1 for i in range(14)))
        face_count = len(cycles(compose(sigma, ALPHA)))
        genus = (9-face_count)//2
        need(9-face_count == 2*genus and 1 <= genus <= 4, 'all orientable Heawood rotation genera')
        census[genus] += 1
        # Two group generators suffice for the known affine Paley group.
        if all(compose(g,sigma) == compose(sigma,g) for g in (paley_actions[1],paley_actions[7])):
            invariant.append(mask)
    need(len(invariant) == 4, 'exactly four Paley-equivariant local rotations')
    quotient_controls = [psl_triangle(5), psl_triangle(7)]
    need(quotient_controls[0]['genus'] == 0 and quotient_controls[1]['genus'] == 3, 'spherical/hyperbolic quotient controls')
    excesses = {q:str(Fraction(1,2)+Fraction(1,3)+Fraction(1,q)-1) for q in (5,6,7)}
    report = dict(status='PROVED rotation-bit mechanism; FINITE-EXACT full Heawood rotation census',
                  graph=dict(vertices=14, edges=21, darts=42, flags=84, abstract_automorphisms=len(aut), paley_subgroup=21),
                  all_rotation_systems=1<<14, genus_counts=dict(sorted(census.items())),
                  paley_invariant_rotation_masks=invariant, maps=reports,
                  tournament_induced_rotation=dict(signs=[-1,1], local_arc_checks=42,
                      genus=3, rule='At both vertex types follow the induced Paley cycle on neighbor labels.'),
                  triangle_excesses=excesses, regular_group_quotients=quotient_controls,
                  loss='The graph, its full flow space, and the Paley subgroup do not determine faces, genus or curvature.',
                  restoration='Retain the local cyclic rotation sigma alongside the intrinsic dart reversal alpha.',
                  boundary='The genus3 Heawood map is type(14,3), not the type(7,3) PSL(2,7) map.')
    print(json.dumps(report, indent=2))
    print('PASS: all checks active under -O; no Collatz transfer asserted')


if __name__ == '__main__':
    main()
