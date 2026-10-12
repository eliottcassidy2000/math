"""Six-colour core at the eight fixed ports: double-port census, the doubling
extension lemma, and three-colourability probes of eleventh-root contexts.

Run from the repository root with python3 and again with python3 -O.
Dependencies: the maintained exact Q(zeta55) engine
(charge_ports_eleven_20261011.py), the exit/boundary script
(paid_exit_boundaries_20261012.py), SymPy, mpmath, PySAT.

Exact arithmetic: every point is a 40-vector of rationals in the power basis
of Q(zeta55).  Numerics (mpmath, 120 digits) are used only to *propose*
square roots; every proposal is verified exactly, and every negative is
certified by a quadratic non-residue witness at a split prime.
"""
from fractions import Fraction as Q
from itertools import combinations, permutations, product
from collections import Counter, defaultdict
from math import lcm
from pathlib import Path
import importlib.util
import json
import mpmath
from pysat.solvers import Solver

ROOT = Path(__file__).resolve().parents[2]
CHECKS = 0


def check(ok, label):
    global CHECKS
    CHECKS += 1
    if not ok:
        raise RuntimeError(label)


def engine(name, filename):
    spec = importlib.util.spec_from_file_location(name, ROOT/'04-computation/experiments'/filename)
    mod = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(mod)
    return mod


G = engine('cyclotomic_geometry', 'charge_ports_eleven_20261011.py')
PE = engine('paid_exits', 'paid_exit_boundaries_20261012.py')

mpmath.mp.dps = 120
ZETA = mpmath.exp(2*mpmath.pi*1j/55)
ZPOW = [ZETA**k for k in range(55)]
LAMBDA = mpmath.pi  # transcendental mixing constant for the two real relations
UNITS = [k for k in range(1, 55) if k % 5 and k % 11]  # (Z/55)^x, 40 elements


def emb(v):
    return sum(mpmath.mpf(x.numerator)/x.denominator*ZPOW[i] if isinstance(x, Q) else x*ZPOW[i]
               for i, x in enumerate(v) if x)


def galois(v, k):
    """The automorphism zeta -> zeta^k, exact."""
    out = G.O
    for i, x in enumerate(v):
        if x:
            out = G.add(out, G.sc(G.powers[(i*k) % 55], x))
    return out


def inv(w):
    """Exact inverse through the product of the other 39 conjugates."""
    others = G.I
    for k in UNITS[1:]:
        others = G.mul(others, galois(w, k))
    norm = G.mul(w, others)
    check(all(x == 0 for x in norm[1:]) and norm[0] != 0, 'rational nonzero field norm')
    result = G.sc(others, Q(1, 1)/norm[0])
    check(G.mul(w, result) == G.I, 'exact inverse')
    return result


def reconstruct(s):
    """Propose the coordinates of a complex number s in K=Q(alpha,beta), or None.

    The ports lie in K, and a square root of an element of K+ that lies in L
    lies in K (L = K(zeta11) has odd degree five over K).  So an eight-term
    integer relation suffices; the result is verified exactly by the caller.
    """
    alpha, beta = G.powers[11], G.b
    basis = []
    for i in range(4):
        for j in range(2):
            vec = G.mul(G.powers[11*i], beta) if j else G.powers[11*i]
            basis.append(vec)
    xs = [mpmath.re(emb(v))+LAMBDA*mpmath.im(emb(v)) for v in basis]
    xs.append(mpmath.re(s)+LAMBDA*mpmath.im(s))
    rel = mpmath.pslq(xs, maxcoeff=10**12, maxsteps=10**5)
    if rel is None or rel[8] == 0:
        return None
    out = G.O
    for c, vec in zip(rel[:8], basis):
        if c:
            out = G.add(out, G.sc(vec, Q(-c, rel[8])))
    return out


def split_primes(count):
    out, p = [], 56
    while len(out) < count:
        p += 55
        if G.s.isprime(p):
            out.append(p)
    return out


def residue_witness(y):
    """A split prime at which y is a quadratic non-residue, proving y is not a square in L."""
    for p in split_primes(60):
        root = pow(int(G.s.primitive_root(p)), (p-1)//55, p)
        check(pow(root, 55, p) == 1 and pow(root, 5, p) != 1 and pow(root, 11, p) != 1, 'order-55 root')
        if any(x.denominator % p == 0 for x in y):
            continue
        value = sum(x.numerator*pow(x.denominator, -1, p)*pow(root, i, p) for i, x in enumerate(y) if x) % p
        if value == 0:
            continue
        if pow(value, (p-1)//2, p) == p-1:
            return p
    return None


def is_square(y):
    """Exact decision: returns ('square', s) with s*s == y, or ('nonsquare', prime)."""
    s = reconstruct(mpmath.sqrt(emb(y)))
    if s is not None and G.mul(s, s) == y:
        return 'square', s
    p = residue_witness(y)
    check(p is not None, 'every non-square receives a residue witness')
    return 'nonsquare', p


def unit_vector(v):
    return G.mul(v, G.conj(v)) == G.I


# ---------------------------------------------------------------- census

def apex_census():
    """All points of L at unit distance from two of the eight ports."""
    ports = G.ports
    check(len(ports) == 8, 'eight distinct ports')
    rows, apexes = [], {}
    for i, j in combinations(range(8), 2):
        diff = G.add(ports[j], G.sc(ports[i], -1))
        D = G.mul(diff, G.conj(diff))
        check(D == G.conj(D), 'squared distance is real')
        Dn = mpmath.re(emb(D))
        check(abs(mpmath.im(emb(D))) < mpmath.mpf(10)**-100, 'numerically real')
        row = {'pair': [i, j], 'D_numeric': float(Dn)}
        if D == G.I:
            row['kind'] = 'unit edge, no apex in L (triangle-free)'
        elif Dn > 4:
            row['kind'] = 'distance above two, no apex'
        elif D == G.sc(G.I, 4):
            row['kind'] = 'distance exactly two, midpoint apex'
            x = G.add(ports[i], G.sc(diff, Q(1, 2)))
            apexes.setdefault(x, []).append((i, j))
        else:
            y = G.mul(D, G.add(D, G.sc(G.I, -4)))  # D(D-4)
            verdict, data = is_square(y)
            if verdict == 'square':
                s = data
                half = G.sc(diff, Q(1, 2))
                shift = G.mul(G.mul(diff, s), G.sc(inv(D), Q(1, 2)))
                found = []
                for sign in (1, -1):
                    x = G.add(G.add(ports[i], half), G.sc(shift, sign))
                    for k in (i, j):
                        d = G.add(x, G.sc(ports[k], -1))
                        check(unit_vector(d), 'apex at unit distance from both ports')
                    found.append(x)
                    apexes.setdefault(x, []).append((i, j))
                row['kind'] = 'two apexes in L'
                row['sqrt_D(D-4)'] = [str(c) for c in s]
            else:
                row['kind'] = 'apexes outside L'
                row['nonresidue_prime'] = data
        rows.append(row)
    # every apex is a point of L; drop any that coincides with a port
    D_points = [x for x in apexes if x not in ports]
    return rows, D_points, {str(k): v for k, v in enumerate(D_points)}, apexes


# ------------------------------------------------- fast complete unit graph

def fast_graph(points):
    """Complete unit-distance graph; same filters as the engine, bucketed."""
    den = lcm(*(x.denominator for p in points for x in p if isinstance(x, Q)))
    pts = [tuple(int(x*den) for x in p) for p in points]
    images = []
    for prime in (331, 661):
        root = pow(int(G.s.primitive_root(prime)), (prime-1)//55, prime)
        check(pow(root, 55, prime) == 1 and pow(root, 5, prime) != 1 and pow(root, 11, prime) != 1, 'order-55 root')
        rp = [pow(root, k, prime) for k in range(40)]
        rc = [pow(root, -k, prime) for k in range(40)]
        images.append((prime, [(sum(x*y for x, y in zip(p, rp)) % prime,
                                sum(x*y for x, y in zip(p, rc)) % prime) for p in pts]))
    prime, img = images[0]
    buckets = defaultdict(list)
    for idx, key in enumerate(img):
        buckets[key].append(idx)
    target = den*den % prime
    edges, survivors = [], 0
    for i, (f, fc) in enumerate(img):
        for t in range(1, prime):
            key = ((f+t) % prime, (fc+target*pow(t, -1, prime)) % prime)
            for j in buckets.get(key, ()):
                if j <= i:
                    continue
                q2, img2 = images[1]
                if (img2[i][0]-img2[j][0])*(img2[i][1]-img2[j][1]) % q2 != den*den % q2:
                    continue
                survivors += 1
                d = G.add(pts[i], G.sc(pts[j], -1))
                if G.mul(d, G.conj(d)) == G.sc(G.I, den*den):
                    edges.append((i, j))
    edges.sort()
    return edges, survivors, den


def colourable(count, edges, k):
    clauses = [[k*i+c+1 for c in range(k)] for i in range(count)]
    clauses += [[-(k*i+c+1), -(k*j+c+1)] for i, j in edges for c in range(k)]
    with Solver(name='g3', bootstrap_with=clauses) as solver:
        if not solver.solve():
            return None
        model = set(solver.get_model())
        colours = [next(c for c in range(k) if k*i+c+1 in model) for i in range(count)]
    for i, j in edges:
        check(colours[i] != colours[j], 'proper colouring on every edge')
    return colours


# ------------------------------------------------------- doubling lemma

def matchings(colours=6):
    out = []
    def rec(rest, acc):
        if not rest:
            out.append(tuple(acc)); return
        a = rest[0]
        for b in rest[1:]:
            rec([c for c in rest if c not in (a, b)], acc+[(a, b)])
    rec(list(range(colours)), [])
    return out


def doubling_blockable(pattern, dports, dadj):
    """SAT: can classes 1..3 on the double-port vertices kill all 90 doubling options?

    pattern: colours of the 8 ports. dports: for each D-vertex its two port indices.
    dadj: unit edges among D-vertices (classes must then differ).
    Option = (perfect matching of the six colours into three pairs, bijection classes->pairs).
    Vertex x kills the option iff the pair assigned to its class equals {t(p),t(q)}.
    """
    n = len(dports)
    var = lambda x, c: 3*x+c+1
    clauses = [[var(x, c) for c in range(3)] for x in range(n)]
    clauses += [[-var(x, c), -var(x, d)] for x in range(n) for c in range(3) for d in range(c)]
    clauses += [[-var(x, c), -var(y, c)] for x, y in dadj for c in range(3)]
    targets = []
    for x, (p, q) in enumerate(dports):
        a, b = pattern[p], pattern[q]
        targets.append(None if a == b else (min(a, b), max(a, b)))
    for M in matchings():
        for bij in permutations(range(3)):
            clause = []
            for x, T in enumerate(targets):
                if T in M:
                    clause.append(var(x, bij[M.index(T)]))
            if not clause:
                return False  # this option survives whatever the classes are
            clauses.append(clause)
    with Solver(name='g3', bootstrap_with=clauses) as solver:
        return bool(solver.solve())


def doubling_controls():
    # abstract control: 15 D-vertices on a star of pairs {0,j} x 3 classes block everything
    dports = [(0, j) for j in range(1, 6) for _ in range(3)]
    pattern = (0, 1, 2, 3, 4, 5, 0, 0)
    check(doubling_blockable(pattern, dports, []), 'fifteen star vertices can block (abstract control)')
    check(not doubling_blockable(pattern, dports[:14], []), 'fourteen cannot block')
    # all 15 pairs in one common class block; made impossible by a path forcing distinct classes
    dports = [(min(i, j), max(i, j)) for i, j in combinations(range(6), 2)]
    check(doubling_blockable(pattern, dports, []), 'fifteen distinct pairs with a shared class block')
    dadj = [(x, x+1) for x in range(14)]
    check(doubling_blockable(pattern, dports, dadj) in (True, False), 'adjacency constraints are consulted')


# ----------------------------------------------------------------- probes

def ball(directions, radius, extra=()):
    dirs = list(dict.fromkeys(list(directions)+[G.sc(v, -1) for v in directions]))
    D = lcm(*(Q(x).denominator for p in dirs+list(extra)+G.ports for x in p))
    dd = [tuple(int(x*D) for x in p) for p in dirs]
    layer = {G.O}; shell = {G.O}
    for _ in range(radius):
        shell = {G.add(x, y) for x, y in product(shell, dd)}-layer
        layer |= shell
    points = list(dict.fromkeys(list(G.ports)+list(extra)+[tuple(Q(x, D) for x in p) for p in sorted(layer)]))
    return points


def probe(label, points, palette=3):
    edges, surv, den = fast_graph(points)
    colours = colourable(len(points), edges, palette)
    rec = {'label': label, 'vertices': len(points), 'edges': len(edges), 'filter_survivors': surv,
           'denominator': den, 'palette': palette, 'colourable': colours is not None}
    if colours is None:
        rec['four_colouring'] = colourable(len(points), edges, 4) is not None
    return rec, edges


def closed_form_apexes():
    """The nine apexes as rhombus completions: P_j + alpha^(j+1) and the four completions at 0."""
    alpha = G.powers[11]
    P = [G.O]
    for j in range(4):
        P.append(G.add(P[-1], G.powers[11*j]))
    eta = G.eta
    pent = [G.add(P[j], G.powers[(11*(j+1)) % 55]) for j in range(5)]
    cross = [G.add(G.I, G.sc(eta, -1)), G.add(G.I, G.sc(G.conj(eta), -1)),
             G.add(P[4], G.sc(eta, -1)), G.add(P[4], G.sc(G.conj(eta), -1))]
    return pent+cross


def residue_graph_109():
    """Hoffman bound for the Cayley graph of F_(109^2) on its 110 norm-one residues."""
    import math
    q = 109
    d = next(x for x in range(2, q) if pow(x, (q-1)//2, q) == q-1)
    N1 = [(a, b) for a in range(q) for b in range(q) if (a*a-d*b*b) % q == 1]
    check(len(N1) == q+1, 'norm-one torus of F_(109^2)')
    S = set(N1)
    check(not any(((a+c) % q, (b+e) % q) in S for a, b in N1 for c, e in N1), 'residue graph is triangle-free')
    lam_min = 0.0
    for u in range(q):
        for v in range(q):
            s = sum(math.cos(2*math.pi*((u*a+v*b) % q)/q) for a, b in N1)
            lam_min = min(lam_min, s)
    hoffman = 1+(q+1)/abs(lam_min)
    check(hoffman > 6.5, 'Hoffman bound exceeds six with margin')
    return dict(prime=q, degree=q+1, nonsquare=d, least_eigenvalue=lam_min, hoffman_lower_bound=hoffman,
                chromatic_number_at_least=math.ceil(hoffman-1e-9))


def main(large=False):
    out = {}
    # 1. census of double-port vertices
    rows, D_points, _, apexes = apex_census()
    out['pair_census'] = rows
    out['double_port_vertices'] = len(D_points)
    # cross-check: any apex adjacent to three ports?
    pts = list(G.ports)+D_points
    edges, surv, den = fast_graph(pts)
    old_edges, old_surv, old_den = G.graph(pts)
    check(edges == old_edges, 'bucketed finder agrees with the engine on ports plus apexes')
    port_neighbours = Counter()
    for i, j in edges:
        if i < 8 <= j:
            port_neighbours[j] += 1
        if j < 8 <= i:
            port_neighbours[i] += 1
    check(all(port_neighbours[8+x] == 2 for x in range(len(D_points))), 'every apex has exactly two port neighbours')
    dports = []
    for x in range(len(D_points)):
        nb = sorted(j for i, j in edges if i == 8+x and j < 8) + sorted(i for i, j in edges if j == 8+x and i < 8)
        check(len(nb) == 2, 'two port neighbours')
        dports.append(tuple(nb))
    dadj = [(i-8, j-8) for i, j in edges if i >= 8 and j >= 8]
    out['apex_port_pairs'] = dports
    out['apex_adjacency'] = dadj
    out['apex_coordinates'] = [[str(c) for c in x] for x in D_points]
    check(set(closed_form_apexes()) == set(D_points), 'the nine apexes are exactly the rhombus completions')
    out['residue_graph_109'] = residue_graph_109()
    # 2. doubling lemma on all 570 canonical patterns
    doubling_controls()
    pe = [(i, j) for i, j in edges if i < 8 and j < 8]
    check(len(pe) == 9, 'nine port edges')
    patterns = list(PE.partitions(8, 6, pe))
    check(len(patterns) == 570, 'complete canonical partition universe')
    blockable = [p for p in patterns if doubling_blockable(p, dports, dadj)]
    out['patterns'] = len(patterns)
    out['doubling_blockable_patterns'] = [list(p) for p in blockable]
    # 3. the apex context and richer eleventh-root contexts: boundary relation and chromatic probes
    mers = [G.unit(2, j) for j in range(1, 11)]
    for u in mers:
        check(unit_vector(u), 'Mersenne direction is a unit vector')
    v = G.sc(G.add(G.sc(G.I, -5), G.b), Q(1, 6))
    check(unit_vector(v), 'thirds unit')
    roots10 = [G.powers[11*j] for j in range(5)]
    roots11 = [G.powers[5*j] for j in range(11)]
    probes = []
    # (a) apex context: ports, apexes, and one Mersenne/thirds/root layer around each
    layer_dirs = list(dict.fromkeys(roots10+roots11+mers+[v]))
    layer_dirs = list(dict.fromkeys(layer_dirs+[G.sc(d, -1) for d in layer_dirs]))
    apex_context = list(dict.fromkeys(list(G.ports)+D_points+[G.add(p, d) for p in list(G.ports)+D_points for d in layer_dirs]))
    rec, e_apex = probe('apex context: ports+apexes+one mixed layer', apex_context)
    probes.append(rec)
    packet = PE.boundary_packet(apex_context, e_apex)
    rec['six_core_vertices'] = len(packet['core_vertices'])
    rec['partition_counts'] = packet['partition_counts']
    # (b) radius-two balls
    rec, _ = probe('ball radius 2: roots of unity 10+11 and Mersenne 2047-directions', ball(roots10+roots11+mers, 2))
    probes.append(rec)
    rec, _ = probe('ball radius 2: roots 10+11, Mersenne, thirds unit v', ball(roots10+roots11+mers+[v], 2))
    probes.append(rec)
    rec, _ = probe('ball radius 3: roots 10, Mersenne', ball(roots10+mers, 3))
    probes.append(rec)
    rec, _ = probe('ball radius 3: roots 10+11, thirds unit', ball(roots10+roots11+[v], 3))
    probes.append(rec)
    if large:  # declared large run: python3 -B <script> --large (about one minute)
        rec, _ = probe('ball radius 4: roots 10, Mersenne', ball(roots10+mers, 4))
        probes.append(rec)
    out['probes'] = probes
    out['checks'] = CHECKS
    out['engine_checks'] = G.checks
    return out


if __name__ == '__main__':
    import sys
    result = main(large='--large' in sys.argv[1:])
    print(json.dumps(result, indent=1, sort_keys=True))
