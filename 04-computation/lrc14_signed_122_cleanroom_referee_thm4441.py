#!/usr/bin/env python3
"""Independent exact referee for the proposed signed-(1,2,2) closure."""

from fractions import Fraction as F
from itertools import combinations, permutations, product
from math import gcd
from collections import defaultdict

R = F(3, 14)
TARGET = F(6, 77)
CHECKS = 0


def require(ok, msg):
    global CHECKS
    CHECKS += 1
    if not ok:
        raise AssertionError(msg)


U = {
    1: (2, 2, -1),
    2: (2, -2, 1),
    3: (1, 2, -2),
    4: (2, 1, -2),
}


def family(w):
    a, b, c = w
    fs = []
    if 2*a + 2*b == c: fs.append(1)
    if 2*a + c == 2*b: fs.append(2)
    if a + 2*b == 2*c: fs.append(3)
    if 2*a + b == 2*c: fs.append(4)
    return tuple(fs)


def all_122_relations(w):
    coeffs = set()
    for p in set(permutations((1, 2, 2))):
        for s in product((-1, 1), repeat=3):
            v = tuple(p[i]*s[i] for i in range(3))
            if v[0] < 0:
                v = tuple(-x for x in v)
            if sum(v[i]*w[i] for i in range(3)) == 0:
                coeffs.add(v)
    return coeffs


def eligible(w):
    return (w[0] < w[1] < w[2] and all(x > 0 and x % 3 for x in w)
            and gcd(gcd(w[0], w[1]), w[2]) == 1)


def rows_parameterized(H):
    rows = set()
    for a in range(1, H+1):
        if a % 3 == 0: continue
        for b in range(a+1, H+1):
            if b % 3 == 0: continue
            candidates = ((1, 2*(a+b)), (2, 2*(b-a)))
            for f, c in candidates:
                w = (a, b, c)
                if c <= H and eligible(w): rows.add((w, f))
    # F3/F4 are also generated directly to avoid assuming the first formulas.
    for c in range(3, H+1):
        if c % 3 == 0: continue
        for a in range(1, c):
            if a % 3 == 0: continue
            b3_num = 2*c-a
            if b3_num % 2 == 0:
                w = (a, b3_num//2, c)
                if eligible(w): rows.add((w, 3))
            b4 = 2*c-2*a
            w = (a, b4, c)
            if eligible(w): rows.add((w, 4))
    return sorted(rows)


def strict_bounds(w):
    # Largest integer n satisfying 14*n < 3*(other-speed sum).
    return tuple((3*(sum(w)-w[i])-1)//14 for i in range(3))


def term(w, C, i):
    jk = [j for j in range(3) if j != i]
    return min(F(3, 7*w[2]),
               F(3*(w[jk[0]]+w[jk[1]])-14*abs(C[i]),
                 14*w[jk[0]]*w[jk[1]]))


def ray_packet(w, f):
    u = U[f]
    bd = strict_bounds(w)
    K = min(bd[i]//abs(u[i]) for i in range(3))
    ks = tuple(k for k in range(1, K+1) if k % 3)
    carriers = {tuple(s*k*u[i] for i in range(3))
                for k in ks for s in (-1, 1)}
    E = tuple(sum(term(w, C, i) for C in carriers) for i in range(3))
    mass = sum(min(term(w, C, i) for i in range(3)) for C in carriers)
    return carriers, E, mass, ks


def raw_packet(w):
    bd = strict_bounds(w)
    live = set()
    # Deliberately brute-force two coordinates; no ray or modular solver used.
    for x in range(-bd[0], bd[0]+1):
        for y in range(-bd[1], bd[1]+1):
            n = -(w[0]*x+w[1]*y)
            if n % w[2]: continue
            z = n//w[2]
            C = (x, y, z)
            if abs(z) > bd[2] or any(q % 3 == 0 for q in C): continue
            require(all(14*abs(C[i]) < 3*(sum(w)-w[i]) for i in range(3)),
                    ("strict roof", w, C))
            live.add(C)
    E = tuple(sum(term(w, C, i) for C in live) for i in range(3))
    mass = sum(min(term(w, C, i) for i in range(3)) for C in live)
    return live, E, mass


def profile_lines(w, f):
    z = tuple(F(q, w[2]) for q in w)
    lines = []
    for i in range(3):
        j, k = [h for h in range(3) if h != i]
        lines.append((R*(z[j]+z[k])/(z[j]*z[k]), F(abs(U[f][i]), 1)/(z[j]*z[k])))
    B = min(R*(z[j]+z[k])/abs(U[f][i])
            for i in range(3) for j, k in [[h for h in range(3) if h != i]])
    return lines, B


def integrate_profiles(w, f):
    lines, B = profile_lines(w, f)
    pts = {F(0), B}
    for a, b in lines:
        x = (a-2*R)/b
        if 0 < x < B: pts.add(x)
    for i in range(3):
        for j in range(i):
            ai, bi = lines[i]; aj, bj = lines[j]
            if bi != bj:
                x = (ai-aj)/(bi-bj)
                if 0 < x < B: pts.add(x)
    pts = sorted(pts)
    I = [F(0), F(0), F(0)]
    P = F(0)
    def value(i, x):
        a, b = lines[i]
        return min(2*R, a-b*x)
    for lo, hi in zip(pts, pts[1:]):
        mid = (lo+hi)/2
        for i in range(3):
            I[i] += (hi-lo)*(value(i, lo)+value(i, hi))/2
        owner = min(range(3), key=lambda i: value(i, mid))
        P += (hi-lo)*(value(owner, lo)+value(owner, hi))/2
    return tuple(I), P, lines, B


def claimed_integrals(w, f):
    t = F(w[0], w[2]) if f != 4 else F(w[1], w[2])
    if f == 1:
        return (R*R*(F(7, 8)+t/4), None, R*R*(1-t+2*t*t))
    if f == 2:
        return (R*R*F(7+16*t, 8*(1+2*t)), None, None)
    if f == 3:
        return (R*R*(F(7, 8)+9*t/16), None, R*R*(1-t/2+t*t/2))
    return (R*R*(1-t/16), None, R*R*(1-t/2+t*t/2))


def envelope_sidecar(w, f, lines, B):
    t = F(w[0], w[2]) if f != 4 else F(w[1], w[2])
    if f == 1: pair, sw = (0, 2), R*(1-t)/2
    elif f == 2: pair, sw = (0, 1), R/2
    elif f == 3 and t <= F(1, 2): pair, sw = (0, 2), R*(1-t)/2
    elif f == 3: pair, sw = (1, 2), R*t/2
    else: pair, sw = (0, 2), R*t/2
    def v(i, x):
        a, b = lines[i]
        return min(2*R, a-b*x)
    require(0 < sw < B, ("switch range", w, f, sw, B))
    require(v(pair[0], sw) == v(pair[1], sw), ("switch equality", w, f))
    probes = {F(0), sw, B}
    for a, b in lines:
        x = (a-2*R)/b
        if 0 < x < B: probes.add(x)
    probes = sorted(probes)
    for lo, hi in zip(probes, probes[1:]):
        x = (lo+hi)/2
        require(min(v(i, x) for i in range(3)) == min(v(pair[0], x), v(pair[1], x)),
                ("omitted coordinate lower", w, f, x))


def update(d, key, value, w, E, mass, ks):
    payload = (w, E, mass, ks)
    if key not in d or value > d[key][0]: d[key] = (value, [payload])
    elif value == d[key][0]: d[key][1].append(payload)


def first_strict_c(base, leader):
    c = 1
    while not base + F(4, 7*c) < leader: c += 1
    return c


def main():
    # All normalized sign vectors, modulo global reversal.
    sign_vectors = sorted({v for p in set(permutations((1, 2, 2)))
                           for s in product((-1, 1), repeat=3)
                           for v0 in [tuple(p[i]*s[i] for i in range(3))]
                           for v in [v0 if v0[0] > 0 else tuple(-x for x in v0)]})
    require(len(sign_vectors) == 12, ("sign vector count", sign_vectors))

    head = rows_parameterized(170)
    require(len(head) == 1951, ("H170 rows", len(head)))
    counts = defaultdict(int)
    leaders = {}
    boundary_hits = 0
    for w, f in head:
        counts[f] += 1
        require(family(w) == (f,), ("cone uniqueness", w, family(w), f))
        rels = all_122_relations(w)
        require(rels == {tuple(U[f][i] if U[f][0] > 0 else -U[f][i] for i in range(3))},
                ("exhaustive signed classification", w, f, rels))
        # Parameter/congruence/parity claims.
        a, b, c = w
        if f == 1:
            require(c == 2*(a+b) and gcd(a,b) == 1 and a % 3 == b % 3, ("F1 params", w))
        elif f == 2:
            require(c == 2*(b-a) and b > 2*a and gcd(a,b) == 1 and a % 3 != b % 3, ("F2 params", w))
        elif f == 3:
            require(a % 2 == 0, ("F3 a even", w)); d = a//2
            require(c == b+d and b > 2*d and gcd(b,d) == 1 and b % 3 == d % 3, ("F3 params", w))
        else:
            require(b % 2 == 0, ("F4 b even", w)); d = b//2
            require(c == a+d and d < a < 2*d and gcd(a,d) == 1 and a % 3 == d % 3, ("F4 params", w))
        require(sum(x % 2 == 0 for x in w) in (1, 2), ("parity", w))

        ray, E, mass, ks = ray_packet(w, f)
        raw, Er, mr = raw_packet(w)
        require((ray, E, mass) == (raw, Er, mr), ("complete ray", w, f))
        update(leaders, ("N", f), min(E), w, E, mass, ks)
        update(leaders, ("P", f), mass, w, E, mass, ks)

        ints, pint, lines, B = integrate_profiles(w, f)
        for i, q in enumerate(claimed_integrals(w, f)):
            if q is not None: require(ints[i] == q, ("integral", w, f, i, ints[i], q))
        require(pint == F(7, 8)*R*R, ("physical integral", w, f, pint))
        envelope_sidecar(w, f, lines, B)

        # A roof equality is genuinely excluded at the strict endpoint.
        u = U[f]
        for k in range(1, max(ks, default=0)+5):
            if k % 3 and any(14*k*abs(u[i]) == 3*(sum(w)-w[i]) for i in range(3)):
                boundary_hits += 1
                C = tuple(k*u[i] for i in range(3))
                require(C not in raw and -k not in ks and k not in ks,
                        ("strict endpoint admitted", w, f, k))

    require(dict(counts) == {1:280, 2:559, 3:744, 4:368}, ("H170 family counts", counts))

    expected = {
        ("N",1):(F(5,77), [(1,10,22)]),
        ("P",1):(F(5,77), [(1,10,22)]),
        ("N",2):(F(51,770), [(1,11,20)]),
        ("P",2):(F(51,770), [(1,11,20)]),
        ("N",3):(F(46,665), [(2,19,20)]),
        ("P",3):(F(173,2660), [(2,19,20)]),
        ("N",4):(F(3,49), [(10,14,17),(13,14,20)]),
        ("P",4):(F(731,12740), [(13,14,20)]),
    }
    for key, (value, ws) in expected.items():
        require(leaders[key][0] == value and [p[0] for p in leaders[key][1]] == ws,
                ("leader/equality locus", key, leaders[key], value, ws))

    # Independent ray-only extension: catches any post-cutoff resurrection through H611.
    ext = rows_parameterized(611)
    ext_counts = defaultdict(int)
    ext_leaders = {}
    for w, f in ext:
        ext_counts[f] += 1
        ray, E, mass, ks = ray_packet(w, f)
        update(ext_leaders, ("N", f), min(E), w, E, mass, ks)
        update(ext_leaders, ("P", f), mass, w, E, mass, ks)
        require(min(E) < TARGET and mass < TARGET, ("6/77 failure", w, f, E, mass))
    require(dict(ext_counts) == {1:3553, 2:7103, 3:9483, 4:4737}, ("H611 counts", ext_counts))
    for key in expected:
        require(ext_leaders[key] == leaders[key], ("leader changed by H611", key, ext_leaders[key], leaders[key]))

    bases = {1:F(87,1568), 2:F(45,784), 3:F(363,6272), 4:F(363,6272)}
    net_cut = {f:first_strict_c(bases[f], expected[("N",f)][0]) for f in range(1,5)}
    phys_base = F(3,56)
    phys_cut = {f:first_strict_c(phys_base, expected[("P",f)][0]) for f in range(1,5)}
    require(net_cut == {1:61,2:65,3:51,4:171}, ("network cutoffs", net_cut))
    require(phys_cut == {1:51,2:46,3:50,4:151}, ("physical cutoffs", phys_cut))
    require(max(v[0] for k,v in expected.items() if k[0] == "N") == F(46,665), "global network")
    require(max(v[0] for k,v in expected.items() if k[0] == "P") == F(51,770), "global physical")
    require(F(46,665) < TARGET and F(51,770) < TARGET, "strict target")

    print("PASS clean-room signed-(1,2,2) referee")
    print("sign_vectors", len(sign_vectors))
    print("H170_counts", dict(sorted(counts.items())), "boundary_hits", boundary_hits)
    print("H611_counts", dict(sorted(ext_counts.items())))
    for key in sorted(expected):
        print(key, leaders[key][0], [p[0] for p in leaders[key][1]])
    print("network_cutoffs", net_cut, "physical_cutoffs", phys_cut)
    print("global", F(46,665), F(51,770), "target", TARGET)
    print("checks", CHECKS)


if __name__ == "__main__":
    main()
