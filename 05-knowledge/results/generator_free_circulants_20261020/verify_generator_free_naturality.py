#!/usr/bin/env python3
"""Independent exact controls. No third-party dependencies; checks survive -O."""
import argparse
import json
from collections import Counter
from math import gcd

CHECKS = 0

def check(test, label):
    global CHECKS
    CHECKS += 1
    if not test:
        raise RuntimeError(label)


def canon(w):
    return min(w[i:] + w[:i] for i in range(len(w)))


def shift(w, k):
    k %= len(w)
    return w[k:] + w[:k]


def decimate(w, a):
    return tuple(w[(a * i) % len(w)] for i in range(len(w)))


def primitive(w):
    for d in divisors(len(w)):
        if w == w[:d] * (len(w) // d):
            return w[:d]
    raise RuntimeError("missing primitive root")


def divisors(n):
    return [d for d in range(1, n + 1) if n % d == 0]


def words(m):
    for bits in range(1 << m):
        yield tuple((bits >> i) & 1 for i in range(m))


def antiperiodic(w):
    return w + tuple(1 - x for x in w)


def projection(f):
    """H consists of even exponents; assign the same result to antipodes."""
    n = len(f)
    m = n // 2
    if m % 2 != 1:
        raise ValueError("odd m required")
    t = sum(f[2 * j] for j in range(m)) % 2
    out = [None] * n
    for j in range(m):
        h = 2 * j
        out[h] = out[(h + m) % n] = f[h] ^ t
    return tuple(out)


def transitions(w):
    return sum(w[i] != w[(i + 1) % len(w)] for i in range(len(w)))


def signed_components(n, a, k):
    """Number of skew words fixed by i -> a*i+k on Z/n."""
    adj = [[] for _ in range(n)]
    for i in range(n):
        for j, b in (((a * i + k) % n, 0), ((i + n // 2) % n, 1)):
            adj[i].append((j, b))
            adj[j].append((i, b))
    seen = {}
    components = 0
    for i in range(n):
        if i in seen:
            continue
        seen[i] = 0
        components += 1
        stack = [i]
        while stack:
            x = stack.pop()
            for y, b in adj[x]:
                v = seen[x] ^ b
                if y in seen:
                    if seen[y] != v:
                        return 0
                else:
                    seen[y] = v
                    stack.append(y)
    return 1 << components


def cyclic_counts(d, a, k):
    """All and even-weight words fixed by affine permutation on Z/d."""
    seen = set()
    sizes = []
    for i in range(d):
        if i in seen:
            continue
        c = 0
        j = i
        while j not in seen:
            seen.add(j)
            c += 1
            j = (a * j + k) % d
        sizes.append(c)
    return 1 << len(sizes), 1 << (len(sizes) - int(any(c % 2 for c in sizes)))


def evenroot_fixed(m, a, k):
    """Separate exact primitive period and parity by divisor recursion."""
    total_exact = {}
    even_exact = {}
    for d in divisors(m):
        total, even = cyclic_counts(d, a, k)
        smaller = [e for e in total_exact if d % e == 0]
        total_exact[d] = total - sum(total_exact[e] for e in smaller)
        even_exact[d] = even - sum(
            even_exact[e] if (d // e) % 2 else total_exact[e]
            for e in smaller
        )
    return sum(even_exact.values())


def type_fixed(m, a):
    sx = sum(signed_components(2 * m, a, k) for k in range(2 * m))
    sy = sum(evenroot_fixed(m, a, k) for k in range(m))
    check(sx % (2 * m) == 0 and sy % m == 0, "Burnside integrality")
    return sx // (2 * m), sy // m


def verify_small(max_m):
    rows = []
    for m in range(1, max_m + 1):
        xs = [antiperiodic(w) for w in words(m)]
        ys = [w for w in words(m) if sum(primitive(w)) % 2 == 0]
        xtypes = {canon(f) for f in xs}
        ytypes = {canon(w) for w in ys}
        check(len(xtypes) == len(ytypes), f"type equinumerosity m={m}")
        units = [a for a in range(1, 2 * m) if gcd(a, 2 * m) == 1]
        fixed_rows = []
        for a in units:
            brute = (
                sum(canon(decimate(f, a)) == f for f in xtypes),
                sum(canon(decimate(w, a)) == w for w in ytypes),
            )
            check(brute == type_fixed(m, a), f"independent fixed counts m={m},a={a}")
            fixed_rows.append([a, *brute])
        if m % 2:
            outputs = Counter()
            orbit_map = {}
            for f in xs:
                out = projection(f)
                w = out[:m]
                check(out == w + w, "output symmetric")
                check(sum(primitive(w)) % 2 == 0, "output target")
                outputs[w] += 1
                key = canon(f)
                value = canon(w)
                check(key not in orbit_map or orbit_map[key] == value, "well-defined orbit map")
                orbit_map[key] = value
                # Generator equivariance implies all shifts; explicit z control too.
                check(projection(shift(f, 1)) == shift(out, 1), "multiplier equivariance")
                check(projection(shift(f, m)) == out, "antipodal fiber")
                for a in units:
                    check(projection(decimate(f, a)) == decimate(out, a), "Aut(C) naturality")
            check(set(outputs) == set(ys), "labelled onto")
            check(set(outputs.values()) == {2}, "exact two-element fibers")
            check(set(orbit_map.values()) == ytypes, "type onto")
            check(len(set(orbit_map.values())) == len(orbit_map), "type injective")
        elif m >= 4:
            a = m + 1
            f = (0,) * m + (1,) * m
            other = decimate(f, a)
            check(transitions(f) == 2, "source transition count")
            check(transitions(other) == 2 * m - 2, "decimated transition count")
            check(canon(f) != canon(other), "moved source type")
            for w in words(m):
                check(decimate(w + w, a) == w + w, "every graph fixed")
            if m % 4 == 0:
                check(type_fixed(m, a)[0] == 0, "stronger fixed-point-free control")
        else:
            check(m == 2 and len(xtypes) == len(ytypes) == 1, "exception m=2")
        rows.append({"m": m, "source_types": len(xtypes), "target_types": len(ytypes), "fixed": fixed_rows})
    return rows


def verify_frobenius():
    expected = {
        (3, 2): ((2, 2), (0, 2)),
        (5, 2): ((172, 172), (0, 46)),
        (3, 3): ((316, 316), (16, 16)),
        (7, 2): ((349536, 349536), (0, 3372)),
        (3, 4): ((13743895360, 13743895360), (0, 2522)),
    }
    rows = []
    for (p, r), want in expected.items():
        q = p ** r
        m = (q - 1) // 2
        total = type_fixed(m, 1)
        fixed = type_fixed(m, p)
        check((total, fixed) == want, f"Frobenius table q={q}")
        if r % 2 == 0:
            check(fixed[0] == 0 and fixed[1] >= 1, "even-degree obstruction")
        rows.append({"p": p, "r": r, "q": q, "types": total, "frobenius_fixed": fixed})
    return rows


def verify_rook():
    # F9 element a+b*t, t^2=-1=2. Independent coordinate implementation.
    def add(x, y): return ((x[0] + y[0]) % 3, (x[1] + y[1]) % 3)
    def neg(x): return ((-x[0]) % 3, (-x[1]) % 3)
    def sub(x, y): return add(x, neg(y))
    def mul(x, y): return ((x[0]*y[0] + 2*x[1]*y[1]) % 3, (x[0]*y[1] + x[1]*y[0]) % 3)
    g = (1, 1)
    powers = []
    x = (1, 0)
    for _ in range(8):
        powers.append(x)
        x = mul(x, g)
    check(len(set(powers)) == 8 and x == (1, 0), "F9 primitive element")
    conn = {powers[i] for i in (2, 3, 6, 7)}
    e, f = (0, 1), (1, 2)
    check(conn == {e, neg(e), f, neg(f)}, "rook connection set")
    basis = {}
    for u in range(3):
        for v in range(3):
            basis[(u, v)] = add(mul((u, 0), e), mul((v, 0), f))
    check(len(set(basis.values())) == 9, "rook basis")
    vertices = list(basis)
    edges = set()
    for i, u in enumerate(vertices):
        for v in vertices[i+1:]:
            if sub(basis[v], basis[u]) in conn:
                edges.add(frozenset((u, v)))
    flip = lambda v: ((-v[0]) % 3, v[1])
    check({frozenset(flip(v) for v in edge) for edge in edges} == edges, "rook involution automorphism")
    swapped = sum(flip(u) == v and flip(v) == u for u, v in (tuple(edge) for edge in edges))
    check(swapped == 3, "odd-sign rook involution")
    check(sum(primitive((0, 0, 1, 1))) % 2 == 0, "even primitive root hostile")
    return {"vertices": 9, "edges": len(edges), "involution_swapped_edges": swapped, "full_sign_even": False}


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("--max-m", type=int, default=13)
    parser.add_argument("--output")
    args = parser.parse_args()
    result = {
        "status": "PASS",
        "exhaustive": verify_small(args.max_m),
        "frobenius": verify_frobenius(),
        "rook_hostile": verify_rook(),
        "checks": CHECKS,
    }
    encoded = json.dumps(result, indent=2, sort_keys=True) + "\n"
    if args.output:
        with open(args.output, "w") as f:
            f.write(encoded)
    print(encoded, end="")


if __name__ == "__main__":
    main()
