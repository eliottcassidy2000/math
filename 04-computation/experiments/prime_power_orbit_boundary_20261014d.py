"""Exact primitive-orbit inversion and a first-hit refinement boundary.

Universe: all labelled functional maps on 1..5 vertices, periods 1..12;
the four explicitly supplied nonzero signed Collatz cycles; and stated
finite controls of an analytically defined countable family. No search
for a new Collatz terminal, nor universal coverage inference, occurs.
Run from the repository root; stdout is deterministic JSON. Checks stay
active with python -O.
"""
from fractions import Fraction
from itertools import product
import json

CHECKS = 0


def require(condition, message):
    global CHECKS
    CHECKS += 1
    if not condition:
        raise ValueError(message)


def mobius(n):
    sign = 1
    p = 2
    while p * p <= n:
        if n % p == 0:
            n //= p
            sign = -sign
            if n % p == 0:
                return 0
        p += 1
    return -sign if n > 1 else sign


def cycles_of(f):
    cycles = set()
    for s in range(len(f)):
        path, seen = [], {}
        x = s
        while x not in seen:
            seen[x] = len(path)
            path.append(x)
            x = f[x]
        cycle = path[seen[x]:]
        k = cycle.index(min(cycle))
        cycles.add(tuple(cycle[k:] + cycle[:k]))
    return sorted(cycles)


def first_hit(f, source, terminals):
    seen, x, steps = set(), source, 0
    while x not in terminals:
        if x in seen:
            return None
        seen.add(x)
        x = f[x]
        steps += 1
    return steps


def trace_coefficients(fixed, degree):
    # zeta'/zeta = sum_{m>=1} fixed[m] z^(m-1).
    out = [Fraction(1)]
    for n in range(1, degree + 1):
        out.append(sum(Fraction(fixed[k]) * out[n-k]
                       for k in range(1, n+1)) / n)
    return out


def product_coefficients(lengths, degree):
    out = [1] + [0] * degree
    for length in lengths:
        for n in range(length, degree + 1):
            out[n] += out[n-length]
    return out


def finite_maps():
    degree = 12
    rows = []
    for n in range(1, 6):
        maps = 0
        for f in product(range(n), repeat=n):
            maps += 1
            cycles = cycles_of(f)
            lengths = [len(c) for c in cycles]
            actual = {m: lengths.count(m) for m in range(1, degree+1)}
            fixed = [0] * (degree+1)
            powers = list(range(n))
            for m in range(1, degree+1):
                powers = [f[x] for x in powers]
                fixed[m] = sum(x == y for x, y in enumerate(powers))
                require(fixed[m] == sum(k for k in lengths if m % k == 0),
                        "primitive-cycle trace")
                inverse = Fraction(sum(mobius(d) * fixed[m//d]
                                       for d in range(1, m+1) if m % d == 0), m)
                require(inverse == actual[m], "cycle Mobius inversion")
            require(trace_coefficients(fixed, degree) ==
                    product_coefficients(lengths, degree), "Euler product")
            terminals = set().union(*(set(c) for c in cycles))
            for s in range(n):
                hit = first_hit(f, s, terminals)
                require(hit is not None and hit < n, "finite all-cycle coverage")
            marked = set(cycles[0])
            h = [Fraction(0) if (t := first_hit(f, s, marked)) is None
                 else Fraction(1, 2**t) for s in range(n)]
            for s in range(n):
                require(h[s] == (1 if s in marked else h[f[s]] / 2),
                        "discounted first-hit equation")
        rows.append({"vertices": n, "maps": maps})
    return {"universe": rows, "periods_and_series_degree": degree}


def shortcut(n):
    return n // 2 if n % 2 == 0 else (3*n+1)//2


def odd_map(n):
    x = 3*n+1
    while x % 2 == 0:
        x //= 2
    return x


def supplied_cycles():
    cycles = [(1, 2), (-1,), (-5, -7, -10),
              (-17, -25, -37, -55, -82, -41, -61, -91, -136, -68, -34)]
    rows = []
    for c in cycles:
        require(len(set(c)) == len(c), "primitive supplied shortcut cycle")
        for a, b in zip(c, c[1:]+c[:1]):
            require(shortcut(a) == b, "supplied shortcut edge")
        odd = tuple(x for x in c if x % 2)
        for a, b in zip(odd, odd[1:]+odd[:1]):
            require(odd_map(a) == b, "supplied odd edge")
        rows.append({"shortcut_cycle": c, "odd_cycle": odd,
                     "shortcut_length": len(c), "odd_length": len(odd)})
    for lengths in ([len(c) for c in cycles],
                    [sum(x % 2 for x in c) for c in cycles]):
        fixed = [0] + [sum(k for k in lengths if m % k == 0)
                       for m in range(1, 133)]
        for m in range(1, 133):
            inverse = Fraction(sum(mobius(d)*fixed[m//d]
                                   for d in range(1, m+1) if m % d == 0), m)
            require(inverse == lengths.count(m), "catalogue inversion")
        require(trace_coefficients(fixed, 132) == product_coefficients(lengths, 132),
                "catalogue product independent route")
    return {"scope": "the four supplied nonzero cycles only, completeness OPEN",
            "cycles": rows}


def truncation_step(n, j):
    # -1 denotes c, which is fixed. All j >= 0 denote r_j.
    if j == -1 or j == n:
        return -1
    return j + 1 if j < n else j - 1


def refinement_controls():
    for n in range(65):
        for j in range(97):
            t = abs(j-n) + 1
            x = j
            for k in range(t):
                require(x != -1, "no premature terminal")
                x = truncation_step(n, x)
            require(x == -1, "exact terminal clock")
            h = Fraction(1, 2**t)
            child = truncation_step(n, j)
            child_h = (Fraction(1) if child == -1
                       else Fraction(1, 2**(abs(child-n)+1)))
            require(h == child_h / 2, "moving terminal first-hit equation")
            if n > j:
                require(truncation_step(n, j) == j+1, "pointwise local agreement")
    # Any prescribed finite initial segment of the escaping ray is replicated
    # by rooted approximants; the first exit lies immediately beyond the segment.
    for depth in range(129):
        for n in (depth, depth+1, depth+17):
            x = 0
            for k in range(depth):
                require(truncation_step(n, x) == x+1, "finite-prefix hostile")
                x += 1
            require(Fraction(1, 2**(n+1)) <= Fraction(1, 2**(depth+1)),
                    "vanishing discounted mass")
    samples = [{"N": n, "source": "r_0", "first_hit": n+1,
                "H_at_z_half": f"1/2^{n+1}"}
               for n in (0, 1, 2, 4, 8, 16, 32, 64, 128, 1024)]
    return {"exact_control_N": [0, 64], "exact_control_j": [0, 96],
            "prefix_depths": [0, 128], "samples": samples,
            "all_N_theorem": "same sole fixed cycle; rooted for every finite N; pointwise limit has an escaping ray",
            "warning": "countermodel to an inference, not a counterexample to Collatz"}


def main():
    data = {"status": "PROVED elementary identities and refinement criterion; FINITE-EXACT controls; universal Collatz OPEN",
            "finite_maps": finite_maps(), "supplied_collatz_cycles": supplied_cycles(),
            "refinement_boundary": refinement_controls()}
    data["checks"] = CHECKS
    print(json.dumps(data, indent=2, sort_keys=True))


if __name__ == "__main__":
    main()
