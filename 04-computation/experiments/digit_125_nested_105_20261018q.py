#!/usr/bin/env python3
"""Exact, deterministic controls for the 125 exponent compiler and 105 archaeology.

Run from any directory; writes a JSON receipt beside the research note.
No floating point or external packages. Assertions are deliberately not used.
"""
from fractions import Fraction
from itertools import combinations, permutations, product
from pathlib import Path
from hashlib import sha256
import json

checks = 0


def check(ok, context):
    global checks
    checks += 1
    if not ok:
        raise RuntimeError(context)


def v2(n):
    if n == 0:
        raise ValueError("valuation of zero")
    return (n & -n).bit_length() - 1


def T(n):
    return (3*n + 1)//2 if n & 1 else n//2


def parity(n, length):
    word = 0
    for i in range(length):
        word |= (n & 1) << i
        n = T(n)
    return word


def post(j, bits):
    # Modulo 2^(bits+2) before division retains bits bits afterward.
    return ((3*pow(125, j, 1 << (bits+2)) + 1)//4) % (1 << bits)


def exponent_for(word, length):
    # Invert one parity bit at a time, without constructing 125**j.
    j = 0
    for k in range(length):
        candidates = [j, j+(1 << k)]
        good = [a for a in candidates
                if parity(post(a, k+1), k+1) == word % (1 << (k+1))]
        check(len(good) == 1, ("unique exponent lift", k, word))
        j = good[0]
    return j


def residue_for(word, length):
    # Independent inverse: ordinary Collatz parity-residue lift first.
    n = 0
    for k in range(length):
        candidates = [n, n+(1 << k)]
        good = [a for a in candidates if parity(a, k+1) == word % (1 << (k+1))]
        check(len(good) == 1, ("unique source lift", k, word))
        n = good[0]
    return n


def discrete_log_for(z, length):
    modulus = 1 << (length+2)
    target = (4*z-1)*pow(3, -1, modulus) % modulus
    j = 0
    for k in range(length):
        mod = 1 << (k+3)
        good = [a for a in (j, j+(1 << k)) if pow(125, a, mod) == target % mod]
        check(len(good) == 1, ("unique power lift", k, target))
        j = good[0]
    return j


def cyclic_corr(a, b):
    n = len(a)
    return tuple(sum(a[x]*b[(x+t) % n] for x in range(n)) for t in range(n))


def layers(chain):
    n = chain[-1]
    ans = [[0]*n for _ in chain[1:]]
    for x in range(1, n):
        for i, (low, high) in enumerate(zip(chain, chain[1:])):
            if x % (n//high) == 0:
                quotient = high//low
                residue = (x//(n//high)) % quotient
                check(residue != 0, ("first occupied subgroup", chain, x))
                ans[i][x] = 1 if residue <= quotient//2 else -1
                break
    return ans


def fold_counts(f):
    # Named runner pair clocks, independent of signed-sequence correlations.
    n = len(f)
    counts = [0]*((n+1)//2)
    for a, b in combinations(range(1, (n+1)//2), 2):
        d = (f[a]*a-f[b]*b) % n
        counts[min(d, n-d)] += 1
    return tuple(counts)


def marked_counts(f, source):
    n = len(f)
    counts = [0]*((n+1)//2)
    for other in range(1, (n+1)//2):
        if other != source:
            d = (f[source]*source-f[other]*other) % n
            counts[min(d, n-d)] += 1
    return tuple(counts)


def matchings(labels):
    if not labels:
        yield ()
        return
    a = labels[0]
    for b in labels[1:]:
        rest = tuple(x for x in labels if x not in (a, b))
        for matching in matchings(rest):
            yield tuple(sorted(((a, b),)+matching))


def half_descent_cover(max_depth=18):
    # (exponent residue, parity word, numerator multiplier, offset, odd steps).
    active = [(0, 0, 3, 1, 1)]  # T^2(125^j) = (3*125^j+1)/4.
    accepted = {}
    census = []
    receipt = sha256()
    examples = []
    for depth in range(1, max_depth+1):
        next_active = []
        new = 0
        denominator = 1 << (depth+2)
        for j, word, a, b, odd in active:
            # The two lifts have identical old bits and opposite next bits.
            branches = []
            for lift in (j, j+(1 << (depth-1))):
                observed = parity(post(lift, depth), depth)
                check(observed % (1 << (depth-1)) == word, "prefix retained")
                bit = observed >> (depth-1)
                branches.append(bit)
                aa = 3*a if bit else a
                bb = 3*b+(denominator//2) if bit else b
                oo = odd+bit
                if 2*(125*aa+bb) <= 125*denominator:
                    check(2*aa < denominator, "strict slope margin")
                    record = (depth, lift, observed, aa, bb, oo)
                    accepted[(depth, lift)] = (aa, bb, oo)
                    receipt.update((str(record)+"\n").encode())
                    if len(examples) < 12:
                        examples.append(record)
                    new += 1
                else:
                    next_active.append((lift, observed, aa, bb, oo))
            check(sorted(branches) == [0, 1], "exact binary refinement")
        active = next_active
        census.append({"tail_depth": depth, "shortcut_clock": depth+2,
                       "new_accepted_cylinders": new, "uncovered_addresses": len(active),
                       "address_modulus": 1 << depth,
                       "covered_fraction": str(1-Fraction(len(active), 1 << depth))})
    replayed = 0
    for j in range(1, 4097):
        n = pow(125, j)
        x = T(T(n))
        odd = 1
        for depth in range(1, max_depth+1):
            odd += x & 1
            x = T(x)
            key = (depth, j % (1 << depth))
            if key in accepted:
                a, b, stored_odd = accepted[key]
                check(x*(1 << (depth+2)) == a*n+b, ("literal affine certificate", j))
                check(2*x <= n and odd == stored_odd, ("actual half descent and clock", j))
                replayed += 1
                break
    return {"census": census, "accepted_cylinders": len(accepted),
            "ordered_certificate_sha256": receipt.hexdigest(), "first_certificates": examples,
            "literal_exponents_tested": 4096, "literal_certificates_accepted": replayed}


def exterior_values(n, d, v):
    delta = 1-4*d*v
    alpha = lambda i: 2*(i+1)*(2*i+1)*v*d
    beta = lambda i: i*(i+1)*delta*d
    u = lambda i: -i*(i+1)*delta
    c = lambda i: d-i-1
    cs = [Fraction(1)]
    for r in range(1, 4):
        cs.append(alpha(n+r)*cs[-1]+(d*beta(n+r)*cs[-2] if r >= 2 else 0))
    es = [Fraction(0), Fraction(0)]
    for length in range(2, 6):
        es.append(u(n+length-1)*es[-1]+c(n+length-1)*cs[length-2])
        direct = 0
        for j in range(1, length):
            term = c(n+j)*cs[j-1]
            for h in range(j+1, length):
                term *= u(n+h)
            direct += term
        check(direct == es[-1], "exterior recurrence/direct agreement")
    return es


def run():
    out = {"status": "FINITE-EXACT controls; elementary proofs in companion note"}
    census = []
    for length in range(1, 15):
        size = 1 << length
        states = [post(j, length) for j in range(size)]
        words = [parity(z, length) for z in states]
        check(len(set(states)) == size, ("exponential bijection", length))
        check(len(set(words)) == size, ("itinerary bijection", length))
        for j in range(size):
            check(post(j+size, length) == states[j], ("period", length, j))
        census.append({"length": length, "states": size, "words": len(set(words))})
    for length in range(1, 10):
        for word in range(1 << length):
            a = exponent_for(word, length)
            b = discrete_log_for(residue_for(word, length), length)
            check(a == b, ("independent compilers", length, word))
            check(parity(T(T(pow(125, a))), length) == word,
                  ("literal trajectory", length, word))
    for j in range(129):
        for k in range(j):
            check(v2(pow(125, j)-pow(125, k)) == 2+v2(j-k), ("LTE", j, k))
    hard = []
    for length in (1, 2, 3, 4, 8, 16, 32, 64, 128, 256):
        word = (1 << length)-1
        j = exponent_for(word, length)
        check(j == discrete_log_for((1 << length)-1, length), ("ones compiler", length))
        check((pow(125, j, 1 << (length+2))*3 + 5) % (1 << (length+2)) == 0,
              ("ones congruence", length))
        hard.append({"shortcut_ones_after_T2": length, "least_exponent_residue": str(j),
                     "exponent_modulus": str(1 << length)})
    run_controls = []
    for j in range(1, 1025):
        n = pow(125, j)
        s = (1+v2(j))//2
        x = n
        for _ in range(s):
            check(v2(3*x+1) == 2, ("a2 run", j))
            x = (3*x+1)//4
        check(x == 1 + pow(3, s)*(n-1)//pow(4, s), ("run endpoint", j))
        check(v2(x-1) == 2+v2(j)-2*s, ("exit valuation", j))
        # First post-run valuation is not 2; may be 1 or >=3.
        check(v2(3*x+1) != 2, ("maximal run", j))
        if j in (1, 2, 4, 8, 16, 64, 256, 1024):
            run_controls.append({"j": j, "initial_a2_run": s,
                                 "terminal_minus_one_valuation": v2(x-1)})
    check([pow(125, j, 105) for j in range(1, 11)] == [20, 85]*5, "nonunit CRT")
    out["exponent_compiler"] = {"census": census, "all_ones": hard,
                                 "run_controls": run_controls}
    out["half_descent_cover"] = half_descent_cover()
    innovations = []
    for k in range(13):
        values = [parity(post(j, k+1), k+1) >> k for j in range(1 << k)]
        for j in range(1 << k):
            check((parity(post(j+(1 << k), k+1), k+1) >> k) == 1-values[j],
                  "triangular innovation coefficient")
        coeffs = values[:]
        for bit in range(k):
            for mask in range(1 << k):
                if mask & (1 << bit):
                    coeffs[mask] ^= coeffs[mask ^ (1 << bit)]
        support = [i for i, c in enumerate(coeffs) if c]
        innovations.append({"bit": k, "ANF_terms": len(support),
                            "ANF_degree": max((i.bit_count() for i in support), default=-1)})
    check(all(((105 >> ((a << 2) | (b << 1) | c)) & 1) == (1 ^ a ^ b ^ c)
              for a, b, c in product((0, 1), repeat=3)), "Rule 105 truth table")
    out["triangular_innovation"] = innovations

    chain_results = []
    chains = [(1, 3, 9, 27), (1, 3, 9, 27, 81), (1, 5, 25), (1, 3, 15, 45)]
    for order in permutations((3, 5, 7)):
        ds = [1]
        for p in order:
            ds.append(ds[-1]*p)
        chains.append(tuple(ds))
    chosen = None
    for chain in chains:
        ls = layers(chain)
        n = chain[-1]
        for a, b in combinations(ls, 2):
            check(not any(cyclic_corr(a, b)), ("cross-layer correlation zero", chain))
        vectors = []
        for signs in product((-1, 1), repeat=len(ls)-1):
            f = [sum(s*layer[x] for s, layer in zip(signs+(1,), ls)) for x in range(n)]
            check(f[0] == 0 and all(f[x] == -f[-x] and abs(f[x]) == 1 for x in range(1, n)),
                  ("half-system", chain))
            vectors.append(f)
        corr = [cyclic_corr([int(x == 1) for x in f], [int(x == 1) for x in f]) for f in vectors]
        clocks = [fold_counts(f) for f in vectors]
        check(len(set(corr)) == len(set(clocks)) == 1, ("homometry", chain))
        check(len({tuple(f) for f in vectors}) == 2**(len(chain)-2), ("cube injection", chain))
        chain_results.append({"chain": chain, "distinct_mod_global": len(vectors),
                              "clock_profile": clocks[0]})
        if chain == (1, 5, 35, 105):
            base = vectors[-1]
            h5 = {21, 42}
            h35 = set(range(3, 53, 3))
            for moved in (h5, h35, h5 ^ h35):
                f = list(base)
                for i in moved:
                    f[i] *= -1
                    f[-i] *= -1
                check(fold_counts(f) == fold_counts(base), ("missed June move", sorted(moved)))
            pair = (base, vectors[1])
            witness = next((i for i in range(1, 53)
                            if marked_counts(pair[0], i) != marked_counts(pair[1], i)), None)
            check(witness is not None, "marked source detects silent pair")
            chosen = {"positive_runners": [[i for i in range(1, 53) if f[i] == 1] for f in vectors],
                      "marked_source": witness,
                      "marked_profiles": [marked_counts(f, witness) for f in pair]}
    out["nested_cosets"] = {"chains": chain_results, "C105_explicit": chosen}

    six = tuple(range(6))
    synthemes = list(matchings(six))
    pentads = [p for p in combinations(synthemes, 5)
               if len(set(edge for m in p for edge in m)) == 15]
    check(len(pentads) == 6, "six pentads")
    dual = {}
    for a, b in combinations(six, 2):
        common = set(pentads[a]) & set(pentads[b])
        check(len(common) == 1, "pentad intersection")
        dual[(a, b)] = next(iter(common))
    check(len(set(dual.values())) == 15, "duad-syntheme bijection")
    for a, b in combinations(dual, 2):
        check((bool(set(a) & set(b))) == (not (set(dual[a]) & set(dual[b]))),
              ("compatibility transport", a, b))
    marked_partitions = [(i, tuple(j)) for i in range(7)
                         for j in combinations([x for x in range(7) if x != i], 2)]
    atlas = {}
    for singleton, pair in marked_partitions:
        remaining = [x for x in range(7) if x != singleton]
        positions = tuple(remaining.index(x) for x in pair)
        matching = tuple(sorted(((singleton, 7),)+tuple(
            (remaining[a], remaining[b]) for a, b in dual[positions])))
        atlas[(singleton, pair)] = matching
    check(len(atlas) == len(set(atlas.values())) == 105, "105 dictionary")
    trans = lambda x: 1-x if x < 2 else x
    fixed_parts = sum(i == trans(i) and set(j) == {trans(x) for x in j}
                      for i, j in marked_partitions)
    eight_matchings = list(matchings(tuple(range(8))))
    fixed_matchings = sum(m == tuple(sorted(tuple(sorted((trans(a), trans(b)))) for a, b in m))
                          for m in eight_matchings)
    check((fixed_parts, fixed_matchings) == (35, 15), "non-equivariance hostile")
    sigma = 0
    c0_scaled = 0
    for p in permutations(range(7)):
        arcs = [1 if p[i] < p[i+1] else -1 for i in range(6)]
        sigma += all(arcs[i] == 1 for i in (0, 1, 2, 4))
        prod = 1
        for a in arcs:
            prod *= a
        c0_scaled += prod
    check(sigma == 105, "transitive sigma31")
    check(Fraction(c0_scaled, 64) == Fraction(-17, 4), "transitive coefficient repair")
    out["tournament_matching"] = {"partitions": len(atlas), "transitive_sigma31": sigma,
        "transposition_fixed_partitions": fixed_parts, "transposition_fixed_matchings": fixed_matchings,
        "transitive_c6": 5040, "transitive_c0": str(Fraction(c0_scaled, 64))}
    # Triangle is Eulerian and admits a balanced orientation, but no balanced edge 2-colouring.
    triangle = ((0, 1), (1, 2), (2, 0))
    balanced_colourings = sum(all(sum(c for c, e in zip(colours, triangle) if v in e) == 0
                                    for v in range(3)) for colours in product((-1, 1), repeat=3))
    check(balanced_colourings == 0, "Euler edge-colouring hostile")
    value_edges = {tuple(sorted((abs(3-y), min((3+y) % 9, 9-(3+y) % 9))))
                   for y in (1, 2, 4)}
    check(value_edges == {(1, 2), (1, 4), (2, 4)}, "native C9 value-graph triangle")
    f9 = [0, 1, -1, 1, 1, -1, -1, 1, -1]
    flipped9 = f9[:]
    flipped9[3], flipped9[6] = -1, 1
    check(fold_counts(f9) == fold_counts(flipped9), "Euler orientation realizes native silent flip")
    exterior_count = 0
    for n in range(1, 101):
        for d in range(1, n+11):
            slope = d*(4*(d-n-2)*(n+3)+2*(2*n+3)*(d-n-3))
            constant = (d-n-2)*(n+3)
            if slope:
                v = Fraction(constant, slope)
                es = exterior_values(n, d, v)
                check(es[3] == 0, "exact E3 root")
                check(es[4] != 0 or es[5] != 0, "no triple cancellation")
                exterior_count += 1
        es = exterior_values(n, n+4, Fraction(n+3, 3*(n+4)*(2*n+5)))
        check(es[2] == 2 and es[3] == es[4] == 0 and es[5] != 0, "double cancellation ray")
        es = exterior_values(n, n+3, Fraction(1, 4*(n+3)))
        check(es[3] == 0 and es[4] != 0, "degenerate Delta control")
    out["exterior_three_observation"] = {"E3_zero_specializations": exterior_count,
                                        "double_ray_indices": 100,
                                        "Delta_zero_controls": 100}
    rational_floor_controls = 0
    for n in range(1, 8):
        for d in (-3, -1, 1, 2, n+2, n+3, n+4, n+5):
            for q in range(1, 10):
                for numerator in range(-2, 4):
                    v = Fraction(numerator, q)
                    es = exterior_values(n, d, v)
                    check((es[3] == es[4] == 0) ==
                          (d == n+4 and v == Fraction(n+3, 3*(n+4)*(2*n+5))),
                          "real double-locus control")
                    check(max(abs(x) for x in es[3:6]) >= Fraction(1, v.denominator**3),
                          "rational E floor")
                    delta = 1-4*d*v
                    if delta:
                        beta_n = n*(n+1)*d*delta
                        check(max(abs(beta_n*x) for x in es[3:6]) >=
                              Fraction(n*(n+1)*abs(d), v.denominator**4), "rational V floor")
                    rational_floor_controls += 1
    out["exterior_three_observation"]["rational_floor_controls"] = rational_floor_controls
    out["checks"] = checks
    path = Path(__file__).resolve().parents[2]/"05-knowledge/results/digit_125_nested_105_20261018q.json"
    path.write_text(json.dumps(out, indent=2, sort_keys=True)+"\n")
    print(json.dumps({"checks": checks, "output": path.name, "C105_marked_source": chosen["marked_source"]}))


if __name__ == "__main__":
    run()
