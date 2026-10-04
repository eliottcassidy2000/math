"""Exact CRT, golden-phase, coset, and guarded Collatz comparisons.

Standard library only; explicit checks survive -O. No inherited code imports.
Run from any directory; deterministic JSON and stdout are written to results.
"""
from fractions import Fraction
from itertools import permutations
from math import gcd, isqrt
from pathlib import Path
import json

CHECKS = 0


def check(condition, label):
    global CHECKS
    CHECKS += 1
    if not condition:
        raise ArithmeticError(label)


def prime(p):
    return p >= 2 and all(p % d for d in range(2, isqrt(p) + 1))


def order(a, p):
    value = 1
    for k in range(1, p):
        value = value * a % p
        if value == 1:
            return k
    raise ArithmeticError("nonunit")


def qmul(z, w, p):
    a, b = z
    c, d = w
    return ((a*c + b*d) % p, (a*d + b*c + b*d) % p)


def qpow(z, k, p):
    out = (1, 0)
    while k:
        if k & 1:
            out = qmul(out, z, p)
        z = qmul(z, z, p)
        k //= 2
    return out


def qorder(z, p):
    out = (1, 0)
    for k in range(1, p*p):
        out = qmul(out, z, p)
        if out == (1, 0):
            return k
    raise ArithmeticError("nonunit in quadratic algebra")


def fib(n):
    a, b = 0, 1
    for _ in range(n):
        a, b = b, a+b
    return a


def step(n):
    z, a = 3*n+1, 0
    while z % 2 == 0:
        z //= 2
        a += 1
    return z, a


def replay(n, word):
    nodes = [n]
    for a in word:
        n, actual = step(n)
        check(a == actual, "valuation guard")
        nodes.append(n)
    return nodes


def word_data(w):
    A = B = 0
    for a in w:
        B = 3*B + 2**A
        A += a
    return len(w), A, B


def compositions(S):
    if S == 0:
        yield ()
    else:
        for a in range(1, S+1):
            for w in compositions(S-a):
                yield (a,) + w


def cuts_compositions(S):
    for mask in range(1 << (S-1)):
        cuts = [0] + [i for i in range(1, S) if mask >> (i-1) & 1] + [S]
        yield tuple(y-x for x, y in zip(cuts, cuts[1:]))


def prime_multiple_source(w, p):
    r, A, B = word_data(w)
    modulus = 2**(A+1)
    residue = (2**A-B)*pow(3**r, -1, modulus) % modulus
    source = residue + modulus*((-residue*pow(modulus, -1, p)) % p)
    return source or modulus*p


def crt_and_phases():
    fiber = [(x, y) for x in range(18) for y in range(12) if (x-y) % 6 == 0]
    check(set(fiber) == {(n % 18, n % 12) for n in range(36)}, "CRT pullback")
    check(len(fiber) == 36, "CRT uniqueness")
    H = [1+6*t for t in range(6)]
    for s in range(6):
        for t in range(6):
            check((1+6*s)*(1+6*t) % 36 == (1+6*((s+t) % 6)), "square-zero law")
    check(1 % 18 == 19 % 18 and 1 % 12 != 19 % 12, "no quotient 18 to 12")
    check(7*7 % 216 != (1+6*2) % 216, "linearization loses quadratic term")
    check({x*x % 18 for x in range(18) if gcd(x, 18) == 1} == {1, 7, 13}, "squares18")
    check({x*x % 9 for x in range(9) if gcd(x, 9) == 1} == {1, 4, 7}, "squares9")
    check(order(2, 13) == 12 and order(2, 19) == 18, "phase periods")
    check(order(2, 247) == 36, "minimal joined period")
    for a in range(1, 37):
        for n in range(247):
            y = (3*n+1)*pow(2, -a, 247) % 247
            for p in (13, 19):
                check(y % p == (3*n+1)*pow(2, -a, p) % p, "phase projection")
    rows = []
    for a in H:
        rows.append({"a": a, "mod18": a % 18, "mod9": a % 9, "mod12": a % 12,
                     "C19": [3*pow(2, -a, 19) % 19, pow(2, -a, 19)],
                     "C13": [3*pow(2, -a, 13) % 13, pow(2, -a, 13)]})
    # CRT guards realize every chosen odd residue with every valuation.
    for p in (13, 19):
        for a in range(1, 37):
            M = 2**(a+1)
            r = (2**a-1)*pow(3, -1, M) % M
            for x in range(p):
                n = r + M*((x-r)*pow(M, -1, p) % p)
                y, actual = step(n)
                check(actual == a and n % p == x, "CRT actual source")
                check(y % p == (3*x+1)*pow(2, -a, p) % p, "actual phase law")
    return rows


def golden():
    periods = {}
    for m in (6, 9, 12, 18, 36):
        a, b = 0, 1
        for k in range(1, 1000):
            a, b = b, (a+b) % m
            if (a, b) == (0, 1):
                break
        check(k == 24, "Pisano24")
        periods[m] = k
    check([fib(1+6*t) % 36 for t in range(5)] == [1, 13, 17, 5, 1], "six-step four-cycle")
    rows = []
    for p, expected_phi, expected_u in ((223, 448, 224), (233, 52, 13)):
        check(prime(p), "prime modulus")
        check(not any((x*x-x-1) % p == 0 for x in range(p)), "golden irreducible")
        phi, u = (0, 1), (p-1, p-1)
        for k in range(1, 2*(p+1)+1):
            check(qpow(phi, k, p) == (fib(k-1) % p, fib(k) % p), "independent Fibonacci power")
        check(qorder(phi, p) == expected_phi, "phi order")
        check(qorder(u, p) == expected_u, "golden norm-one order")
        norm_one = {(a, b) for a in range(p) for b in range(p) if (a*a+a*b-b*b) % p == 1}
        check(len(norm_one) == p+1, "norm-one census")
        rows.append({"p": p, "phi_order": expected_phi, "u_order": expected_u,
                     "norm_one_size": len(norm_one)})
    p = 233
    u, generator, tau = (232, 232), (4, 189), (42, 77)
    check(qpow((0, 1), 13, p) == (144, 0), "phi13")
    check(144**2 % p == p-1, "phi26")
    check(qorder(generator, p) == 234, "norm-one generator")
    check(qpow(generator, 13, p) == tau and qorder(tau, p) == 18, "order18 factor")
    products = {qmul(qpow(tau, i, p), qpow(u, j, p), p) for i in range(18) for j in range(13)}
    check(products == norm_one and len(products) == 234, "direct product covers norm-one circle")
    for n in range(13):
        for a in range(1, 37):
            y = (3*n+1)*pow(2, -a, 13) % 13
            check(qpow(u, y, p) == qpow(qmul(u, qpow(qpow(u, n, p), 3, p), p), pow(2, -a, 13), p), "golden Collatz phase")
    return {"Pisano": periods, "quadratic_fields": rows,
            "order18_generator": tau, "order13_generator": u,
            "Fibonacci_sample": [fib(1+6*t) for t in range(5)]}


def cosets(p, d, h):
    check(prime(p) and order(2, p) == h and order(3, p) == p-1, "coset orders")
    H = {pow(2, i, p) for i in range(h)}
    check(H == {pow(x, d, p) for x in range(1, p)}, "power subgroup")
    classes = [{pow(3, j, p)*x % p for x in H} for j in range(d)]
    label = {x: j for j, C in enumerate(classes) for x in C}
    check(len(label) == p-1, "coset partition")
    counts = [[0]*d for _ in range(d)]
    for x in range(1, p):
        y = (3*x+1) % p
        if y:
            counts[label[x]][label[y]] += 1
    missing = [(i, j) for i in range(d) for j in range(d) if counts[i][j] == 0]
    check(missing == [(d-1, d//2)], "unique forbidden transition")
    check(not any((u+v+1) % p == 0 for u in H for v in H), "subgroup Fermat obstruction")
    # Independent projective-chart enumeration, including coordinate zeros.
    affine = sum((pow(x, d, p)+pow(y, d, p)+1) % p == 0 for x in range(p) for y in range(p))
    infinity = sum((pow(x, d, p)+1) % p == 0 for x in range(p))
    check(affine+infinity == 0, "independent empty Fermat curve")
    # Hostile control: quadratic curves are not empty at these primes.
    check(any((x*x+y*y+1) % p == 0 for x in range(p) for y in range(p)), "nonempty positive control")
    return {"p": p, "index": d, "ord2": h, "ord3": p-1,
            "missing": missing, "transition_counts": counts, "Fermat_projective_points": 0}


def returns(p, budget):
    hits = []
    for S in range(1, budget+1):
        words = list(compositions(S))
        check(set(words) == set(cuts_compositions(S)), "independent word universe")
        for w in words:
            r, A, B = word_data(w)
            B_direct = sum(3**(r-1-i)*2**sum(w[:i]) for i in range(r))
            check(B == B_direct, "independent ordered carry")
            if B % p == 0:
                gap = 2**A-3**r
                check(gap > 0 and p*gap > B, "all positive p-multiple returns descend")
                n = prime_multiple_source(w, p)
                ns = replay(n, w)
                check(ns[-1] < n and ns[-1] % p == 0, "guarded descending return")
                hits.append({"word": w, "r": r, "A": A, "B": B,
                             "threshold": str(Fraction(B, gap)),
                             "least_guarded_source": n, "endpoint": ns[-1],
                             "source_period": p*2**(A+1)})
    exact_prefixes = []
    for S in range(8):
        for prefix in compositions(S):
            if word_data(prefix+(1,))[2] == p:
                exact_prefixes.append(prefix)
    check(exact_prefixes == ([(2, 3, 1)] if p == 223 else [(5, 2)]), "complete exact-carry prefix search")
    hostile_word = (1, 1, 1, 1, 2, 1, 1, 1) if p == 223 else (1, 3, 1, 1, 1, 3, 1)
    n = prime_multiple_source(hostile_word, p)
    ns = replay(n, hostile_word)
    check(sum(hostile_word) == budget+1, "sharp budget")
    check(ns[-1] > n and ns[-1] % p == 0 and all(x % p for x in ns[1:-1]), "expanding first return")
    return {"p": p, "budget": budget, "compositions_checked": 2**budget-1,
            "return_words": hits, "exact_carry_prefix": exact_prefixes[0],
            "hostile_word": hostile_word, "hostile_nodes": ns}


def return_entry_census():
    """One complete source period, independent of the word-carry census."""
    eligible = []
    for k in range(1, 2048, 2):
        n = 233*k
        x, A = n, 0
        while A <= 10:
            x, a = step(x)
            A += a
            if A > 10:
                break
            if x % 233 == 0:
                check(x < n, "source census descent")
                eligible.append(n)
                break
    check(len(eligible) == 8 and min(eligible) == 55221, "complete entry period")
    check(Fraction(len(eligible), 1024) == Fraction(1, 128), "entry density")
    return {"source_period": 233*2048, "odd_233_multiples": 1024,
            "eligible": eligible, "relative_density": "1/128",
            "least_eligible": 55221, "eligible_at_most_10000": 0}


def source_routes_and_graph():
    routes = {}
    for start in (111, 223, 233):
        n, ns, ks = start, [start], []
        while n != 1:
            check(len(ns) < 100, "bounded trajectory control")
            n, a = step(n)
            ns.append(n)
            ks.append(a)
        routes[start] = {"nodes": ns, "word": ks}
    check(routes[223]["nodes"][5] == routes[233]["nodes"][10] == routes[111]["nodes"][5] == 425, "common future")
    # Independently replay the inherited 233 -> 223 tournament response square.
    base = [252, 169, 242, 164, 202, 144, 42, 64]
    counts = []
    for flips in ((), ((0, 3),), ((2, 3),), ((0, 3), (2, 3))):
        adj = base.copy()
        for u, v in flips:
            adj[u] ^= 1 << v
            adj[v] ^= 1 << u
        H = sum(all(adj[u] >> v & 1 for u, v in zip(w, w[1:])) for w in permutations(range(8)))
        counts.append(H)
    check(counts == [233, 291, 123, 223], "inherited tournament response")
    check(counts[3]-counts[1]-counts[2]+counts[0] == 42, "two-flip interaction")
    check(2**13-3**8 == 7*233, "golden shape gap")
    return {"routes": routes, "tournament_response": counts}


def main():
    data = {"crt_phase_table": crt_and_phases(), "golden": golden(),
            "cosets": [cosets(223, 6, 37), cosets(233, 8, 29)],
            "return_budgets": [returns(223, 8), returns(233, 10)],
            "return_entry_census": return_entry_census(),
            "sources_and_graph": source_routes_and_graph()}
    data["checks"] = CHECKS
    directory = Path(__file__).resolve().parents[2] / "05-knowledge/results"
    stem = "triplet_crt_223_233_20261004"
    (directory / (stem+".json")).write_text(json.dumps(data, indent=2)+"\n")
    lines = ["CRT bridge: Z/36 = Z/18 x_(Z/6) Z/12 = Z/9 x Z/4.",
             "Valuation phases: order_19(2)=18, order_13(2)=12, order_247(2)=36.",
             "Golden 233: phi order52, -phi^2 order13; norm-one group order234=18*13.",
             "Fibonacci pair periods at moduli6,9,12,18,36: all24.",
             "223: index6, unique forbidden coset5->3, empty Fermat sextic.",
             "233: index8, unique forbidden coset7->4, empty Fermat octic.",
             "Sharp descending return budgets: 223 cost<=8; 233 cost<=10.",
             "233 rule entry: least source55221; relative density1/128 among odd233-multiples.",
             "Common future: U^5(223)=U^5(111)=U^10(233)=425.",
             "Inherited tournament response independently replayed: 233,291,123,223.",
             f"PASS: {CHECKS} explicit checks; standard library; no assert statements."]
    out = "\n".join(lines)+"\n"
    (directory / (stem+".out")).write_text(out)
    print(out, end="")


if __name__ == "__main__":
    main()
