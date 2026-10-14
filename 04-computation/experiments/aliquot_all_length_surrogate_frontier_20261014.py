#!/usr/bin/env python3
"""Bounded falsification and exact restrictions for an open all-length question.

No result here proves either all-length conjecture. Graph checks concern
cycles wholly inside the finite window, not convergence of all starts.
"""
from collections import Counter
from itertools import product
import json


def require(condition, message):
    if not condition:
        raise RuntimeError(message)


def core(n):
    return n & -n


def f(n):
    t = core(n)
    p = n // t
    return (t-1)*p+2*t-1


def canonical(word):
    return min(tuple(word[i:]+word[:i]) for i in range(len(word)))


def cycle_packet(cycle):
    k = len(cycle)
    require(all(f(n) == cycle[(i+1)%k] for i, n in enumerate(cycle)), "cycle edge")
    ts = [core(n) for n in cycle]
    ps = [n//t for n, t in zip(cycle, ts)]
    x = min(ts)
    require(sum(ps) == 2*sum(ts)-k, "global cofactor identity")
    require(all((p+1) % x == 0 for p in ps), "minimum-core divisibility")
    hs, us = [(p+1)//x for p in ps], [t//x for t in ts]
    require(sum(hs) == 2*sum(us), "normalized sum identity")
    require(all((hs[i]-us[i]-us[(i+1)%k]) % x == 0 for i in range(k)),
            "successor-indexed congruence")
    if 1 in ps:
        require(x == 2, "cofactor-one minimum-core boundary")
    return {"cycle": cycle, "cores": ts, "cofactors": ps, "H": hs, "u": us,
            "minimum_core": x, "primitive_length": k}


def finite_graph(bound):
    seen = bytearray(bound//2+1)
    cycles = []
    local_counts = Counter()
    for n in range(2, bound+1, 2):
        t = core(n)
        p = n//t
        m = f(n)
        t2 = core(m)
        require(m > 0 and m % 2 == 0, "evenness")
        if m > n:
            local_counts["increasing"] += 1
            if t2 >= t:
                require(p == t-1 and m == t*t and m//t2 == 1,
                        "exceptional increase must land on cofactor one")
                local_counts["exceptional_quadratic_increase"] += 1
            else:
                local_counts["increase_with_core_drop"] += 1
        elif m < n:
            local_counts["decreasing"] += 1
        else:
            require(p == 2*t-1, "fixed-point formula")
            local_counts["fixed"] += 1
        if t2 >= t and m != n and m//t2 != 1:
            require(m < n, "core-size tradeoff")
        a = t.bit_length()-1
        if a % 2:
            require((m+n) % 3 == 0, "odd-exponent mod-three rule")
        if p == 3 and a % 2 and a >= 3:
            require(t2 == 4, "odd-exponent prime-three escape")
    for start in range(2, bound+1, 2):
        if seen[start//2]:
            continue
        path, positions = [], {}
        n = start
        while 2 <= n <= bound and not seen[n//2] and n not in positions:
            positions[n] = len(path)
            path.append(n)
            n = f(n)
        if n in positions:
            cycles.append(canonical(path[positions[n]:]))
        for m in path:
            seen[m//2] = 1
    cycles.sort()
    nonconstant = [c for c in cycles if len(c)>1]
    require(nonconstant == [(8,22,14,10), (16,46,26)], "finite nonconstant cycle list")
    fixed = [c[0] for c in cycles if len(c)==1]
    require(fixed == [6,28,120,496,2016,8128,32640,130816,523776], "finite fixed-point list")
    return {"bound_all_cycle_members": bound, "vertices": bound//2,
            "claim": "Only cycles wholly contained in the finite window are classified",
            "histogram": dict(Counter(map(len, cycles))),
            "cycles": [cycle_packet(list(c)) for c in cycles],
            "local_restriction_inputs_checked": bound//2,
            "local_counts": dict(local_counts)}


def template_checks():
    results = []
    for length in range(4,8):
        checked, closures = 0, []
        for exponents in product(range(1,7), repeat=length):
            if len(set(exponents)) == 1:
                continue
            if exponents != min(exponents[j:]+exponents[:j] for j in range(length)):
                continue
            checked += 1
            ts = [1 << a for a in exponents]
            slope_num, intercept_num, denominator = 1, 0, 1
            for i, t in enumerate(ts):
                slope_num, intercept_num, denominator = (
                    (t-1)*slope_num,
                    (t-1)*intercept_num+(2*t-1)*denominator,
                    ts[(i+1)%length]*denominator)
            p, remainder = divmod(intercept_num, denominator-slope_num)
            if remainder:
                continue
            ps = [p]
            for i in range(length-1):
                p, remainder = divmod((ts[i]-1)*p+2*ts[i]-1, ts[i+1])
                require(remainder == 0, "integral closure propagation")
                ps.append(p)
            numbers = [t*p for t,p in zip(ts,ps)]
            require(all(f(n) == numbers[(i+1)%length] for i,n in enumerate(numbers)),
                    "template/deterministic-map equivalence")
            require(1 in ps, "bounded counterexample to cofactor-one conjecture")
            closures.append({"exponents": exponents, "cofactors": ps, "numbers": numbers})
        results.append({"word_length": length, "nonconstant_cyclic_representatives": checked,
                        "integral_closures": closures})
    require([r["nonconstant_cyclic_representatives"] for r in results]
            == [330,1554,7820,39990], "template representative counts")
    require([len(r["integral_closures"]) for r in results] == [1,0,1,0], "bounded closures")
    return {"exponent_range": [1,6], "word_length_range": [4,7], "results": results,
            "warning": "Word length can be a multiple of primitive cycle period"}


def main():
    result = {
        "status": "open all-length conjectures; proved local and modular restrictions",
        "C1": "Every nonconstant finite F-cycle contains cofactor one (OPEN)",
        "C2": "The only nonconstant primitive F-cycles are [8,22,14,10] and [16,46,26] (OPEN)",
        "finite_graph": finite_graph(1_000_000),
        "template_checks": template_checks(),
        "proved_prime_specific_restriction": "Odd cycle length and all dyadic exponents odd force the fixed point 6",
        "not_claimed": ["convergence of every start in the window", "global boundedness of F orbits",
                        "all-length prime-chamber cycle classification", "general aliquot cycle classification"],
    }
    print(json.dumps(result, indent=2, sort_keys=True))


if __name__ == "__main__":
    main()
