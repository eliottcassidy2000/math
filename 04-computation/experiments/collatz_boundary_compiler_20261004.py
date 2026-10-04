#!/usr/bin/env python3
"""Exact inputs for the all-length carry bound and a source-aware compiler.

The infinite proof uses the explicitly cited Matveev theorem; this script
checks its constants and an exact rational logarithm bridge, not that theorem.
All decisions use integers/Fractions, including when run with python -O.
"""
from fractions import Fraction as F
from math import gcd
import argparse
import json


def check(test, message="check failed"):
    if not test:
        raise AssertionError(message)


def compose(word):
    P = Q = 1
    B = 0
    for a in word:
        check(a >= 1)
        B = 3 * B + Q
        P *= 3
        Q <<= a
    return P, Q, B


def step(n):
    z = 3 * n + 1
    check(z != 0)
    a = (abs(z) & -abs(z)).bit_length() - 1
    return z // (1 << a), a


def replay(n, length):
    path, word = [n], []
    for _ in range(length):
        n, a = step(n)
        path.append(n)
        word.append(a)
    return path, tuple(word)


def log_interval(z, terms=100):
    """log(z), z>1: rational atanh series plus rigorous positive tail."""
    u = (z - 1) / (z + 1)
    check(0 < u < 1)
    power = u
    lower = F(0)
    for k in range(terms):
        lower += 2 * power / (2 * k + 1)
        power *= u * u
    tail = 2 * power / ((2 * terms + 1) * (1 - u * u))
    return lower, lower + tail


def carry_bound_inputs():
    rho = F(1121, 3328)
    maxima = []
    P = 1
    for j in range(1, 5001):
        P *= 3
        A = P.bit_length()
        Q = 1 << A
        B = 3 ** (j - 1) + (1 << (A - j + 1)) * (3 ** (j - 1) - (1 << (j - 1)))
        ratio = F(B, Q * (Q - P))
        check(ratio <= rho, f"finite carry bound j={j}")
        if j <= 64:
            maxima.append((ratio, j, A, B))
    maxima.sort(reverse=True)
    check(maxima[0] == (rho, 5, 8, 1121))
    check(maxima[1][0] < rho < F(1, 2))
    print("Sharp finite maximum j=1..64:", maxima[0])
    print("Independent extended integer check j=1..5000: PASS")

    p, q = 9115015689657667, 5750934602875680
    approx, error = F(p, q), F(1, q * q)
    l2, u2 = log_interval(F(2))
    l3, u3 = log_interval(F(3))
    alpha_lower, alpha_upper = l3 / u2, u3 / l2
    check(approx - error < alpha_lower < alpha_upper < approx + error)
    # Independent series: log(3)=log(2)+log(3/2), with parameter 1/5.
    l32, u32 = log_interval(F(3, 2), terms=80)
    alpha_lower2, alpha_upper2 = 1 + l32 / u2, 1 + u32 / l2
    check(approx - error < alpha_lower2 < alpha_upper2 < approx + error)
    check(gcd(p, q) == 1 and (1 << 43) < q < (1 << 53))
    check(F(1, 2) < l2 < u2 < 1 and u3 < 2)
    # At 65<=j<=2^42: Lambda>1/(4q)>2^-55, rho_j<2^(55-j).
    check(F(1, 4 * q) > F(1, 1 << 55))
    check(F(1, 1024) < rho)
    cutoff, constant = 1 << 42, 10 ** 10
    # 1.4*30^5*2^4.5*(1)*(2) < the following integer upper bound.
    check(2 * 30 ** 5 * 32 * 2 < constant)
    check(44 * constant < cutoff // 4 and 4 * constant < cutoff)
    print("Two exact log-series certificates: |log_2(3)-p/q|<1/q^2")
    print("p,q:", p, q, "; bridge 65..2^42; Matveev C=10^10 thereafter")
    return {"sharp_ratio": str(rho), "equality_word": [4, 1, 1, 1, 1],
            "log_approximant": [p, q], "matveev_C": constant, "cutoff": cutoff}


def compositions(total, length):
    if length == 1:
        yield (total,)
    else:
        for a in range(1, total - length + 2):
            for tail in compositions(total - a, length - 1):
                yield (a,) + tail


def independent_small_words():
    count = 0
    for j in range(1, 8):
        amin = (3 ** j).bit_length()
        for A in range(amin, amin + 3):
            Q = 1 << A
            maximum = -1
            maximizers = []
            for w in compositions(A, j):
                P, q, B = compose(w)
                check(q == Q and F(B, Q * (Q - P)) <= F(1121, 3328))
                if B > maximum:
                    maximum, maximizers = B, [w]
                elif B == maximum:
                    maximizers.append(w)
                r = (-B * pow(P, -1, Q)) % Q
                for k in (0, 1, 3):
                    n = r + k * Q
                    path, actual = replay(n, j)
                    check(actual[:-1] == w[:-1] and actual[-1] >= w[-1])
                    formal = (P * n + B) // Q
                    check(P * n + B == Q * formal)
                    check(path[-1] <= formal)
                    if k:
                        check(path[-1] < n, "tail did not descend")
                for k in (1, 3):
                    n = r - k * Q
                    path, actual = replay(n, j)
                    check(actual[:-1] == w[:-1] and actual[-1] >= w[-1])
                    check(n < path[-1] < 0, "negative contraction did not reduce magnitude")
                count += 1
            expected = (A - j + 1,) + (1,) * (j - 1)
            check(maximizers == [expected], "carry maximum not uniquely front-loaded")
    print("Independent exhaustive words: j<=7, three least contracting totals:", count)
    return count


def growth_words(bound):
    stack = [((), 0)]
    while stack:
        word, A = stack.pop()
        j = len(word) + 1
        if j > bound:
            continue
        maxA = (3 ** j).bit_length() - 1
        for a in range(1, maxA - A + 1):
            child = word + (a,)
            yield child
            stack.append((child, A + a))


def compile_growth(word, m=1):
    repeated = word * m
    P, Q, B = compose(repeated)
    check(P > Q)
    t = 1
    while (Q << t) <= P:
        t += 1
    modulus = Q << t
    coarse = repeated[:-1] + (repeated[-1] + t,)
    r = (-B * pow(P, -1, modulus)) % modulus
    path, actual = replay(r, len(repeated))
    check(actual[:-1] == coarse[:-1] and actual[-1] >= coarse[-1])
    boundary = "descent" if path[-1] < r else "root" if r == 1 else "unresolved"
    # Universal theorem: every positive lift k>=1 descends, whatever r does.
    for k in (1, 7, 29):
        n = r + k * modulus
        lifted, aw = replay(n, len(repeated))
        check(aw[:-1] == coarse[:-1] and aw[-1] >= coarse[-1])
        check(all(x > n for x in lifted[1:-1]))
        check(lifted[-1] < n)
    if boundary == "descent":
        check(all(x > r for x in path[1:-1]))
    return {"nominal_word": list(word), "repetitions": m, "extra_guard": t,
            "coarse_word": list(coarse), "residue": r, "modulus": modulus,
            "boundary_endpoint": path[-1], "boundary_status": boundary}


def atlas_controls():
    words = list(growth_words(8))
    classes, improved, previously_excluded, unresolved = set(), [], 0, []
    counts = {"descent": 0, "root": 0, "unresolved": 0}
    for w in words:
        P, Q, B = compose(w)
        a = F(B, P - Q)
        h = a.numerator
        for m in range(1, 7):
            row = compile_growth(w, m)
            counts[row["boundary_status"]] += 1
            if row["boundary_status"] == "unresolved":
                unresolved.append(row)
            classes.add((row["residue"], row["modulus"]))
            if Q ** m <= h:
                previously_excluded += 1
            else:
                oldt = 1
                while (1 << oldt) * (Q ** m - h) <= P ** m - h:
                    oldt += 1
                check(oldt >= row["extra_guard"])
                if oldt > row["extra_guard"]:
                    improved.append((oldt - row["extra_guard"], w, m, oldt, row))
    print("All marked prefix-growth words length<=8:", len(words), "; m=1..6")
    print("Compiled rows / unique coarse classes:", len(words) * 6, len(classes))
    print("Boundary statuses:", counts)
    print("Old Q^m>h restriction excludes:", previously_excluded, "rows")
    print("Less restrictive guard on old-admissible rows:", len(improved))
    if improved:
        best = max(improved, key=lambda row: row[0])
        print("Largest guard saving:", best[:4], "new t", best[4]["extra_guard"])
    seven = compile_growth((1, 1, 2, 2))
    check((seven["residue"], seven["modulus"], seven["boundary_endpoint"]) == (7, 128, 5))
    print("7 repair:", seven)
    print("Repeat constraint removed; no claim of new coverage against every older bank")
    return {"word_count": len(words), "rows": len(words) * 6,
            "unique_classes": len(classes), "boundary_statuses": counts,
            "old_excluded": previously_excluded, "old_guard_improved": len(improved),
            "seven": seven, "unresolved": unresolved}


def first_coefficient_census(depth=16):
    """One coarse exit cylinder per growing prefix; arbitrary final valuation.

    The final valuation is only bounded below. Therefore no finite cap on it
    is imposed or omitted, unlike an exact-word enumeration.
    """
    stack = [(0, 0, 0)]  # length, valuation sum, carry of growth prefix
    counts = [0] * (depth + 1)
    exceptions = []
    mass = F(0)
    while stack:
        length, A, B = stack.pop()
        j = length + 1
        if j > depth:
            continue
        P = 3 ** j
        exitA = P.bit_length()
        Q = 1 << exitA
        Bj = 3 * B + (1 << A)
        check(exitA - A >= 2)
        r = (-Bj * pow(P, -1, Q)) % Q
        counts[j] += 1
        path, actual = replay(r, j)
        check(P * r + Bj == (1 << sum(actual)) * path[-1])
        actualQ = 1
        for i, a in enumerate(actual, 1):
            actualQ <<= a
            check((actualQ > 3 ** i) == (i == j))
        check(path[-1] < r or r == 1)
        if r * (Q - P) <= Bj:
            exceptions.append({"source": r, "length": j,
                               "endpoint": path[-1], "word": list(actual)})
        for a in range(1, exitA - A):
            stack.append((j, A + a, Bj))
    for j in range(1, depth + 1):
        mass += F(counts[j], 1 << (3 ** j).bit_length())
    check(exceptions == [{"source": 1, "length": 1, "endpoint": 1, "word": [2]}])
    # An independent residue sweep for depths<=8 tests guards and completeness.
    M = 1 << (3 ** 8).bit_length()
    direct = [0] * 9
    for n in range(1, M, 2):
        x, Q = n, 1
        for j in range(1, 9):
            x, a = step(x)
            Q <<= a
            if Q > 3 ** j:
                direct[j] += 1
                break
    for j in range(1, 9):
        check(direct[j] == counts[j] * M // (1 << (3 ** j).bit_length()))
    print("First-coefficient exit cylinders per length 1..16:", counts[1:])
    print("Census total:", sum(counts), "; threshold exceptions:", exceptions)
    print("Exact density among all positive integers:", mass)
    print("Independent residue sweep modulo", M, "through depth8: PASS")
    return {"depth": depth, "counts": counts[1:], "density": str(mass),
            "exceptions": exceptions}


def boundary_switch_census():
    """Actual positive sources 3..9999, odd prefixes<=96, stopping at root.

    Failed long-word endpoint tests are retained as boundary obligations,
    then checked against the earlier source-relative descent certificate.
    No conclusion is drawn about unvisited sources or unbounded words.
    """
    def sibling_min(n):
        while n % 8 == 5:
            n = (n - 1) // 4
        return n

    bad, common_future_rows = [], []
    prefixes = 0
    for source in range(3, 10000, 2):
        n, P, Q, B = source, 1, 1, 0
        first_descent = None
        a0 = sibling_min(source)
        common = (0, a0) if a0 < source else None
        for j in range(1, 97):
            n, a = step(n)
            B = 3 * B + Q
            P *= 3
            Q <<= a
            prefixes += 1
            if n < source and first_descent is None:
                first_descent = (j, n)
            a0 = sibling_min(n)
            if common is None and a0 < source:
                check(step(n)[0] == step(a0)[0], "sibling common future lost")
                common = (j, a0)
            if Q > P and n >= source:
                check((-B * pow(P, -1, Q)) % Q == source)
                check(first_descent is not None and first_descent[0] < j)
                bad.append((source, j, n, *first_descent))
                check(common is not None and common[0] <= first_descent[0])
                check(replay(source, common[0] + 1)[0][-1] == step(common[1])[0])
                common_future_rows.append((source, j, *common, first_descent[0]))
            if n == 1:
                break
    check(bad and bad[0] == (165, 17, 167, 1, 31))
    print("Boundary-switch census: odd sources3..9999, <=96 odd steps, stop at1")
    print("Actual prefixes:", prefixes, "; failed contracting endpoint tests:", len(bad))
    print("Distinct sources:", len({x[0] for x in bad}), "; all rescued by earlier source-relative descent")
    improvements = [x for x in common_future_rows if x[2] < x[4]]
    print("Sibling common-future obligations appear earlier than actual descent:", len(improvements), "rows")
    return {"source_max": 9999, "odd_step_cap": 96, "prefixes": prefixes,
            "boundary_failures": len(bad), "distinct_sources": len({x[0] for x in bad}),
            "rows": bad, "common_future_rows": common_future_rows,
            "earlier_common_future_rows": len(improvements)}


def hostile_controls():
    w = (4, 1, 1, 1, 1, 2, 2, 1, 2, 1, 1, 2, 1, 1, 1, 2, 3)
    P, Q, B = compose(w)
    path, actual = replay(165, len(w))
    check(actual == w and P < Q and path[-1] == 167)
    check((-B * pow(P, -1, Q)) % Q == 165)
    check(path[1] == 31 and 165 * (Q - P) < B)
    for k in (1, 3):
        n = 165 + k * Q
        check(replay(n, len(w))[0][-1] < n)
    print("Generic contracting-word hostile: 165 ->167 in17 odd steps; first step31")
    print("Hostile pole:", F(B, Q - P), "; class modulus:", Q)
    # Two exact finite paths can certify a common future without orbit descent.
    p7, _ = replay(7, 4)
    p3, _ = replay(3, 1)
    check(p7[-1] == p3[-1] == 5)
    check(replay(5, 1)[0][-1] == 1)
    # Source must remain 27, not the smaller moving checkpoint 41.
    p27, _ = replay(27, 2)
    check(p27 == [27, 41, 31] and p27[-1] > p27[0])
    # Signed cycle tags are separate certified terminals, not one numeric order.
    terminal_cycles = {1: (2,), -1: (1,), -5: (1, 2), -17: (1, 1, 1, 2, 1, 1, 4)}
    for root, word in terminal_cycles.items():
        path, aw = replay(root, len(word))
        check(path[-1] == root and aw == word)
    # The same compiler interface must allow a nontrivial positive terminal
    # for 5n+1: a small boundary is not automatically the designated root.
    n, p5, q5, b5 = 13, 1, 1, 0
    for a in (1, 1, 5):
        z = 5 * n + 1
        check((z & -z).bit_length() - 1 == a)
        n = z >> a
        b5 = 5 * b5 + q5
        p5 *= 5
        q5 <<= a
    check(n == 13 and (p5, q5, b5) == (125, 128, 39))
    check(F(b5, q5 - p5) == 13)
    print("5n+1 hostile: 13 ->33 ->83 ->13, word115, positive pole13")
    return {"source": 165, "endpoint": 167, "word": list(w),
            "modulus": Q, "pole": str(F(B, Q - P)), "smaller_first_step": 31}


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("--json", help="optional exact certificate summary")
    args = parser.parse_args()
    result = {"carry_bound": carry_bound_inputs()}
    result["small_word_controls"] = independent_small_words()
    result["atlas"] = atlas_controls()
    result["census"] = first_coefficient_census()
    result["hostile"] = hostile_controls()
    result["boundary_switch"] = boundary_switch_census()
    print("PASS: exact finite checks; all-length argument is in the companion proof")
    if args.json:
        with open(args.json, "w") as handle:
            json.dump(result, handle, indent=2)
            handle.write("\n")


if __name__ == "__main__":
    main()
