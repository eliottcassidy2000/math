#!/usr/bin/env python3
"""Exact audits for signed repairs, section returns, six-letter guards and height.

Standard library only. The Graver theorem is proved/cited in the note; this
program does not claim to compute unrestricted Graver bases or prove Collatz.
"""
from collections import Counter
from fractions import Fraction
from itertools import product
from math import comb
import json

CHECKS = 0


def check(condition, message):
    global CHECKS
    CHECKS += 1
    if not condition:
        raise ValueError(message)


def v2(n):
    n = abs(n)
    return (n & -n).bit_length()-1


def U(n):
    return (3*n+1)//2**v2(3*n+1)


def path(n, limit=10000):
    result, seen = [], set()
    while n != 1:
        check(n not in seen and len(result) < limit, "finite root control")
        seen.add(n)
        result.append(n)
        n = U(n)
    return result


def tidy(x):
    return {k: v for k, v in sorted(x.items()) if v}


def add(*vectors):
    out = Counter()
    for vector in vectors:
        for k, v in vector.items():
            out[k] += v
    return tidy(out)


def scale(vector, scalar):
    return tidy({k: scalar*v for k, v in vector.items()})


def boundary(c):
    out = Counter()
    for n, a in c.items():
        out[n] += a
        out[U(n)] -= a
    return tidy(out)


def value(c):
    x = Fraction(1)
    for n, a in c.items():
        x *= Fraction(n, U(n))**a
    return x


def D(n, c):
    return add(boundary(c), {n: -1}, {1: 1})


def norm(v):
    return sum(abs(x) for x in v.values())


def factors(n):
    out = Counter()
    p = 2
    while p*p <= n:
        while n % p == 0:
            out[p] += 1
            n //= p
        p += 1
    if n > 1:
        out[n] += 1
    return dict(out)


def prime_image(v):
    return add(*(scale(factors(n), a) for n, a in v.items()))


def K(n):
    f = factors(n)
    return add({n: 1}, scale(f, -1), {1: sum(f.values())-1})


def push(v):
    return add(*({U(n): a} for n, a in v.items()))


def conformal(x, y):
    return all(a*b >= 0 and abs(a) <= abs(b)
               for k in set(x) | set(y) for a, b in [(x.get(k, 0), y.get(k, 0))])


def section(v):
    check(v > 0 and v % 2 and v % 3, "section domain")
    a = next(a for a in range(1, 7) if (2**a*v-1) % 9 == 0)
    return (2**a*v-1)//3, a


def word_data(word):
    p, q, b = 1, 1, 0
    for a in word:
        p, q, b = 3*p, q*2**a, 3*b+q
    return p, q, b


def execute(n, word):
    out = [n]
    for a in word:
        check(v2(3*n+1) == a, "exact word valuation")
        n = U(n)
        out.append(n)
    return out


# r, binary denominator cost, carry, native source residue, source modulus.
LETTERS = {
    "H": (6, 10, 669, 155, 2048), "G": (2, 3, 5, 11, 16),
    "A": (4, 7, 85, 187, 256), "B": (4, 7, 73, 7, 256),
    "L": (2, 4, -3, 219, 256), "J": (5, 8, 147, 799, 1024),
}


def controller(word):
    p, q, b = 1, 1, 0
    for letter in word:
        r, a, carry, _, _ = LETTERS[letter]
        p, q, b = 3**r*p, 2**a*q, 3**r*b+carry*q
    return p, q, b


def independent_guard(word):
    """Intersect each primitive's source guard, retaining its actual prefix."""
    residue, modulus, p, q, b = 0, 1, 1, 1, 0
    for letter in word:
        r, a, carry, guard, mod = LETTERS[letter]
        required_mod = q*mod
        required = ((guard*q-b)*pow(p, -1, required_mod)) % required_mod
        small = min(modulus, required_mod)
        if (residue-required) % small:
            return None
        if required_mod > modulus:
            residue, modulus = required, required_mod
        p, q, b = 3**r*p, 2**a*q, 3**r*b+carry*q
    return residue, modulus


def decoder(beta):
    reverse = []
    while beta:
        check(beta.denominator & (beta.denominator-1) == 0, "dyadic decoder")
        cost = beta.denominator.bit_length()-1
        tag = beta.numerator*pow(beta.denominator, -1, 9) % 9
        if tag == 3:
            t = beta.numerator*pow(beta.denominator, -1, 243) % 243
            possible = [l for l in "HJ" if
                LETTERS[l][2]*pow(2**LETTERS[l][1], -1, 243) % 243 == t]
            check(len(possible) == 1, "five-digit H/J tag")
            letter = possible[0]
        else:
            letter = {4: "G", 2: "A", 5: "B", 6: "L"}.get(tag)
            check(letter is not None, "native tag")
        r, a, b, _, _ = LETTERS[letter]
        check(cost >= a, "decoder cost decreases")
        beta = (2**a*beta-b)/3**r
        reverse.append(letter)
    return ''.join(reversed(reverse))


def safe_moves():
    c0 = {3: 2, 5: 2}
    target = Counter(path(9))
    z = add(target, scale(c0, -1))
    check(value(z) == 1, "signed neutral repair")
    check(add(c0, z) == tidy(target), "nonnegative repair endpoint")
    check(norm(D(9, c0)) == 4 and not D(9, target), "9 defect repaired")
    # Exhaust every conformal subvector: no nonzero proper prime-balanced one.
    keys = list(z)
    proper = 0
    for amounts in product(*(range(abs(z[u])+1) for u in keys)):
        small = {u: (1 if z[u] > 0 else -1)*a for u, a in zip(keys, amounts) if a}
        if small and small != z and value(small) == 1:
            proper += 1
    check(proper == 0, "primitive signed move in this seven-edge pool")
    # Every fibre of 3^k in this pool is parameterized by this same move.
    fibre_rows = []
    for k in range(1, 13):
        states = []
        for t in range(k//2+1):
            c = add({3: k, 5: k}, scale(z, t))
            check(all(a >= 0 for a in c.values()), "fibre stock")
            check(value(c) == 3**k, "fibre prime value")
            states.append(norm(D(3**k, c)))
        check(all(b < a for a, b in zip(states, states[1:])), "finite-pool descent")
        fibre_rows.append({"k": k, "defect_norms": states})
    neutral = {7: 2, 11: 1, 17: 1, 55: 1, 65: 1, 83: 1, 5: 1}
    check(value(neutral) == 1, "neutral stock witness")
    check(all(n % 3 for n in neutral), "3-source obstruction")
    requests = [{7: 5, 55: 2}, {11: 4, 83: 3}, {5: 1, 65: 2}]
    buffer = scale(neutral, 4)
    check(all(all(buffer.get(u, 0) >= a for u, a in req.items()) for req in requests),
          "one neutral context for finite unit-stock requests")
    root5 = {5: 1}
    c5, d5 = add(root5, neutral), add(root5, {1: 1})
    common = add(root5, neutral, {1: 1})
    check(value(c5) == value(d5) == value(common) == 5, "common extension value")
    check(value(add(common, scale(c5, -1))) == value(add(common, scale(d5, -1))) == 1,
          "both extension differences neutral")
    check({u: a for u, a in c0.items() if u % 3 == 0}
          != {u: a for u, a in target.items() if u % 3 == 0},
          "same-value9 receipts have incompatible neutral profiles")
    return {"nine_move": z, "nine_root_path": list(target),
            "proper_conformal_subrelations": proper, "power3_fibres": fibre_rows,
            "common_neutral_context_multiple": 4}


def section_returns():
    table = {r: section(r)[1] for r in (1, 5, 7, 11, 13, 17)}
    for v in range(1, 10002, 2):
        if v % 3 == 0:
            continue
        n, a = section(v)
        check(n % 3 == 0 and U(n) == v and v2(3*n+1) == a, "canonical predecessor")
        for t in range(4):
            large = (2**(a+6*t)*v-1)//3
            check(section(U(large))[0] == n and large >= n, "section normalization")
    levels = Counter()
    samples = []
    # First section-contracting exit among words with no forward contraction.
    def walk(word, p, q, b):
        depth = len(word)
        if depth >= 2 and 2*p < 3*q and b*pow(q, -1, 9) % 9 == 5:
            residue = (q-b)*pow(p, -1, 2*q) % (2*q)
            d = Fraction(2*(p*residue+b)-q, 3*q)
            check(d.denominator == 1 and 0 < d < residue, "whole-cell payment at least source")
            check(d % 3 == 0, "section return is a multiple of3")
            levels[depth] += 1
            if len(samples) < 12:
                samples.append({"word": list(word), "source": residue,
                                "modulus": 2*q, "child": int(d)})
            return Fraction(1, q)
        if depth == 12:
            return Fraction()
        total = Fraction()
        # 2^(new A) < 3^(depth+1), hence exact integer threshold.
        for a in range(1, (3*p).bit_length()-(q.bit_length()-1)):
            total += walk(word+(a,), 3*p, q*2**a, 3*b+q)
        return total
    mass = walk((), 1, 1, 0)
    check(min(levels) == 6 and levels[6] == 1, "first non-forward paid section word")
    w = (1, 1, 1, 1, 2, 3)
    for t in range(1000):
        n = 799+1024*t
        route = execute(n, w)
        child = (243*n+147)//256
        check(all(x > n for x in route[1:]), "no actual descent in six steps")
        check(0 < child < n and child % 3 == 0 and U(child) == route[-1], "new join")
        check(n % 16 == 15 and all(n % mod != g for _, _, _, g, mod in list(LETTERS.values())[:5]),
              "outside old five native first letters")
        check(Fraction(child+5, n+5) <= Fraction(191, 201), "sharp energy bound")
    check(9*191**3 < 8*201**3, "one third growth-credit payment")
    return {"residue_mod9_to_exponent": {str(v % 9): a for v, a in table.items()},
            "first_exit_counts_by_depth": dict(sorted(levels.items())),
            "first_exit_odd_density_through12": str(mass), "samples": samples,
            "J": {"guard": [799, 1024], "child": [243, 147, 256],
                  "left_word": list(w), "right_word": [1],
                  "sharp_energy_ratio": "191/201", "third_credit_ratio": str(Fraction(9*191**3, 8*201**3))}}


def controllers():
    total = valid = forbidden = 0
    for length in range(1, 6):
        for letters in product(LETTERS, repeat=length):
            w = ''.join(letters)
            p, q, b = controller(w)
            check(decoder(Fraction(b, q)) == w, "six-letter exact decoding")
            g = independent_guard(w)
            expected = not ('LB' in w or 'LJ' in w)
            check((g is not None) == expected, "complete forbidden-pair language")
            if not expected:
                forbidden += 1
                continue
            s, t = (16, 11) if w[-1] == 'L' else ((4, 3) if w[-1] == 'J' else (2, 1))
            m = s*q
            residue = (t*q-b)*pow(p, -1, m) % m
            check(g == (residue, m), "compositional terminal guard")
            for k in (0, 1, 3):
                x = residue+k*m
                if x == 0:
                    x += m
                source = x
                for l in w:
                    r, a, carry, native, mod = LETTERS[l]
                    check(x % mod == native, "literal primitive source")
                    check((3**r*x+carry) % 2**a == 0, "literal primitive integrality")
                    x = (3**r*x+carry)//2**a
                    check(x > 0 and x % 2, "positive primitive output")
                check(q*x == p*source+b, "compiled exact endpoint")
            valid += 1
            total += 1
    repeat_examples = []
    for q in range(1, 7):
        g = independent_guard('J'*q)
        for t in (0, 1, 2):
            n = g[0]+t*g[1]
            predicted = max(0, (v2(13*n-147)-2)//8)
            x, copies = n, 0
            while x % 1024 == 799:
                x = (243*x+147)//256
                copies += 1
            check(copies == predicted and copies >= q, "exact J repeat fuel")
        repeat_examples.append({"q": q, "least_source": g[0], "modulus": g[1]})
    return {"formal_words": valid+forbidden, "native_words": valid,
            "forbidden_words": forbidden, "max_word_length": 5,
            "HJ_translation_difference": str(Fraction(669, 1024)-Fraction(147, 256)),
            "JJJG_guard": independent_guard('JJJG'), "J_repeat_examples": repeat_examples}


def transport():
    count = 0
    for n in range(1, 1002, 2):
        lhs = push(K(n))
        charge = add(factors(U(n)), *(scale(factors(U(p)), -e) for p, e in factors(n).items()))
        rhs = add(K(U(n)), *(scale(K(U(p)), -e) for p, e in factors(n).items()),
                  charge, {1: -sum(charge.values())})
        check(lhs == rhs, "prime-plus-composite coupled transport")
        count += 1
    check(push(K(9)) == {1: 1, 5: -2, 7: 1}, "K9 leaves factorization kernel")
    check(prime_image(push(K(9))) == {5: -2, 7: 1}, "lost charge witness")
    neutral = {7: 2, 11: 1, 17: 1, 55: 1, 65: 1, 83: 1, 5: 1}
    for t in range(32):
        n = 799+1024*t
        child = (243*n+147)//256
        left = Counter(execute(n, (1, 1, 1, 1, 2, 3))[:-1])
        right = {child: 1}
        child_receipt = add(Counter(path(child)), neutral)
        source_receipt = add(child_receipt, left, scale(right, -1))
        check(all(v >= 0 for v in source_receipt.values()), "replacement retains stock")
        check(value(source_receipt) == n and value(child_receipt) == child, "replacement values")
        check(D(n, source_receipt) == D(child, child_receipt) == boundary(neutral),
              "common-future replacement preserves the full nonzero defect")
        check(add(source_receipt, right, scale(left, -1)) == child_receipt,
              "inverse replacement with retained prefix")
    return {"cases": count, "K9_after_U": push(K(9)), "prime_charge": prime_image(push(K(9)))}


def portfolios():
    seeds = [n for n in range(1, 100, 2) if n % 3][:16]
    for s in seeds:
        path(s)
    rows = []
    for h in range(1, 7):
        modulus = 3**h
        best = [None]*modulus
        for i, s in enumerate(seeds, 1):
            for a in range(1, 2*modulus+1):
                if (2**a*s-1) % 3:
                    continue
                n = (2**a*s-1)//3
                r = n % modulus
                if best[r] is None or n < best[r]:
                    best[r] = n
            check(all(n is not None for n in best), "one-step portfolio local coverage")
            if i in (1, 2, 4, 8, 16):
                worst = max(best)
                ell = worst.bit_length()-1
                check(i*(ell+2) >= modulus, "depth-one universal height bound")
                rows.append({"h": h, "seeds": i, "largest_least_source_bits": worst.bit_length(),
                             "worst_residue": best.index(worst)})
    # Independent direct census versus the bound for finite seed sets/depths.
    bounded = []
    for X in (63, 255, 1023, 4095):
        ell = X.bit_length()-1
        for bank in ([1], seeds[:4]):
            for depth in range(1, 6):
                count = 0
                for n in range(1, X+1, 2):
                    x, cost, hit = n, 0, n in bank
                    for d in range(1, depth+1):
                        cost += v2(3*x+1)
                        x = U(x)
                        check(cost <= ell+2*d, "ordinary-height halving-cost bound")
                        hit = hit or x in bank
                    count += hit
                bound = len(bank)*sum(comb(ell+2*d, d) for d in range(depth+1))
                check(count <= bound, "finite-depth source capacity")
                bounded.append({"X": X, "seeds": len(bank), "depth": depth,
                                "sources": count, "bound": bound})
    return {"seed_bank": seeds, "local_portfolio_rows": rows, "source_capacity_rows": bounded}


def main():
    result = {"signed_repairs": safe_moves(), "section_returns": section_returns(),
              "controllers": controllers(), "coupled_transport": transport(),
              "height_portfolios": portfolios()}
    result["checks"] = CHECKS
    result["status"] = "PASS; exact finite audits; universal Collatz entry remains OPEN"
    print(json.dumps(result, indent=2, sort_keys=True))


if __name__ == '__main__':
    main()
