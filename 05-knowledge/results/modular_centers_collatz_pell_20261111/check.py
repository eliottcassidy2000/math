"""Exact finite controls for the accompanying all-integer proofs; stdlib only.

Run from anywhere: python check.py (or python -O check.py).
No assertion is disabled by -O. Output is deterministic, with explicit bounds.
"""
from collections import Counter
from fractions import Fraction
from itertools import product
from math import gcd, isqrt, prod
from pathlib import Path
import hashlib
import json
import tournament_checks

HERE = Path(__file__).resolve().parent


def require(ok, label, witness=None):
    if not ok:
        raise ValueError((label, witness))


def factors(n):
    out = {}
    p = 2
    while p * p <= n:
        while n % p == 0:
            out[p] = out.get(p, 0) + 1
            n //= p
        p += 1
    if n > 1:
        out[n] = out.get(n, 0) + 1
    return out


def phi(n):
    return prod((p - 1) * p ** (e - 1) for p, e in factors(n).items())


def divisors(n):
    return [d for d in range(1, n + 1) if n % d == 0]


def raw(x, sign):
    return 3 * x + sign if x % 2 else x // 2


def shortcut(x, sign):
    return (3 * x + sign) // 2 if x % 2 else x // 2


def modulus_step(m, sign):
    return 3 * m + 2 * sign if m % 8 == (4 + sign) % 8 else (m + sign) // 2


def central(m):
    k = m // 2
    indices = [k, k + 1] if m % 2 else [k - 1, k, k + 1]
    return [[a * b % m for b in indices] for a in indices]


def square_parents(n):
    # Literal graph definition, independent of the floor formula.
    return [v for v in range(1, n) if isqrt(v + n) ** 2 == v + n]


def square_degree(n):
    return isqrt(2 * n - 1) - isqrt(n)


def word_affine(word):
    a, b, denominator = 1, 0, 1
    for bit in word:
        if bit:
            a, b = 3 * a, 3 * b + denominator
        denominator *= 2
    return a, b, denominator


def literal_word(x, sign, length):
    word = []
    for _ in range(length):
        word.append(x % 2)
        x = shortcut(x, sign)
    return tuple(word), x


def check_centers_and_edges():
    for m in range(3, 10002, 2):
        x, sign = ((m - 1) // 4, 1) if m % 4 == 1 else ((m + 1) // 4, -1)
        b = 3 * x + sign
        expected = [[b, x], [x, b]] if sign == 1 else [[x, b], [b, x]]
        require(central(m) == expected, "odd center", m)
        require(all(sum(row) == m for row in expected), "row sums", m)
        require(all(gcd(z, m) == 1 for row in expected for z in row), "central units", m)
    for m in range(4, 10001, 2):
        k = m // 2
        d = k if k % 2 else 0
        e = k - d
        expected = [[d + 1, e, d - 1], [e, d, e], [d - 1, e, d + 1]]
        require(central(m) == [[z % m for z in row] for row in expected], "even center", m)
    for sign in (-1, 1):
        for x in range(1, 100001):
            m = 4 * x + sign
            b = 3 * x + sign
            decoded = b if x % 2 else (b - sign) // 6
            require(decoded == raw(x, sign), "center edge", (x, sign))
            require(modulus_step(m, sign) == 4 * raw(x, sign) + sign, "conjugacy", (x, sign))
        # Build incoming edges independently from all possible raw sources.
        parents = {x: [] for x in range(1, 10001)}
        for y in range(1, 20001):
            target = raw(y, sign)
            if target in parents:
                parents[target].append(4 * y + sign)
        gate = 17 if sign == 1 else 7
        for x, literal in parents.items():
            m = 4 * x + sign
            predicted = [2 * m - sign]
            if m % 24 == gate:
                predicted.append((m - 2 * sign) // 3)
            require(sorted(literal) == sorted(predicted), "inverse gate", (m, sign))
    return {"odd_modulus_max": 10001, "even_modulus_max": 10000,
            "raw_inputs_per_sign": 100000, "inverse_targets_per_sign": 10000,
            "24_25_26_centers": {str(m): central(m) for m in (24, 25, 26)}}


def check_words():
    total = 0
    examples = []
    for length in range(1, 13):
        q = 2 ** length
        words = list(product((0, 1), repeat=length))
        for sign in (-1, 1):
            # First enumerate literal trajectories on a complete positive residue system.
            observed = {}
            for x in range(1, q + 1):
                word, endpoint = literal_word(x, sign, length)
                require(word not in observed, "unique word", (length, sign, x))
                observed[word] = (x % q, endpoint)
            for word in words:
                a, b, denominator = word_affine(word)
                residue = (-sign * b * pow(a, -1, denominator)) % denominator
                require(observed[word][0] == residue, "word residue", (word, sign))
                m_residue = (4 * residue + sign) % (4 * denominator)
                require(gcd(m_residue, 4 * denominator) == 1, "Dirichlet coprimality", word)
                k = 4 * b + denominator - a
                # Large lifts ensure this is not just a representative coincidence.
                for lift in (1, 7, 100):
                    x = residue + lift * denominator
                    actual_word, endpoint = literal_word(x, sign, length)
                    require(actual_word == word, "lifted guard", (word, sign, lift))
                    require(a * x + sign * b == denominator * endpoint, "affine endpoint", word)
                    m = 4 * x + sign
                    require(a * m + sign * k == denominator * (4 * endpoint + sign), "modulus endpoint", word)
                total += 1
                if length <= 3:
                    examples.append({"sign": sign, "word": "".join(map(str, word)),
                                     "x_residue": residue, "x_period": denominator,
                                     "m_residue": m_residue, "m_period": 4 * denominator,
                                     "A": a, "B": b, "K": k})
    return {"maximum_length": 12, "signed_word_guards_checked": total,
            "positive_lifts_per_guard": [1, 7, 100], "small_guards": examples}


def check_tables():
    cells = 0
    for m in range(2, 202):
        output_counts = Counter()
        row_collisions = 0
        for a in range(m):
            counts = Counter(a * b % m for b in range(m))
            g = gcd(a, m)
            require(len(counts) == m // g and set(counts.values()) == {g}, "row fibres", (m, a))
            output_counts.update(counts)
            row_collisions += sum(c * c for c in counts.values())
        ds = divisors(m)
        for t in range(m):
            formula = sum(d * phi(m // d) for d in ds if t % d == 0)
            require(output_counts[t] == formula, "product fibres", (m, t))
        z = sum(d * phi(m // d) for d in ds)
        require(row_collisions == m * z, "row collisions", m)
        for p, e in factors(m).items():
            observed = Fraction(sum(sum(gcd(a, m) % (p ** j) == 0 for j in range(1, e + 1))
                                    for a in range(m)), m)
            expected = sum((Fraction(1, p ** j) for j in range(1, e + 1)), Fraction(0))
            require(observed == expected, "entropy coefficient", (m, p))
        cells += m * m
    rows = []
    for m in range(3, 1002, 2):
        fs = factors(m)
        prime = len(fs) == 1 and next(iter(fs.values())) == 1
        # Scan actual nonzero row indices outward from the middle.
        first = None
        for radius in range(1, (m - 1) // 2 + 1):
            u = 2 * radius - 1
            if gcd((m + u) // 2, m) > 1 or gcd((m - u) // 2, m) > 1:
                first = radius
                break
        require(first == (None if prime else (min(fs) + 1) // 2), "first singular window", m)
        z = sum(gcd(a, m) for a in range(m))
        rows.append({"m": m, "unit_rows": phi(m), "row_types": len(divisors(m)),
                     "zero_product_pairs": z, "row_collision_triples": m * z,
                     "first_singular_radius": first})
    # CRT preserves the actual multiplication, not just a row histogram.
    for a_mod, b_mod in ((3, 5), (5, 7), (9, 25), (8, 9)):
        m = a_mod * b_mod
        images = {(z % a_mod, z % b_mod) for z in range(m)}
        require(len(images) == m, "CRT bijection", m)
        for a in range(m):
            for b in range(m):
                z = a * b % m
                require((z % a_mod, z % b_mod) == (a * b % a_mod, a * b % b_mod), "CRT product", m)
    require(Counter(a * 5 % 25 for a in range(25)) == Counter({0: 5, 5: 5, 10: 5, 15: 5, 20: 5}), "25 hostile")
    return {"full_moduli": [2, 201], "literal_cells": cells, "odd_profile_max": 1001,
            "CRT_pairs": [[3, 5], [5, 7], [9, 25], [8, 9]], "profiles": rows}


def check_pell():
    families = []
    for epsilon in (1, 2):
        r, s = (5, 7) if epsilon == 1 else (3, 4)
        records = []
        for j in range(1, 101):
            r2 = r * r
            require(s * s == 2 * r2 - epsilon, "Pell norm", (epsilon, j))
            require(r2 % 8 == 1, "odd square", r2)
            require(len({square_degree(r2 - 1), square_degree(r2), square_degree(r2 + 1)}) == 1, "plateau", r2)
            h = (r - 1) // 2
            x = h * (h + 1)
            require(r2 == 4 * x + 1 and x % 2 == 0, "pronic source", r2)
            incoming = 2 * r2 - 1
            require(modulus_step(incoming, 1) == r2, "incoming table", r2)
            require(r2 % 24 in (1, 9), "unique incoming gate", r2)
            require(modulus_step(r2, 1) == (r2 + 1) // 2, "outgoing table", r2)
            if epsilon == 1:
                require(incoming == s * s, "square incoming", r2)
            else:
                require(incoming == s * s + 1, "quadratic incoming", r2)
            if r2 <= 1000000:
                middle = list(range(r + 1, s))
                predicted = [[1] + [q * q - r2 + 1 for q in middle],
                             [q * q - r2 for q in middle] + [r2 - epsilon],
                             [q * q - r2 - 1 for q in middle] + [r2 - epsilon - 1]]
                for n, portset in zip((r2 - 1, r2, r2 + 1), predicted):
                    require(square_parents(n) == portset, "literal port template", n)
            if j <= 8:
                records.append({"j": j, "center": r2, "r": r, "s": s,
                                "epsilon": epsilon, "incoming": incoming,
                                "outgoing": (r2 + 1) // 2, "birth_degree": square_degree(r2)})
            r, s = 3 * r + 2 * s, 4 * r + 3 * s
        families.extend(records)
    # Independent scalar recurrence for the epsilon=2 incoming moduli.
    m_values = [1, 17]
    for j in range(2, 129):
        m_values.append(34 * m_values[-1] - m_values[-2])
    r, t = 1, 0
    for j in range(129):
        require(m_values[j] == 2 * r * r - 1, "Pell trace bridge", j)
        for d in range(1, j + 1):
            if j % d == 0 and (j // d) % 2 == 1:
                require(m_values[j] % m_values[d] == 0, "odd index divisibility", (j, d))
        r, t = 3 * r + 4 * t, 2 * r + 3 * t
    dyadic = [m_values[2 ** a] for a in range(8)]
    for i, x in enumerate(dyadic):
        if i:
            require(x == 2 * dyadic[i - 1] ** 2 - 1, "dyadic recurrence", i)
        for y in dyadic[:i]:
            require(gcd(x, y) == 1, "coprime dyadic terms", (x, y))
    certificates = {19601: [17, 1153], 886731088897: [257, 1409, 2448769]}
    for n, primes in certificates.items():
        require(prod(primes) == n and all(factors(p) == {p: 1} for p in primes), "factor certificate", n)
    require(m_values[8] == 886731088897, "prime-index hostile")
    return {"terms_per_Pell_family": 100, "literal_parent_sets_center_max": 1000000,
            "24_25_26_parents": {str(n): square_parents(n) for n in (24, 25, 26)},
            "sample_plateaus": sorted(families, key=lambda z: z["center"]),
            "trace_divisibility_index_max": 128,
            "first_incoming_traces": m_values[:9], "dyadic_terms": list(map(str, dyadic)),
            "factor_certificates": certificates}


def main():
    results = {"status": "FINITE-EXACT; universal claims use PROOFS.md",
               "centers_and_edges": check_centers_and_edges(), "parity_words": check_words(),
               "multiplication_tables": check_tables(), "Pell": check_pell(),
               "tournaments": tournament_checks.run()}
    output = HERE / "verification.json"
    output.write_text(json.dumps(results, indent=2, sort_keys=True) + "\n", encoding="utf-8", newline="\n")
    manifest = {p.name: hashlib.sha256(p.read_text(encoding="utf-8").replace("\r\n", "\n").encode()).hexdigest()
                for p in (HERE / "check.py", HERE / "tournament_checks.py", output)}
    (HERE / "manifest.json").write_text(json.dumps(manifest, indent=2, sort_keys=True) + "\n", encoding="utf-8", newline="\n")
    print(json.dumps({"status": "PASS", "signed_guards": results["parity_words"]["signed_word_guards_checked"],
                      "literal_table_cells": results["multiplication_tables"]["literal_cells"],
                      "Pell_terms": 200, "output": output.name, "hashes": manifest}, indent=2))


if __name__ == "__main__":
    main()
