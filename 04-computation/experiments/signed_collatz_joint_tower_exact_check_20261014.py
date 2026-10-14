#!/usr/bin/env python3
"""Exact signed Collatz joint laws modulo 3^j; standard library only."""
import hashlib
import json
from fractions import Fraction

checks = 0


def require(ok, message):
    global checks
    checks += 1
    if not ok:
        raise ArithmeticError(message)


def valuation(n):
    return (n & -n).bit_length() - 1


def odd_step(n, sign):
    x = 3 * n + sign
    return x // (1 << valuation(x))


def bareiss(a):
    a = [row[:] for row in a]
    n, previous, sign = len(a), 1, 1
    if n == 1:
        return a[0][0]
    for k in range(n - 1):
        if not a[k][k]:
            p = next((i for i in range(k + 1, n) if a[i][k]), None)
            if p is None:
                return 0
            a[k], a[p] = a[p], a[k]
            sign = -sign
        pivot = a[k][k]
        for i in range(k + 1, n):
            for j in range(k + 1, n):
                value = a[i][j] * pivot - a[i][k] * a[k][j]
                require(value % previous == 0, "Bareiss division exact")
                a[i][j] = value // previous
            a[i][k] = 0
        previous = pivot
    return sign * a[-1][-1]


def rank(a):
    a = [[Fraction(x) for x in row] for row in a]
    r = 0
    for c in range(len(a[0])):
        p = next((i for i in range(r, len(a)) if a[i][c]), None)
        if p is None:
            continue
        a[r], a[p] = a[p], a[r]
        pivot = a[r][c]
        a[r] = [x / pivot for x in a[r]]
        for i in range(r + 1, len(a)):
            if a[i][c]:
                scale = a[i][c]
                a[i] = [x - scale * y for x, y in zip(a[i], a[r])]
        r += 1
        if r == len(a):
            break
    return r


def conditional_matrix(row, h, source_class):
    q, d, t = row["modulus"], row["matrix_dimension"], row["modulus"] // 3
    p, strip = 3 ** h, 3 ** (h + 1)
    units = row["units_order"]
    ix = {u: i for i, u in enumerate(units)}
    M, crt = [[0] * d for _ in range(d)], [[0] * d for _ in range(d)]
    ustrip = (3 * source_class + 1) * pow(2, -1, strip) % strip
    vstrip = (3 * source_class - 1) * pow(2, -1, strip) % strip
    for k in range(2, d + 2):
        factor, weight = pow(2, k - 1, q), p * (1 << (d + 1 - k))
        for u in units:
            if u % 3 == ((-1) ** k) % 3:
                v = (factor * u - 1) % q
                if v % strip == vstrip:
                    M[ix[u]][ix[v]] += weight
        for v in units:
            if v % 3 == -((-1) ** k) % 3:
                u = (factor * v + 1) % q
                if u % strip == ustrip:
                    M[ix[u]][ix[v]] += weight
        for sign in (-1, 1):
            period = 1 << (k + 1)
            root = (((1 << k) - sign) * pow(3, -1, period)) % period
            base = root + period * (((source_class - root) * pow(period, -1, p)) % p)
            for z in range(t // p):
                n = base + period * p * z
                u, v = odd_step(n, 1) % q, odd_step(n, -1) % q
                require(n % p == source_class, "conditioned CRT source class")
                require(valuation(3 * n + sign) == k, "conditioned CRT valuation")
                crt[ix[u]][ix[v]] += weight
    require(M == crt, "conditioned line and independent CRT matrices agree")
    require(sum(map(sum, M)) == row["denominator"], "conditioned total mass")
    require(rank(M) == 2 * t // p, "conditioned exact rank")
    return M


def construct(level):
    q, t = 3 ** level, 3 ** (level - 1)
    d, D = 2 * t, 2 * t * (4 ** t - 1)
    units = [u for u in range(q) if u % 3]
    ix = {u: i for i, u in enumerate(units)}
    M = [[0] * d for _ in range(d)]
    crt = [[0] * d for _ in range(d)]
    require(pow(2, d, q) == 1, "order divides 2t")
    require(all(pow(2, r, q) != 1 for r in range(1, d)), "exact order 2t")

    for k in range(2, d + 2):
        h, weight = pow(2, k - 1, q), 1 << (d + 1 - k)
        for u in units:
            if u % 3 == ((-1) ** k) % 3:
                v = (h * u - 1) % q
                require(v % 3 == 1, "minus low valuation residue")
                M[ix[u]][ix[v]] += weight
        for v in units:
            if v % 3 == -((-1) ** k) % 3:
                u = (h * v + 1) % q
                require(u % 3 == 2, "plus low valuation residue")
                M[ix[u]][ix[v]] += weight

        # Independent construction from exact source cylinders, not output lines.
        for sign in (-1, 1):
            period = 1 << (k + 1)
            root = (((1 << k) - sign) * pow(3, -1, period)) % period
            pairs = []
            for z in range(t):
                n = root + z * period
                require(valuation(3 * n + sign) == k, "source cylinder high valuation")
                require(valuation(3 * n - sign) == 1, "source cylinder other valuation")
                u, v = odd_step(n, 1) % q, odd_step(n, -1) % q
                pairs.append((u, v))
                crt[ix[u]][ix[v]] += weight
            # Check that shifting k by its claimed period preserves the pair set.
            k2 = k + d
            period2 = 1 << (k2 + 1)
            root2 = (((1 << k2) - sign) * pow(3, -1, period2)) % period2
            pairs2 = [(odd_step(root2 + z * period2, 1) % q,
                       odd_step(root2 + z * period2, -1) % q) for z in range(t)]
            require(sorted(pairs) == sorted(pairs2), "valuation-period fibre identity")

    require(M == crt, "affine-line compiler equals independent CRT compiler")
    require(sum(map(sum, M)) == D, "total mass one")
    for i, u in enumerate(units):
        require(3 ** level * sum(M[i]) == D * (1 if u % 3 == 1 else 2), "plus marginal")
        for j, v in enumerate(units):
            require((M[i][j] == 0) == (u % 3 == 1 and v % 3 == 2), "exact support")
    for j, v in enumerate(units):
        require(q * sum(M[i][j] for i in range(d)) == D * (2 if v % 3 == 1 else 1),
                "minus marginal")

    H = [ix[u] for u in units if u % 3 == 1]
    K = [ix[u] for u in units if u % 3 == 2]
    A = [[M[i][j] for j in H] for i in H]
    C = [[M[i][j] for j in K] for i in K]
    determinant, detA, detC = bareiss(M), bareiss(A), bareiss(C)
    expected_block = 2 ** t * (4 ** t - 1) ** (t - 1)
    require(abs(detA) == expected_block, "first circulant block determinant")
    require(abs(detC) == expected_block, "second circulant block determinant")
    require(abs(determinant) == expected_block ** 2, "full determinant formula")
    require(determinant != 0, "full rank")
    return {"level": level, "modulus": q, "units_order": units, "denominator": D,
            "integer_matrix": M, "support_size": 3 * t * t, "matrix_dimension": d,
            "integer_determinant": determinant,
            "block_determinants": [detA, detC],
            "absolute_determinant_formula": f"2^({2*t}) * ({4**t-1})^({2*t-2})"}


rows = [construct(j) for j in range(1, 4)]
for low, high in zip(rows, rows[1:]):
    q, units, M, D = low["modulus"], low["units_order"], low["integer_matrix"], low["denominator"]
    q2, units2, M2, D2 = high["modulus"], high["units_order"], high["integer_matrix"], high["denominator"]
    ix = {u: i for i, u in enumerate(units)}
    projected = [[0] * len(units) for _ in units]
    for i, u in enumerate(units2):
        for j, v in enumerate(units2):
            projected[ix[u % q]][ix[v % q]] += M2[i][j]
    for i in range(len(units)):
        for j in range(len(units)):
            require(projected[i][j] * D == M[i][j] * D2, "tower projection coherence")

conditional_rows, conditional_cache = [], {}
for row in rows:
    j, q = row["level"], row["modulus"]
    for h in range(j):
        matrices = []
        for source_class in range(3 ** h):
            M = conditional_matrix(row, h, source_class)
            matrices.append(M)
            conditional_cache[(j, h, source_class)] = M
        for i in range(len(M)):
            for k in range(len(M)):
                require(sum(A[i][k] for A in matrices) == 3 ** h * row["integer_matrix"][i][k],
                        "conditional decomposition recovers full matrix")
        # One proper union, when available, tests the partial-sum rank theorem.
        if len(matrices) > 1:
            union = [[matrices[0][i][k] + matrices[-1][i][k] for k in range(len(M))]
                     for i in range(len(M))]
            require(rank(union) == 4 * 3 ** (j - h - 1), "two-class union exact rank")
        conditional_rows.append({"level": j, "conditioned_digits": h,
                                 "classes_checked": 3 ** h, "rank_each": 2 * 3 ** (j - h - 1)})

native_rows = []
words = [(1,), (2,), (1, 1), (1, 2), (2, 1), (1, 1, 1), (2, 1, 3)]
for sign in (-1, 1):
    for word in words:
        L, A, B, used = len(word), sum(word), 0, 0
        for a in word:
            B = 3 * B + (1 << used)
            used += a
        modulus = 1 << (A + 1)
        root = (((1 << A) - sign * B) * pow(3 ** L, -1, modulus)) % modulus
        endpoint = None
        for z in range(16):
            n = root + modulus * z
            x = n
            for a in word:
                require(valuation(3 * x + sign) == a, "native prefix valuation")
                x = odd_step(x, sign)
            require((1 << A) * x == 3 ** L * n + sign * B, "native affine endpoint")
            if endpoint is None:
                endpoint = x
            require(x == endpoint + 2 * 3 ** L * z, "native endpoint progression")
        ranks = []
        for row in rows:
            j, h = row["level"], min(L, row["level"] - 1)
            source_class = sign * B * pow(1 << A, -1, 3 ** h) % (3 ** h)
            require(endpoint % (3 ** h) == source_class, "native ternary endpoint class")
            actual = rank(conditional_cache[(j, h, source_class)])
            require(actual == 2 * 3 ** max(j - L - 1, 0), "native path conditional rank")
            ranks.append(actual)
        native_rows.append({"sign": sign, "word": word, "source_root": root,
                            "source_modulus": modulus, "endpoint_root": endpoint,
                            "paired_ranks_at_levels_1_2_3": ranks})

result = {"status": "FINITE-EXACT, with general proofs in companion note", "checks": checks,
          "levels": rows, "conditional_rank_checks": conditional_rows,
          "native_path_rank_checks": native_rows,
          "methods": ["Affine-line geometric compiler", "Independent exact source CRT compiler",
                      "Fraction-free Bareiss determinant", "Projection to preceding level"],
          "scope": "One-step same-source coupling; no temporal Markov or convergence claim."}
payload = json.dumps(result, sort_keys=True, separators=(",", ":")).encode()
result["payload_sha256"] = hashlib.sha256(payload).hexdigest()
print(json.dumps(result, sort_keys=True, indent=2))
