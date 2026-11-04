"""Independent finite audit of realized path-colour rank and affine probes.

No producer or inherited checker is imported. All tournament profiles in the
finite census and witness family come from literal vertex permutations.
"""

from argparse import ArgumentParser
from fractions import Fraction
from hashlib import sha256
from itertools import combinations, permutations
from math import comb, factorial
from pathlib import Path
import json


SOURCE = Path(__file__)
SOURCE_HASH = sha256(SOURCE.read_bytes()).hexdigest()
CHECKS = 0


def require(condition, message):
    global CHECKS
    CHECKS += 1
    if not condition:
        raise ValueError(message)


def tournament(n, code):
    out = [0] * n
    for bit, (i, j) in enumerate(combinations(range(n), 2)):
        u, v = (i, j) if code >> bit & 1 else (j, i)
        out[u] |= 1 << v
    return tuple(out)


def transitive(n):
    return tuple(sum(1 << j for j in range(i + 1, n)) for i in range(n))


def witness(n, j):
    out = list(transitive(n))
    if j:
        end = 2 * j
        out[0] ^= 1 << end
        out[end] ^= 1
    return tuple(out)


def profile(out):
    """Cut every backward adjacency; choose other cuts independently."""
    n = len(out)
    hist = [0] * n
    for word in permutations(range(n)):
        b = sum(not (out[u] >> v & 1) for u, v in zip(word, word[1:]))
        hist[b] += 1
    require(hist == hist[::-1], "permutation-reversal sidecar")
    q = [0] * (n + 1)
    for k in range(1, n + 1):
        q[k] = sum(hist[b] * comb(n - 1 - b, k - 1 - b)
                   for b in range(k))
    return q


def power_coefficients(q):
    n = len(q) - 1
    result = [Fraction(0)] * (n + 1)
    falling = [Fraction(1)]
    for k in range(1, n + 1):
        updated = [Fraction(0)] * (len(falling) + 1)
        for i, c in enumerate(falling):
            updated[i] -= (k - 1) * c
            updated[i + 1] += c
        falling = updated
        for i, c in enumerate(falling):
            result[i] += q[k] * c / factorial(k)
    return result


def rank(rows):
    if not rows:
        return 0
    echelon = {}
    for row in rows:
        v = list(map(Fraction, row))
        for i, basis in sorted(echelon.items()):
            if v[i]:
                c = v[i]
                v = [a - c * b for a, b in zip(v, basis)]
        pivot = next((i for i, c in enumerate(v) if c), None)
        if pivot is not None:
            c = v[pivot]
            echelon[pivot] = [a / c for a in v]
    return len(echelon)


def solve(matrix, values):
    n = len(matrix)
    if n == 0:
        return []
    a = [list(map(Fraction, row)) + [Fraction(v)]
         for row, v in zip(matrix, values)]
    for i in range(n):
        pivot = next(j for j in range(i, n) if a[j][i])
        a[i], a[pivot] = a[pivot], a[i]
        scale = a[i][i]
        a[i] = [x / scale for x in a[i]]
        for j in range(n):
            if j != i:
                scale = a[j][i]
                a[j] = [x - scale * y for x, y in zip(a[j], a[i])]
    return [row[-1] for row in a]


def interpolate(q, monic):
    n = len(q) - 1
    powers = list(range(1 if n % 2 else 2, n + 1, 2))
    if monic:
        powers.pop()
    values = [sum(q[k] * comb(t, k) for k in range(1, t + 1))
              - (t ** n if monic else 0)
              for t in range(1, len(powers) + 1)]
    matrix = [[t ** e for e in powers] for t in range(1, len(powers) + 1)]
    result = [Fraction(0)] * (n + 1)
    for e, c in zip(powers, solve(matrix, values)):
        result[e] = c
    if monic:
        result[n] = 1
    return result


def validate_profile(q):
    n = len(q) - 1
    coeff = power_coefficients(q)
    require(q[n] == factorial(n), "top ordered-cover count")
    require(coeff[n] == 1 and coeff[0] == 0, "monic/no constant")
    require(all(not c for e, c in enumerate(coeff) if (e - n) % 2), "parity")
    require(interpolate(q, False) == coeff, "linear half-grid recovery")
    require(interpolate(q, True) == coeff, "affine shortened-grid recovery")
    return coeff


def skeleton(r):
    """Enumerate exterior permutations and chosen hole slots literally."""
    m = 2 * r
    row = [0] * (r + 2)
    for order in permutations(range(m)):
        for k in range(1, r + 2):
            for slots in combinations(range(m + 1), k):
                slots = set(slots)
                word = []
                for i in range(m + 1):
                    if i in slots:
                        word.append(-1)
                    if i < m:
                        word.append(order[i])
                legal = True
                for u, v in zip(word, word[1:]):
                    if u == -1:
                        legal = v % 2 == 0
                    elif v == -1:
                        legal = u % 2 == 1
                    else:
                        legal = u < v
                    if not legal:
                        break
                row[k] += legal
    require(row[-1] == factorial(r), "actual-probe diagonal")
    return row


def wrap(out, r):
    n = len(out)
    extended = list(out) + [0] * (2 * r)
    for i in range(2 * r):
        for j in range(i + 1, 2 * r):
            extended[n + i] |= 1 << (n + j)
        for v in range(n):
            if i % 2 == 0:
                extended[v] |= 1 << (n + i)
            else:
                extended[n + i] |= 1 << v
    return tuple(extended)


def hamilton(out):
    """Different mechanism from profile: endpoint subset dynamic program."""
    n = len(out)
    paths = [[0] * n for _ in range(1 << n)]
    for v in range(n):
        paths[1 << v][v] = 1
    for mask in range(1, 1 << n):
        for end in range(n):
            count = paths[mask][end]
            if count:
                todo = out[end] & ~mask
                while todo:
                    bit = todo & -todo
                    todo -= bit
                    paths[mask | bit][bit.bit_length() - 1] += count
    return sum(paths[-1])


def probe_recovery(out, q, rows):
    n = len(out)
    d = (n + 1) // 2
    recovered = [0] * (n + 1)
    observations = []
    for r in range(d - 1):
        h = hamilton(wrap(out, r))
        observations.append(h)
        row = rows[r]
        expected = sum(row[k] * q[k] for k in range(1, len(row)))
        require(h == expected, "actual-context count")
        numerator = h - sum(row[k] * recovered[k] for k in range(1, r + 1))
        require(numerator % row[-1] == 0, "probe division exact")
        recovered[r + 1] = numerator // row[-1]
    require(interpolate(recovered, True) == power_coefficients(q), "actual affine probes recover all")
    return observations


def expected_witness(n, j):
    coeff = [0] * (n + 1)
    coeff[n] = 1
    if j:
        for k in range(1, 2 * j, 2):
            coeff[n - 1 - k] = 2 * comb(2 * j - 1, k)
    return coeff


def run():
    rows = [skeleton(r) for r in range(4)]
    census = []
    for n in range(1, 6):
        profiles = []
        count = 1 << (n * (n - 1) // 2)
        for code in range(count):
            out = tournament(n, code)
            q = profile(out)
            coeff = validate_profile(q)
            require(hamilton(out) == q[1], "independent Hamilton agreement")
            probe_recovery(out, q, rows)
            profiles.append(coeff)
        linear = rank(profiles)
        affine = rank([[a - b for a, b in zip(row, profiles[0])] for row in profiles[1:]])
        require(linear == (n + 1) // 2, "census linear rank")
        require(affine == (n + 1) // 2 - 1, "census affine rank")
        census.append({"n": n, "labelled_tournaments": count, "linear_rank": linear,
                       "affine_rank": affine})
    family = []
    for n in range(1, 9):
        vectors = []
        probes = []
        for j in range((n + 1) // 2):
            out = witness(n, j)
            q = profile(out)
            coeff = validate_profile(q)
            require(coeff == expected_witness(n, j), "one-back-edge colour formula")
            vectors.append(coeff)
            probes.append(probe_recovery(out, q, rows))
        require(rank(vectors) == (n + 1) // 2, "actual family rank")
        family.append({"n": n, "family_size": len(vectors), "affine_probe_values": probes})
    symbolic = []
    for n in range(1, 129):
        vectors = [expected_witness(n, j) for j in range((n + 1) // 2)]
        # Exact triangular support check; no claim of graph enumeration here.
        for j, coeff in enumerate(vectors):
            lowest = next(e for e, c in enumerate(coeff) if c)
            require(lowest == n - 2 * j, "triangular symbolic support")
            require(coeff[lowest] == (2 if j else 1), "triangular symbolic coefficient")
        symbolic.append((n + 1) // 2)
    a = (16, 9, 3, 5, 14)
    b = (24, 17, 3, 6, 12)
    qa, qb = profile(a), profile(b)
    require(qa[1] == qb[1] == 15 and qa != qb, "scalar hostile")
    actual = [hamilton(wrap(t, 1)) for t in (a, b)]
    require(actual == [123, 135], "actual scalar hostile wrapper")
    two = [transitive(2), (0, 1)]
    require(profile(two[0]) == profile(two[1]), "nonuniform hostile same profile")
    # Add x with x->0 and 1->x in both labelled bases.
    nonuniform = []
    for out in two:
        extended = (out[0], out[1] | 4, 1)
        nonuniform.append(hamilton(extended))
    require(nonuniform == [3, 1], "uniformity hypothesis is essential")
    require(hamilton(two[1]) == 1, "m=2 boundary: no 1+2^0 formula")
    require(sha256(SOURCE.read_bytes()).hexdigest() == SOURCE_HASH, "source unchanged during run")
    return {
        "status": "FINITE-EXACT independent audit passed",
        "source_sha256": SOURCE_HASH,
        "checks": CHECKS,
        "census": census,
        "actual_family": family,
        "skeleton_rows": rows,
        "symbolic_triangular_orders": [1, 128],
        "scalar_hostile": {"q_a": qa[1:], "q_b": qb[1:], "wrapped_H": actual},
        "nonuniform_hostile_H": nonuniform,
        "scope": "All labelled tournaments n=1..5; explicit family n=1..8; symbolic triangular support n=1..128. Universal proof is separate."
    }


if __name__ == "__main__":
    parser = ArgumentParser()
    parser.add_argument("--output", type=Path)
    args = parser.parse_args()
    data = json.dumps(run(), indent=2, sort_keys=True) + "\n"
    if args.output:
        args.output.write_text(data, encoding="utf-8", newline="\n")
    else:
        print(data, end="")
