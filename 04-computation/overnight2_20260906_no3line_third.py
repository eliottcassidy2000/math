"""Exact new-cycle coefficients of E[(X_3)_3] at n=8.

A compiled, exhaustive unordered-event census retains whole union incidence.
Python independently generates the grid triples and all skeleton edge-copy
profiles, applies rational label weights, and restores the ordered factor 6.
Only >=8-edge unions are needed for c8 and c4^2 by the proved cycle budget.
Literal n=4 permutation labels provide an independent normalization control.
All validity gates remain active under -O; no Monte Carlo is used.
"""
from collections import Counter, defaultdict
from fractions import Fraction as Q
from hashlib import sha256
from itertools import combinations, permutations
from math import comb, factorial
from pathlib import Path
import json
import subprocess
import tempfile

ROOT = Path(__file__).resolve().parents[1]
SOURCE = Path(__file__).with_suffix(".cpp")
DEST = ROOT / "05-knowledge/results/overnight2_20260906_no3line_third_certificates.json"
CHECKS = 0


def need(ok, message):
    global CHECKS
    CHECKS += 1
    if not ok:
        raise RuntimeError(message)


def literal_triples(n):
    points = list((r, c) for r in range(n) for c in range(n))
    result = []
    for a, b, c in combinations(points, 3):
        if len({a[0], b[0], c[0]}) < 3 or len({a[1], b[1], c[1]}) < 3:
            continue
        if (b[0]-a[0])*(c[1]-a[1]) == (c[0]-a[0])*(b[1]-a[1]):
            result.append(frozenset((a, b, c)))
    return result


def incidence(points):
    adjacency = defaultdict(set)
    for r, c in points:
        adjacency[(0, r)].add((1, c))
        adjacency[(1, c)].add((0, r))
    if any(len(v) > 2 for v in adjacency.values()):
        return None
    unseen = set(adjacency)
    components = []
    while unseen:
        seed = next(iter(unseen))
        stack, seen = [seed], set()
        while stack:
            vertex = stack.pop()
            if vertex not in seen:
                seen.add(vertex)
                stack.extend(adjacency[vertex]-seen)
        unseen.difference_update(seen)
        nr = sum(v[0] == 0 for v in seen)
        nc = len(seen)-nr
        edges = sum(len(adjacency[v]) for v in seen)//2
        cycle = all(len(adjacency[v]) == 2 for v in seen)
        components.append((edges, nr, nc, cycle))
    return tuple(sorted(components))


def encode(sig):
    tokens = sorted(3*m + (1 if r > c else 2 if c > r else 0)
                    for m, r, c, cycle in sig)
    return sum(token << (5*i) for i, token in enumerate(tokens))


def decode(key):
    components = []
    while key:
        token, key = key & 31, key >> 5
        m, variant = divmod(token, 3)
        if variant:
            r, c = m//2, m//2
            if variant == 1:
                r += 1
            else:
                c += 1
            cycle = False
        else:
            cycle = m % 2 == 0
            r = c = m//2 if cycle else (m+1)//2
        components.append((m, r, c, cycle))
    return tuple(sorted(components))


def skeleton(parts):
    points, start = [], 0
    for half_length in parts:
        need(half_length >= 2, "simple cycle length")
        for r in range(half_length):
            points.extend(((start+r, start+r),
                           (start+r, start+(r+1) % half_length)))
        start += half_length
    return points


def profiles(parts, maximum):
    points = skeleton(parts)
    counts = Counter()
    for m in range(1, min(maximum, len(points))+1):
        before = sum(counts.values())
        for subset in combinations(points, m):
            sig = incidence(subset)
            need(sig is not None, "degree-two skeleton subset")
            key = encode(sig)
            counts[key] += 1
        need(sum(counts.values())-before == comb(len(points), m), "complete edge-subset universe")
    return counts


def grid_copies(n, key):
    sig = decode(key)
    nr = sum(c[1] for c in sig)
    nc = sum(c[2] for c in sig)
    if nr > n or nc > n:
        return 0
    automorphisms = 1
    for (m, r, c, cycle), multiplicity in Counter(sig).items():
        one = m if cycle else 2 if m % 2 == 0 else 1
        automorphisms *= one**multiplicity * factorial(multiplicity)
    ordered = factorial(n)//factorial(n-nr) * (factorial(n)//factorial(n-nc))
    need(ordered % automorphisms == 0, "shore-copy orbit denominator")
    return ordered//automorphisms


def native_census(executable, n, minimum, literal=False):
    command = [str(executable), str(n), str(minimum)] + (["literal"] if literal else [])
    raw = subprocess.run(command, check=True,
                         capture_output=True, text=True).stdout
    rows = raw.splitlines()
    metadata = list(map(int, rows[0].split()[1:]))
    masks = list(map(int, rows[1].split()[1:]))
    geometry = {int(row.split()[1]): int(row.split()[2]) for row in rows[2:]}
    triples = literal_triples(n)
    reference = sorted(sum(1 << (8*r+c) for r, c in t) for t in triples)
    need(masks == reference, "independent determinant versus slope grid triples")
    need(metadata[:3] == [n, minimum, len(triples)], "native metadata")
    need(metadata[3] == comb(len(triples), 3), "all distinct unordered triple events")
    need(metadata[3] == sum(metadata[4:7]), "degree/size/accepted partition")
    need(sum(geometry.values()) == metadata[6] and len(geometry) == metadata[7], "native profile totals")
    for key in geometry:
        need(encode(decode(key)) == key, "component signature roundtrip")
        need(minimum <= sum(c[0] for c in decode(key)) <= 9, "union edge budget")
    return metadata, geometry, triples


def literal_label_moments(parts):
    points = skeleton(parts)
    n = sum(parts)
    totals = [0, 0, 0]
    labels = list(permutations(range(n)))
    for rows in labels:
        for cols in labels:
            board = [(rows[r], cols[c]) for r, c in points]
            x = 0
            for a, b, c in combinations(board, 3):
                if len({a[0], b[0], c[0]}) < 3 or len({a[1], b[1], c[1]}) < 3:
                    continue
                x += (b[0]-a[0])*(c[1]-a[1]) == (c[0]-a[0])*(b[1]-a[1])
            totals[0] += x
            totals[1] += x*(x-1)
            totals[2] += x*(x-1)*(x-2)
    return tuple(Q(t, len(labels)**2) for t in totals)


def main():
    with tempfile.TemporaryDirectory(prefix="overnight2_no3line_third_") as temporary:
        executable = Path(temporary)/"exact_census.exe"
        flags = [] if __debug__ else ["-DNDEBUG"]
        subprocess.run(["g++", "-std=c++17", "-O3", "-Wall", "-Wextra"]
                       + flags + [str(SOURCE), "-o", str(executable)], check=True)
        meta4, geom4, triples4 = native_census(executable, 4, 3)
        direct = Counter()
        rejected = 0
        for a, b, c in combinations(triples4, 3):
            sig = incidence(a | b | c)
            if sig is None:
                rejected += 1
            else:
                direct[encode(sig)] += 1
        need(dict(direct) == geom4 and rejected == meta4[4], "full Python/C++ n4 incidence census")
        controls = {}
        for parts in [(4,), (2, 2)]:
            bank = profiles(parts, 8)
            calculated = sum((Q(6*count*bank[key], grid_copies(4, key))
                              for key, count in geom4.items()), Q())
            literal = literal_label_moments(parts)
            need(calculated == literal[2], "full-label independent third factorial moment")
            need(literal[0] == 2, "inherited n4 mean control")
            expected_var = Q(25, 18) if parts == (4,) else Q(56, 9)
            need(literal[1]+literal[0]-literal[0]**2 == expected_var, "inherited n4 variance control")
            controls[str(parts)] = [str(x) for x in literal]
            print("n4 literal labels", parts, "factorial moments 1,2,3 =", controls[str(parts)])

        meta8, geometry, _ = native_census(executable, 8, 8)
        literal_meta8, literal_geometry, _ = native_census(executable, 8, 8, literal=True)
        need((meta8, geometry) == (literal_meta8, literal_geometry),
             "complete literal-degree replay versus saturated-vertex filter")
    print("n8 metadata n,min_edges,triples,unordered_total,degree_rejected,small_rejected,accepted,types:", meta8)
    need(meta8 == [8, 8, 648, 45139896, 20781112, 2939296, 21419488, 150], "frozen exhaustive universe")

    cycle_types = [(8,), (4, 4), (3, 5), (2, 6), (2, 3, 3), (2, 2, 4), (2, 2, 2, 2)]
    banks = {parts: profiles(parts, 9) for parts in cycle_types}
    universe = sorted(set(geometry).union(*(set(bank) for bank in banks.values())))
    coefficients = {}
    for key in universe:
        a = Q(banks[(8,)][key])
        square = Q(banks[(2, 2, 2, 2)][key] - 4*banks[(2, 6)][key] + 3*a, 12)
        b = banks[(2, 6)][key]-a-square
        d = banks[(3, 5)][key]-a
        e = Q(banks[(4, 4)][key]-a, 2)
        coefficients[key] = (a, b, d, e, square)
        for parts in cycle_types:
            c4, c6, c8 = parts.count(2), parts.count(3), parts.count(4)
            need(banks[parts][key] == a+b*c4+d*c6+e*c8+square*c4*c4,
                 "all cycle types obey exact weighted-degree profile law")
        if sum(c[0] for c in decode(key)) <= 7:
            need(e == square == 0, "short-edge annihilation control")
    need(banks[(4, 4)][24] == 2 and banks[(8,)][24] == 0, "eight-cycle positive witness")
    need(banks[(2, 2, 2, 2)][396] == 6, "two-four-cycle automorphism witness")

    totals = [Q(), Q()]
    signed = [[Q(), Q()], [Q(), Q()]]
    by_size = {8: [Q(), Q()], 9: [Q(), Q()]}
    rows = []
    for key, count in sorted(geometry.items()):
        denominator = grid_copies(8, key)
        need(denominator > 0, "geometric profile occurs in grid")
        e, square = coefficients[key][3:]
        contributions = [Q(6*count)*x/denominator for x in (e, square)]
        size = sum(c[0] for c in decode(key))
        for j, contribution in enumerate(contributions):
            totals[j] += contribution
            signed[j][0 if contribution >= 0 else 1] += contribution
            by_size[size][j] += contribution
        rows.append({"key": key, "signature": decode(key), "unordered_event_triples": count,
                     "grid_copies": denominator,
                     "skeleton_copies": {str(parts): banks[parts][key] for parts in cycle_types},
                     "profile_coefficients_A_B_D_E_F": list(map(str, coefficients[key])),
                     "ordered_weighted_contributions_E_F": list(map(str, contributions))})
    need(totals == [Q(172483, 529200), Q(11881, 50400)], "frozen new-cycle coefficients")
    need(by_size == {8: [Q(768463, 2116800), Q(456371, 2116800)],
                     9: [Q(-26177, 705600), Q(42631, 2116800)]}, "frozen union-size contributions")
    need(all(value > 0 for value in totals), "both new-cycle terms are actually detected")
    print("All 45,139,896 event triples replayed with literal degrees; profile counts match.")
    print("Full third factorial moment coefficients: c8 =", totals[0], "; c4^2 =", totals[1])
    for size in (8, 9):
        print("Contribution from", size, "union edges:", list(map(str, by_size[size])))
    for label, (positive, negative), total in zip(("c8", "c4^2"), signed, totals):
        print(label, "positive total", positive, "; negative total", negative, "; net", total)
    print("Raw contrasts: 2C8 minus C16 =", 2*totals[0],
          "; 4C4 minus 4(C4+C12) plus 3C16 =", 12*totals[1])
    manifest = {"scope": "n=8, full third factorial moment; coefficients c8 and c4^2 only",
                "native_metadata_n4": meta4, "native_metadata_n8": meta8,
                "n4_literal_factorial_moments": controls,
                "ordered_factor": 6,
                "coefficients": {"c8": str(totals[0]), "c4_squared": str(totals[1])},
                "contributions_by_union_size": {str(k): list(map(str, v)) for k, v in by_size.items()},
                "signed_contributions": {label: list(map(str, vals)) for label, vals in zip(("c8", "c4_squared"), signed)},
                "profiles": rows}
    encoded = (json.dumps(manifest, indent=2, sort_keys=True)+"\n").encode()
    DEST.write_bytes(encoded)
    print("Certificate SHA256:", sha256(encoded).hexdigest())
    print("Exact gates:", CHECKS)
    print("PASS: exhaustive geometric multiplicities; no independence, sampling, or density conclusion.")


if __name__ == "__main__":
    main()
