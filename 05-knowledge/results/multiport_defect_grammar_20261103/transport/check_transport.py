"""Integral endpoint-forest transport, including flips and growing ports.

Standard library only. Explicit require() checks survive python -O.
"""
from __future__ import annotations
import argparse
from collections import Counter
from functools import lru_cache
import hashlib
from itertools import combinations, permutations
import json
from pathlib import Path

HERE = Path(__file__).resolve().parent
RESULTS = HERE.parent.parent
STAR = -1


def require(value, message):
    if not value:
        raise RuntimeError(message)


def clean(vector):
    return Counter({s: c for s, c in vector.items() if c})


def canonical(k, pairs):
    pairs = list(pairs)
    k += pairs.count((STAR, STAR))
    return k, tuple(sorted(p for p in pairs if p != (STAR, STAR)))


@lru_cache(None)
def forests(rows):
    """Independent literal permutation-and-cut enumeration, unordered covers."""
    n = len(rows)
    if not n:
        return ((),)
    found = set()
    for word in permutations(range(n)):
        mandatory = sum(1 << i for i in range(n - 1)
                        if not rows[word[i]] >> word[i + 1] & 1)
        for cuts in range(1 << (n - 1)):
            if mandatory & ~cuts:
                continue
            paths, start = [], 0
            for i in range(n - 1):
                if cuts >> i & 1:
                    paths.append(word[start:i + 1])
                    start = i + 1
            paths.append(word[start:])
            found.add(tuple(sorted(paths)))
    return tuple(sorted(found))


def project(rows, ports):
    ports = frozenset(ports)
    result = Counter()
    for cover in forests(tuple(rows)):
        pairs = [(p[0] if p[0] in ports else STAR,
                  p[-1] if p[-1] in ports else STAR) for p in cover]
        result[canonical(0, pairs)] += 1
    return result


def join_signature(signature, a, b):
    k, pairs = signature
    ia = next((i for i, p in enumerate(pairs) if p[1] == a), None)
    ib = next((i for i, p in enumerate(pairs) if p[0] == b), None)
    if ia is None or ib is None or ia == ib:
        return None
    merged = (pairs[ia][0], pairs[ib][1])
    return canonical(k, [p for i, p in enumerate(pairs) if i not in (ia, ib)] + [merged])


def join(vector, a, b):
    result = Counter()
    for state, count in vector.items():
        target = join_signature(state, a, b)
        if target is not None:
            result[target] += count
    return clean(result)


def flip(vector, a, b):
    result = Counter(vector)
    for state, count in join(vector, b, a).items():
        result[state] += count
    for state, count in join(vector, a, b).items():
        result[state] -= count
    return clean(result)


def forget(vector, port):
    result = Counter()
    for (k, pairs), count in vector.items():
        target = canonical(k, [(STAR if a == port else a, STAR if b == port else b)
                               for a, b in pairs])
        result[target] += count
    return clean(result)


def converse(vector):
    result = Counter()
    for (k, pairs), count in vector.items():
        result[canonical(k, [(b, a) for a, b in pairs])] += count
    return clean(result)


def birth(vector, new, incoming_ports, anonymous_incoming):
    """Old active cut is arbitrary; ALL forgotten vertices have one orientation.

    ANON is one formal anonymous component, selected with multiplicity k.
    A single uniform birth cannot choose anonymous components for both roles.
    """
    incoming_ports = frozenset(incoming_ports)
    ANON = -2
    result = Counter()
    for (k, pairs), count in vector.items():
        incoming = [None] + [i for i, (_, b) in enumerate(pairs)
                             if (anonymous_incoming if b == STAR else b in incoming_ports)]
        outgoing = [None] + [i for i, (a, _) in enumerate(pairs)
                             if not (anonymous_incoming if a == STAR else a in incoming_ports)]
        if k:
            (incoming if anonymous_incoming else outgoing).append(ANON)
        for left in incoming:
            for right in outgoing:
                if left is not None and left == right:
                    continue
                lp = (STAR, STAR) if left == ANON else (pairs[left] if left is not None else None)
                rp = (STAR, STAR) if right == ANON else (pairs[right] if right is not None else None)
                selected_anon = int(left == ANON or right == ANON)
                remaining = [p for i, p in enumerate(pairs) if i not in (left, right)]
                remaining.append((lp[0] if lp else new, rp[1] if rp else new))
                target = canonical(k - selected_anon, remaining)
                result[target] += count * (k if selected_anon else 1)
    return clean(result)


def readout(vector, root=None):
    h = f = 0
    for (k, pairs), count in vector.items():
        d = k + len(pairs)
        if d == 1:
            h += count
        if root is not None and d in (1, 2) and any(b == root for _, b in pairs):
            f += count
    return {"H": h, "F": f, "defect": f - h} if root is not None else {"H": h}


def literal_birth(rows, ports, incoming, anon_incoming):
    n = len(rows)
    result = list(rows) + [0]
    for v in range(n):
        into = v in incoming if v in ports else anon_incoming
        if into:
            result[v] |= 1 << n
        else:
            result[n] |= 1 << v
    return tuple(result)


def literal_flip(rows, a, b):
    result = list(rows)
    result[a] ^= 1 << b
    result[b] ^= 1 << a
    return tuple(result)


def strong(rows):
    n = len(rows)
    full = (1 << n) - 1
    reach = [row | (1 << v) for v, row in enumerate(rows)]
    for k in range(n):
        for i in range(n):
            if reach[i] >> k & 1:
                reach[i] |= reach[k]
    return all(r == full for r in reach)


def tournaments(n):
    edges = list(combinations(range(n), 2))
    for bits in range(1 << len(edges)):
        rows = [0] * n
        for i, (a, b) in enumerate(edges):
            u, v = (a, b) if bits >> i & 1 else (b, a)
            rows[u] |= 1 << v
        yield tuple(rows)


def small_controls():
    counts = Counter()
    for n in range(1, 5):
        for rows in tournaments(n):
            for pmask in range(1 << n):
                ports = {i for i in range(n) if pmask >> i & 1}
                old = project(rows, ports)
                reversed_rows = tuple(((1 << n) - 1) ^ row ^ (1 << v) for v, row in enumerate(rows))
                require(converse(old) == project(reversed_rows, ports), "converse")
                counts["converse"] += 1
                for r in ports:
                    require(forget(old, r) == project(rows, ports - {r}), "forget")
                    counts["forget"] += 1
                for a, b in combinations(sorted(ports), 2):
                    if not rows[a] >> b & 1:
                        a, b = b, a
                    changed = flip(old, a, b)
                    require(all(c >= 0 for c in changed.values()), "actual flip negativity")
                    require(changed == project(literal_flip(rows, a, b), ports), "flip")
                    require(flip(changed, b, a) == old, "inverse flip")
                    require(not join(join(old, a, b), b, a), "opposite joins")
                    counts["flip"] += 1
                pp = sorted(ports)
                for imask in range(1 << len(pp)):
                    incoming = {p for i, p in enumerate(pp) if imask >> i & 1}
                    for eps in (False, True):
                        actual = literal_birth(rows, ports, incoming, eps)
                        predicted = birth(old, n, incoming, eps)
                        require(predicted == project(actual, ports | {n}), "birth")
                        counts["birth"] += 1
            # All directed endpoint joins commute, including inadmissible pairs.
            old = project(rows, set(range(n)))
            arcs = [(a, b) for a in range(n) for b in range(n) if a != b]
            endpoint = EndpointEvaluator(rows)
            for a, b in arcs:
                joined = join(old, a, b)
                require(readout(joined)["H"] == endpoint.J(a, b), "two-slot join readout")
                for root in range(n):
                    require(readout(joined, root)["F"] == endpoint.F_after_join(root, a, b),
                            "three-slot join readout")
                    counts["join_readout"] += 1
            for e, f in combinations(arcs, 2):
                require(join(join(old, *e), *f) == join(join(old, *f), *e), "join commutation")
                counts["commutation"] += 1
            # C3(T,1,1) is an all-order positive two-birth+forget control.
            plain = project(rows, set())
            first = birth(plain, n, set(), True)
            second = birth(first, n + 1, {n}, False)
            collapsed = forget(forget(second, n), n + 1)
            polynomial = Counter()
            for (k, pairs), amount in plain.items():
                require(not pairs, "anonymous profile")
                for degree, coefficient in ((k + 2, 1), (k + 1, 2 * k + 1),
                                             (k, k * (k + 2)), (k - 1, k * (k - 1))):
                    if coefficient:
                        polynomial[(degree, ())] += amount * coefficient
            require(collapsed == polynomial, "two-birth differential operator")
            counts["C3_operator"] += 1
            # The inherited adaptive reset factors through the C3 macro and a
            # full star correction. All old physical vertices are retained.
            for w in range(n):
                xx = birth(old, n, set(range(n)), True)
                yy = birth(xx, n + 1, {n}, False)
                for v in range(n):
                    if v != w:
                        yy = flip(yy, v, n)
                expected_h = endpoint.h[-1] + 2 * endpoint.F(w)
                require(readout(yy, n) == {"H": expected_h, "F": expected_h, "defect": 0},
                        "zero-ready reset factorization")
                direct = literal_birth(rows, set(range(n)), {w}, False)
                direct = literal_birth(direct, set(range(n + 1)), {n}, False)
                require(yy == project(direct, set(range(n + 2))), "full reset profile")
                counts["reset_factorization"] += 1
    # Forgetting is exact, but arbitrary reactivation is not a legal operation.
    a, b = (2, 0), (0, 1)
    require(project(a, set()) == project(b, set()), "reactivation source equality")
    xa = literal_birth(a, {0, 1}, {0}, False)
    xb = literal_birth(b, {0, 1}, {0}, False)
    pair = [readout(project(t, set()))["H"] for t in (xa, xb)]
    require(pair == [1, 3], "reactivation hostile")
    # End sets agree, but only the first pairing permits the new1->2 join.
    s1 = canonical(0, [(0, 1), (2, 3)])
    s2 = canonical(0, [(0, 3), (2, 1)])
    require(join_signature(s1, 1, 2) is not None and join_signature(s2, 1, 2) is None,
            "component pairing hostile")
    return dict(counts), {"rows_pair": [a, b], "same_forgotten_polynomial": "z+z^2",
                         "reactivate_0_birth_H_pair": pair,
                         "pairing_hostile": [s1, s2]}


def endpoint_dp(rows, first=None):
    n = len(rows)
    full = (1 << n) - 1
    ends = [[0] * n for _ in range(full + 1)]
    h = [0] * (full + 1)
    h[0] = int(first is None)
    for v in range(n):
        if first is None or v == first:
            ends[1 << v][v] = 1
    for mask in range(1, full + 1):
        for v, count in enumerate(ends[mask]):
            if not count:
                continue
            h[mask] += count
            choices = rows[v] & (full ^ mask)
            while choices:
                bit = choices & -choices
                choices -= bit
                ends[mask | bit][bit.bit_length() - 1] += count
    return h, ends


def submasks(mask):
    subset = mask
    while True:
        yield subset
        if not subset:
            break
        subset = (subset - 1) & mask


def triple(a, b, c, fa, fb, fc):
    if fa & fb or fa & fc or fb & fc:
        return 0
    full = len(a) - 1
    free = full ^ (fa | fb | fc)
    total = 0
    for aa in submasks(free):
        am = aa | fa
        av = a[am]
        if not av:
            continue
        rest = free ^ aa
        for bb in submasks(rest):
            bm = bb | fb
            bv = b[bm]
            if bv:
                total += av * bv * c[(rest ^ bb) | fc]
    return total


class EndpointEvaluator:
    def __init__(self, rows):
        self.rows = tuple(rows)
        n = len(rows)
        self.full = (1 << n) - 1
        self.h, ends = endpoint_dp(rows)
        reverse = [self.full ^ row ^ (1 << i) for i, row in enumerate(rows)]
        _, starts = endpoint_dp(reverse)
        self.e = [[row[i] for row in ends] for i in range(n)]
        self.s = [[row[i] for row in starts] for i in range(n)]
        self.fixed = {}

    def g(self, first, last):
        if first not in self.fixed:
            self.fixed[first] = endpoint_dp(self.rows, first)[1]
        return [row[last] for row in self.fixed[first]]

    def F(self, root):
        return sum(self.e[root][m] * self.h[self.full ^ m] for m in range(self.full + 1))

    def J(self, a, b):
        return sum(self.e[a][m] * self.s[b][self.full ^ m] for m in range(self.full + 1))

    def F_after_join(self, root, a, b):
        first = triple(self.e[a], self.g(b, root), self.h,
                       1 << a, (1 << b) | (1 << root), 0)
        second = triple(self.e[a], self.s[b], self.e[root],
                        1 << a, 1 << b, 1 << root)
        return first + second


def collar_control():
    path = RESULTS / "zero_ready_strong_switch_20261103/history_phase/collar.json"
    data = json.loads(path.read_text())
    source, target = data["records"]
    names, rows = source["names"], tuple(source["rows"])
    root = names.index(32)
    transcript = []
    for left, right in ((4, 33), (22, 32), (31, 33)):
        a, b = names.index(left), names.index(right)
        if not rows[a] >> b & 1:
            a, b = b, a
        before = EndpointEvaluator(rows)
        plus_h, minus_h = before.J(b, a), before.J(a, b)
        plus_f = before.F_after_join(root, b, a)
        minus_f = before.F_after_join(root, a, b)
        changed = literal_flip(rows, a, b)
        after = EndpointEvaluator(changed)
        old_h, old_f = before.h[-1], before.F(root)
        new_h, new_f = after.h[-1], after.F(root)
        require(new_h == old_h + plus_h - minus_h, "collar H join transport")
        require(new_f == old_f + plus_f - minus_f, "collar F three-slot transport")
        if root in (a, b):
            other = b if root == a else a
            vv = [i for i in range(len(rows)) if i != root]
            sub = [sum(1 << j for j, v in enumerate(vv) if rows[u] >> v & 1) for u in vv]
            deletion_f = EndpointEvaluator(sub).F(vv.index(other))
            require(new_f - old_f == (1 if root == a else -1) * deletion_f,
                    "root-incident deletion identity")
        else:
            deletion_f = None
        incoming = [i for i in range(len(rows)) if changed[i] >> root & 1]
        defect_parts = {str(names[i]): after.J(root, i) for i in incoming}
        require(sum(defect_parts.values()) == new_f - new_h, "positive defect partition")
        transcript.append({"old_arc": [names[a], names[b]],
                           "before": {"H": old_h, "F32": old_f, "b32": old_f - old_h},
                           "H_join_plus_minus": [plus_h, minus_h],
                           "F_join_plus_minus": [plus_f, minus_f],
                           "incident_deleted_F": deletion_f,
                           "after": {"H": new_h, "F32": new_f, "b32": new_f - new_h},
                           "positive_defect_parts": defect_parts})
        rows = changed
    require(list(rows) == target["rows"], "exact collar target adjacency")
    require(transcript[-1]["after"] == {"H": 80319, "F32": 83279, "b32": 2960}, "collar target")
    require(transcript[-1]["positive_defect_parts"] == {"22": 2070, "33": 890}, "collar split")
    # Independent Boolean cube: the highest mixed term really samples five
    # initial components in the two-component F readout.
    edge_names = [(4, 33), (22, 32), (31, 33)]
    corners = []
    for mask in range(8):
        rr = tuple(source["rows"])
        for i, (a, b) in enumerate(edge_names):
            if mask >> i & 1:
                rr = literal_flip(rr, names.index(a), names.index(b))
        e = EndpointEvaluator(rr)
        h, f, er = e.h[-1], e.F(root), e.e[root][-1]
        corners.append([h, f, er, f - er, f - h])
    mobius = [r.copy() for r in corners]
    for i in range(3):
        for mask in range(8):
            if mask >> i & 1:
                mobius[mask] = [a - b for a, b in zip(mobius[mask], mobius[mask ^ (1 << i)])]
    require(mobius[7] == [396, 396, 0, 396, 0], "five-component mixed response")

    # Peel port33 from D=B14-32. Q is fixed throughout arbitrary33-star edits.
    qnames = [v for v in names if v not in (32, 33)]
    qindices = [names.index(v) for v in qnames]
    qrows = [sum(1 << j for j, v in enumerate(qindices) if source["rows"][u] >> v & 1)
             for u in qindices]
    qe = EndpointEvaluator(qrows)
    u = qnames.index(22)
    base = qe.e[u][-1] + qe.s[u][-1]
    arms = []
    for j, name in enumerate(qnames):
        if name == 22:
            continue
        out_value, in_value = qe.J(u, j), qe.J(j, u)
        old_out = bool(source["rows"][names.index(33)] >> names.index(name) & 1)
        arms.append({"name": name, "33_out_value": out_value, "33_in_value": in_value,
                     "source_33_out": old_out, "protected": name == 40})
    current = base + sum(a["33_out_value"] if a["source_33_out"] else a["33_in_value"] for a in arms)
    minimum = base + sum(a["33_in_value"] if a["protected"] else min(a["33_out_value"], a["33_in_value"])
                         for a in arms)
    maximum = base + sum(a["33_in_value"] if a["protected"] else max(a["33_out_value"], a["33_in_value"])
                         for a in arms)
    free_arms = [a for a in arms if not a["protected"]]
    spectrum = Counter()
    for mask in range(1 << len(free_arms)):
        value = base + sum(a["33_in_value"] for a in arms if a["protected"])
        value += sum(a["33_out_value"] if mask >> i & 1 else a["33_in_value"]
                     for i, a in enumerate(free_arms))
        spectrum[value] += 1
    require((current, minimum, maximum) == (2900, 614, 7838), "star bank bounds")
    require(min(spectrum) == minimum and max(spectrum) == maximum, "sharp star bounds")
    require(sum(spectrum.values()) == 1024 and not spectrum.get(0), "all star arm choices")
    arm4 = next(a for a in arms if a["name"] == 4)
    arm31 = next(a for a in arms if a["name"] == 31)
    delta4 = arm4["33_out_value"] - arm4["33_in_value"]
    delta31 = arm31["33_in_value"] - arm31["33_out_value"]
    require((delta4, delta31, current + delta4 + delta31) == (-198, 258, 2960),
            "exact star defect transport")
    extremizers = []
    for minimum_case in (True, False):
        rr = literal_flip(source["rows"], names.index(32), names.index(22))
        for a in arms:
            if a["protected"]:
                desired_out = False
            else:
                desired_out = (a["33_out_value"] <= a["33_in_value"]) if minimum_case else \
                              (a["33_out_value"] >= a["33_in_value"])
            w, x = names.index(33), names.index(a["name"])
            if bool(rr[w] >> x & 1) != desired_out:
                rr = literal_flip(rr, w, x)
        evaluator = EndpointEvaluator(rr)
        defect = evaluator.F(root) - evaluator.h[-1]
        require(defect == (minimum if minimum_case else maximum), "literal star extremizer")
        extremizers.append({"kind": "minimum" if minimum_case else "maximum", "rows": rr,
                            "strong": strong(rr), "H": evaluator.h[-1], "F32": evaluator.F(root),
                            "b32": defect})
    star = {"fixed_Q_names": qnames, "fixed_Q_rows": qrows, "H_Q": qe.h[-1],
            "end22_Q": qe.e[u][-1], "start22_Q": qe.s[u][-1], "arms": arms,
            "initial_pair_kernel": current, "actual_arm_increments": [delta4, delta31],
            "protected_prefix": 40, "inport_pair": [22, 33],
            "sharp_minimum": minimum, "sharp_maximum": maximum,
            "free_arm_patterns": 1024, "distinct_defects": len(spectrum),
            "defect_histogram": sorted(spectrum.items()), "literal_extremizers": extremizers}
    cube = {"edge_names": edge_names, "column_names": ["H", "F32", "end32", "two_components_ending32", "b32"],
            "corners": corners, "mobius_coefficients": mobius}
    return transcript, cube, star, {str(path.relative_to(RESULTS)): hashlib.sha256(path.read_bytes()).hexdigest()}


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("--output", type=Path, default=HERE / "verification.json")
    args = parser.parse_args()
    controls, hostile = small_controls()
    collar, cube, star, pins = collar_control()
    g = [1, 4]
    for p in range(2, 9):
        g.append(4 * g[-1] + 2 * (p - 1) * g[-2])
    result = {"status": "PROVED typed transport; FINITE-EXACT controls; two-branch closure OPEN",
              "small_universe": "Every labeled tournament n=1..4, every active subset, every legal active flip, every uniform-anonymous birth cut, every forget; independent literal cover enumeration.",
              "controls": controls, "hostile_boundaries": hostile,
              "signature_dimensions_p0_to8": g, "collar_transcript": collar,
              "collar_flip_cube": cube, "closed_star_defect_bank": star,
              "source_sha256": pins,
              "checker_sha256": hashlib.sha256(Path(__file__).read_bytes()).hexdigest()}
    args.output.write_text(json.dumps(result, indent=2) + "\n")
    print(json.dumps({"controls": controls, "signature_dimensions": g,
                      "collar_transcript": collar}, indent=2))


if __name__ == "__main__":
    main()
