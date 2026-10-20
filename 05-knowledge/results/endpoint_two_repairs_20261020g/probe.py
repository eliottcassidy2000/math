"""Exact endpoint-two blue presentations and reusable labelled K-ledger compiler.

Universe: all standard closed blue paths with endpoint scores (n-3,2),
n=4..8, filtered explicitly for strong, mixed endpoints. No isomorphism
quotient and no broad tournament census. Standard library only.
"""
from collections import Counter
from hashlib import sha256
from itertools import combinations
from pathlib import Path
import json


def need(condition, message):
    if not condition:
        raise RuntimeError(message)


def arc(t, u, v):
    return (t[u] >> v) & 1


def edge(u, v):
    return tuple(sorted((u, v)))


def edges(p):
    return {edge(u, v) for u, v in zip(p, p[1:])}


def path(t, p):
    return len(p) == len(t) and len(set(p)) == len(t) and all(arc(t, u, v) for u, v in zip(p, p[1:]))


def blue(t, p):
    n = len(t)
    return all(arc(t, p[i], p[j]) == arc(t, p[n-1-j], p[n-1-i])
               for i, j in combinations(range(n), 2))


def anti_witness(t, p):
    n = len(t)
    for i, j in combinations(range(n), 2):
        if arc(t, p[i], p[j]) != arc(t, p[n-1-j], p[n-1-i]):
            return [p[i], p[j], p[n-1-j], p[n-1-i]]
    return None


def complement(t, p):
    need(path(t, p), "K requires a directed spanning word")
    full = (1 << len(t)) - 1
    rows = [full ^ row ^ (1 << u) for u, row in enumerate(t)]
    for u, v in zip(p, p[1:]):
        rows[u] ^= 1 << v
        rows[v] ^= 1 << u
    return tuple(rows)


def reverse_edges(t, p):
    rows = list(t)
    for u, v in zip(p, p[1:]):
        rows[u] ^= 1 << v
        rows[v] ^= 1 << u
    return tuple(rows)


def strong(t):
    n = len(t)
    for start in range(n):
        seen = 1 << start
        while True:
            expanded = seen
            for v in range(n):
                if seen >> v & 1:
                    expanded |= t[v]
            if expanded == seen:
                break
            seen = expanded
        if seen != (1 << n) - 1:
            return False
    return True


def paths(t, vertices=None, start=None, end=None):
    vertices = tuple(range(len(t))) if vertices is None else tuple(vertices)
    full = sum(1 << v for v in vertices)
    def visit(word, unused):
        if not unused:
            if end is None or word[-1] == end:
                yield tuple(word)
            return
        choices = unused if not word else unused & t[word[-1]]
        if end is not None and unused != 1 << end:
            choices &= ~(1 << end)
        while choices:
            bit = choices & -choices
            choices -= bit
            yield from visit(word + [bit.bit_length()-1], unused ^ bit)
    if start is None:
        yield from visit([], full)
    else:
        yield from visit([start], full ^ (1 << start))


def blue_presentations(n):
    free = {(u, v) for u in range(n) for v in range(u+2, n)}
    orbits = []
    while free:
        u, v = min(free)
        orbit = {(u, v), (n-1-v, n-1-u)}
        orbits.append(orbit)
        free -= orbit
    for mask in range(1 << len(orbits)):
        rows = [0] * n
        for u in range(n-1):
            rows[u] |= 1 << (u+1)
        for i, orbit in enumerate(orbits):
            for u, v in orbit:
                source, target = (u, v) if mask >> i & 1 else (v, u)
                rows[source] |= 1 << target
        yield tuple(rows)


def endpoint_two(n):
    """Construct from an arbitrary blue middle and one marked exceptional port."""
    for middle in blue_presentations(n-2):
        for x in range(2, n-1):
            y = n-1-x
            rows = [sum(1 << v for v in range(1, n-1))]
            rows += [(row << 1) | (1 << (n-1)) for row in middle]
            rows += [1]
            for u, v in [(0, x), (y, n-1)]:
                rows[u] ^= 1 << v
                rows[v] ^= 1 << u
            yield tuple(rows), x, y


def compile_ledger(n, words):
    """Exact original-arc direction guards, XOR target, affine score charges.

Guard [u,v,bit] means the original sorted-pair arc u->v equals bit.
After j K moves its orientation has parity j plus the prior path XOR.
This function uses no tournament and is reusable for proposed word templates.
"""
    used = set()
    requirements = {}
    conflicts = []
    sign, constants = 1, [0] * n
    charges = []
    for j, word in enumerate(words):
        need(len(word) == n and set(word) == set(range(n)), "ledger word is a permutation")
        s, t = word[0], word[-1]
        charges.append({"score_sign": sign, "score_vertices": [s, t],
                        "constant": constants[s] + constants[t] - (n-1)})
        for u, v in zip(word, word[1:]):
            pair = edge(u, v)
            desired = int(u < v) ^ (j & 1) ^ int(pair in used)
            if pair in requirements and requirements[pair] != desired:
                conflicts.append([*pair, requirements[pair], desired])
            requirements[pair] = desired
        used ^= edges(word)
        sign = -sign
        constants = [n-1-value for value in constants]
        constants[s] += 1
        constants[t] -= 1
    return {"guards": [[*pair, bit] for pair, bit in sorted(requirements.items())],
            "guard_conflicts": conflicts,
            "odd_K_parity": bool(len(words) & 1),
            "path_XOR_is_standard": used == edges(tuple(range(n))),
            "endpoint_charges": charges}


def evaluate(t, words):
    now = t
    statuses = []
    for j, word in enumerate(words):
        missing = [[u, v] for u, v in zip(word, word[1:]) if not arc(now, u, v)]
        if missing:
            return {"success": False, "invalid_step": j+1, "missing_arcs": missing,
                    "steps": statuses}
        witness = anti_witness(now, word)
        statuses.append({"black": witness is not None, "anti_failure": witness,
                         "endpoint_excess": now[word[0]].bit_count() + now[word[-1]].bit_count() - len(t) + 1,
                         "positional_reversal": [word[len(t)-1-word.index(v)] for v in range(len(t))]})
        now = complement(now, word)
    target = complement(t, tuple(range(len(t))))
    need(now == target, "valid ledger did not reach the labelled target")
    return {"success": all(item["black"] for item in statuses), "steps": statuses,
            "blue_steps": [i+1 for i, item in enumerate(statuses) if not item["black"]]}


def cyclic_words(n):
    a, b, c, z = 0, 1, n-2, n-1
    w = tuple(range(2, n-1))
    yield "release", ((b,)+w+(z,a), (z,)+w+(b,a), (a,z)+w+(b,))
    w = tuple(range(2, n-2))
    if len(w) >= 2:
        yield "closed_no_chord", ((c,)+w+(b,z,a), (z,c)+w+(b,a), (a,z,b)+w+(c,))
        yield "closed_chord", ((c,b,z,a)+w, (z,c,b,a)+w, (a,z,b)+w+(c,))


def reflect_words(n, words):
    return tuple(tuple(n-1-v for v in word[::-1]) for word in words)


def two_block(t):
    n = len(t)
    for pivot in range(1, n-1):
        left, right = tuple(range(pivot)), tuple(range(pivot+1, n))
        q = left + right + (pivot,)
        if not path(t, q):
            continue
        flipped = reverse_edges(reverse_edges(t, left), right)
        left_paths = list(paths(flipped, left, left[0], left[-1]))
        right_paths = list(paths(flipped, right, right[0], right[-1]))
        for hl in left_paths:
            for hr in right_paths:
                words = (q, tuple(reversed(hl+(pivot,)+hr)), (pivot,)+hr[::-1]+hl[::-1])
                result = evaluate(t, words)
                if result["success"]:
                    return {"kind": "two_block", "pivot": pivot, "words": words, "evaluation": result}
    return None


def square(t, all_paths):
    p = tuple(range(len(t)))
    pe = edges(p)
    for q in all_paths:
        if not (edges(q) & pe) and not blue(t, q):
            words = (q, p[::-1], q[::-1])
            result = evaluate(t, words)
            if result["success"]:
                return {"kind": "square", "words": words, "evaluation": result}
    return None


def open_hub_from_partner(t, partner_paths):
    n = len(t)
    p = tuple(range(n))
    partner = complement(t, p)
    for q in partner_paths:
        if (q[0], q[-1]) != (0, n-1) or blue(partner, q):
            continue
        for i, (u, v) in enumerate(zip(q, q[1:])):
            if v != u+1:
                continue
            rotation_q, rotation_p = q[i+1:]+q[:i+1], p[v:]+p[:v]
            words = (rotation_p, rotation_q, q)
            result = evaluate(t, words)
            if result["success"]:
                return {"kind": "open_partner_hub", "words": words, "evaluation": result}
    return None


def main():
    summaries, residuals, records = [], [], []
    compiled = {}
    for n in range(4, 9):
        stats = Counter()
        positions = {}
        direct_universe = set()
        for t, x, y in endpoint_two(n):
            need(t not in direct_universe, "two-port presentation duplicated")
            direct_universe.add(t)
            p = tuple(range(n))
            need(path(t, p) and blue(t, p) and arc(t, n-1, 0), "source blue closed path")
            need(t[0].bit_count() == n-3 and t[-1].bit_count() == 2, "endpoint scores")
            need([v for v in range(1, n-1) if arc(t, v, 0)] == [x], "marked incoming port")
            need([v for v in range(1, n-1) if arc(t, n-1, v)] == [y], "marked outgoing port")
            need(strong(t), "closed Hamilton path must be strong")
            stats["all_endpoint_two_presentations"] += 1
            partner = complement(t, p)
            if not strong(partner):
                stats["nonstrong_partner"] += 1
                continue
            hp, hp_partner = list(paths(t)), list(paths(partner))
            black = sum(not blue(t, word) for word in hp)
            black_partner = sum(not blue(partner, word) for word in hp_partner)
            if not black or not black_partner:
                stats["pure_endpoint"] += 1
                continue
            stats["strong_mixed"] += 1
            branch = "coincident" if x == y else ("distinct_transitive" if arc(t, x, y) else "distinct_R4")
            stats[branch] += 1
            position = positions.setdefault(str(x), Counter())
            position["total"] += 1
            evaluations, certs = {}, []
            for name, words in cyclic_words(n):
                for reflected in [False, True]:
                    label = name + ("_reflected" if reflected else "")
                    actual = reflect_words(n, words) if reflected else words
                    key = f"n{n}:{label}"
                    compiled.setdefault(key, {"words": actual, **compile_ledger(n, actual)})
                    need(not compiled[key]["guard_conflicts"] and compiled[key]["path_XOR_is_standard"], "inherited symbolic XOR")
                    result = evaluate(t, actual)
                    evaluations[label] = result
                    if result["success"]:
                        certs.append({"kind": label, "words": actual, "evaluation": result})
                        stats[label] += 1
                        position[label] += 1
            for finder, arguments in [(two_block, (t,)), (square, (t, hp)),
                                      (open_hub_from_partner, (t, hp_partner))]:
                cert = finder(*arguments)
                if cert:
                    certs.append(cert)
                    stats[cert["kind"]] += 1
                    position[cert["kind"]] += 1
            record = {"n": n, "rows": t, "partner_rows": partner, "x": x, "y": y,
                      "branch": branch, "scores": [row.bit_count() for row in t],
                      "path_counts": {"source": len(hp), "source_black": black,
                                      "partner": len(hp_partner), "partner_black": black_partner},
                      "cyclic_evaluations": evaluations,
                      "certificates": certs}
            records.append(record)
            if certs:
                stats["covered"] += 1
            else:
                stats["uncovered"] += 1
                position["uncovered"] += 1
                residuals.append(record)
        # Independent generator comparison: filter all standard blue presentations.
        filtered = {t for t in blue_presentations(n) if arc(t, n-1, 0) and t[0].bit_count() == n-3}
        need(direct_universe == filtered, "two-port model does not exhaust endpoint layer")
        summary = {"n": n, "counts": dict(stats), "by_x_position": {key: dict(value) for key, value in positions.items()}}
        summaries.append(summary)
        print(json.dumps(summary), flush=True)
    receipt = {"status": "FINITE-EXACT labelled presentation diagnostic; no all-order coverage claim",
               "source_sha256": sha256(Path(__file__).read_bytes()).hexdigest(),
               "universe": "All closed standard-blue presentations with d(first)=n-3,d(last)=2, n4..8; retain strong mixed endpoints; no isomorphism quotient",
               "summary": summaries, "residuals": residuals,
               "compiled_inherited_templates": compiled, "records": records}
    Path(__file__).with_suffix(".json").write_text(json.dumps(receipt, indent=2)+"\n", encoding="utf-8")
    print("PASS; residual count", len(residuals), flush=True)
    return receipt


if __name__ == "__main__":
    main()
