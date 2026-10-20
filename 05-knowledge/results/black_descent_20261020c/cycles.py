"""Exact controls for the cycle-changing proof. Run from any directory.

The cross-set census is independent of tournament enumeration. The latter
uses the previously authenticated complete class files and examines every
self-converse strong class at n=8,9; non-self-converse classes have only
black paths and satisfy descent by the cycle-cut averaging proof.
All checks use explicit exceptions, so they survive python -O.
"""
from collections import Counter
import csv
import gzip
from itertools import combinations, product
import json
from pathlib import Path

HERE = Path(__file__).resolve().parent
ROOT = HERE.parents[2]


def need(ok, why):
    if not ok:
        raise RuntimeError(why)


def rows(code, n):
    adj = [0] * n
    bit = n * (n - 1) // 2 - 1
    for u in range(n):
        for v in range(u + 1, n):
            a, b = (u, v) if code >> bit & 1 else (v, u)
            adj[a] |= 1 << b
            bit -= 1
    return adj


def blue(adj, path):
    n = len(adj)
    rev = [0] * n
    for i, v in enumerate(path):
        rev[v] = path[n - 1 - i]
    return all((adj[u] >> v & 1) == (adj[rev[v]] >> rev[u] & 1)
               for u in range(n) for v in range(n) if u != v)


def strongly_connected(adj):
    n = len(adj)
    for backwards in (False, True):
        seen = 1
        while True:
            new = seen
            for u in range(n):
                if seen >> u & 1:
                    if backwards:
                        new |= sum(1 << v for v in range(n)
                                   if adj[v] >> u & 1)
                    else:
                        new |= adj[u]
            if new == seen:
                break
            seen = new
        if seen != (1 << n) - 1:
            return False
    return True


def hamiltonian_cycles(adj):
    n = len(adj)

    def visit(path, used):
        if len(path) == n:
            if adj[path[-1]] & 1:
                yield tuple(path)
            return
        available = adj[path[-1]] & ~used
        while available:
            one = available & -available
            v = one.bit_length() - 1
            yield from visit(path + [v], used | one)
            available ^= one
    yield from visit([0], 1)


def descent(adj, cycle):
    score = [r.bit_count() for r in adj]
    for i in range(len(adj)):
        path = cycle[i:] + cycle[:i]
        if score[path[-1]] >= score[path[0]] and not blue(adj, path):
            return path
    return None


def check_descent_certificate(adj, path):
    n = len(adj)
    need(len(path) == n and len(set(path)) == n, "descent certificate does not span")
    need(all(adj[a] >> b & 1 for a, b in zip(path, path[1:]+path[:1])),
         "descent certificate is not a closed path")
    need(not blue(adj, path), "descent certificate is blue")
    changed = list(adj)
    for a, b in zip(path, path[1:]):
        changed[a] ^= 1 << b
        changed[b] ^= 1 << a

    def triangles(arcs):
        return sum(bool((arcs[a] >> b & 1) and (arcs[b] >> c & 1) and (arcs[c] >> a & 1))
                   or bool((arcs[b] >> a & 1) and (arcs[c] >> b & 1) and (arcs[a] >> c & 1))
                   for a, b, c in combinations(range(n), 3))
    difference = triangles(changed) - triangles(adj)
    expected = adj[path[0]].bit_count() - adj[path[-1]].bit_count() - 1
    need(difference == expected and difference < 0, "literal triangle descent check failed")
    return difference


def terminal_splice(adj, cycle):
    """Literal construction of Section 6 for a bad even cycle with two blue cuts."""
    n = len(adj)
    need(n >= 6 and n % 2 == 0, "terminal splice requires even n>=6")
    need(descent(adj, cycle) is None, "terminal-splice input already has descent")
    score = [r.bit_count() for r in adj]
    start = min(range(n), key=lambda i: score[cycle[i]])
    c = cycle[start:] + cycle[:start]
    m = n // 2
    need([i for i in range(n) if blue(adj, c[i:]+c[:i])] == [0, m],
         "terminal-splice input does not have two opposite blue cuts")
    delta = score[c[0]]
    need(delta <= m // 2, "two increasing blocks violate the score bound")
    yi = max(i for i in range(1, n-2) if adj[c[i]] >> c[0] & 1)
    k = yi + 1
    t = n - k
    need(2 <= t < m and t <= delta + 1, "terminal block length bound")
    x, u = c[k], c[-1]
    need(score[c[yi]] < score[x] and score[u] - score[x] >= t-1,
         "terminal scores do not strictly rise")
    candidates = [j for j in range(k-1) if adj[c[j]] >> x & 1]
    need(len(candidates) >= delta, "not enough prefix vertices beat block start")
    j = next(j for j in candidates if adj[u] >> c[j+1] & 1)
    new_cycle = c[:j+1] + c[k:] + c[j+1:k]
    need(len(set(new_cycle)) == n, "spliced cycle does not span")
    need(all(adj[a] >> b & 1 for a, b in
             zip(new_cycle, new_cycle[1:]+new_cycle[:1])),
         "spliced cycle has missing arc")
    need(not blue(adj, new_cycle), "spliced cut is not black")
    need(score[new_cycle[-1]] >= score[new_cycle[0]], "spliced cut does not descend")
    difference = check_descent_certificate(adj, new_cycle)
    return dict(original_minimum_rooted_cycle=c, terminal_block=c[k:],
                insertion_after=c[j], repaired_cycle=new_cycle,
                triangle_change=difference)


def cross_certificate(p, sets, eps):
    ell = len(sets)
    field = set(range(p))
    order = list(range(ell))
    chosen = sets
    if all(s != field for s in sets) and sum(map(len, sets)) < (p * ell + 1) // 2:
        order = [0] + list(range(ell - 1, 0, -1))
        chosen = [{(-d) % p for d in field - sets[v]}
                  for v in order[1:] + order[:1]]
    # Dynamic programming stores a literal boundary choice for every sum.
    witnesses = {0: []}
    for s in chosen:
        following = {}
        for total, word in witnesses.items():
            for d in sorted(s):
                following[(total + d) % p] = word + [d]
        witnesses = following
    need(len(witnesses) == p, "chosen boundary sumset is not the prime field")
    ds = witnesses[sum(eps) % p]
    starts = [0]
    for i, d in zip(order, ds):
        starts.append((starts[-1] + d - eps[i]) % p)
    need(starts[-1] == starts[0], "boundary certificate does not close")
    for j, d in enumerate(ds):
        need(d in chosen[j], "chosen cross difference missing")
        need((starts[j+1] - starts[j] + eps[order[j]]) % p == d,
             "boundary direction sign error")
    return order, starts[:-1], ds


def make_rotation_tournament(p, sets, eps):
    ell = len(sets)
    n = p * ell
    adj = [0] * n
    for a in range(n):
        i, j = a % ell, a // ell
        for b in range(a + 1, n):
            k, h = b % ell, b // ell
            if i == k:
                connection = set(range(1, (p + 1) // 2))
                if eps[i] == -1:
                    connection.remove(1)
                    connection.add(p - 1)
                forward = (h - j) % p in connection
            elif k == (i + 1) % ell:
                forward = (h - j) % p in sets[i]
            elif i == (k + 1) % ell:
                forward = (j - h) % p not in sets[k]
            else:
                forward = i < k
            u, v = (a, b) if forward else (b, a)
            adj[u] |= 1 << v
    return adj


def cross_controls():
    result = []
    for p, ell in ((3, 3), (3, 4), (5, 3)):
        nonempty = [{d for d in range(p) if mask >> d & 1}
                    for mask in range(1, 1 << p)]
        tested = literal = 0
        for sets in product(nonempty, repeat=ell):
            for eps in product((-1, 1), repeat=ell):
                order, starts, ds = cross_certificate(p, sets, eps)
                if tested % 101 == 0:
                    adj = make_rotation_tournament(p, sets, eps)
                    cycle = tuple(i + ell * ((start + j * eps[i]) % p)
                                  for i, start in zip(order, starts)
                                  for j in range(p))
                    need(len(set(cycle)) == p * ell, "orbit blocks do not span")
                    need(all(adj[a] >> b & 1 for a, b in
                             zip(cycle, cycle[1:] + cycle[:1])),
                         "literal constructed cycle has a missing arc")
                    rotation = lambda x: (x + ell) % (p * ell)
                    need(all((adj[a] >> b & 1) ==
                             (adj[rotation(a)] >> rotation(b) & 1)
                             for a in range(p * ell) for b in range(p * ell)
                             if a != b), "constructed rotation is not an automorphism")
                    literal += 1
                tested += 1
        result.append(dict(p=p, orbits=ell, boundary_systems=len(nonempty)**ell,
                           signed_certificates=tested, literal_tournaments=literal))
    return result


def census_controls():
    results = []
    certificates = []
    for n in (8, 9):
        if n == 8:
            path = ROOT / "04-computation/experiments/tiling_line_metagraph_census_20261020_classes_n8.tsv"
            stream = path.open(encoding="utf-8")
        else:
            path = HERE.parent / "blue_quotient_20261020b/n9_classes.tsv.gz"
            stream = gzip.open(path, "rt", encoding="utf-8")
        classes = strong = 0
        ranks = Counter()
        no_equal = []
        extremal_tested = 0
        extremal_failures = []
        max_rank_example = None
        with stream:
            for record in csv.DictReader(stream, delimiter="\t"):
                code = int(record["canon" if n == 8 else "code"])
                sc = int(record["SC"]) if n == 8 else code == int(record["converse_code"])
                if not sc:
                    continue
                classes += 1
                adj = rows(code, n)
                if not strongly_connected(adj):
                    continue
                strong += 1
                score = [r.bit_count() for r in adj]
                found = None
                equal = False
                regular = min(score) == max(score)
                extremal = regular
                rank = 0
                for cycle in hamiltonian_cycles(adj):
                    rank += 1
                    equal |= any(score[cycle[i-1]] == score[cycle[i]] for i in range(n))
                    extremal |= any(score[cycle[i]] == min(score) and
                                    score[cycle[i-1]] < max(score) for i in range(n))
                    if found is None:
                        good = descent(adj, cycle)
                        if good is not None:
                            found = (rank, good)
                    if found is not None and equal and extremal:
                        break
                need(found is not None, f"no black closed descent at n={n}, code={code}")
                ranks[found[0]] += 1
                if not equal:
                    no_equal.append(code)
                if not regular:
                    extremal_tested += 1
                    if not extremal:
                        extremal_failures.append(code)
                if max_rank_example is None or found[0] > max_rank_example["cycle_rank"]:
                    max_rank_example = dict(code=code, cycle_rank=found[0],
                                            path=found[1], scores=score)
                certificates.append(dict(n=n, code=code, cycle_rank=found[0],
                                         path=found[1],
                                         triangle_change=check_descent_certificate(adj, found[1])))
        results.append(dict(n=n, self_converse_classes=classes,
                            strong_self_converse_classes=strong,
                            descent_cycle_rank_histogram=dict(sorted(ranks.items())),
                            strongest_rank_example=max_rank_example,
                            extremal_predecessor_tested_classes=extremal_tested,
                            extremal_predecessor_failures=extremal_failures,
                            no_equal_endpoint_closed_path_classes=no_equal))
    (HERE / "cycles_certificates.json").write_text(
        json.dumps(certificates, indent=2) + "\n", encoding="utf-8")
    return results


def smaller_extremal_controls():
    results = []
    for n in (4, 6):
        path = ROOT / f"04-computation/experiments/tiling_line_metagraph_census_20261020_classes_n{n}.tsv"
        tested = 0
        failures = []
        with path.open(encoding="utf-8") as stream:
            for record in csv.DictReader(stream, delimiter="\t"):
                if not int(record["SC"]):
                    continue
                code = int(record["canon"])
                adj = rows(code, n)
                if not strongly_connected(adj):
                    continue
                tested += 1
                scores = [r.bit_count() for r in adj]
                found = any(any(scores[c[i]] == min(scores) and
                                scores[c[i-1]] < max(scores) for i in range(n))
                            for c in hamiltonian_cycles(adj))
                if not found:
                    failures.append(code)
        results.append(dict(n=n, strong_self_converse_classes=tested,
                            extremal_predecessor_failures=failures))
    return results


def all_small_terminal_splices():
    results = []
    for n in (6, 8):
        path = ROOT / f"04-computation/experiments/tiling_line_metagraph_census_20261020_classes_n{n}.tsv"
        cycles = 0
        other_bad = Counter()
        bad = []
        with path.open(encoding="utf-8") as stream:
            for record in csv.DictReader(stream, delimiter="\t"):
                if not int(record["SC"]):
                    continue
                code = int(record["canon"])
                adj = rows(code, n)
                if not strongly_connected(adj):
                    continue
                for cycle in hamiltonian_cycles(adj):
                    cycles += 1
                    if descent(adj, cycle) is not None:
                        continue
                    B = [i for i in range(n) if blue(adj, cycle[i:]+cycle[:i])]
                    if len(B) != 2:
                        other_bad[len(B)] += 1
                        continue
                    witness = terminal_splice(adj, cycle)
                    witness["code"] = code
                    bad.append(witness)
        results.append(dict(n=n, all_rooted_cycles_tested=cycles,
                            two_blue_bad_cycles_repaired=len(bad),
                            other_bad_cycles_by_blue_cut_count=dict(other_bad),
                            witnesses=bad))
    return results


def fixed_cycle_hostiles():
    records = []
    for adj in ([2,12,25,17,35,15],
                [2,20,25,115,33,71,151,63],
                [98,92,409,273,225,206,140,267,115]):
        n = len(adj)
        cycle = tuple(range(n))
        need(all(adj[a] >> b & 1 for a,b in zip(cycle, cycle[1:]+cycle[:1])),
             "hostile word not a cycle")
        need(descent(adj, cycle) is None, "fixed cycle hostile has descent")
        blue_cuts = [i for i in range(n) if blue(adj, cycle[i:]+cycle[:i])]
        changed = next((i, descent(adj,c)) for i,c in enumerate(hamiltonian_cycles(adj),1)
                       if descent(adj,c) is not None)
        score = [r.bit_count() for r in adj]
        records.append(dict(n=n, adjacency=adj, scores=score, blue_cuts=blue_cuts,
                            changed_cycle_rank=changed[0], descent_path=changed[1]))
    return records


def main():
    report = dict(cross_set_controls=cross_controls(),
                  fixed_cycle_hostiles=fixed_cycle_hostiles(),
                  complete_small_terminal_splice_controls=all_small_terminal_splices(),
                  smaller_extremal_predecessor_controls=smaller_extremal_controls(),
                  complete_self_converse_class_controls=census_controls())
    text = json.dumps(report, indent=2) + "\n"
    (HERE / "cycles.out").write_text(text, encoding="utf-8")
    print(text, end="")


if __name__ == "__main__":
    main()
