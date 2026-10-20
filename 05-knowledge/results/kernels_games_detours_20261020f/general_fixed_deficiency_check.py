"""Finite hostile controls for the new full-position owner and its bulk law."""
from collections import Counter
from functools import lru_cache
from importlib.util import module_from_spec, spec_from_file_location
from itertools import combinations, product
from pathlib import Path

spec = spec_from_file_location("port_kernel_controls", Path(__file__).resolve().parents[1] /
                              "rotation_owners_spectrum_20261020e/owners_check.py")
K = module_from_spec(spec)
spec.loader.exec_module(K)
O, S, H, need = K.O, K.S, K.H, K.need


def full_record(rows, weights, order, budget):
    edges = tuple((a, b) for a, b in combinations(range(len(order)), 2)
                  if rows[order[b]] >> order[a] & 1)
    if sum(weights[order[a]] * weights[order[b]] for a, b in edges) > budget:
        return None
    return budget - len(edges), tuple(weights[v] for v in order), edges


@lru_cache(None)
def new_owner(rows, weights, budget):
    heavy = sorted((v for v, w in enumerate(weights) if w == H),
                   key=lambda v: -sum(rows[v] >> z & 1 for z, w in enumerate(weights) if w == H))
    need(heavy, "controls require positive heavy count")
    gaps = [[] for _ in range(len(heavy) + 1)]
    for v, w in enumerate(weights):
        if w != H:
            bits = tuple(rows[z] >> v & 1 for z in heavy)
            i = sum(bits)
            need(bits == (1,) * i + (0,) * (len(heavy) - i), "fixed intrinsic light gap")
            gaps[i].append(v)
    cross = sum(weights[a] * weights[b]
                for i, j in combinations(range(len(gaps)), 2)
                for a in gaps[i] for b in gaps[j] if rows[b] >> a & 1)
    choices = [O.feasible_gap_orders(S.G.induced(rows, gap), tuple(weights[v] for v in gap),
                                    budget - cross) for gap in gaps]
    best = None
    for parts in product(*choices):
        if sum(cost for _, cost in parts) > budget - cross:
            continue
        order = []
        for i, (part, _) in enumerate(parts):
            order.extend(gaps[i][v] for v in part)
            if i < len(heavy):
                order.append(heavy[i])
        candidate = full_record(rows, weights, order, budget)
        need(candidate is not None, "local cost decomposition")
        if best is None or candidate < best:
            best = candidate
    need(best is not None, "new owner exists")
    return best


def is_new_owner(record, budget):
    rows, weights, edges = O.ordered_template(record[5], record[3], record[4])
    if len(O.interval_chains(len(rows), edges)) != len(rows):
        return False
    own = new_owner(rows, weights, budget)
    return full_record(rows, weights, tuple(range(len(rows))), budget) == own


def main():
    # A mandatory failed simplification: unguarded edge-list lex is not insertion-stable.
    need((1,) < (1, 3) and (1, 4) > (1, 3, 4), "proper-prefix hostile")
    results, examples, changed, total = {}, [], 0, 0
    for j in (2, 3):
        counts, sample = Counter(), Counter()
        for record in K.kernels(j):
            accepted = is_new_owner(record, j)
            old_accepted = K.is_owner(record, j)
            changed += accepted != old_accepted
            total += 1
            ell = record[4].count(3)
            if accepted:
                counts[ell] += 1
            category = record[:3], ell, accepted
            if sample[category] < 2:
                examples.append((j, record, accepted))
                sample[category] += 1
        results[j] = counts
        print("new-owner D3 j=", j, "coefficients=", tuple(counts[i] for i in range(4)), flush=True)
    n2 = results[2][0]
    n3 = 3 * sum(results[2].values()) + sum(results[3].values())
    need(n2 == 756 and n3 == 36473, "independent full-class small counts")
    print("old/new owner decisions compared=", total, "changed=", changed,
          "same full-class totals=", (n2, n3), flush=True)

    bulk, substituted, mixed = 0, 0, 0
    for j, record, accepted in examples:
        for extra in (1, 2):
            for child in K.completions(record, j, extra):
                need(is_new_owner(child, j + extra) == accepted, "new-owner arbitrary matching insertion")
                bulk += 1
        positions = [i for i, value in enumerate(record[4]) if value == 3]
        for values in product((3, 4, 7), repeat=len(positions)):
            gap = list(record[4])
            for p, value in zip(positions, values):
                gap[p] = value
            changed_record = record[:4] + (tuple(gap), record[5])
            need(is_new_owner(changed_record, j) == accepted, "independent large-light substitutions")
            substituted += 1
            if len(set(values)) > 1:
                mixed += 1
            # Joint operation, retaining actual matching ports rather than only gap types.
            for child in K.completions(changed_record, j, 1):
                need(is_new_owner(child, j + 1) == accepted, "joint L substitution and insertion")
                bulk += 1
    need(mixed > 0, "controls contain independent unequal L values")
    print("bulk matching and joint controls=", bulk, "L substitution controls=", substituted,
          "unequal-L controls=", mixed, flush=True)

    # D=0 empty-endpoint kernel, and its arbitrary two-edge matching completions.
    empty = (0, 0, 0, (), (H,), ())
    need(is_new_owner(empty, 0), "D0 empty kernel")
    zero_count = 0
    for child in K.completions(empty, 0, 2):
        need(is_new_owner(child, 2), "D0 completed owner")
        zero_count += 1
    need(zero_count == 3, "D0 two-edge matching count")
    print("D0 empty kernel and all three two-edge completions pass; prefix hostile retained", flush=True)
    print("PASS full-position owner, arbitrary matching ports, and independent large-light controls")


if __name__ == "__main__":
    main()
