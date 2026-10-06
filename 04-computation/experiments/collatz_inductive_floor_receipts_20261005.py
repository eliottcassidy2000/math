"""Canonical descent obligations and independent floors on a two-step tree.

Production APIs accept no source ROOT word or measured moment values.
Missing or assumed leaves remain explicit. Only ROOT is an unconditional axiom.
"""
from dataclasses import dataclass, replace
from fractions import Fraction as F
import json
import collatz_floor_transport_deadlines_20261005 as floors

CHECKS = 0


def need(ok, message):
    global CHECKS
    CHECKS += 1
    if not ok:
        raise ValueError(message)


def odd(n):
    if type(n) is not int or n <= 0 or n % 2 != 1:
        raise ValueError("positive odd exact integer required")


def natural(n):
    if type(n) is not int or n < 0:
        raise ValueError("exact natural integer required")


@dataclass(frozen=True)
class Descent:
    source: int
    word: tuple
    endpoint: int


def validate(receipt):
    """Verify the unique first strict descent, not merely some later cut."""
    if type(receipt) is not Descent:
        raise ValueError("typed descent receipt required")
    odd(receipt.source)
    odd(receipt.endpoint)
    if receipt.source == 1 or type(receipt.word) is not tuple or not receipt.word:
        raise ValueError("nonempty descent from a nonroot source required")
    current = receipt.source
    for index, requested in enumerate(receipt.word):
        if type(requested) is not int or requested < 1:
            raise ValueError("positive exact valuation required")
        current, actual = floors.step(current)
        if actual != requested:
            raise ValueError("actual valuation guard failed")
        if index + 1 < len(receipt.word) and current < receipt.source:
            raise ValueError("the cut is later than the first descent")
    if current != receipt.endpoint or current >= receipt.source:
        raise ValueError("exact strictly smaller endpoint required")
    return receipt


def discover(source, cap):
    """Bounded local observations only; None means no descent seen in this cap."""
    odd(source)
    natural(cap)
    if source == 1:
        return None
    current, word = source, []
    for _ in range(cap):
        current, a = floors.step(current)
        word.append(a)
        if current < source:
            return validate(Descent(source, tuple(word), current))
    return None


def compile_floor(receipts, source, assumptions=None):
    """Follow only supplied decreasing receipts; never search a missing suffix.

    A numerical assumed floor stays conditional. Aggregate the whole prefix
    before scalarizing, so subdivision does not repeatedly weaken the floor.
    """
    odd(source)
    if type(receipts) is not dict:
        raise ValueError("finite receipt dictionary required")
    assumptions = {} if assumptions is None else assumptions
    if type(assumptions) is not dict:
        raise ValueError("explicit assumption dictionary required")
    for n, row in receipts.items():
        odd(n)
        validate(row)
        if n != row.source:
            raise ValueError("receipt key/source mismatch")
    for n, epsilon in assumptions.items():
        odd(n)
        floors.floor_parameter(epsilon)
        if n == 1 or F(epsilon) > F(1, 3):
            raise ValueError("nonroot point floor assumption required")
    current, total, cuts = source, [], []
    while current != 1 and current in receipts:
        row = receipts[current]
        total.extend(row.word)
        cuts.append(current)
        current = row.endpoint
    word = tuple(total)
    r, k = len(word), sum((a - 1)//2 for a in word)
    result = {"source": source, "anchor": current, "word": word,
              "cuts": tuple(cuts), "r": r, "k": k}
    if current == 1:
        result.update(status="grounded", floor=floors.word_weight(word),
                      obligations=())
    elif current in assumptions:
        result.update(status="conditional",
                      floor=floors.transport_floor(F(assumptions[current]), r, k),
                      obligations=((current, F(assumptions[current])),))
    else:
        result.update(status="pending", floor=None, obligations=(current,))
    return result


def palette(parent):
    """The two unit two-step inverse children; no ROOT suffix is queried."""
    odd(parent)
    if parent <= 1 or parent % 3 == 0:
        raise ValueError("nonroot unit parent required")
    a0 = next(a for a in range(3, 9) if ((1 << (a+1))*parent - 5) % 9 == 0)
    result = []
    for a in (a0, a0+6, a0+12):
        child = ((1 << (a+1))*parent-5)//9
        if child % 3:
            result.append(validate(Descent(child, (1, a), parent)))
    if len(result) != 2:
        raise ValueError("palette theorem failed")
    return tuple(result)


def recognize(source):
    """Total recognizer for the proper binary tree rooted at 5."""
    odd(source)
    current, rows = source, {}
    while current != 5:
        if current == 1 or current % 3 == 0:
            return {"member": False, "frontier": current, "receipts": rows}
        row = discover(current, 2)
        if row is None or len(row.word) != 2:
            return {"member": False, "frontier": current, "receipts": rows}
        parent = row.endpoint
        if parent <= 1 or parent % 3 == 0 or row not in palette(parent):
            return {"member": False, "frontier": current, "receipts": rows}
        rows[current] = row
        current = parent
    rows[5] = validate(Descent(5, (4,), 1))
    return {"member": True, "frontier": 5, "receipts": rows}


def source_floor(source):
    """A source-size formula, not an evaluation of its completed ROOT word."""
    report = recognize(source)
    if not report["member"]:
        return None
    b = source.bit_length()
    return floors.weight(4*b, 18*b+1)


def one_step_palette(parent):
    """The inherited SF5 palette, retained with its actual source guards."""
    odd(parent)
    if parent <= 1 or parent % 3 == 0:
        raise ValueError("nonroot unit parent required")
    exponents = (2, 4, 6) if parent % 3 == 1 else (3, 5, 7)
    result = []
    for a in exponents:
        child = ((1 << a)*parent-1)//3
        if child % 3:
            result.append(validate(Descent(child, (a,), parent)))
    if len(result) != 2:
        raise ValueError("one-step palette theorem failed")
    return tuple(result)


def hybrid_palette(parent):
    return one_step_palette(parent) + palette(parent)


def recognize_hybrid(source, only_one_step=False):
    """Canonical parse by the first valuation; every accepted block decreases."""
    odd(source)
    if type(only_one_step) is not bool:
        raise ValueError("exact Boolean mode required")
    current, rows = source, {}
    while current != 5:
        if current == 1 or current % 3 == 0:
            return {"member": False, "frontier": current, "receipts": rows}
        endpoint, first = floors.step(current)
        if first == 1:
            if only_one_step:
                return {"member": False, "frontier": current, "receipts": rows}
            row = discover(current, 2)
            native = palette
        else:
            row = validate(Descent(current, (first,), endpoint))
            native = one_step_palette
        if row is None:
            return {"member": False, "frontier": current, "receipts": rows}
        parent = row.endpoint
        if parent <= 1 or parent % 3 == 0 or row not in native(parent):
            return {"member": False, "frontier": current, "receipts": rows}
        rows[current] = row
        current = parent
    rows[5] = validate(Descent(5, (4,), 1))
    return {"member": True, "frontier": 5, "receipts": rows}


def least_leaf(source):
    odd(source)
    if source <= 1 or source % 3 == 0:
        raise ValueError("nonroot unit required")
    a = next(a for a in range(1, 7) if ((1 << a)*source-1) % 9 == 0)
    leaf = ((1 << a)*source-1)//3
    if leaf % 6 != 3 or floors.step(leaf) != (source, a):
        raise ValueError("designated leaf guard failed")
    return leaf, a


def localization_degree(epsilon):
    """Least d with 8*(16/25)^d <= epsilon/2, by exact integer arithmetic."""
    epsilon = floors.floor_parameter(epsilon)
    left, right, d = 16*epsilon.denominator, epsilon.numerator, 0
    while left > right:
        left *= 16
        right *= 25
        d += 1
    return d


def measurement_receipt(source):
    """Prove a positive exact localized readout; no moment values are inputs.

    The conclusion concerns the actual lambda measure via the inherited
    universal selector error, not an arbitrary caller-supplied distribution.
    """
    report = recognize(source)
    if not report["member"]:
        return None
    b = source.bit_length()
    leaf, a = least_leaf(source)
    epsilon = floors.weight(4*b+1, 18*b+3)
    d = localization_degree(epsilon)
    return {"source": source, "leaf": leaf, "leaf_valuation": a,
            "index": (leaf-3)//6, "bit_length": b,
            "atom_floor": epsilon, "degree": d,
            "exact_readout_floor": epsilon/2,
            "moment_values_consumed": 0, "target_ROOT_words_consumed": 0}


def hybrid_measurement_receipt(source):
    report = recognize_hybrid(source)
    if not report["member"]:
        return None
    b = source.bit_length()
    leaf, a = least_leaf(source)
    epsilon = floors.weight(6*b+1, 27*b+3)
    return {"source": source, "leaf": leaf, "leaf_valuation": a,
            "index": (leaf-3)//6, "bit_length": b,
            "atom_floor": epsilon, "degree": localization_degree(epsilon),
            "exact_readout_floor": epsilon/2,
            "moment_values_consumed": 0, "target_ROOT_words_consumed": 0}


def literal_root(source):
    """Independent experiment-only replay, never called by production APIs."""
    odd(source)
    n, word = source, []
    while n != 1:
        n, a = floors.step(n)
        word.append(a)
        if len(word) > 4096:
            raise ValueError("audit cap")
    return tuple(word)


def rejected(call):
    try:
        call()
    except (ValueError, TypeError):
        return True
    return False


def main():
    report = {}
    universe = tuple(range(1, 512, 2))
    census = []
    for cap in (1, 2, 4, 8, 16):
        rows = {}
        observed = 0
        for n in universe[1:]:
            row = discover(n, cap)
            if row is not None:
                rows[n] = row
                observed += len(row.word)
            else:
                observed += cap
        tokens = [compile_floor(rows, n) for n in universe]
        grounded = [x for x in tokens if x["status"] == "grounded"]
        pending = sorted({x["anchor"] for x in tokens if x["status"] == "pending"})
        for token in grounded:
            need(floors.replay(token["source"], token["word"]) == 1,
                 "grounded receipt is a strict actual ROOT path")
            need(token["word"] == literal_root(token["source"]),
                 "unique canonical word independently replayed")
        census.append({"cap": cap, "sources": len(universe),
                       "local_receipts": len(rows), "grounded": len(grounded),
                       "pending_anchor_count": len(pending),
                       "first_pending_anchors": pending[:8],
                       "literal_queries": observed})
    report["finite_induction"] = census

    rows = {7: discover(7, 4)}
    pending = compile_floor(rows, 7)
    conditional = compile_floor(rows, 7, {5: F(1, 3)})
    need(pending["status"] == "pending" and pending["anchor"] == 5,
         "unpaid endpoint remains explicit")
    need(conditional["status"] == "conditional" and
         conditional["floor"] > 0 and conditional["obligations"] == ((5, F(1, 3)),),
         "a numerical premise never silently grounds itself")
    rows[5] = discover(5, 1)
    closed = compile_floor(rows, 7)
    need(closed["status"] == "grounded" and closed["floor"] == F(1, 84),
         "ROOT-only discharge produces the exact weight")
    report["seven"] = {"descent": [list(rows[7].word), 5],
                       "conditional_floor": str(conditional["floor"]),
                       "grounded_weight": str(closed["floor"])}

    levels, all_nodes = [5], {5}
    level_counts, max_bits = [], 0
    measurements = []
    for depth in range(9):
        level_counts.append(len(levels))
        for n in levels:
            recognition = recognize(n)
            need(recognition["member"], "generated member recognized")
            token = compile_floor(recognition["receipts"], n)
            audit_word = literal_root(n)
            need(token["word"] == audit_word, "two-step tree canonical proof")
            need(len(audit_word) == 2*depth+1, "tree depth/odd rank")
            k = sum((a-1)//2 for a in audit_word)
            need(k <= 9*depth+1, "uniform sibling-cost bound")
            b = n.bit_length()
            max_bits = max(max_bits, b)
            need(depth <= 2*b, "size controls generation depth")
            need(source_floor(n) <= token["floor"], "analytic source floor")
            if depth:
                need(audit_word[0] == 1 and n % 8 == 3,
                     "every new node rises first; disjoint from SF5")
            if depth <= 3:
                mr = measurement_receipt(n)
                d, epsilon = mr["degree"], mr["atom_floor"]
                need(8*F(16, 25)**d <= epsilon/2, "localization budget")
                need(d == 0 or 8*F(16, 25)**(d-1) > epsilon/2,
                     "least sufficient degree")
                leaf_word = (mr["leaf_valuation"],) + audit_word
                need(floors.replay(mr["leaf"], leaf_word) == 1,
                     "independent audit of designated leaf")
                need(floors.word_weight(leaf_word) >= epsilon,
                     "analytic leaf bound versus literal audit")
                measurements.append({"source": n, "leaf": mr["leaf"],
                                     "degree": d, "floor_denominator_bits":
                                     epsilon.denominator.bit_length()})
        if depth < 8:
            next_level = [row.source for n in levels for row in palette(n)]
            need(len(set(next_level)) == len(next_level), "distinct children")
            need(not all_nodes.intersection(next_level), "no merged generations")
            all_nodes.update(next_level)
            levels = next_level
    need(tuple(row.source for row in palette(5)) == (35, 2275), "first two children")
    need(not recognize(7)["member"] and not recognize(27)["member"],
         "known rooted sources outside the proper tree")
    need(measurement_receipt(7) is None, "family rejection is not a global floor")
    report["binary_family"] = {"levels": level_counts, "sources": len(all_nodes),
                               "max_source_bits": max_bits,
                               "first_children": [35, 2275],
                               "localized_measurement_controls": measurements}

    levels, mixed_nodes = [5], {5}
    mixed_counts = []
    for depth in range(6):
        mixed_counts.append(len(levels))
        for n in levels:
            recognition = recognize_hybrid(n)
            need(recognition["member"], "hybrid generated member recognized")
            token = compile_floor(recognition["receipts"], n)
            need(token["word"] == literal_root(n), "hybrid canonical proof")
            need(len(token["cuts"]) == depth+1, "hybrid macro depth")
            b = n.bit_length()
            need(depth <= 3*b, "mixed growth controls macro depth")
            need(len(token["word"]) <= 2*depth+1, "hybrid odd rank")
            need(token["k"] <= 9*depth+1, "hybrid cost")
            need(floors.weight(6*b, 27*b+1) <= token["floor"],
                 "hybrid analytic source floor")
        if depth < 5:
            next_level = [r.source for n in levels for r in hybrid_palette(n)]
            need(len(set(next_level)) == len(next_level), "four distinct child branches")
            need(not mixed_nodes.intersection(next_level), "hybrid generations disjoint")
            mixed_nodes.update(next_level)
            levels = next_level
    mixed = recognize_hybrid(739)
    need(mixed["member"] and not recognize(739)["member"] and
         not recognize_hybrid(739, True)["member"], "strict enlargement of pure-tree union")
    need(compile_floor(mixed["receipts"], 739)["word"] == (1, 8, 3, 4),
         "mixed witness exact word")
    mm = hybrid_measurement_receipt(739)
    need(8*F(16, 25)**mm["degree"] <= mm["atom_floor"]/2,
         "hybrid positive readout budget")
    need(floors.word_weight(literal_root(mm["leaf"])) >= mm["atom_floor"],
         "independent hybrid leaf bound")
    report["hybrid_family"] = {"levels": mixed_counts, "sources": len(mixed_nodes),
                              "mixed_witness": 739, "word": [1, 8, 3, 4],
                              "leaf": mm["leaf"], "degree": mm["degree"]}

    hostiles = (
        lambda: discover(True, 2), lambda: discover(3.0, 2),
        lambda: discover(3, False), lambda: validate(Descent(1, (2,), 1)),
        lambda: validate(Descent(3, (1,), 5)),
        lambda: validate(Descent(7, (1, 1, 2, 3, 4), 1)),
        lambda: validate(replace(rows[7], endpoint=3)),
        lambda: validate(replace(rows[7], word=(1.0, 1, 2, 3))),
        lambda: compile_floor({7: rows[5]}, 7),
        lambda: compile_floor({}, 7, {7: F(1, 2)}),
        lambda: palette(3), lambda: measurement_receipt(False),
        lambda: least_leaf(1), lambda: localization_degree(F(0)),
        lambda: recognize_hybrid(739, 1), lambda: one_step_palette(True),
    )
    for call in hostiles:
        need(rejected(call), "typed, source, first-hit, rank or floor hostile rejected")
    report["hostiles"] = len(hostiles)
    report["checks"] = CHECKS
    report["scope"] = ("Infinite proper induction domain and positive exact measurement "
                       "guarantees; universal descent-template coverage remains OPEN.")
    print(json.dumps(report, indent=2, sort_keys=True))


if __name__ == "__main__":
    main()
