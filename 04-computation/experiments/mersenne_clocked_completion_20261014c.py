"""Exact fixed-head coincidence compiler and necessary phase-cycle clock bills.

No search for a ROOT certificate occurs in the public APIs. Main has a small
declared control universe; algebraic coincidences are separately authenticated.
"""
from fractions import Fraction
from itertools import product
from math import comb
import json
from collatz_fibre_integration_20261009 import fan_members
import wilson_defect_trees_20261014c as wilson

CHECKS = 0


def require(condition, message):
    global CHECKS
    CHECKS += 1
    if not condition:
        raise ValueError(message)


def integer(value, lower, name):
    require(type(value) is int and value >= lower, name)


def word(value):
    require(type(value) is tuple, "word must be a tuple")
    require(all(type(a) is int and a > 0 for a in value), "positive valuations")
    require(not value or value[0] >= 2, "head starts after the maximal ones run")
    return value


def carrier(w):
    word(w)
    p, q, b = 1, 1, 0
    for a in w:
        p, q, b = 3 * p, q * 2**a, 3 * b + q
    return p, q, b


def power_exponent(value, base):
    if type(value) is not int or value < 1:
        return None
    e = 0
    while value % base == 0:
        value //= base
        e += 1
    return e if value == 1 else None


def classify(drop, source_head, child_head):
    """Complete rational identity classification for E >= drop+2.

    An identity is not a native valuation guard. A singleton candidate is not
    an authenticated join. authenticate() checks that separate obligation.
    """
    integer(drop, 1, "positive exponent drop")
    pu, qu, bu = carrier(source_head)
    pv, qv, bv = carrier(child_head)
    c = qv * pu * 3**drop - qu * pv
    r = qu * (bv - pv) - qv * (bu - pu)
    defect = (sum(source_head) - sum(child_head),
              drop + len(source_head) - len(child_head))
    result = {"kind": "none", "C": c, "R": r,
              "normalized_defect": list(defect), "exponent": None}
    if c == 0:
        if r == 0:
            require(defect == (0, 0), "identity must have zero clock defect")
            result["kind"] = "identity"
        return result
    ratio = Fraction(r, 2 * c)
    k = power_exponent(ratio.numerator, 3) if ratio.denominator == 1 else None
    if k is not None and k >= 1:
        result.update(kind="singleton", exponent=drop + 1 + k)
    return result


def rational_endpoint(exponent, head):
    value = Fraction(2 * 3**(exponent - 1) - 1)
    for a in head:
        value = (3 * value + 1) / 2**a
    return value


def native_endpoint(exponent, head):
    integer(exponent, 2, "nontrivial Mersenne exponent")
    word(head)
    value = 2 * 3**(exponent - 1) - 1
    for a in head:
        require(value > 1, "prefix extends beyond first ROOT")
        numerator = 3 * value + 1
        require(numerator % 2**a == 0, "nonintegral step")
        value = numerator // 2**a
        require(value % 2 == 1, "wrong exact valuation")
    require(value > 1, "join endpoint must precede ROOT")
    return value


def authenticate(exponent, drop, source_head, child_head):
    integer(exponent, 2, "source exponent")
    integer(drop, 1, "drop")
    require(exponent >= drop + 2, "nontrivial child exponent")
    classification = classify(drop, source_head, child_head)
    require(classification["kind"] == "identity" or
            classification["exponent"] == exponent, "wrong algebraic source")
    left = native_endpoint(exponent, source_head)
    right = native_endpoint(exponent - drop, child_head)
    require(left == right, "endpoint mismatch")
    return {"source_exponent": exponent, "child_exponent": exponent - drop,
            "endpoint": left, "normalized_defect": classification["normalized_defect"]}


def cycle_bill(period, time_defect, odd_defect):
    """Necessary only: one repeating grounded phase cycle, constant defects."""
    integer(period, 1, "positive exponent period")
    require(type(time_defect) is int and type(odd_defect) is int, "integer clocks")
    return {"period": period, "time_defect": time_defect, "odd_defect": odd_defect,
            "passes_necessary_bill": odd_defect >= period and
            time_defect >= 0 and 2**time_defect > 3**period}


def odd_budget_witness(base_exponent, base_odd_rank, period, odd_defect):
    """First contradictory lift if a rooted phase claims too few odd steps."""
    integer(base_exponent, 2, "base exponent")
    integer(base_odd_rank, base_exponent, "valid initial odd-rank lower bound")
    integer(period, 1, "period")
    require(type(odd_defect) is int and odd_defect < period, "underpaid odd cycle")
    j = (base_odd_rank - base_exponent) // (period - odd_defect) + 1
    return {"lift_index": j, "source_exponent": base_exponent + j * period,
            "claimed_odd_rank": base_odd_rank + j * odd_defect}


def contracting_inverse_words(exponent, length, node_cap=10000000):
    """All native inverse words of this length and contracting slope.

    Exact residues suffice; never construct 2**exponent. An exhausted budget
    raises an error and is never returned as an empty complete search.
    """
    integer(exponent, 2, "source exponent")
    integer(length, 1, "positive word length")
    integer(node_cap, 1, "explicit node cap")
    powers = [3**j for j in range(length + 1)]
    cap = powers[-1].bit_length() - 1
    nodes = 0
    hits = []

    def visit(residue, remaining, cost, backwards):
        nonlocal nodes
        nodes += 1
        if nodes > node_cap:
            raise ValueError("incomplete inverse search: node budget exhausted")
        if remaining == 0:
            hits.append(tuple(reversed(backwards)))
            return
        if residue % 3 == 0:
            return
        for a in range(1, cap - cost - (remaining - 1) + 1):
            numerator = 2**a * residue - 1
            if numerator % 3 == 0:
                visit((numerator // 3) % powers[remaining - 1], remaining - 1,
                      cost + a, backwards + (a,))

    visit((pow(2, exponent, powers[-1]) - 1) % powers[-1], length, 0, ())
    return {"length": length, "cost_cap": cap, "full_word_count": comb(cap, length),
            "native_prefix_nodes": nodes, "hits": hits}


def compositions(total, length):
    if length == 1:
        yield (total,)
    else:
        for first in range(1, total - length + 2):
            for tail in compositions(total - first, length - 1):
                yield (first,) + tail


def direct_inverse_guards(exponent, length):
    """Independent exhaustive carrier path, with no reduced-head restriction."""
    cap = (3**length).bit_length() - 1
    hits = []
    for total in range(length, cap + 1):
        for w in compositions(total, length):
            p, q, b = 1, 1, 0
            for a in w:
                p, q, b = 3*p, 2**a*q, 3*b+q
            require(p == 3**length and q == 2**total and q < p,
                    "complete contracting carrier universe")
            if pow(2, exponent + total, p) == (q + b) % p:
                hits.append(w)
    return hits


def rejects(action):
    try:
        action()
    except ValueError:
        return
    raise ValueError("hostile unexpectedly accepted")


def main():
    # 36 reduced/empty heads, all pairs, drops1..3; equality independently
    # evaluated as successive Fraction maps for E=drop+2,...,24.
    heads = [()]
    for length in (1, 2):
        for first in range(2, 7):
            for tail in product(range(1, 7), repeat=length - 1):
                heads.append((first,) + tail)
    counts = {"identity": 0, "singleton": 0, "none": 0}
    singletons, native, algebraic_hostiles = [], [], []
    for drop in range(1, 4):
        for u in heads:
            for v in heads:
                result = classify(drop, u, v)
                counts[result["kind"]] += 1
                for e in range(drop + 2, 25):
                    observed = rational_endpoint(e, u) == rational_endpoint(e - drop, v)
                    predicted = result["kind"] == "identity" or result["exponent"] == e
                    require(observed == predicted, "independent affine evaluation")
                if result["kind"] == "singleton":
                    e = result["exponent"]
                    require(rational_endpoint(e, u) == rational_endpoint(e - drop, v),
                            "unique candidate independently evaluated")
                    row = {"drop": drop, "source_head": list(u), "child_head": list(v),
                           "exponent": e, "normalized_defect": result["normalized_defect"]}
                    singletons.append(row)
                    try:
                        receipt = authenticate(e, drop, u, v)
                        native.append(dict(row, endpoint=receipt["endpoint"]))
                    except ValueError:
                        algebraic_hostiles.append(row)
    positive = authenticate(3, 1, (2, 3), ())
    require(positive["normalized_defect"] == [5, 3] and positive["endpoint"] == 5,
            "inherited seven/three nonzero-defect control")
    require(classify(1, (3,), (2, 1))["kind"] == "identity", "empty native phase identity")
    for e in range(3, 49):
        rejects(lambda e=e: authenticate(e, 1, (3,), (2, 1)))
    for e in range(4, 100, 2):
        valuation = power_exponent(e & -e, 2)
        reset = 3 + valuation
        receipt = authenticate(e, 1, (reset,), (2, reset - 2))
        require(receipt["normalized_defect"] == [0, 0], "inherited even reset")
    require(bool(algebraic_hostiles), "native guard hostile present")
    sample = algebraic_hostiles[0]
    rejects(lambda: authenticate(sample["exponent"], sample["drop"],
                                 tuple(sample["source_head"]), tuple(sample["child_head"])))
    hostiles = [
        lambda: classify(True, (), ()), lambda: classify(0, (), ()),
        lambda: classify(1, [2], ()), lambda: classify(1, (True,), ()),
        lambda: classify(1, (1,), ()), lambda: authenticate(True, 1, (2, 3), ()),
        lambda: authenticate(4, 1, (2, 3), ()),
        lambda: native_endpoint(2, (4,)),  # ROOT endpoint is not pre-ROOT.
        lambda: native_endpoint(2, (4, 2)),  # Never pad with ROOT self-loops.
        lambda: cycle_bill(True, 2, 1), lambda: cycle_bill(1, True, 1),
        lambda: odd_budget_witness(3, 2, 1, 0),
        lambda: odd_budget_witness(3, 5, 1, 1),
    ]
    for action in hostiles:
        rejects(action)
    bills = [cycle_bill(1, 0, 0), cycle_bill(1, 1, 1), cycle_bill(1, 2, 0),
             cycle_bill(1, 2, 1), cycle_bill(1, 5, 3), cycle_bill(2, 4, 2)]
    require([b["passes_necessary_bill"] for b in bills] ==
            [False, False, False, True, True, True], "independent two-coordinate budget")
    witness = odd_budget_witness(3, 5, 1, 0)
    require(witness == {"lift_index": 3, "source_exponent": 6, "claimed_odd_rank": 5},
            "concrete impossible zero-clock grounded extension")
    # Inherited labels only: no giant source or long absorption is replayed.
    fan = fan_members()
    require(len(fan) == 1949 and min(fan) == 99708991741, "inherited giant fan")
    fan_controls = []
    for phase, period, expected in ((323, 512, 4), (199, 256, 7), (79, 128, 15)):
        matches = sorted(e for e in fan if e % period == phase)
        escaped = [e - 2 for e in matches if e - 2 not in fan]
        require(len(matches) == expected and not escaped, "new phase stays in old fan")
        fan_controls.append({"phase": phase, "period": period, "source_matches": matches,
                             "children_outside_inherited_fan": escaped})
    marked_giant_sources = [99708991741, 99708993705, 99708993713]
    require(all(not any(e % row["period"] == row["phase"] for row in fan_controls)
                for e in marked_giant_sources), "marked giant sources still miss")
    # Two different trees have a prefix-preserving address bijection. Native
    # arithmetic is deliberately NOT inferred from that formal bijection.
    tree_addresses = 0
    for depth in range(7):
        sections = set()
        for digits in product(range(3), repeat=depth):
            address = tuple((d,) for d in digits)
            section = wilson.encode(3, address)
            require(wilson.decode(section) == address, "exact ternary address round trip")
            sections.add(section.residues)
            if digits:
                parent, last = wilson.unlift(section)
                require(parent == wilson.encode(3, address[:-1]) and last == address[-1],
                        "formal inverse-word prefix retained")
            tree_addresses += 1
        require(len(sections) == 3**depth, "all labelled sections at this depth")
    giant = min(fan)
    inverse_rows = [contracting_inverse_words(giant, length) for length in range(1, 65)]
    require(all(not row["hits"] for row in inverse_rows), "declared inverse search has no hit")
    require(sum(row["native_prefix_nodes"] for row in inverse_rows) == 28087511,
            "independently observed complete search node count")
    require(sum(row["full_word_count"] for row in inverse_rows) ==
            7818118225914629120771112104, "complete formal contracting-word universe")
    for length in range(1, 65):
        require(giant >= 3*length + 2, "noncontracting inverse words cannot pay at this source")
        require((3**length + 1) * (3**length - 2**length) < 2**(4*length + 1),
                "integer version of the sufficient carry bound")
    for row in inverse_rows[:12]:
        require(set(direct_inverse_guards(giant, row["length"])) == set(row["hits"]),
                "independent full-carrier enumeration agrees")
    for e in range(2, 33):
        for length in range(1, 7):
            require(set(contracting_inverse_words(e, length)["hits"]) ==
                    set(direct_inverse_guards(e, length)), "positive/hostile small-source comparison")
    require((1, 2) in contracting_inverse_words(5, 2)["hits"], "31 has the smaller child27")
    rejects(lambda: contracting_inverse_words(giant, 16, node_cap=1))
    rejects(lambda: contracting_inverse_words(True, 4))
    output = {"status": "PROVED scoped compiler and necessary cycle bills; universal ROOT OPEN",
              "checks": CHECKS, "universe": {"heads": len(heads), "drops": [1, 2, 3],
              "equality_control_max_exponent": 24}, "classification_counts": counts,
              "native_singletons": native, "nonnative_singleton_count": len(algebraic_hostiles),
              "first_nonnative_singleton": sample, "inherited_nonzero_control": positive,
              "cycle_bill_controls": bills, "underpaid_cycle_witness": witness,
              "inherited_giant_fan_controls": fan_controls,
              "marked_giant_sources_missing_all_three_new_guards": marked_giant_sources,
              "inverse_search_source_exponent": giant, "complete_inverse_search": inverse_rows,
              "inverse_search_exponent_period": 2*3**63,
              "formal_ternary_addresses_checked": tree_addresses,
              "hostile_api_cases": len(hostiles) + 1}
    print(json.dumps(output, indent=2, sort_keys=True))


if __name__ == "__main__":
    main()
