#!/usr/bin/env python3
"""Scratch-first AP-centred H6 common-scale-30 carrier preflight."""

from collections import Counter
from itertools import combinations, product
from math import gcd, lcm, prod

import numpy as np


P = 13
C = 30
ORDERS = (1, 2, 3, 5, 6, 10, 15, 30)
UNITS = {
    d: ((0,) if d == 1 else tuple(u for u in range(1, d) if gcd(u, d) == 1))
    for d in ORDERS
}
FULL = (1 << C) - 1


def require(condition, message):
    if not condition:
        raise RuntimeError(message)


def centered(value, modulus):
    residue = value % modulus
    return residue - modulus if 2 * residue > modulus else residue


def crt_base(label, order, unit):
    if order == 1:
        return label
    step = ((unit - order * label) * pow(P, -1, order)) % order
    return (order * label + P * step) % (P * order)


def local_mask(label, order, unit, owner):
    base = crt_base(label, order, unit)
    owner_inverse = pow(owner, -1, P)
    mask = 0
    for sheet in range(C):
        value = centered(base * (owner_inverse + P * sheet), P * order)
        if -order < value <= order:
            mask |= 1 << sheet
    return mask


def analytic_cardinality(label, order, owner):
    ratio = label * pow(owner, -1, P) % P
    target = order * ratio % P
    one_period = sum(x % P == target for x in range(-order + 1, order + 1))
    return (C // order) * one_period


def build_tables():
    masks = {}
    cards = np.zeros((P, len(ORDERS), P), dtype=np.int16)
    card_profiles = {}
    for label in range(1, P):
        for order_index, order in enumerate(ORDERS):
            for unit in UNITS[order]:
                base = crt_base(label, order, unit)
                answers = [x for x in range(P * order)
                           if x % P == order * label % P and x % order == unit % order]
                require(answers == [base], "CRT mismatch")
                for owner in range(1, P):
                    mask = local_mask(label, order, unit, owner)
                    card = analytic_cardinality(label, order, owner)
                    require(mask.bit_count() == card, "mask/card mismatch")
                    masks[label, order, unit, owner] = mask
                    cards[label, order_index, owner] = card
                    for sheet in range(C):
                        require(((mask >> sheet) & 1)
                                == ((mask >> ((sheet + order) % C)) & 1),
                                "periodicity mismatch")
    for order in ORDERS:
        card_profiles[order] = tuple(
            analytic_cardinality(ratio, order, 1) for ratio in range(1, P)
        )
    return masks, cards, card_profiles


def hereditary(index_word):
    word = tuple(ORDERS[index] for index in index_word)
    return all(lcm(*(word[j] for j in range(6) if j != omitted)) == C
               for omitted in range(6))


def grammar():
    words = []
    state_words = 0
    multiplicities = Counter()
    for index_word in product(range(len(ORDERS)), repeat=6):
        word = tuple(ORDERS[index] for index in index_word)
        via_lcm = hereditary(index_word)
        via_valuations = all(
            sum(order % prime == 0 for order in word) >= 2
            for prime in (2, 3, 5)
        )
        require(via_lcm == via_valuations, "valuation grammar mismatch")
        if via_lcm:
            words.append(index_word)
            state_words += prod(len(UNITS[d]) for d in word)
            multiplicities[tuple(word.count(d) for d in ORDERS)] += 1
    return np.asarray(words, dtype=np.int8), state_words, multiplicities


def scalar_census(words, cards):
    survivors = []
    owner_hist = np.zeros(7, dtype=np.int64)
    support_hist = Counter()
    multiplicities = Counter()
    capacity_vectors = set()
    for support_tuple in combinations(range(1, P), 6):
        support = np.asarray(support_tuple, dtype=np.int8)
        capacities = np.zeros((len(words), 6), dtype=np.int16)
        for provider, label in enumerate(support):
            capacities += cards[label, words[:, provider, None], support[None, :]]
        live = (capacities >= C).sum(axis=1)
        owner_hist += np.bincount(live, minlength=7)
        indices = np.flatnonzero(live == 6)
        support_hist[len(indices)] += 1
        for index in indices:
            word = tuple(ORDERS[int(i)] for i in words[index])
            capacity = tuple(int(x) for x in capacities[index])
            survivors.append((support_tuple, word, capacity))
            multiplicities[tuple(word.count(d) for d in ORDERS)] += 1
            capacity_vectors.add(capacity)
    return {
        "survivors": tuple(survivors),
        "owner_hist": tuple(int(x) for x in owner_hist),
        "support_hist": support_hist,
        "multiplicities": multiplicities,
        "capacity_vectors": capacity_vectors,
    }


def anchor_union_bank(support, word, owner, anchor_orders, masks):
    bank = frozenset((0,))
    for label, order in zip(support, word):
        if order not in anchor_orders:
            continue
        options = frozenset(masks[label, order, u, owner] for u in UNITS[order])
        bank = frozenset(partial | option for partial in bank for option in options)
    return bank


def relaxation_bound(support, word, owner, anchor_orders, bank, masks):
    best = 0
    for qmask in bank:
        bound = qmask.bit_count()
        for label, order in zip(support, word):
            if order in anchor_orders:
                continue
            bound += max(
                (masks[label, order, unit, owner] & ~qmask).bit_count()
                for unit in UNITS[order]
            )
        best = max(best, bound)
    return best


def flag_census(survivors, masks, flags):
    summaries = {}
    rows = {}
    for name, anchors in flags.items():
        owner_bound_hist = Counter()
        bank_hist = Counter()
        context_live_hist = Counter()
        hard_rows = []
        for support, word, capacity in survivors:
            bounds = []
            for owner_index, owner in enumerate(support):
                bank = anchor_union_bank(support, word, owner, anchors, masks)
                bound = relaxation_bound(support, word, owner, anchors, bank, masks)
                require(bound <= capacity[owner_index], "bound exceeds scalar capacity")
                bounds.append(bound)
                owner_bound_hist[bound] += 1
                bank_hist[len(bank)] += 1
            live = sum(bound >= C for bound in bounds)
            context_live_hist[live] += 1
            if live == 6:
                hard_rows.append((support, word, capacity, tuple(bounds)))
            rows[name, support, word] = tuple(bounds)
        summaries[name] = {
            "owner_bound_hist": owner_bound_hist,
            "bank_hist": bank_hist,
            "context_live_hist": context_live_hist,
            "hard_rows": tuple(hard_rows),
        }
    return rows, summaries


def fmt(counter):
    return " ".join(f"{k}:{counter[k]}" for k in sorted(counter))


def main():
    masks, cards, profiles = build_tables()
    words, state_words, grammar_mult = grammar()
    scalar = scalar_census(words, cards)
    base_flags = {
        "mod2": frozenset((1, 2)),
        "mod3": frozenset((1, 3)),
        "mod5": frozenset((1, 5)),
        "base23": frozenset((1, 2, 3)),
        "base25": frozenset((1, 2, 5)),
        "base35": frozenset((1, 3, 5)),
        "base235": frozenset((1, 2, 3, 5)),
    }
    _rows, flags = flag_census(scalar["survivors"], masks, base_flags)
    print("orders", ORDERS, "units", tuple(len(UNITS[d]) for d in ORDERS))
    print("card profiles by ratio 1..12")
    for order in ORDERS:
        print(order, profiles[order])
    print("hereditary words", len(words), "state words/support", state_words,
          "raw labelled", 924 * state_words, "grammar mult patterns", len(grammar_mult))
    print("scalar contexts", 924 * len(words), "owner hist", scalar["owner_hist"])
    print("support survivor hist", fmt(scalar["support_hist"]))
    print("survivors", len(scalar["survivors"]), "supports",
          len({s for s, _, _ in scalar["survivors"]}), "capacity vectors",
          len(scalar["capacity_vectors"]), "mult patterns", len(scalar["multiplicities"]))
    for name in base_flags:
        item = flags[name]
        print(name, "bank hist", fmt(item["bank_hist"]))
        print(name, "owner bound hist", fmt(item["owner_bound_hist"]))
        print(name, "context live hist", fmt(item["context_live_hist"]),
              "all-live rows", len(item["hard_rows"]),
              "hard mult patterns", len({
                  tuple(word.count(d) for d in ORDERS)
                  for _support, word, _capacity, _bounds in item["hard_rows"]
              }))


if __name__ == "__main__":
    main()
