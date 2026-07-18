#!/usr/bin/env python3
"""Exploratory scalar census for primitive AP-centred H6 scale c=28."""

from collections import Counter
from itertools import combinations, product
from math import gcd, lcm, prod
import time

import numpy as np

P = 13
C = 28
ORDERS = (1, 2, 4, 7, 14, 28)
PHI = tuple(1 if d == 1 else sum(gcd(u, d) == 1 for u in range(1, d)) for d in ORDERS)


def card(ratio, order):
    target = order * ratio % P
    return (C // order) * sum(x % P == target for x in range(-order + 1, order + 1))


def hereditary(word):
    return all(lcm(*(word[j] for j in range(6) if j != i)) == C for i in range(6))


def centered(x, modulus):
    x %= modulus
    return x - modulus if 2 * x > modulus else x


def literal_base(label, order, unit):
    for x in range(P * order):
        if x % P == order * label % P and x % order == unit % order:
            return x
    raise AssertionError


def local_mask(label, order, unit, owner):
    base = literal_base(label, order, unit)
    yi = pow(owner, -1, P)
    mask = 0
    for sheet in range(C):
        x = centered(base * (yi + P * sheet), P * order)
        if -order < x <= order:
            mask |= 1 << sheet
    return mask


def owner_local(support, word, owner, masks):
    reachable = {0}
    for label, order in zip(support, word):
        choices = masks[label, order, owner]
        reachable = {partial | choice for partial in reachable for choice in choices}
    maximum = max(x.bit_count() for x in reachable)
    return maximum == C, maximum, len(reachable)


def projection_capacity_flag(support, word, owner, masks, moduli):
    thresholds = tuple(C // modulus for modulus in moduli for _ in range(modulus))
    states = {tuple(0 for _ in thresholds)}
    for label, order in zip(support, word):
        options = set()
        for mask in masks[label, order, owner]:
            signature = []
            for modulus in moduli:
                signature.extend(
                    sum((mask >> sheet) & 1 for sheet in range(residue, C, modulus))
                    for residue in range(modulus)
                )
            options.add(tuple(signature))
        states = {
            tuple(min(bound, left + right) for bound, left, right in zip(thresholds, state, option))
            for state in states
            for option in options
        }
    return thresholds in states, len(states)


def main():
    t0 = time.time()
    cards = np.array([[card(r, d) for r in range(P)] for d in ORDERS], dtype=np.int8)
    words = []
    weighted = 0
    for iw in product(range(len(ORDERS)), repeat=6):
        word = tuple(ORDERS[i] for i in iw)
        h = hereditary(word)
        flags = sum(d % 4 == 0 for d in word) >= 2 and sum(d % 7 == 0 for d in word) >= 2
        assert h == flags
        if h:
            words.append(iw)
            weighted += prod(PHI[i] for i in iw)
    words = np.array(words, dtype=np.int8)
    print("ratio")
    for i, d in enumerate(ORDERS):
        print(d, tuple(map(int, cards[i, 1:])))
    print("grammar", len(words), weighted)

    feasible_hist = np.zeros(7, dtype=np.int64)
    support_hist = Counter()
    mult = Counter()
    cap_hist = Counter()
    min_slack = Counter()
    max_slack = Counter()
    survivors = []
    for support in combinations(range(1, P), 6):
        # ratio card tensor [provider, order index, owner index]
        rc = np.empty((6, len(ORDERS), 6), dtype=np.int16)
        for p, label in enumerate(support):
            for o in range(len(ORDERS)):
                for y, owner in enumerate(support):
                    rc[p, o, y] = cards[o, label * pow(owner, -1, P) % P]
        caps = np.zeros((len(words), 6), dtype=np.int16)
        for p in range(6):
            caps += rc[p, words[:, p], :]
        fc = (caps >= C).sum(axis=1)
        feasible_hist += np.bincount(fc, minlength=7)
        inds = np.flatnonzero(fc == 6)
        support_hist[len(inds)] += 1
        for k in inds:
            word = tuple(ORDERS[int(i)] for i in words[k])
            cv = tuple(int(x) for x in caps[k])
            survivors.append((support, word, cv))
            mult[tuple(word.count(d) for d in ORDERS)] += 1
            cap_hist[cv] += 1
            min_slack[min(cv) - C] += 1
            max_slack[max(cv) - C] += 1
    print("support/order rows", 924 * len(words))
    print("feasible_hist", dict(enumerate(map(int, feasible_hist))))
    print("survivors", len(survivors), "supports", sum(k > 0 for k in support_hist.elements()), "support_hist", sorted(support_hist.items()))
    print("multiplicity", sorted(mult.items()))
    print("capacity vectors", len(cap_hist), "min_slack", sorted(min_slack.items()), "max_slack", sorted(max_slack.items()))
    print("sample")
    for row in survivors[:30]: print(row)
    masks = {}
    for label in range(1, P):
        for order in ORDERS:
            units = (0,) if order == 1 else tuple(u for u in range(1, order) if gcd(u, order) == 1)
            for owner in range(1, P):
                options = frozenset(local_mask(label, order, u, owner) for u in units)
                assert {m.bit_count() for m in options} == {card(label * pow(owner, -1, P) % P, order)}
                masks[label, order, owner] = options
    owner_hist = Counter()
    max_hist = Counter()
    reach_hist = Counter()
    feasible_rows = []
    feasible_owner_keys = Counter()
    feasible_profile = Counter()
    row_profile_feasible = Counter()
    per_order_owner_max = Counter()
    flag_confusion = {mods: Counter() for mods in ((4,), (7,), (4, 7))}
    flag_state_hist = {mods: Counter() for mods in flag_confusion}
    for row_index, (support, word, cv) in enumerate(survivors):
        local = []
        for owner_index, owner in enumerate(support):
            z = owner_local(support, word, owner, masks)
            local.append(z)
            max_hist[z[1]] += 1
            reach_hist[z[2]] += 1
            per_order_owner_max[word[owner_index], z[1]] += 1
            for mods in flag_confusion:
                flag, states = projection_capacity_flag(support, word, owner, masks, mods)
                flag_confusion[mods][z[0], flag] += 1
                flag_state_hist[mods][states] += 1
        fc = sum(z[0] for z in local)
        owner_hist[fc] += 1
        profile = tuple(word.count(order) for order in ORDERS)
        row_profile_feasible[profile, fc] += 1
        for owner_index, (owner, z) in enumerate(zip(support, local)):
            if z[0]:
                oi = pow(owner, -1, P)
                key = tuple(sorted((label * oi % P, order) for label, order in zip(support, word)))
                feasible_owner_keys[key] += 1
                feasible_profile[profile] += 1
        if fc: feasible_rows.append((support, word, cv, tuple(local)))
        if row_index % 250 == 0: print("owner progress", row_index, "elapsed", time.time() - t0, flush=True)
    print("owner_hist", sorted(owner_hist.items()))
    print("max_hist", sorted(max_hist.items()))
    print("reach range", min(reach_hist), max(reach_hist), "distinct", len(reach_hist), "top", reach_hist.most_common(20))
    print("order-owner-max", sorted(per_order_owner_max.items()))
    print("projection capacity flags", {mods: sorted(hist.items()) for mods, hist in flag_confusion.items()})
    print("projection state ranges", {mods: (min(hist), max(hist), len(hist)) for mods, hist in flag_state_hist.items()})
    print("row-profile-feasible", sorted(row_profile_feasible.items()))
    print("feasible-profile-incidences", sorted(feasible_profile.items()))
    print("feasible-owner-keys", len(feasible_owner_keys))
    for key, count in sorted(feasible_owner_keys.items()): print("FKEY", key, count)
    print("feasible row sample")
    for row in feasible_rows[:30]: print(row)
    print("runtime", time.time() - t0)


if __name__ == "__main__":
    main()
