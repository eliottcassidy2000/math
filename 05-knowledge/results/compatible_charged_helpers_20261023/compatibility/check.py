"""Exact compatibility intervals, neutral-fibre quotient, and native controls."""
from __future__ import annotations
import argparse
from hashlib import sha256
from itertools import combinations, permutations, product
import json
from pathlib import Path
import random

HERE = Path(__file__).resolve().parent
ROOT = HERE.parents[3]
CHECKS = 0


def need(value, message):
    global CHECKS
    CHECKS += 1
    if not value:
        raise RuntimeError(message)


def submasks(mask):
    out = []
    sub = mask
    while True:
        out.append(sub)
        if not sub:
            return out[::-1]
        sub = (sub-1) & mask


class Poset:
    def __init__(self, n, edges):
        self.n = n
        self.all = (1 << n)-1
        self.up_rows = [(1 << i) for i in range(n)]
        for i, j in edges:
            self.up_rows[i] |= 1 << j
        for k in range(n):
            for i in range(n):
                if self.up_rows[i] & (1 << k):
                    self.up_rows[i] |= self.up_rows[k]
        self.down_rows = [sum(1 << i for i in range(n)
                              if self.up_rows[i] & (1 << j)) for j in range(n)]
        self.upsets = [x for x in range(1 << n) if self.up(x) == x]

    def up(self, mask):
        result = 0
        for i in range(self.n):
            if mask & (1 << i):
                result |= self.up_rows[i]
        return result

    def down(self, mask):
        result = 0
        for i in range(self.n):
            if mask & (1 << i):
                result |= self.down_rows[i]
        return result

    def interval(self, forced_one, forced_zero):
        low = self.up(forced_one)
        high = self.all ^ self.down(forced_zero)
        return None if low & ~high else (low, high)


def natural_posets(n):
    pairs = list(combinations(range(n), 2))
    seen = {}
    for mask in range(1 << len(pairs)):
        p = Poset(n, [e for i, e in enumerate(pairs) if mask & (1 << i)])
        seen[tuple(p.up_rows)] = p
    return [seen[k] for k in sorted(seen)]


def weights_table(weights):
    return [sum(w for i, w in enumerate(weights) if mask & (1 << i))
            for mask in range(1 << len(weights))]


def interval_members(p, interval):
    low, high = interval
    return [u for u in p.upsets if low & ~u == 0 and u & ~high == 0]


def distance_and_witness(a, b, wt):
    la, ra = a
    lb, rb = b
    wa = la | (lb & ra)
    wb = lb | (la & rb)
    return wt[la & ~rb] + wt[lb & ~ra], wa, wb


def components(count, edges):
    parent = list(range(count))
    def find(i):
        while parent[i] != i:
            parent[i] = parent[parent[i]]
            i = parent[i]
        return i
    for i, j in edges:
        a, b = find(i), find(j)
        if a != b:
            parent[b] = a
    return [find(i) for i in range(count)]


def interval_tests():
    rows = []
    for n in range(1, 5):
        count_unary = count_pairs = 0
        pp = natural_posets(n)
        for p in pp:
            for forced_one in range(1 << n):
                for forced_zero in range(1 << n):
                    compiled = p.interval(forced_one, forced_zero)
                    literal = [u for u in p.upsets
                               if not forced_one & ~u and not forced_zero & u]
                    expected = [] if compiled is None else interval_members(p, compiled)
                    need(literal == expected, "Unary compatibility closure mismatch")
                    count_unary += 1
            intervals = [(a, b) for a in p.upsets for b in p.upsets if not a & ~b]
            members = {a: interval_members(p, a) for a in intervals}
            # Zero weights specifically test zero distance != nonempty intersection.
            for weights in [(1,)*n, tuple((2*i) % 3 for i in range(n))]:
                wt = weights_table(weights)
                for a in intervals:
                    for b in intervals:
                        bound, x, y = distance_and_witness(a, b, wt)
                        literal = min(wt[u ^ v] for u in members[a] for v in members[b])
                        need(bound == literal, "Interval weighted distance mismatch")
                        need(x in members[a] and y in members[b], "Boundary witnesses outside fibre")
                        need(wt[x ^ y] == bound, "Boundary witnesses not optimal")
                        intersection = bool(set(members[a]) & set(members[b]))
                        la, ra = a; lb, rb = b
                        need(intersection == (((la | lb) & ~(ra & rb)) == 0),
                             "Intersection compatibility mismatch")
                        count_pairs += 1
        rows.append({"elements": n, "naturally_labelled_posets": len(pp),
                     "unary_forcing_controls": count_unary,
                     "weighted_interval_pair_controls": count_pairs})
    return rows


def quotient_case(p, intervals, q, helper_charge, weights, budget, helper_distance=1,
                  group_moduli=(2,)):
    """Compare complete literal graph with independently compiled quotient."""
    n = p.n
    wt = weights_table(weights)
    def norm(v):
        if isinstance(v, int):
            v = (v,)
        if len(v) != len(group_moduli):
            raise ValueError("Charge dimension mismatch")
        return tuple(x % modulus if modulus else x for x, modulus in zip(v, group_moduli))
    def add(a, b):
        return norm(tuple(x+y for x, y in zip(a, b)))
    zero = (0,)*len(group_moduli)
    q = [norm(x) for x in q]
    helper_charge = [norm(x) for x in helper_charge]
    charge = []
    for mask in range(1 << n):
        value = zero
        for i in range(n):
            if mask & (1 << i):
                value = add(value, q[i])
        charge.append(value)
    literal_states = [(h, u) for h, interval in enumerate(intervals)
                      for u in interval_members(p, interval)
                      if add(charge[u], helper_charge[h]) == zero]
    visible = sum(1 << i for i in range(n) if q[i] != zero or weights[i] > budget)
    neutral = p.all ^ visible
    fibres = {}
    for h, (forced_one, high) in enumerate(intervals):
        forced_zero = p.all ^ high
        for t in submasks(visible):
            compiled = p.interval(forced_one | t, forced_zero | (visible ^ t))
            if compiled is None or add(charge[t], helper_charge[h]) != zero:
                continue
            low, upper = compiled
            need(low & visible == t and upper & visible == t, "Visible bit fixation failed")
            fibres[h, t] = (low & neutral, upper & neutral)
    qkeys = sorted(fibres)
    qindex = {key: i for i, key in enumerate(qkeys)}
    need(set(qkeys) == {(h, u & visible) for h, u in literal_states},
         "Quotient state existence mismatch")
    # Check each induced neutral interval against full records, without assuming
    # that restricting to N retains all cross-visible implications.
    for key, (low, high) in fibres.items():
        h, t = key
        predicted = [u for u in submasks(neutral)
                     if low & ~u == 0 and u & ~high == 0
                     and all((p.up_rows[i] & neutral) & ~u == 0
                             for i in range(n) if u & (1 << i))]
        actual = [u & neutral for hh, u in literal_states if hh == h and u & visible == t]
        need(sorted(predicted) == sorted(actual), "Neutral induced-poset fibre mismatch")
        internal = [(i, j) for i, u in enumerate(actual) for j, v in enumerate(actual[:i])
                    if (u ^ v).bit_count() == 1 and wt[u ^ v] <= budget]
        cc = components(len(actual), internal)
        need(len(set(cc)) == 1, "Light neutral fibre is disconnected")
    literal_edges = []
    for i, (h, u) in enumerate(literal_states):
        for j, (g, v) in enumerate(literal_states[:i]):
            cost = wt[u ^ v] + helper_distance*(h != g)
            if cost <= budget:
                literal_edges.append((i, j))
    quotient_edges = []
    for i, a in enumerate(qkeys):
        for j, b in enumerate(qkeys[:i]):
            boundary, x, y = distance_and_witness(fibres[a], fibres[b], wt)
            cost = wt[a[1] ^ b[1]] + helper_distance*(a[0] != b[0]) + boundary
            xx = (a[0], a[1] | x)
            yy = (b[0], b[1] | y)
            need(xx in literal_states and yy in literal_states, "Full boundary reconstruction failed")
            need(wt[xx[1] ^ yy[1]] + helper_distance*(xx[0] != yy[0]) == cost,
                 "Full boundary metric mismatch")
            if cost <= budget:
                quotient_edges.append((i, j))
    lc = components(len(literal_states), literal_edges)
    qc = components(len(qkeys), quotient_edges)
    for i, (h, u) in enumerate(literal_states):
        for j, (g, v) in enumerate(literal_states[:i]):
            need((lc[i] == lc[j]) ==
                 (qc[qindex[h, u & visible]] == qc[qindex[g, v & visible]]),
                 "Full/quotient bounded connectivity disagreement")
    return len(literal_states), len(qkeys), len(set(lc))


def quotient_tests():
    cases = 0
    maximum_compression = 1
    for n in (1, 2):
        for p in natural_posets(n):
            intervals = [(a, b) for a in p.upsets for b in p.upsets if not a & ~b]
            for pair in product(intervals, repeat=2):
                for q in product((0, 1), repeat=n):
                    for a in product((0, 1), repeat=2):
                        for weights in [(1,)*n, tuple(i+1 for i in range(n))]:
                            for b in (0, 1, 2, 3):
                                full, compressed, _ = quotient_case(p, pair, q, a, weights, b)
                                if compressed:
                                    maximum_compression = max(maximum_compression, full/compressed)
                                cases += 1
    rng = random.Random(2026102307)
    random_cases = 1200
    for case in range(random_cases):
        n = rng.randrange(3, 8)
        p = Poset(n, [e for e in combinations(range(n), 2) if rng.randrange(4) == 0])
        intervals = []
        for h in range(rng.randrange(2, 5)):
            low = rng.choice(p.upsets)
            high = rng.choice([u for u in p.upsets if not low & ~u])
            intervals.append((low, high))
        q = [int(rng.randrange(4) == 0) for _ in range(n)]
        a = [rng.randrange(2) for _ in intervals]
        weights = [rng.randrange(0, 5) for _ in range(n)]
        b = rng.randrange(0, 7)
        full, compressed, _ = quotient_case(p, intervals, q, a, weights, b,
                                          helper_distance=rng.randrange(1, 4))
        if compressed:
            maximum_compression = max(maximum_compression, full/compressed)
    # Explicit heavy-neutral hostile: omitting the heavy coordinate would
    # merge disconnected full records into one purportedly connected fibre.
    p = Poset(1, [])
    full, compressed, comps = quotient_case(p, [(0, 1)], [0], [0], [2], 1)
    need((full, compressed, comps) == (2, 2, 2), "Heavy neutral control failed")
    # Empty intersection is not necessarily a forbidden coupled move.
    for b in (1, 2):
        full, compressed, comps = quotient_case(p, [(0, 0), (1, 1)], [0], [0, 0], [1], b)
        need(comps == (2 if b == 1 else 1), "Coupled boundary threshold failed")
    # Exact Z and a finite coordinate remain separate: erasing the former
    # changes the accepted records, not just their presentation.
    p3 = Poset(3, [])
    for b in (1, 2):
        full, compressed, comps = quotient_case(
            p3, [(0, 7)], [(1, 0), (-1, 0), (0, 0)], [(0, 0)], [1, 1, 1], b,
            group_moduli=(0, 4))
        need((full, compressed, comps) == (4, 2, 2 if b == 1 else 1),
             "Exact integer charge/product-group control failed")
    full, compressed, comps = quotient_case(p3, [(0, 7)], [0, 0, 0], [0], [1, 1, 1], 1,
                                           group_moduli=(4,))
    need((full, compressed, comps) == (8, 1, 1), "Charge-erasure hostile not distinguished")
    full, compressed, comps = quotient_case(Poset(2, []), [(0, 3)], [2, 2], [0], [1, 1], 1,
                                           group_moduli=(4,))
    need((full, compressed, comps) == (2, 2, 2), "Finite carry component was erased")
    return {"exhaustive_n1_n2_cases": cases, "seeded_larger_cases": random_cases,
            "seed": 2026102307, "maximum_observed_state_compression": maximum_compression,
            "explicit_heavy_and_coupled_controls": 3,
            "explicit_exact_integer_product_group_controls": 4,
            "scope": "All natural posets on1/2 elements, every pair of nonempty helper intervals, every Z2 charge, stated weights and budgets; seeded n3..7."}


P = ((3,0,6,7,1,5,4,2), (3,1,6,7,0,5,4,2))
G = (4,0,2,3,1,5,6,7)


def native_source(r):
    core = [[0]*8 for _ in range(8)]
    for bit, (u, v) in enumerate(combinations(range(8), 2)):
        core[u][v] = (29558785 >> (27-bit)) & 1
        core[v][u] = 1-core[u][v]
    p = 8*r
    return tuple(tuple(0 if u == v else
                       core[u % 8][v % 8] if u < p and v < p and u//8 == v//8 else
                       int(u < v) if u < p and v < p else
                       int(u < p) if (u < p) != (v < p) else
                       int((v-u) % 3 == 1)
                       for v in range(p+3)) for u in range(p+3))


def native_path(r, b, k):
    return tuple(8*j+u for j in range(r) for u in P[(b >> j) & 1]) + tuple(8*r+(k+j)%3 for j in range(3))


def partner(a, p):
    clocks = {u: j for j, u in enumerate(p)}
    return tuple(tuple(0 if u == v else a[u][v] if abs(clocks[u]-clocks[v]) == 1 else a[v][u]
                       for v in range(len(a))) for u in range(len(a)))


def word(a, p):
    return tuple(a[p[i]][p[j]] for i in range(len(p)) for j in range(i+2, len(p)))


def native_controls():
    rows = []
    for r in (1, 2, 3):
        a = native_source(r)
        n = len(a)
        paths = {(b, k): native_path(r, b, k) for b in range(1 << r) for k in range(3)}
        words = {z: word(a, p) for z, p in paths.items()}
        partners = {z: partner(a, p) for z, p in paths.items()}
        rev = tuple(zip(*a))
        raw = set(words.values()) | {word(rev, p[::-1]) for p in paths.values()}
        reservoir = raw | {tuple(1-v for v in w) for w in raw}
        need(len(reservoir) == 4*(1 << r), "Native normalized resource reservoir wrong")
        for b in range(1 << r):
            need(len({words[b, k] for k in range(3)}) == 1, "Kernel lifts created false capacity")
            prefix = paths[b, 0][:8*r]
            complete = [prefix+suffix for suffix in permutations(range(8*r, n))
                        if all(a[u][v] for u, v in zip(prefix+suffix, (prefix+suffix)[1:]))]
            need(set(complete) == {paths[b, k] for k in range(3)},
                 "Native prefix-pinned kernel menu is not complete")
        for (b, k), p in paths.items():
            need(sorted(p) == list(range(n)), "Native path not permutation")
            need(all(a[u][v] for u, v in zip(p, p[1:])), "Native path illegal")
            need(partner(rev, p[::-1]) == tuple(zip(*partners[b, k])), "Native converse square fails")
            need(p[::-1][2] == p[8*r], "Actual shared-owner converse clock wrong")
            for kk in range(3):
                delta = (kk-k) % 3
                gg = tuple(u if u < 8*r else 8*r+(u-8*r+delta)%3 for u in range(n))
                need(tuple(gg[u] for u in p) == paths[b, kk], "Native kernel lift mismatch")
                need(all(partners[b, k][u][v] == partners[b, kk][gg[u]][gg[v]]
                         for u in range(n) for v in range(n)), "Native target isomorphism fails")
            for j in range(r):
                bb = b ^ (1 << j)
                perm = G if not (b >> j) & 1 else tuple(G.index(u) for u in range(8))
                gg = tuple(8*j+perm[u % 8] if u//8 == j else u for u in range(n))
                need(all(partners[b, k][u][v] == partners[bb, k][gg[u]][gg[v]]
                         for u in range(n) for v in range(n)), "Native prefix target isomorphism fails")
        for slots in range(2, min(5, 1 << r)+1):
            accepted = []
            for phase in product(range(3), repeat=slots):
                records = [paths[b, phase[b]] for b in range(slots)]
                need(len({words[b, phase[b]] for b in range(slots)}) == slots,
                     "Private native prefix resources collided")
                if len({p[8*r] for p in records}) == 1:
                    accepted.append(phase)
            need(accepted == [(k,)*slots for k in range(3)], "Native shared-owner fibre not diagonal")
            for aa, bb in combinations(accepted, 2):
                changed_slots = sum(x != y for x, y in zip(aa, bb))
                clock_cost = sum(sum(x != y for x, y in zip(paths[i, aa[i]], paths[i, bb[i]]))
                                 for i in range(slots))
                need(changed_slots == slots and clock_cost == 3*slots, "Native sharp kernel budget wrong")
            # One helper, one old SCC; boundary is not a missing charge.
            m = slots-1
            pp = Poset(1, [])
            for budget in (m, m+1):
                _, _, cc = quotient_case(pp, [(0, 0), (1, 1)], [0], [0, 0], [m], budget)
                need(cc == (2 if budget == m else 1), "Native quotient bridge threshold wrong")
            rows.append({"prefix_blocks": r, "source_order": n, "slots": slots,
                         "full_kernel_menu_choices": 3**slots, "accepted_states": 3,
                         "sharp_changed_slot_budget": slots, "sharp_owner_clock_budget": 3*slots,
                         "normalized_resources": slots, "shared_forward_clock": 8*r,
                         "shared_converse_clock": 2})
    return rows


def input_hashes():
    paths = [
        "05-knowledge/results/native_packet_decomposition_20261023/interface/ENDPOINT_FILTERS.md",
        "05-knowledge/results/native_packet_decomposition_20261023/interface/COMPENSATOR_LAYERS.md",
        "05-knowledge/results/native_packet_decomposition_20261023/workspace/ONE_BIT_WORKSPACE.md",
        "05-knowledge/results/native_packet_decomposition_20261023/workspace/check.py",
        "05-knowledge/results/balanced_native_20261023/SIX_TRADE_LATTICE.md",
    ]
    return {p: sha256((ROOT/p).read_bytes().replace(b"\r\n", b"\n")).hexdigest() for p in paths}


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("--write", action="store_true")
    args = parser.parse_args()
    result = {"status": "FINITE-EXACT supporting scoped PROVED lifting theorem",
              "input_sha256_lf_normalized": input_hashes(),
              "intervals": interval_tests(), "quotients": quotient_tests(),
              "native_shared_kernel_controls": native_controls()}
    result["checks_before_replay"] = CHECKS
    result = json.loads(json.dumps(result))
    output = HERE/"receipt.json"
    if args.write:
        output.write_text(json.dumps(result, indent=2, sort_keys=True)+"\n", encoding="utf-8")
    else:
        need(result == json.loads(output.read_text(encoding="utf-8")), "Frozen receipt changed")
    print(f"PASS {CHECKS} checks; {'wrote' if args.write else 'replayed'} {output.name}")


if __name__ == "__main__":
    main()
