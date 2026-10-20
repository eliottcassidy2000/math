"""Exact, Burnside-free tournament/tiling census and compression hostile probes.

Run from repository root:
  python 05-knowledge/results/metagraph_growth_20261020/enumeration.py

No imported research implementation. The algorithm reuses the documented
class-BFS idea, not its code. Canonical codes use invariant score refinement;
independent full-permutation minima audit every class through order six.
Bit convention: 1 means i -> j for i < j, including fixed forward path arcs.
"""
from collections import Counter, defaultdict, deque
from itertools import permutations, product
from math import comb, factorial


def need(test, message):
    if not test:
        raise RuntimeError(message)


def pairs(n):
    return tuple((i, j) for i in range(n) for j in range(i + 1, n))


def tile_pairs(n):
    return tuple((i, j) for i, j in pairs(n) if j > i + 1)


def matrix(n, bits):
    rows = [0] * n
    for k, (i, j) in enumerate(pairs(n)):
        if bits >> k & 1:
            rows[i] |= 1 << j
        else:
            rows[j] |= 1 << i
    return tuple(rows)


def rowcode(rows):
    n = len(rows)
    z = 0
    for r in rows:
        z = (z << n) | r
    return z


def decode(n, code):
    mask = (1 << n) - 1
    return tuple((code >> (n * (n - i - 1))) & mask for i in range(n))


def ordercode(rows, order):
    return rowcode(tuple(sum(((rows[v] >> w) & 1) << j
                             for j, w in enumerate(order)) for v in order))


def canonical(rows, aut=False):
    n = len(rows)
    colour = tuple(r.bit_count() for r in rows)
    while True:
        labels = sorted(set(colour))
        signature = tuple((colour[v], tuple(sum((rows[v] >> w) & 1
                                                for w in range(n) if colour[w] == c)
                                           for c in labels)) for v in range(n))
        ranks = {s: k for k, s in enumerate(sorted(set(signature)))}
        new = tuple(ranks[s] for s in signature)
        if new == colour:
            break
        colour = new
    cells = [tuple(v for v in range(n) if colour[v] == c)
             for c in sorted(set(colour))]
    best, count = None, 0
    for blocks in product(*(permutations(c) for c in cells)):
        order = tuple(v for b in blocks for v in b)
        code = ordercode(rows, order)
        if best is None or code < best:
            best, count = code, 1
        elif code == best:
            count += 1
    return (best, count) if aut else best


def slow_canonical(rows):
    return min(ordercode(rows, p) for p in permutations(range(len(rows))))


def converse(rows):
    full = (1 << len(rows)) - 1
    return tuple(full ^ r ^ (1 << v) for v, r in enumerate(rows))


def flip(rows, i, j):
    new = list(rows)
    new[i] ^= 1 << j
    new[j] ^= 1 << i
    return tuple(new)


def ham_count(rows):
    n = len(rows)
    dp = {}
    for v in range(n):
        dp[(1 << v, v)] = 1
    for mask in range(1, 1 << n):
        for v in range(n):
            value = dp.get((mask, v), 0)
            for w in range(n):
                if value and not (mask >> w & 1) and rows[v] >> w & 1:
                    key = (mask | (1 << w), w)
                    dp[key] = dp.get(key, 0) + value
    return sum(dp.get(((1 << n) - 1, v), 0) for v in range(n))


def tiling(n, t):
    rows = [0] * n
    for i in range(n - 1):
        rows[i] |= 1 << (i + 1)
    for k, (i, j) in enumerate(tile_pairs(n)):
        if t >> k & 1:
            rows[i] |= 1 << j
        else:
            rows[j] |= 1 << i
    return tuple(rows)


def sigma(n, t):
    ts = tile_pairs(n)
    index = {p: i for i, p in enumerate(ts)}
    return sum(((t >> k) & 1) << index[(n - 1 - j, n - 1 - i)]
               for k, (i, j) in enumerate(ts))


def extension(rows, mask):
    n = len(rows)
    r = list(rows) + [0]
    for v in range(n):
        if mask >> v & 1:
            r[v] |= 1 << n
        else:
            r[n] |= 1 << v
    return tuple(r)


def census(n):
    start = canonical(matrix(n, (1 << comb(n, 2)) - 1))
    todo = deque([start])
    reps = {start: decode(n, start)}
    flipadj = defaultdict(set)
    calls = 0
    while todo:
        c = todo.popleft()
        for i, j in pairs(n):
            nxt = canonical(flip(reps[c], i, j))
            calls += 1
            if nxt != c:
                flipadj[c].add(nxt)
            if nxt not in reps:
                reps[nxt] = decode(n, nxt)
                todo.append(nxt)
    dual = {c: canonical(converse(r)) for c, r in reps.items()}
    need(all(dual[dual[c]] == c for c in reps), "converse involution")
    merged = {c: min(c, dual[c]) for c in reps}
    fibre = Counter()
    blue = Counter()
    tc = []
    m = comb(n - 1, 2)
    for t in range(1 << m):
        c = canonical(tiling(n, t))
        need(c in reps, "path cube class missing from BFS")
        tc.append(c)
        fibre[c] += 1
        if sigma(n, t) == t:
            blue[c] += 1
        need(canonical(tiling(n, sigma(n, t))) == dual[c],
             "sigma is converse on classes")
    need(set(fibre) == set(reps), "BFS/path cube class mismatch")
    ha = {}
    for c, r in reps.items():
        _, aut = canonical(r, aut=True)
        H = ham_count(r)
        ha[c] = (H, aut)
        need(H == aut * fibre[c], "H/Aut fibre law")
        need(fibre[c] % 2 == 1 and aut % 2 == 1, "odd fibre and automorphism order")
    need(sum(blue.values()) == 1 << (((n - 1) ** 2) // 4), "blue global count")
    if n <= 6:
        slow = {slow_canonical(r) for r in reps.values()}
        need(len(slow) == len(reps), "refinement canonicalization false split")
        # Exhaustive labelled orbit mass checks completeness independently
        # of the path cube and of connectivity.
        mass = sum(factorial(n) // canonical(r, aut=True)[1] for r in reps.values())
        need(mass == 1 << comb(n, 2), "labelled orbit mass")
    stat = {u: [0] * 6 for u in set(merged.values())}
    # Fields F,B,blue_other,black_other,blue_self,black_self.
    for c in reps:
        u = merged[c]
        stat[u][0] += fibre[c]
        stat[u][1] += blue[c]
    partners = {u: Counter() for u in stat}
    if m:
        complement_mask = (1 << m) - 1
        for t in range(1 << (m - 1)):
            q = t ^ complement_mask
            u, v = merged[tc[t]], merged[tc[q]]
            colour = 0 if sigma(n, t) == t else 1
            partners[u][(colour, v)] += 1
            partners[v][(colour, u)] += 1
            if u == v:
                stat[u][4 + colour] += 1
            else:
                stat[u][2 + colour] += 1
                stat[v][2 + colour] += 1
        for u, s in stat.items():
            F, B, bo, ko, bs, ks = s
            need(F == bo + ko + 2 * bs + 2 * ks, "degree/fibre equality")
            need(B == bo + 2 * bs, "blue incidence equality")
            need(ko % 2 == 0, "black external incidence parity")
            need((F % 2 == 1) == (dual[u] == u), "merged fibre parity")
    tau_witness = None
    targetsets = defaultdict(dict)
    for t, c in enumerate(tc):
        partner = merged[tc[t ^ ((1 << m) - 1)]]
        targetsets[c].setdefault(partner, t)
        if len(targetsets[c]) > 1 and tau_witness is None:
            entries = list(targetsets[c].items())[:2]
            tau_witness = (c, entries)
    return dict(n=n, reps=reps, dual=dual, merged=merged, tc=tc,
                fibre=fibre, blue=blue, stat=stat, partners=partners,
                flipadj=flipadj, calls=calls, tau_witness=tau_witness, ha=ha)


def print_counter(counter):
    return str(sorted(counter.items()))


def local_signature(data, u):
    # The six incidence coordinates, plus simple partner degree and
    # merged single-arc-flip degree. Self loops do not count in these degrees.
    cache = data.setdefault("signature_cache", {})
    if u in cache:
        return cache[u]
    merged_flip = {data["merged"][v]
                   for c in data["reps"] if data["merged"][c] == u
                   for v in data["flipadj"][c]}
    merged_flip.discard(u)
    partner_degree = len({v for _, v in data["partners"][u] if v != u})
    result = tuple(data["stat"][u]) + (partner_degree, len(merged_flip))
    cache[u] = result
    return result


def ha_signature(data, u):
    return local_signature(data, u) + data["ha"][u]


def basic_ha_signature(data, u):
    return tuple(data["stat"][u][:2]) + data["ha"][u]


def profile_collision(data, profile, signature=local_signature):
    groups = defaultdict(list)
    for u in sorted(data["stat"]):
        groups[signature(data, u)].append(u)
    collisions = []
    for sig, us in groups.items():
        for i, u in enumerate(us):
            for v in us[i + 1:]:
                if profile[u] != profile[v]:
                    collisions.append((sig, u, v))
    return collisions


def main():
    data = {}
    expected = {1: (1, 1), 2: (1, 1), 3: (2, 2), 4: (4, 2),
                5: (12, 8), 6: (56, 12), 7: (456, 88)}
    print("EXACT UNIVERSE: all tournament classes and all fixed-path tilings, n=1..7.")
    print("BFS uses single-arc flips; no Burnside sum. Canonicalization retains all arcs.")
    for n in range(1, 8):
        d = data[n] = census(n)
        a = len(d["reps"])
        s = sum(c == q for c, q in d["dual"].items())
        M = len(d["stat"])
        need((a, s) == expected[n], "known count control")
        need(a == 2 * M - s, "merged count formula")
        print(f"\nn={n}: classes={a}, SC={s}, merged={M}, merged_pairs={M-s}, "
              f"tilings={len(d['tc'])}, BFS_arc_calls={d['calls']}")
        print("unmerged fibre histogram:", print_counter(Counter(d["fibre"].values())))
        print("merged (fibre,blue) histogram:",
              print_counter(Counter((v[0], v[1]) for v in d["stat"].values())))
        if n >= 3:
            print("blue/black lines:", sum(v[1] for v in d["stat"].values()) // 2,
                  (len(d["tc"]) - sum(v[1] for v in d["stat"].values())) // 2)
        witness = d["tau_witness"]
        if witness:
            c, entries = witness
            print("tau-not-functional merged-target witness:",
                  dict(source_class=c, rows=d["reps"][c],
                       tiling_target_pairs=[(t, u) for u, t in entries]))
        else:
            print("tau-not-functional merged-target witness: none")
        profiles = {}
        for u in d["stat"]:
            profiles[u] = Counter()
            for (colour, v), count in d["partners"][u].items():
                profiles[u][(colour, local_signature(d, v))] += count
        hostile = profile_collision(d, profiles)
        print("same eight local coordinates/different neighbour-profile pairs:", len(hostile))
        rich_hostile = profile_collision(d, profiles, ha_signature)
        print("same eight coordinates plus H,Aut/different neighbour-profile pairs:",
              len(rich_hostile))
        if hostile:
            sig, u, v = hostile[0]
            print("line-profile witness:", dict(signature=sig, u=u, v=v,
                                               rows_u=d["reps"][u], rows_v=d["reps"][v],
                                               HA_u=d["ha"][u], HA_v=d["ha"][v]))
            print("profile u:", print_counter(profiles[u]))
            print("profile v:", print_counter(profiles[v]))
    print("\nONE-VERTEX EXTENSION PROBE: every 2^n incidence pattern, n=1..6.")
    for n in range(1, 7):
        d, child = data[n], data[n + 1]
        full, coarse = {}, {}
        for u in sorted(d["stat"]):
            full[u], coarse[u] = Counter(), Counter()
            for mask in range(1 << n):
                c = canonical(extension(d["reps"][u], mask))
                need(c in child["reps"], "extension absent from child census")
                v = child["merged"][c]
                full[u][v] += 1
                coarse[u][local_signature(child, v)] += 1
        hostile = profile_collision(d, coarse)
        fullhostile = profile_collision(d, full)
        ha_hostile = profile_collision(d, coarse, ha_signature)
        basic_hostile = profile_collision(d, coarse, basic_ha_signature)
        print(f"n={n}: coarse extension-profile collision pairs={len(hostile)}, "
              f"exact child-orbit profile collision pairs={len(fullhostile)}, "
              f"plus H,Aut collisions={len(ha_hostile)}, "
              f"fibre,blue,H,Aut collisions={len(basic_hostile)}")
        if hostile:
            sig, u, v = hostile[0]
            print("extension-profile witness:", dict(signature=sig, u=u, v=v,
                                                    rows_u=d["reps"][u], rows_v=d["reps"][v]))
            print("child signatures u:", print_counter(coarse[u]))
            print("child signatures v:", print_counter(coarse[v]))
        if basic_hostile:
            sig, u, v = basic_hostile[0]
            print("basic H/Aut extension witness:",
                  dict(signature=sig, rows_u=d["reps"][u], rows_v=d["reps"][v],
                       child_distinct_u=len(full[u]), child_distinct_v=len(full[v]),
                       local_u=local_signature(d, u), local_v=local_signature(d, v)))
    print("\nMINIMAL EXPLICIT HOSTILES:")
    a, b = tiling(4, 0), tiling(4, 1)
    iso = (1, 2, 3, 0)
    need(all(((a[i] >> j) & 1) == ((b[iso[i]] >> iso[j]) & 1)
             for i in range(4) for j in range(4)), "explicit tau source isomorphism")
    score_a = sorted(r.bit_count() for r in tiling(4, 7))
    score_b = sorted(r.bit_count() for r in tiling(4, 6))
    need(score_a == [0, 1, 2, 3] and score_b == [0, 2, 2, 2],
         "explicit tau target distinction")
    print("n4 tile order:", tile_pairs(4))
    print("tau hostile t=0,1: rows", a, b, "isomorphism", iso,
          "complement scores", score_a, score_b)
    for t in (0, 1):
        profile = Counter(tuple(data[4]["stat"][data[4]["merged"][data[4]["tc"][q]]])
                          for q in range(8) if q % 2 == t)
        print("fixed-path append source tiling", t, "scores",
              sorted(r.bit_count() for r in tiling(3, t)),
              "four children merged incidence profiles", print_counter(profile))
    print("\nALL CHECKS PASSED (explicit exceptions; active under python -O).")


if __name__ == "__main__":
    main()
