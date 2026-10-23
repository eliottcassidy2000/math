"""Exact founder, phase, and coupled-coloring checks; normal/-O invariant."""
from collections import Counter, deque
from itertools import product
from math import gcd, isqrt
from pathlib import Path
import argparse
import json

HERE = Path(__file__).resolve().parent
CHECKS = 0


def require(statement, message):
    global CHECKS
    CHECKS += 1
    if not statement:
        raise RuntimeError(message)


def older(n):
    return [k*k-n for k in range(isqrt(n)+1, isqrt(2*n-1)+1)]


def edges(n):
    return [(a, b) for b in range(2, n+1) for a in older(b)]


def components(n, removed=()):
    adjacency = {a: set() for a in range(1, n+1) if a not in removed}
    for a,b in edges(n):
        if a in adjacency and b in adjacency:
            adjacency[a].add(b)
            adjacency[b].add(a)
    unseen = set(adjacency)
    result = []
    while unseen:
        seen = {min(unseen)}
        stack = list(seen)
        while stack:
            x = stack.pop()
            for y in adjacency[x]-seen:
                seen.add(y)
                stack.append(y)
        result.append(sorted(seen))
        unseen -= seen
    return sorted(result)


def founder_data(limit):
    seeds = {1: 1, 2: 2, 4: 4}
    ancestry = {}
    count = {}
    parent = {}
    for n in range(1, limit+1):
        neighbors = older(n)
        if not neighbors:
            require(n in seeds, "unexpected founder")
            ancestry[n] = seeds[n]
            count[n] = tuple(int(seeds[n] == 1 << j) for j in range(3))
            parent[n] = {j: None for j in range(3) if seeds[n] & 1 << j}
        else:
            ancestry[n] = 0
            for x in neighbors:
                ancestry[n] |= ancestry[x]
            count[n] = tuple(sum(count[x][j] for x in neighbors) for j in range(3))
            parent[n] = {j: min(x for x in neighbors if ancestry[x] & 1 << j)
                         for j in range(3) if ancestry[n] & 1 << j}
        require(ancestry[n] == sum(1 << j for j in range(3) if count[n][j]),
                "reachability and multiplicity disagree")
    return ancestry, count, parent


def founder_path(n, j, parent):
    path = [n]
    while parent[path[-1]][j] is not None:
        path.append(parent[path[-1]][j])
    return path


def proper(colors, n):
    return all(colors[a] != colors[b] for a,b in edges(n))


def first_fit(q, limit):
    coloring = {}
    for n in range(1, limit+1):
        used = {coloring[x] for x in older(n)}
        choices = set(range(q))-used
        if not choices:
            return coloring, n
        coloring[n] = min(choices)
    return coloring, None


def berggren_word(m, n):
    reverse = []
    while (m,n) != (1,2):
        if n < 2*m:
            reverse.append("L")
            m,n = 2*m-n,m
        elif n < 3*m:
            reverse.append("M")
            m,n = n-2*m,m
        else:
            reverse.append("R")
            m,n = m,n-2*m
        require(0 < m < n, "Berggren inverse escaped chamber")
    return "".join(reversed(reverse))


def run():
    data = {}
    # Independent adjacency path for all small n, including no-loop cases.
    for n in range(1, 101):
        slow = {(a,b) for b in range(2,n+1) for a in range(1,b)
                if isqrt(a+b)**2 == a+b}
        require(set(edges(n)) == slow, "square edge generation")
    data["islands_12"] = components(12)
    data["islands_18_delete_13_14_15"] = components(18, (13,14,15))
    require(data["islands_12"] == [[1,3,6,8,10],[2,7,9],[4,5,11,12]], "initial islands")
    require(data["islands_18_delete_13_14_15"] ==
            [[1,3,6,8,10,17],[2,7,9,16,18],[4,5,11,12]], "supplied islands")
    data["component_counts_1_25"] = [len(components(n)) for n in range(1,26)]

    ancestry, counts, parent = founder_data(10000)
    require(ancestry[22] == 7 and all(ancestry[n] != 7 for n in range(1,22)), "first complete ancestor set")
    require(ancestry[70] == 5 and all(ancestry[n] == 7 for n in range(71,10001)), "least complete tail")
    for n in range(100,10001):
        k = isqrt(n+71)
        if k*k < n+71:
            k += 1
        require(71 <= k*k-n < n, "induction predecessor interval")
        require(n*n-150*n+5045 > 0, "square-gap inequality")
    # Independent directed DFS determines the same three reachability bits.
    for n in range(1,121):
        seen = {n}
        stack = [n]
        while stack:
            x = stack.pop()
            for y in older(x):
                if y not in seen:
                    seen.add(y)
                    stack.append(y)
        require(ancestry[n] == sum(bit for seed,bit in ((1,1),(2,2),(4,4)) if seed in seen), "directed DFS")
    data["ancestry_1_100"] = [ancestry[n] for n in range(1,101)]
    data["induction_base_71_99"] = [ancestry[n] for n in range(71,100)]
    data["last_incomplete"] = {"vertex":70,"mask":ancestry[70],"older":older(70)}
    data["founder_paths_71"] = [founder_path(71,j,parent) for j in range(3)]
    data["path_counts"] = {str(n):list(counts[n]) for n in (22,25,70,71,100,1000,10000)}

    # Founder-copy colorings, not proper graph colorings. Shared parents retained.
    states = [()]
    prefix_counts = []
    first_unary_failure = None
    marginal_product = 1
    for n in range(1,36):
        neighbors = older(n)
        if not neighbors:
            fixed = {1:0,2:1,4:2}[n]
            states = [s+(fixed,) for s in states]
        else:
            states = [s+(c,) for s in states for c in sorted({s[x-1] for x in neighbors})]
        marginal_product *= ancestry[n].bit_count()
        prefix_counts.append(len(states))
        if first_unary_failure is None and len(states) != marginal_product:
            first_unary_failure = n
        for c in range(3):
            require(any(s[-1] == c for s in states) == bool(ancestry[n] & (1 << c)), "marginal founder color")
        if n == 22:
            require(not any(s[13] == 1 and s[21] == 2 for s in states), "shared-parent hostile")
            require(any(s[13] == 1 for s in states) and any(s[21] == 2 for s in states), "hostile marginals live")
    require(first_unary_failure == 22, "earliest unary-factorization failure")
    require(prefix_counts[24] == 192 and prefix_counts[-1] == 66780, "founder-color counts")
    data["founder_copy_counts_1_35"] = prefix_counts
    data["first_unary_failure"] = {"n":22,"forbidden": {"14":1,"22":2},"actual":32,"unary_product":48}

    # A finite, explicitly defined candidate family rather than an unspecified sweep.
    candidates = []
    for q in (2,3):
        for coeff in product(range(q),repeat=4):
            if not any(coeff):
                continue
            colors = {n:sum(a*b for a,b in zip(coeff, (n,(n&-n).bit_length()-1,n.bit_count(),isqrt(n)))) % q
                      for n in range(1,121)}
            witness = next(((a,b) for a,b in edges(120) if colors[a] == colors[b]),None)
            require(witness is not None, "candidate survived finite hostile scan")
            candidates.append({"q":q,"coefficients":list(coeff),"first_bad_edge":list(witness)})
    data["candidate_family"] = candidates
    data["candidate_worst_first_failure"] = {str(q):max(r["first_bad_edge"][1] for r in candidates if r["q"]==q) for q in (2,3)}

    # Actual three-vertex simultaneous recoloring, with every outside vertex pinned.
    two, fail2 = first_fit(2,100)
    original, fail3 = first_fit(3,100)
    require(fail2 == 15 and fail3 == 35, "first-fit boundaries")
    require(proper(original,34), "initial proper coloring")
    block = {1,15,24}
    for x in block:
        for color in set(range(3))-{original[x]}:
            single = dict(original)
            single[x] = color
            require(not proper(single,34), "single-site unexpectedly movable")
    repaired = dict(original)
    for x in block:
        repaired[x] = 2-original[x]
    require(proper(repaired,34), "simultaneous block repair")
    require(all(repaired[x] == original[x] for x in range(1,35) if x not in block), "outside pin change")
    repaired[35] = 0
    require(proper(repaired,35), "arrival extension")
    boundary = sorted({a if b in block else b for a,b in edges(34) if (a in block) != (b in block)})
    require(all(original[x] == 1 for x in boundary), "boundary monochromatic guard")
    legal_block_states = []
    for values in product(range(3), repeat=3):
        candidate = dict(original)
        candidate.update(zip(sorted(block),values))
        if proper(candidate,34):
            legal_block_states.append(list(values))
    require(legal_block_states == [[0,2,2],[2,0,0]], "exact pinned phase fiber")
    data["coupled_move"] = {"block":sorted(block),"pinned_boundary":boundary,
       "old": [original[n] for n in range(1,35)],"new":[repaired[n] for n in range(1,36)],
       "new_vertex_older":older(35),"legal_block_states":legal_block_states,
       "color_zero_counts_on_block":[1,2]}

    # Label transport on reflections: finite modes and the three-state residue action.
    for modulus in range(2,33):
        for t1 in range(0,12):
            for t2 in range(0,12):
                for x in range(modulus):
                    require((t2-(t1-x)) % modulus == (x+t2-t1) % modulus, "reflection composition")
        require((modulus*modulus+3*modulus*modulus) == (2*modulus)**2,
                "periodic-color hostile")
    for a,b in edges(1000):
        k = isqrt(a+b)
        require((a % 2) ^ (b % 2) == k % 2, "parity coboundary")
        require(b % 3 == (k*k-a) % 3, "three-state reflection")
    for n in range(2,101):
        require(all(isqrt(4*a+4*b)**2 == 4*a+4*b for a,b in edges(n)), "square-scale lift")
        require(all((isqrt(a+b)**2 == a+b) == (isqrt(4*a+4*b)**2 == 4*a+4*b)
                    for a in range(1,n) for b in range(a+1,n+1)), "induced square-scale lift")
    data["residue_reflections_mod3"] = {"sum0":[0,2,1],"sum1":[1,0,2]}

    # Inherited Pythagorean construction: three labelled matchings form a path.
    pyth = []
    for m,n,label in ((1,2,"root"),(2,3,"L"),(2,5,"M"),(1,4,"R")):
        a,b,c = n*n-m*m,2*m*n,n*n+m*m
        s,t = sorted((a,b))
        limit = t*t-1
        targets = {s*s,t*t,c*c}
        adjacency = {x:[] for x in range(1,limit+1)}
        edge_count = 0
        for x in range(1,limit+1):
            for target in targets:
                y=target-x
                if x < y <= limit:
                    adjacency[x].append(y);adjacency[y].append(x);edge_count += 1
        seen={1};stack=[1]
        while stack:
            for y in adjacency[stack.pop()]:
                if y not in seen:seen.add(y);stack.append(y)
        require(len(seen)==limit and edge_count==limit-1 and max(map(len,adjacency.values()))==2,
                "three-square Pythagorean chain")
        pyth.append({"branch":label,"triple":[a,b,c],"n":limit,"endpoints":[x for x in adjacency if len(adjacency[x])==1]})
    data["pythagorean_branch_chains"] = pyth
    collision = []
    for m,n in ((14,15),(10,21)):
        triple = sorted((n*n-m*m,2*m*n))+[n*n+m*m]
        require(triple[1]==420 and gcd(triple[0],triple[1])==1, "equal-order Pythagorean hostile")
        collision.append({"parameters":[m,n],"triple":triple,"word":berggren_word(m,n)})
    require(collision[0]["word"][0] != collision[1]["word"][0], "first-branch collision")
    data["unlabelled_chain_hostile"] = collision
    data["checks"] = CHECKS
    return data


if __name__ == "__main__":
    parser = argparse.ArgumentParser()
    parser.add_argument("--write", action="store_true")
    args = parser.parse_args()
    result = run()
    target = HERE / "receipt.json"
    if args.write:
        target.write_text(json.dumps(result,indent=2,sort_keys=True)+"\n",encoding="utf-8")
    else:
        require(result == json.loads(target.read_text(encoding="utf-8")), "stored receipt")
    print(json.dumps({"status":"PASS","checks":CHECKS,"last_incomplete_founder":70,
                      "colorings_at_35":66780,"coupled_block":[1,15,24]},sort_keys=True))
