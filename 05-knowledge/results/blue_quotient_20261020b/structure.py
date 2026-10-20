"""Exact controls for the structural black-path lemmas; no output writes."""
from pathlib import Path
from collections import Counter, defaultdict
import itertools
import json
import runpy

ROOT = Path(__file__).resolve().parent
CANON = runpy.run_path(str(ROOT.parent / "merged_metagraph_distance_moment_constraints_20261020" / "probe.py"))["canon_tilings"]


def require(ok, message):
    if not ok:
        raise RuntimeError(message)


def path_is_blue(out, path):
    n = len(out)
    r = [0] * n
    for i, v in enumerate(path):
        r[v] = path[n - 1 - i]
    return all(((out[u] >> v) & 1) == ((out[r[v]] >> r[u]) & 1)
               for u in range(n) for v in range(n) if u != v)


def is_path(out, p, closed=False):
    return len(set(p)) == len(out) and all((out[a] >> b) & 1 for a, b in zip(p, p[1:] + (p[:1] if closed else [])))


def reverse_path(out, path):
    result = list(out)
    for a, b in zip(path, path[1:]):
        require((result[a] >> b) & 1, "reversing a nonexistent directed path arc")
        result[a] ^= 1 << b
        result[b] ^= 1 << a
    return result


def fixed_path_checks(n):
    free, codes, fibre = CANON(n)
    index = {p: i for i, p in enumerate(free)}
    mirror = [index[(n-1-j, n-1-i)] for i, j in free]
    hist = defaultdict(Counter)
    coset_checks = switch_checks = 0
    for x, u in enumerate(fibre):
        out = [0] * n
        for i in range(n-1):
            out[i] |= 1 << (i+1)
        for k, (i, j) in enumerate(free):
            a, b = (i, j) if (x >> k) & 1 else (j, i)
            out[a] |= 1 << b
        scores = [row.bit_count() for row in out]
        blue = all(((x >> i) & 1) == ((x >> j) & 1) for i, j in enumerate(mirror))
        closed = not ((x >> index[(0, n-1)]) & 1)
        h = hist[u]
        h["black"] += not blue
        h["closed"] += closed
        if closed:
            delta = scores[-1] - scores[0]
            h["blue_current" if blue else "black_current"] += delta
            h["total_current"] += delta
            h["equal"] += delta == 0
            if not blue:
                h["black_closed"] += 1
                h["good"] += delta >= 0
                h["black_equal"] += delta == 0
        if n <= 6 and closed:
            cycle = list(range(n))
            B = {i for i in range(n) if path_is_blue(out, cycle[i:] + cycle[:i])}
            if B:
                i = min(B)
                H = {(j-i) % n for j in B}
                require(all((a+b) % n in H for a in H for b in H), "blue cuts not a coset")
                require(all((j in H) == all(((out[a] >> b) & 1) == ((out[(a+2*j) % n] >> ((b+2*j) % n)) & 1)
                                           for a in range(n) for b in range(n) if a != b) for j in range(n)),
                        "blue coset differs from rotational automorphisms")
            require(len(B) == n or n-len(B) >= 2, "one isolated black cut")
            coset_checks += 1
        if n <= 6 and not closed and not blue:
            middle = reverse_path(out, list(range(n)))
            cycle = list(reversed(range(n)))
            require(is_path(middle, cycle, closed=True), "open path did not become cycle")
            choices = [cycle[i:]+cycle[:i] for i in range(1,n)
                       if not path_is_blue(middle, cycle[i:]+cycle[:i])]
            require(choices, "no second black cut")
            for q in choices:
                target = reverse_path(middle, q)
                changed = [(a,b) for a in range(n) for b in range(a+1,n)
                           if ((out[a] ^ target[a]) >> b) & 1]
                require(len(changed) == 2 and (0,n-1) in changed,
                        "two-step operation is not shortcut plus one arc")
                require(any(b == a+1 for a,b in changed), "other changed arc is not a path arc")
                switch_checks += 1
    active = {u for u,h in hist.items() if h["closed"] and h["black"]}
    require(all(hist[u]["black_closed"] for u in active), "black closed path missing")
    require(all(hist[u]["good"] for u in active), "FINITE DESCENT CANDIDATE FAILED")
    require(all(h["total_current"] == 0 for h in hist.values()), "cycle-cut current not zero")
    positive = [(codes[u],h["blue_current"]) for u,h in hist.items() if h["blue_current"] > 0]
    if n == 5:
        u = codes.index(10)
        require(hist[u]["blue_current"] == 2 and hist[u]["black_current"] == -2 and hist[u]["good"] == 2,
                "mean-current hostile changed")
        require(hist[codes.index(12)]["black_equal"] == 0, "equal-endpoint hostile changed")
    if n == 7:
        require(hist[codes.index(310)]["equal"] == 0, "order-seven equal-endpoint hostile changed")
    return {"n":n,"tilings":len(fibre),"strong_black_active":len(active),
            "black_closed_descent_failures":0,"positive_blue_current_classes":len(positive),
            "positive_blue_current_first":positive[:2],"cycle_cosets_checked":coset_checks,
            "literal_two_black_switches_checked":switch_checks}


def even_constructive_checks():
    total = 0
    for k in (3,5,7,9,11):
        # All inner circulants and cross words for k=3,5; structured samples later.
        inner_masks = range(1 << ((k-1)//2)) if k <= 5 else (0,1,(1 << ((k-1)//2))-1)
        cross_masks = range(1 << k) if k <= 5 else (1,3,(1 << (k-1))-1,((1 << k)-1)^2)
        for im in inner_masks:
            S = {d if (im >> (d-1)) & 1 else k-d for d in range(1,(k+1)//2)}
            epsilon = 1 if 1 in S else -1
            for am in cross_masks:
                if not (am & 1) or ((am >> (k-1)) & 1):
                    continue
                a = lambda d: (am >> (d % k)) & 1
                out = [0] * (2*k)
                for i,j in itertools.permutations(range(k),2):
                    if (j-i) % k in S:
                        out[2*i] |= 1 << (2*j)
                        out[2*i+1] |= 1 << (2*j+1)
                for i,j in itertools.product(range(k),repeat=2):
                    u,w = 2*i,2*j+1
                    if a(j-i):out[u] |= 1 << w
                    else:out[w] |= 1 << u
                original = list(range(2*k))
                require(is_path(out,original,closed=True),"initial alternating cycle invalid")
                require(all(path_is_blue(out,original[j:]+original[:j]) for j in range(2*k)),
                        "two-orbit model not all-blue along original cycle")
                t = next(t for t in range(k) if a(t+epsilon) and not a(t-epsilon))
                cycle = [2*((j*epsilon)%k) for j in range(k)] + [2*((t+j*epsilon)%k)+1 for j in range(k)]
                require(is_path(out,cycle,closed=True),"concatenated cycle invalid")
                path = cycle[1:]+cycle[:1]  # remove first internal U arc
                require(not path_is_blue(out,path),"constructed path unexpectedly blue")
                require(out[path[0]].bit_count() == out[path[-1]].bit_count(),"endpoint scores unequal")
                total += 1
    return {"even_two_orbit_constructions":total,"orders":[6,10,14,18,22],"failures":0}


def same_cycle_hostile():
    out = [98,92,409,273,225,206,140,267,115]
    n = len(out)
    require(all(((out[u] >> v) & 1) + ((out[v] >> u) & 1) == 1
                for u in range(n) for v in range(u+1,n)), "not a tournament")
    scores = [r.bit_count() for r in out]
    require(scores == [3,4,5]*3, "hostile scores changed")
    cycle = list(range(n))
    require(is_path(out,cycle,True), "hostile cycle invalid")
    blue = [i for i in range(n) if path_is_blue(out,cycle[i:]+cycle[:i])]
    require(blue == [0,3,6], "hostile blue cuts changed")
    require(all(scores[(i-1)%n]-scores[i] == -1 for i in range(n) if i not in blue),
            "a black cut descends on the hostile cycle")
    free = [(i,j) for i in range(n) for j in range(i+2,n)]
    mask = sum(((out[i] >> j) & 1) << k for k,(i,j) in enumerate(free))
    require(mask == 48440728, "hostile tiling mask changed")
    path = [0,1,2,3,8,4,5,6,7]
    require(is_path(out,path,True) and not path_is_blue(out,path), "replacement path invalid or blue")
    require(scores[path[-1]]-scores[path[0]] == 1, "replacement does not descend")
    return {"same_cycle_hostile_n":n,"out_bitsets":out,"tiling_mask":mask,
            "scores":scores,"blue_cuts":blue,"black_cut_endpoint_difference":-1,
            "changed_cycle_black_closed_path":path,"triangle_count_change":-2}


if __name__ == "__main__":
    for n in range(3,8):
        print(json.dumps(fixed_path_checks(n)),flush=True)
    print(json.dumps(even_constructive_checks()),flush=True)
    print(json.dumps(same_cycle_hostile()),flush=True)
    print("ALL STRUCTURAL CONTROLS PASSED; all-order black connectivity remains OPEN.")
