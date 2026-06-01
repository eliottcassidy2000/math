"""Investigations 2,3,4: smallest Galois-lonely prime, resonance structure,
adversarial search for the largest smallest-Galois-lonely-prime.
"""
from lrc import *
from itertools import combinations, product

# ---------- Investigation 3 helper: small integer relations sum c_i v_i = 0 ----------
def small_relations(speeds, cmax=3):
    """Find nonzero integer vectors c with |c_i|<=cmax and sum c_i v_i = 0.
    Returns list of (c-tuple) reduced (primitive, first nonzero positive)."""
    m = len(speeds)
    found = set()
    rng = range(-cmax, cmax+1)
    for c in product(rng, repeat=m):
        if all(x == 0 for x in c):
            continue
        if sum(ci*vi for ci, vi in zip(c, speeds)) == 0:
            # normalize sign: first nonzero positive
            cc = c
            for x in cc:
                if x != 0:
                    if x < 0:
                        cc = tuple(-y for y in cc)
                    break
            # primitivity
            g = 0
            for x in cc:
                g = gcd(g, abs(x))
            if g > 1:
                cc = tuple(x//g for x in cc)
            found.add(cc)
    return sorted(found, key=lambda t: (sum(abs(x) for x in t), t))


def min_relation_weight(speeds, cmax=4):
    rels = small_relations(speeds, cmax)
    if not rels:
        return None, None
    best = min(rels, key=lambda t: sum(abs(x) for x in t))
    return sum(abs(x) for x in best), best


# ---------- Investigation 2 & 4 ----------
def report_smallest(label, n, sets, pmax=4000):
    m = n-1
    print(f"\n### {label}  (smallest prime q with O_q cap B_q^m nonempty)")
    print(f"    {'speeds':<22} {'q_min':>6}   min-relation (|c|<=4)")
    rows = []
    for sp in sets:
        q = smallest_galois_lonely_prime(sp, n, pmax)
        w, rel = min_relation_weight(sp, 4)
        relstr = f"{rel} (wt {w})" if rel else "none"
        print(f"    {str(sp):<22} {str(q):>6}   {relstr}")
        rows.append((sp, q, w, rel))
    return rows


if __name__ == "__main__":
    print("="*78)
    print("INVESTIGATION 2: SMALLEST GALOIS-LONELY PRIME")
    print("="*78)

    # AP and structured m=3
    report_smallest("m=3,n=4: APs and structured", 4, [
        [1,2,3],[1,2,4],[1,2,5],[1,3,5],[1,2,6],[2,3,4],[3,4,5],
        [1,4,6],[2,4,6],[1,5,7],[2,3,5],[1,3,9],[1,2,8],
    ])

    report_smallest("m=4,n=5: APs and structured", 5, [
        [1,2,3,4],[1,2,4,8],[1,3,5,7],[2,3,5,7],[1,2,3,5],
        [1,2,4,6],[1,4,9,16],[2,4,6,8],[1,5,9,13],[1,2,5,11],
    ])

    report_smallest("m=5,n=6: APs and structured", 6, [
        [1,2,3,4,5],[1,2,4,8,16],[1,3,5,7,9],[2,3,5,7,11],
        [1,2,3,5,8],[1,2,4,6,8],[2,4,6,8,10],[1,5,9,13,17],
    ])

    print("\n" + "="*78)
    print("INVESTIGATION 3 + 4: ADVERSARIAL SEARCH for LARGEST q_min")
    print("="*78)

    # Search many sets, find those with largest smallest-galois-lonely-prime.
    def adversarial(n, gen, topk=12, pmax=6000):
        m = n-1
        results = []
        for sp in gen:
            q = smallest_galois_lonely_prime(sp, n, pmax)
            qv = q if q is not None else 10**9
            w, rel = min_relation_weight(sp, 4)
            results.append((qv, sp, q, w, rel))
        results.sort(reverse=True)
        print(f"\n--- {n=}, top {topk} resistant sets (largest q_min) ---")
        print(f"    {'speeds':<24} {'q_min':>7}  min-rel")
        for qv, sp, q, w, rel in results[:topk]:
            relstr = f"{rel} (wt {w})" if rel else "none"
            qs = str(q) if q is not None else f">{pmax}"
            print(f"    {str(sp):<24} {qs:>7}  {relstr}")
        return results

    # m=3 brute over speed sets up to 20
    gen3 = (list(c) for c in combinations(range(1,21), 3))
    adversarial(4, gen3, topk=15, pmax=4000)

    # m=4 over speeds up to 14
    gen4 = (list(c) for c in combinations(range(1,15), 4))
    adversarial(5, gen4, topk=15, pmax=4000)

    # m=5 over speeds up to 12
    gen5 = (list(c) for c in combinations(range(1,13), 5))
    adversarial(6, gen5, topk=15, pmax=4000)
