#!/usr/bin/env python3
"""Small fresh audit of the bounded-defect mechanism, with no producer imports."""
import itertools as it
import argparse
import gzip
import json
import math
from collections import Counter
from pathlib import Path

HERE = Path(__file__).resolve().parent


def matrix(n, mask):
    A = [[0]*n for _ in range(n)]
    for bit,(i,j) in enumerate(it.combinations(range(n),2)):
        A[i][j] = (mask >> bit) & 1
        A[j][i] = 1-A[i][j]
    return A


def strong(A):
    n = len(A)
    reach = [sum(value << j for j,value in enumerate(row)) | (1 << i)
             for i,row in enumerate(A)]
    for k in range(n):
        for i in range(n):
            if reach[i] & (1 << k):
                reach[i] |= reach[k]
    return all(row == (1 << n)-1 for row in reach)


def check(n, mask):
    A = matrix(n,mask)
    assert strong(A)
    perms = list(it.permutations(range(n)))
    aut = [p for p in perms if all(A[i][j] == A[p[i]][p[j]]
                                  for i in range(n) for j in range(i+1,n))]
    paths = [p for p in perms if all(A[i][j] for i,j in zip(p,p[1:]))]
    blue = set()
    for p in paths:
        r = dict(zip(p,p[::-1]))
        if all(A[i][j] == A[r[j]][r[i]]
               for i in range(n) for j in range(i+1,n)):
            blue.add(p)
    assert len(paths) % len(aut) == len(blue) % len(aut) == 0
    Delta = (len(paths)-len(blue)) // len(aut)
    cycles = [p for p in paths if p[0] == 0 and A[p[-1]][0]]
    pure = black = 0
    minout = min(map(sum,A))
    cycle_profiles = Counter()
    for p in cycles:
        rotations = []
        blue_cuts = 0
        for j in range(n):
            cut = p[j:] + p[:j]
            rotation = dict(zip(p,cut))
            if tuple(rotation[v] for v in range(n)) in aut:
                rotations.append(rotation)
            blue_cuts += cut in blue
        group_order = len(rotations)
        q = n//group_order
        assert n % group_order == 0 and group_order % 2 == 1
        assert blue_cuts % group_order == 0
        blue_orbits = blue_cuts // group_order
        assert blue_orbits in (0, math.gcd(2,q))
        assert q-blue_orbits <= Delta and q <= Delta+2
        assert (2*minout+1)*q >= n
        for v in range(n):
            orbit = {r[v] for r in rotations}
            assert len(orbit) == group_order
            assert sum(A[v][w] for w in orbit) == (group_order-1)//2
        if blue_cuts == n:
            pure += 1
            square = dict(zip(p,p[2:] + p[:2]))
            assert tuple(square[v] for v in range(n)) in aut
        else:
            black += 1
        cycle_profiles[str((q,blue_orbits))] += 1
    assert cycles and pure+black == len(cycles)
    assert black <= len(paths)-len(blue) == Delta*len(aut)
    assert pure <= (1 if n%2 else n//2)*len(aut)
    assert len(cycles) <= (n+Delta)*len(aut)
    assert (2*minout+1)*(Delta+2) >= n
    return dict(n=n,mask=mask,H=len(paths),aut=len(aut),b=len(blue)//len(aut),
                Delta=Delta,HC=len(cycles),all_blue_cycles=pure,
                black_cycles=black,minout=minout,cycle_profiles=dict(cycle_profiles))


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--source", type=Path,
        default=HERE.parent / "independent/verification_results.json.gz",
        help="Existing independent census JSON or .json.gz; used only to select NS6 masks.")
    args = parser.parse_args()
    totals = []
    for n in (3,4,5):
        records = [check(n,mask) for mask in range(1 << math.comb(n,2))
                   if strong(matrix(n,mask))]
        totals.append(dict(n=n, strong=len(records),
                           cycles=sum(r["HC"] for r in records)))
    # One n=6 non-self-converse representative from each reported invariant
    # profile. Source values choose masks only; all invariants are recomputed.
    source = args.source
    contents = gzip.decompress(source.read_bytes()).decode() if source.suffix == ".gz" else source.read_text()
    records = json.loads(contents)["tournaments"]["individual_cases"]
    samples = {}
    for r in records:
        if r["n"] == 6 and not r["self_converse"]:
            samples.setdefault((r["H"],r["A"],r["HC"]),r["mask"])
    checked_ns = [check(6,mask) for mask in samples.values()]
    assert all(r["b"] == 0 for r in checked_ns)
    # Fresh square-root fiber enumeration via rooted cycle words.
    roots = []
    for n in range(2,10):
        counts = Counter()
        for tail in it.permutations(range(1,n)):
            p = (0,)+tail
            sq = dict(zip(p,p[2:]+p[:2]))
            counts[tuple(sq[v] for v in range(n))] += 1
        wanted = 1 if n%2 else n//2
        assert set(counts.values()) == {wanted}
        roots.append(dict(n=n,full_cycles=sum(counts.values()),
                          square_images=len(counts),fiber_size=wanted))
    out = dict(all_strong_through5=totals, ns6_profile_controls=checked_ns,
               full_cycle_square_roots=roots)
    (HERE/"verification.json").write_text(json.dumps(out,indent=2)+"\n")
    print(json.dumps(dict(all_strong_through5=totals,
                          ns6_profiles=len(checked_ns),
                          full_cycle_square_roots=roots),indent=2))
    print("ALL FRESH BOUNDED-DEFECT CHECKS PASSED")

if __name__ == "__main__":
    main()
