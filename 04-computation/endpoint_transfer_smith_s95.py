#!/usr/bin/env python3
"""
endpoint_transfer_smith_s95.py

Smith normal form probe for endpoint-transfer matrices.

The GF(2) rank data asks whether row information survives mod 2.  Smith
normal form asks the stronger integral question: which invariant factors carry
2-adic torsion?

This script computes tournament and merged Smith forms through transition
5->6, where those matrices are still small.  It also computes the even-graph
Smith form through 6->7, where the rank defect is largest but the matrix is
still only 16x54.
"""

import os
import sys
from collections import Counter

import sympy as sp
from sympy.matrices.normalforms import smith_normal_form
from sympy.polys.domains import ZZ

sys.path.insert(0, os.path.dirname(__file__))

import endpoint_transfer_bucket_recursion_s95 as tourn
import even_graph_endpoint_transfer_s95 as even


def matrix_from_rows(rows, ncols):
    return sp.Matrix([[row.get(j, 0) for j in range(ncols)] for row in rows])


def v2(x):
    x = abs(int(x))
    if x == 0:
        return None
    c = 0
    while x % 2 == 0:
        c += 1
        x //= 2
    return c


def smith_summary(M):
    S = smith_normal_form(M, domain=ZZ)
    diag = []
    for i in range(min(S.rows, S.cols)):
        val = int(S[i, i])
        if val != 0:
            diag.append(abs(val))
    v2s = [v2(x) for x in diag]
    return {
        "rank_Z": len(diag),
        "diag": diag,
        "v2": v2s,
        "v2_spectrum": Counter(v2s),
        "odd_factors": sum(1 for x in diag if x % 2 == 1),
    }


def print_summary(label, M):
    s = smith_summary(M)
    print(f"  {label}: shape={M.rows}x{M.cols}, rank_Z={s['rank_Z']}, odd_factors={s['odd_factors']}")
    print(f"    invariant factors={s['diag']}")
    print(f"    v2 spectrum={dict(sorted(s['v2_spectrum'].items()))}")
    return s


def main():
    print("=" * 78)
    print("ENDPOINT TRANSFER SMITH S95")
    print("=" * 78)
    print("Exact Smith normal forms for small endpoint-transfer matrices.")

    t_levels = {n: tourn.build_level(n) for n in range(2, 7)}
    e_levels = {n: even.build_even_level(n) for n in range(2, 8)}

    sequences = {
        "t_odd_factors": [],
        "m_odd_factors": [],
        "e_odd_factors": [],
        "e_even_factors": [],
    }

    for n in range(2, 6):
        print("\n" + "-" * 78)
        print(f"n={n}->{n+1}")
        td = tourn.analyze_transfer(t_levels[n], t_levels[n + 1])
        ed = even.analyze_transfer(e_levels[n], e_levels[n + 1])

        tM = matrix_from_rows(td["class_rows"], len(t_levels[n + 1]["cans"]))
        mM = matrix_from_rows(td["merged_rows"], len(t_levels[n + 1]["merged"]))
        eM = matrix_from_rows(ed["rows"], len(e_levels[n + 1]["classes"]))

        ts = print_summary("tournament", tM)
        ms = print_summary("merged", mM)
        es = print_summary("even graph", eM)

        sequences["t_odd_factors"].append(ts["odd_factors"])
        sequences["m_odd_factors"].append(ms["odd_factors"])
        sequences["e_odd_factors"].append(es["odd_factors"])
        sequences["e_even_factors"].append(es["rank_Z"] - es["odd_factors"])

    print("\n" + "-" * 78)
    print("n=6->7, even graph only")
    ed = even.analyze_transfer(e_levels[6], e_levels[7])
    eM = matrix_from_rows(ed["rows"], len(e_levels[7]["classes"]))
    es = print_summary("even graph", eM)
    sequences["e_odd_factors"].append(es["odd_factors"])
    sequences["e_even_factors"].append(es["rank_Z"] - es["odd_factors"])

    print("\n" + "=" * 78)
    print("SEQUENCES")
    print("=" * 78)
    for key, value in sequences.items():
        print(f"  {key}: {value}")

    print("\n" + "=" * 78)
    print("INTERPRETATION")
    print("=" * 78)
    print("""
1. The number of odd Smith factors equals the GF(2) row rank.
   Full GF(2) row rank is equivalent to all row-rank invariant factors being
   odd.

2. The tournament matrices through 5->6 have only odd Smith factors.  This is
   stronger than mod-2 rank: the endpoint transfer has no visible 2-primary
   cokernel obstruction in these small cases.

3. The even-graph matrices have even Smith factors exactly where parity rank is
   lost.  The rank defect is therefore an integral 2-adic torsion phenomenon,
   not just an artifact of reducing entries modulo 2.
""")


if __name__ == "__main__":
    main()
