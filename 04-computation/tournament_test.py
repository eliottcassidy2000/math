#!/usr/bin/env python3
"""
tournament_test.py — Command-line tool for tournament ranking significance.

Usage:
  python3 tournament_test.py --n 7 --h 189
  python3 tournament_test.py --n 20
  python3 tournament_test.py --matrix "1,1,1,0;0,1,1,0;0,0,1,0;1,1,1,1"

Outputs:
  Expected H, CV, sigma, Z-score, and interpretation.
"""

import argparse
from math import sqrt, factorial
from itertools import permutations

def E_H(n):
    return factorial(n) / 2**(n-1)

def cv(n):
    return sqrt(2.0/n)

def H_from_matrix(T, n):
    count = 0
    for perm in permutations(range(n)):
        ok = True
        for i in range(n-1):
            if T[perm[i]][perm[i+1]] != 1:
                ok = False
                break
        if ok:
            count += 1
    return count

def interpret(z):
    if z > 3: return "VERY SIGNIFICANT — clear skill hierarchy"
    if z > 2: return "SIGNIFICANT — ranking is real"
    if z > 1: return "MODERATE — some evidence of ranking"
    if z > -1: return "NOT SIGNIFICANT — consistent with random"
    if z > -2: return "BELOW AVERAGE — fewer rankings than expected"
    return "VERY LOW — highly contradictory results"

def main():
    parser = argparse.ArgumentParser(description="Tournament ranking significance test")
    parser.add_argument("--n", type=int, help="Number of teams/items")
    parser.add_argument("--h", type=int, help="Hamiltonian path count (if known)")
    parser.add_argument("--matrix", type=str, help="Adjacency matrix (semicolon-separated rows, comma-separated entries)")
    parser.add_argument("--info", action="store_true", help="Print formula information")
    args = parser.parse_args()

    if args.info:
        print("Tournament Ranking Significance Test")
        print("Based on: CV^2(H) = 2/n (Cayley-Delannoy formula)")
        print()
        print("H(T) = number of directed Hamiltonian paths in tournament T")
        print("E[H] = n! / 2^{n-1} for a random tournament")
        print("CV = sqrt(2/n) = coefficient of variation")
        print("Z = (H - E[H]) / (E[H] * CV) = significance score")
        print()
        print("Z > 2: ranking is statistically significant")
        print("Z ~ 0: consistent with random coin flips")
        print("Z < -2: results are more contradictory than random")
        return

    if args.matrix:
        rows = args.matrix.split(";")
        n = len(rows)
        T = []
        for row in rows:
            T.append([int(x) for x in row.split(",")])
        H = H_from_matrix(T, n)
        args.n = n
        args.h = H
        print(f"Tournament matrix ({n}x{n}): H(T) = {H}")

    if not args.n:
        parser.print_help()
        return

    n = args.n
    mu = E_H(n)
    c = cv(n)
    sigma = mu * c

    print(f"Tournament Ranking Test (n={n} teams)")
    print("=" * 40)
    print(f"E[H]  = {mu:.2f}")
    print(f"CV    = sqrt(2/{n}) = {c:.4f}")
    print(f"sigma = {sigma:.2f}")

    if args.h is not None:
        H = args.h
        z = (H - mu) / sigma
        print(f"H(T)  = {H}")
        print(f"H/E[H]= {H/mu:.3f}")
        print(f"Z     = {z:+.2f}")
        print(f"Result: {interpret(z)}")
    else:
        print()
        print("Significance thresholds:")
        for sig_level, z_val in [("95% (Z>1.65)", 1.65), ("99% (Z>2.33)", 2.33), ("99.9% (Z>3.09)", 3.09)]:
            h_threshold = mu + z_val * sigma
            print(f"  {sig_level}: H > {h_threshold:.0f} (H/E[H] > {h_threshold/mu:.3f})")

if __name__ == "__main__":
    main()
