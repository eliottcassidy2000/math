#!/usr/bin/env python3
"""Quick test of delta_I = 2*(gained-lost) at n=6."""
import sys
sys.path.insert(0, r"C:\Users\Eliott\Documents\GitHub\math\03-artifacts\code")
from tournament_lib import *
import random

def flip_arc(T, i, j):
    Tp = [row[:] for row in T]
    Tp[i][j] = 1 - T[i][j]
    Tp[j][i] = 1 - T[j][i]
    return Tp

def compute_I(T):
    cycles = find_odd_cycles(T)
    if not cycles:
        return 1
    cg = conflict_graph(cycles)
    return independence_poly_at(cg, 2)

def uses_arc(cyc, a, b):
    L = len(cyc)
    return any(cyc[k] == a and cyc[(k+1) % L] == b for k in range(L))

rng = random.Random(42)

for n in [6]:
    print(f"=== n={n} ===")
    ok_I = ok_H = 0
    fail_I = fail_H = 0
    total = 0
    fail_examples = []

    for trial in range(100):
        T = random_tournament(n, rng)
        ht = hamiltonian_path_count(T)
        it = compute_I(T)
        cyc_T = find_odd_cycles(T)

        # Pick one random arc to flip
        i, j = rng.sample(range(n), 2)
        if not T[i][j]:
            i, j = j, i

        Tp = flip_arc(T, i, j)
        htp = hamiltonian_path_count(Tp)
        itp = compute_I(Tp)
        cyc_Tp = find_odd_cycles(Tp)

        lost = sum(1 for c in cyc_T if uses_arc(c, i, j))
        gained = sum(1 for c in cyc_Tp if uses_arc(c, j, i))
        predicted = 2 * (gained - lost)

        delta_I = itp - it
        delta_H = htp - ht

        total += 1
        if delta_I == predicted:
            ok_I += 1
        else:
            fail_I += 1
            if len(fail_examples) < 3:
                fail_examples.append(f"  dI={delta_I} pred={predicted} lost={lost} gained={gained}")

        if delta_H == predicted:
            ok_H += 1
        else:
            fail_H += 1

        if trial % 20 == 0:
            print(f"  {trial}/100: I formula {ok_I}/{total} pass, H formula {ok_H}/{total} pass")

    print(f"FINAL n={n}: {total} tests")
    print(f"  delta_I = 2*(gained-lost): {ok_I}/{total} ({fail_I} fail)")
    print(f"  delta_H = 2*(gained-lost): {ok_H}/{total} ({fail_H} fail)")
    if fail_examples:
        print("  Failure examples (delta_I):")
        for ex in fail_examples:
            print(ex)
